# ============================================================================
#  歪対称行列に基づく非対称MDS：APM / ADD / TPD
#
#  前提：data/trade_2023_comtrade.csv（行=輸出国、列=輸出先、単位=10億米ドル）
#  必要：install.packages(c("smacof", "vegan", "ggplot2", "ggrepel"))
#
#  構成：前半が関数（パッケージ化するときは R/ に移す）、後半が実行部。
#        入力は非対称行列 M（対角は無視）。ここでは輸出額の対数をそのまま渡す。
# ============================================================================

library(smacof)
library(vegan)
library(ggplot2)
library(ggrepel)


# ----------------------------------------------------------------------------
#  関数
# ----------------------------------------------------------------------------

# 対称成分 S と歪対称成分 A。対角は無視する（何が入っていても A の対角は 0）。
skew_decompose <- function(M) {
  M <- as.matrix(M)
  diag(M) <- 0
  list(S = (M + t(M)) / 2, A = (M - t(M)) / 2)
}

# A から APM・ADD・TPD の3行列と、d_APM^2 に占める 2*delta_ADD^2 の割合を返す。
apm_matrices <- function(M) {
  A   <- skew_decompose(M)$A
  APM <- as.matrix(dist(t(A)))                 # 列ベクトル間ユークリッド距離
  ADD <- abs(A)
  # TPD は定義どおり k ≠ i, j の和として独立に計算する（恒等式から逆算すると検算にならない）
  n   <- nrow(A)
  TPD <- matrix(0, n, n, dimnames = dimnames(A))
  for (i in 1:n) for (j in 1:n) if (i != j) {
    k <- setdiff(1:n, c(i, j))
    TPD[i, j] <- sqrt(sum((A[k, i] - A[k, j])^2))
  }
  stopifnot(all(abs(APM^2 - (2 * ADD^2 + TPD^2)) < 1e-10))   # 式(7)の検算
  share <- 2 * A^2 / APM^2
  diag(share) <- NA
  list(A = A, APM = APM, ADD = ADD, TPD = TPD, share = share)
}

# 3つの布置。APM は cmdscale、ADD・TPD は smacof。
# ADD・TPD の布置は vegan::procrustes（回転・平行移動のみ、拡大縮小なし）で APM に揃える。
apm_mds <- function(M, ndim = 2, seed = 123) {
  mats <- apm_matrices(M)
  nm   <- rownames(mats$A)

  fit_apm <- cmdscale(as.dist(mats$APM), k = ndim, eig = TRUE)
  set.seed(seed); fit_add <- mds(as.dist(mats$ADD), ndim = ndim, type = "ratio")
  set.seed(seed); fit_tpd <- mds(as.dist(mats$TPD), ndim = ndim, type = "ratio")

  X_apm <- scale(fit_apm$points, scale = FALSE)
  X_add <- procrustes(X_apm, fit_add$conf, scale = FALSE)$Yrot
  X_tpd <- procrustes(X_apm, fit_tpd$conf, scale = FALSE)$Yrot

  as_df <- function(X, method)
    data.frame(method = method, country = nm, Dim1 = X[, 1], Dim2 = X[, 2],
               row.names = NULL)
  conf <- rbind(as_df(X_apm, "APM"), as_df(X_add, "ADD"), as_df(X_tpd, "TPD"))
  conf$method <- factor(conf$method, levels = c("APM", "ADD", "TPD"))

  ev <- fit_apm$eig
  list(matrices = mats, conf = conf, eig = ev,
       gof = sum(ev[1:ndim]) / sum(abs(ev)),
       stress = c(ADD = fit_add$stress, TPD = fit_tpd$stress),
       fits = list(apm = fit_apm, add = fit_add, tpd = fit_tpd))
}

# 布置の図。methods で1枚ずつも、3枚並べても出せる。
plot_apm <- function(fit, methods = c("APM", "ADD", "TPD"), seed = 123) {
  df <- subset(fit$conf, method %in% methods)
  ggplot(df, aes(Dim1, Dim2, label = country)) +
    geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_point(size = 3.2, colour = "steelblue") +
    geom_text_repel(size = 4.2, seed = seed) +
    facet_wrap(~ method, nrow = 1) +
    coord_equal() +
    labs(x = "Dimension 1", y = "Dimension 2") +
    theme_minimal(base_size = 12)
}

# 表3：全ペアの APM・ADD・TPD と割合（%）。割合の大きい順。
# a は符号つきの a_ij（pair の左が i、右が j）。larger は m_ij と m_ji の大きい向きで、
# "X → Y" は m_XY > m_YX を表す。ペアを書く順に依らず向きが読める。
pair_table <- function(mats) {
  nm <- rownames(mats$A)
  u  <- which(upper.tri(mats$APM), arr.ind = TRUE)
  a  <- mats$A[u]
  tbl <- data.frame(pair   = paste(nm[u[, 1]], nm[u[, 2]], sep = "-"),
                    a      = a,
                    larger = ifelse(a > 0, paste(nm[u[, 1]], "→", nm[u[, 2]]),
                             ifelse(a < 0, paste(nm[u[, 2]], "→", nm[u[, 1]]), "均衡")),
                    APM    = mats$APM[u],
                    ADD    = mats$ADD[u],
                    TPD    = mats$TPD[u],
                    share  = 100 * mats$share[u])
  tbl[order(-tbl$share), ]
}

# 三角不等式の違反（順序三つ組）と、最悪の組。
tri_violations <- function(M) {
  n <- nrow(M); cnt <- 0; gmax <- 0; wi <- wj <- wl <- NA
  for (i in 1:n) for (j in 1:n) for (l in 1:n) {
    if (i != j && j != l && i != l) {
      gap <- M[i, j] - (M[i, l] + M[l, j])
      if (gap > 1e-12) {
        cnt <- cnt + 1
        if (gap > gmax) { gmax <- gap; wi <- i; wj <- j; wl <- l }
      }
    }
  }
  list(count = cnt, i = wi, j = wj, l = wl)
}

# 二重中心化行列の固有値。最小値が負ならユークリッド距離行列でない。
dc_eigen <- function(M) {
  n <- nrow(M); J <- diag(n) - matrix(1 / n, n, n)
  eigen(-0.5 * J %*% (M^2) %*% J, symmetric = TRUE, only.values = TRUE)$values
}


# ----------------------------------------------------------------------------
#  実行部
# ----------------------------------------------------------------------------

# 1. データ
trade <- read.csv("data/trade_2023_comtrade.csv", fileEncoding = "UTF-8-BOM")
rownames(trade) <- trade$country
T  <- data.matrix(trade[, -1])
nm <- rownames(T)
n  <- nrow(T)
print(T)

# 2. 輸出額の対数をそのまま非対称行列として渡す（非類似度への反転はしない）。
#    a_ij = (1/2) log(T_ij / T_ji)。符号を反転しても A の符号が変わるだけで、
#    APM・ADD・TPD はいずれも不変。
M <- log(T)
diag(M) <- 0
print(round(M, 3))

# 3. 分析
fit  <- apm_mds(M, ndim = 2)
mats <- fit$matrices

print(round(mats$A, 3))
print(round(mats$APM, 3))
print(round(mats$ADD, 3))
print(round(mats$TPD, 3))
print(round(100 * mats$share, 1))

# 4. 適合度
cat("APM 2次元説明率 :", round(100 * fit$gof, 1), "%\n")
cat("APM 最小固有値  :", round(min(fit$eig), 6), "\n")
cat("ADD stress-1    :", round(fit$stress["ADD"], 4), "\n")
cat("TPD stress-1    :", round(fit$stress["TPD"], 4), "\n")

# 5. 図1〜3
print(plot_apm(fit, "APM"))
print(plot_apm(fit, "ADD"))
print(plot_apm(fit, "TPD"))
print(plot_apm(fit))                       # 3枚並べ

# 6. 表3と要約
tbl <- pair_table(mats)
print(format(tbl, digits = 3, nsmall = 1), row.names = FALSE)

u <- upper.tri(mats$APM)
sh <- 100 * mats$share
cat("2δ²ADD/d²APM の割合 : 最小", round(min(sh[u]), 1), "%  最大", round(max(sh[u]), 1),
    "%  平均", round(mean(sh[u]), 1), "%  中央値", round(median(sh[u]), 1), "%\n")
cat("相関 ADD-APM :", round(cor(mats$ADD[u], mats$APM[u]), 3), "\n")
cat("相関 TPD-APM :", round(cor(mats$TPD[u], mats$APM[u]), 3), "\n")
cat("相関 ADD-TPD :", round(cor(mats$ADD[u], mats$TPD[u]), 3), "\n")

# 7. 三角不等式とユークリッド性（別の性質）
cat("\n--- 三角不等式とユークリッド性 ---\n")
for (lbl in c("APM", "ADD", "TPD")) {
  D <- mats[[lbl]]; v <- tri_violations(D); ev <- dc_eigen(D)
  cat(sprintf("%s : 三角不等式違反 %d 組 / 最小固有値 %.4f（最大固有値比 %.1f%%）\n",
              lbl, v$count, min(ev), 100 * abs(min(ev)) / max(ev)))
  if (v$count > 0)
    cat(sprintf("     最悪 : %s-%s = %.3f > %s-%s + %s-%s = %.3f\n",
                nm[v$i], nm[v$j], D[v$i, v$j], nm[v$i], nm[v$l], nm[v$l], nm[v$j],
                D[v$i, v$l] + D[v$l, v$j]))
}

# 8. ADD の適合の悪さは次元不足ではない
cat("\n--- ADD の次元別 stress-1 ---\n")
for (k in 2:(n - 1)) {
  set.seed(123)
  cat(sprintf("  %d次元 : %.4f\n", k, mds(as.dist(mats$ADD), ndim = k, type = "ratio")$stress))
}

# 9. APM の CMDS は A の列の中心化主成分分析（B = J A'A J）
A <- mats$A
J <- diag(n) - matrix(1 / n, n, n)
cat("\n--- APM の構造 ---\n")
cat("  B と J(A'A)J の最大差 :",
    format(max(abs(-0.5 * J %*% (mats$APM^2) %*% J - J %*% (t(A) %*% A) %*% J)), digits = 3), "\n")
cat("  A の特異値（対で現れる）:", paste(round(svd(A)$d, 3), collapse = " "), "\n")
cat("  A の行平均（正なら相手全体に対して輸出超過）\n")
print(round(rowMeans(A), 3))

# 10. TPD 布置の重心からの距離
X_tpd <- as.matrix(subset(fit$conf, method == "TPD")[, c("Dim1", "Dim2")])
cat("\n--- TPD 布置の重心からの距離 ---\n")
print(round(sort(setNames(sqrt(rowSums(scale(X_tpd, scale = FALSE)^2)), nm)), 3))
