# ============================================================================
#  歪対称行列に基づく非対称MDS：APM / ADD / TPD の3手法
#
#  前提：data/trade_2023_comtrade.csv（このリポジトリに同梱）
#        （UN Comtrade 2023年、行=輸出国、列=輸出先、単位=10億米ドル）
#  必要：install.packages(c("ggplot2", "ggrepel", "smacof"))
# ============================================================================

library(ggplot2)
library(ggrepel)
library(smacof)


# ============================================================================
#  1. データを行列にする
# ============================================================================
# シートは1列目が国名、2列目以降が各輸出先への輸出額。
# data.matrix で数値部分だけを行列にし、行名を国名にする。

trade <- read.csv("data/trade_2023_comtrade.csv", fileEncoding = "UTF-8-BOM")
rownames(trade) <- trade$country
T <- data.matrix(trade[, -1])   # 1列目（国名）を除く

diag(T) <- NA                   # 自国から自国はデータが無い

nm <- rownames(T)               # 国名。図のラベルに使う
n  <- nrow(T)                   # 対象数。ここでは 8

print(T)


# ============================================================================
#  2. 輸出額を非類似度に変換する
# ============================================================================
# 輸出額が大きいほど「近い」ので、対数をとって符号を反転する。
#     δ_ij = -log(T_ij)
#
# こうすると歪対称成分の要素が
#     a_ij = (1/2) * log(T_ji / T_ij)
# になる。二方向の貿易額の対数比の半分で、解釈しやすい。
#
# 対角は 0 にしておく。ただし a_ii = (δ_ii - δ_ii)/2 = 0 なので、
# 実は対角に何を入れても A の対角は 0 になる。

Delta <- -log(T)
diag(Delta) <- 0

print(round(Delta, 3))


# ============================================================================
#  3. 対称成分 S と歪対称成分 A に分ける
# ============================================================================
#     S = (Δ + Δ') / 2   … 二方向の平均。今回は使わないが確認用に作る
#     A = (Δ - Δ') / 2   … 二方向の差。これが本研究の対象

S <- (Delta + t(Delta)) / 2
A <- (Delta - t(Delta)) / 2

# 歪対称になっているか目で確かめる
cat("A + t(A) がゼロか :", all(abs(A + t(A)) < 1e-12), "\n")
cat("A の対角がゼロか   :", all(abs(diag(A)) < 1e-12), "\n")

print(round(A, 3))


# ============================================================================
#  4. APM を計算する
# ============================================================================
#     d_APM^2(i,j) = sum_{k=1}^{n} (a_ki - a_kj)^2
#
# A の第 i 列と第 j 列のユークリッド距離。
# 「対象 i が他の全対象に対してもつ方向差」を並べたベクトル同士を比べている。
# dist(t(A)) 一行でも書けるが、何をしているか見えるように二重ループで書く。

d_APM2 <- matrix(0, n, n, dimnames = list(nm, nm))

for (i in 1:n) {
  for (j in 1:n) {
    if (i != j) {
      d_APM2[i, j] <- sum((A[, i] - A[, j])^2)   # A[, i] は第 i 列
    }
  }
}

APM <- sqrt(d_APM2)

print(round(APM, 3))


# ============================================================================
#  5. ADD を計算する
# ============================================================================
#     δ_ADD^2(i,j) = a_ij^2
#
# 対象ペア自身の非対称量。A の各要素を二乗するだけ。第三者は関与しない。

delta_ADD2 <- A^2
ADD  <- sqrt(delta_ADD2)      # abs(A) と同じ

print(round(ADD, 3))


# ============================================================================
#  6. TPD を計算する
# ============================================================================
#     δ_TPD^2(i,j) = sum_{k != i,j} (a_ki - a_kj)^2
#
# APM の和から k = i と k = j の項を除いたもの。
# setdiff(1:n, c(i,j)) が「i でも j でもない k」の並びを作る。

delta_TPD2 <- matrix(0, n, n, dimnames = list(nm, nm))

for (i in 1:n) {
  for (j in 1:n) {
    if (i != j) {
      k <- setdiff(1:n, c(i, j))          # 第三者のインデックス（n-2 個）
      delta_TPD2[i, j] <- sum((A[k, i] - A[k, j])^2)
    }
  }
}

TPD <- sqrt(delta_TPD2)

print(round(TPD, 3))


# ============================================================================
#  7. 分解の検算
# ============================================================================
# 命題：d_APM^2 = 2 * δ_ADD^2 + δ_TPD^2
# ここがずれたらどこかで間違えている。

cat("分解定理が成立するか :", all(abs(d_APM2 - (2 * delta_ADD2 + delta_TPD2)) < 1e-10), "\n")

# 各ペアで直接項が占める割合（本文の表3になる）
share <- 2 * delta_ADD2 / d_APM2
diag(share) <- NA
print(round(share * 100, 1))


# ============================================================================
#  8. APM の布置（古典的MDS）
# ============================================================================
# APM は R^n 上のベクトル間のユークリッド距離なので、必ずユークリッド距離行列。
# したがって cmdscale で問題ない。eig = TRUE で固有値が返る。

fit_apm <- cmdscale(as.dist(APM), k = 2, eig = TRUE)
X_apm   <- fit_apm$points          # 8行2列の座標

ev <- fit_apm$eig
cat("APM 2次元説明率 :", round(sum(ev[1:2]) / sum(abs(ev)) * 100, 1), "%\n")
cat("APM 最小固有値  :", round(min(ev), 6), " （0以上ならユークリッド）\n")


# ============================================================================
#  9. ADD の布置（SMACOF）
# ============================================================================
# ADD はユークリッド距離行列である保証が無い。cmdscale だと負の固有値が出て解が歪むので、
# ストレスを直接下げる SMACOF を使う。
# type = "ratio" は「非類似度と距離が比例する」という一番素直な当てはめ。
# 初期値で結果が変わりうるので set.seed で固定しておく。

set.seed(123)
fit_add <- mds(as.dist(ADD), ndim = 2, type = "ratio")
X_add   <- fit_add$conf

cat("ADD stress-1 :", round(fit_add$stress, 4), "\n")


# ============================================================================
#  10. TPD の布置（SMACOF）
# ============================================================================
# ADD と同じ手順。TPD もユークリッド距離行列である保証は無い。

set.seed(123)
fit_tpd <- mds(as.dist(TPD), ndim = 2, type = "ratio")
X_tpd   <- fit_tpd$conf

cat("TPD stress-1 :", round(fit_tpd$stress, 4), "\n")


# ============================================================================
#  11. 3枚の向きを揃える（プロクラステス回転）
# ============================================================================
# MDS の布置は回転・鏡映しても距離が変わらないので、向き自体に意味は無い。
# そのままだと3枚がバラバラの向きで出て見比べられないので、
# APM を基準に ADD と TPD を回して重ねる。
#
#   1) 両方を中心化する（重心を原点へ）
#   2) t(X) %*% ref を特異値分解する
#   3) u %*% t(v) が最適な回転行列になる
#   4) X にそれを掛ける
#
# 回転と鏡映だけで伸縮はしない。だから布置内の距離関係は変わらない。

X_apm <- scale(X_apm, center = TRUE, scale = FALSE)   # 基準を中心化

# --- ADD を APM に合わせる ---
X_add <- scale(X_add, center = TRUE, scale = FALSE)
sv    <- svd(t(X_add) %*% X_apm)
R     <- sv$u %*% t(sv$v)
X_add <- X_add %*% R

# --- TPD を APM に合わせる（同じことをもう一度書く）---
X_tpd <- scale(X_tpd, center = TRUE, scale = FALSE)
sv    <- svd(t(X_tpd) %*% X_apm)
R     <- sv$u %*% t(sv$v)
X_tpd <- X_tpd %*% R


# ============================================================================
#  12. 作図
# ============================================================================
# 3手法とも同じ作りなので3回書く。coord_equal() は MDS 図では必須。
# 縦横の縮尺が違うと距離が正しく見えなくなるため。

# --- 図2：APM-MDS ---
df_apm <- data.frame(country = nm, Dim1 = X_apm[, 1], Dim2 = X_apm[, 2])

p_apm <- ggplot(df_apm, aes(Dim1, Dim2, label = country)) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_point(size = 3.2, colour = "steelblue") +
  geom_text_repel(size = 4.2, seed = 123) +
  coord_equal() +
  labs(title = "APM-MDS", x = "Dimension 1", y = "Dimension 2") +
  theme_minimal(base_size = 12)

print(p_apm)

# --- 図3：ADD-MDS ---
df_add <- data.frame(country = nm, Dim1 = X_add[, 1], Dim2 = X_add[, 2])

p_add <- ggplot(df_add, aes(Dim1, Dim2, label = country)) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_point(size = 3.2, colour = "steelblue") +
  geom_text_repel(size = 4.2, seed = 123) +
  coord_equal() +
  labs(title = "ADD-MDS", x = "Dimension 1", y = "Dimension 2") +
  theme_minimal(base_size = 12)

print(p_add)

# --- 図4：TPD-MDS ---
df_tpd <- data.frame(country = nm, Dim1 = X_tpd[, 1], Dim2 = X_tpd[, 2])

p_tpd <- ggplot(df_tpd, aes(Dim1, Dim2, label = country)) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_point(size = 3.2, colour = "steelblue") +
  geom_text_repel(size = 4.2, seed = 123) +
  coord_equal() +
  labs(title = "TPD-MDS", x = "Dimension 1", y = "Dimension 2") +
  theme_minimal(base_size = 12)

print(p_tpd)


# ============================================================================
#  13. 表3のもとになる一覧
# ============================================================================
# 28ペア（8国から2つ選ぶ組み合わせ）を1行ずつ並べる。i < j だけ拾えばよい。

tbl <- data.frame()

for (i in 1:(n - 1)) {
  for (j in (i + 1):n) {
    tbl <- rbind(tbl, data.frame(
      pair  = paste(nm[i], nm[j], sep = "-"),
      ADD   = round(ADD[i, j], 3),
      TPD   = round(TPD[i, j], 3),
      APM   = round(APM[i, j], 3),
      share = round(2 * delta_ADD2[i, j] / d_APM2[i, j] * 100, 1)
    ))
  }
}

tbl <- tbl[order(-tbl$share), ]     # 直接項の寄与が大きい順
print(tbl, row.names = FALSE)

# 要約統計と相関は、表示用に丸めた tbl の列ではなく、丸め前の行列から計算する
u <- upper.tri(APM)                  # 28ペア（上三角）だけを使う
sh <- 100 * 2 * delta_ADD2 / d_APM2  # 直接項の寄与（%）

cat("直接項の寄与 : 最小", round(min(sh[u]), 1), "%  最大", round(max(sh[u]), 1),
    "%  平均", round(mean(sh[u]), 1), "%  中央値", round(median(sh[u]), 1), "%\n")

cat("相関 ADD-APM :", round(cor(ADD[u], APM[u]), 3), "\n")
cat("相関 TPD-APM :", round(cor(TPD[u], APM[u]), 3), "\n")
cat("相関 ADD-TPD :", round(cor(ADD[u], TPD[u]), 3), "\n")


# ============================================================================
#  14. 検証：三角不等式とユークリッド性は別の性質である
# ============================================================================
# SMACOF を使う理由は「三角不等式が破れるから」ではなく、「ユークリッド距離
# 行列である保証がないから」である。この2つは別の性質で、本データの TPD は
# 三角不等式をすべて満たしながら非ユークリッドになる。

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

# 二重中心化行列の最小固有値。負ならユークリッド空間に埋め込めない。
min_eig <- function(M) {
  n <- nrow(M); J <- diag(n) - matrix(1 / n, n, n)
  min(eigen(-0.5 * J %*% (M^2) %*% J, symmetric = TRUE, only.values = TRUE)$values)
}

cat("\n--- 三角不等式とユークリッド性 ---\n")
for (lbl in c("APM", "ADD", "TPD")) {
  M <- get(lbl); v <- tri_violations(M)
  cat(sprintf("%s : 三角不等式違反 %d 組（順序三つ組） / 最小固有値 %.4f\n",
              lbl, v$count, min_eig(M)))
  if (v$count > 0)
    cat(sprintf("     最悪 : %s-%s = %.3f > %s-%s + %s-%s = %.3f\n",
                nm[v$i], nm[v$j], M[v$i, v$j], nm[v$i], nm[v$l], nm[v$l], nm[v$j],
                M[v$i, v$l] + M[v$l, v$j]))
}


# ============================================================================
#  15. 検証：ADD の適合の悪さは次元の不足ではない
# ============================================================================
# 三角不等式が破れているため、どの次元のユークリッド空間にも計量的には
# 埋め込めない。実際 stress-1 は次元を上げても下がらない。

cat("\n--- ADD の次元別 stress-1 ---\n")
for (k in 2:(n - 1)) {
  set.seed(123)
  cat(sprintf("  %d次元 : %.4f\n", k,
              mds(as.dist(ADD), ndim = k, type = "ratio")$stress))
}


# ============================================================================
#  16. 検証：APM の古典的MDSは A の列の中心化主成分分析である
# ============================================================================
# APM の二乗距離行列を二重中心化したグラム行列は J(A'A)J と厳密に一致する。
# よって「負の固有値が出ない」は経験則ではなく定理であり、固有値が対で
# 現れるのも歪対称行列の特異値が対で現れることの帰結。
# なお中心化で差し引かれる列ベクトルの重心は勾配 φ_i = (1/n) Σ_j a_ij に一致する。

J <- diag(n) - matrix(1 / n, n, n)
cat("\n--- APM の構造 ---\n")
cat("  B と J(A'A)J の最大差 :",
    format(max(abs(-0.5 * J %*% (APM^2) %*% J - J %*% (t(A) %*% A) %*% J)), digits = 3), "\n")
cat("  A の特異値（対で現れる）:", paste(round(svd(A)$d, 3), collapse = " "), "\n")
cat("  列ベクトルの重心（= 勾配 φ）\n")
print(round(rowMeans(A), 3))


# ============================================================================
#  17. 検証：TPD 布置（図4）の原点からの距離
# ============================================================================

cat("\n--- TPD 布置の原点（重心）からの距離 ---\n")
print(round(sort(sqrt(rowSums(scale(X_tpd, scale = FALSE)^2))), 3))
