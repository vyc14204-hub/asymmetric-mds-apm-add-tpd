# ============================================================================
#  歪対称行列に基づく非対称MDS：APM / ADD / TPD の3手法
#
#  前提：R Canvas のシート trade_2023_comtrade
#        （UN Comtrade 2023年、行=輸出国、列=輸出先、単位=米ドル）
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

trade <- trade_2023_comtrade
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
#     d_ij = -log(T_ij)
#
# こうすると歪対称成分の要素が
#     a_ij = (1/2) * log(T_ji / T_ij)
# になる。二方向の貿易額の対数比の半分で、解釈しやすい。
#
# 対角は 0 にしておく。ただし a_ii = (d_ii - d_ii)/2 = 0 なので、
# 実は対角に何を入れても A の対角は 0 になる。

D <- -log(T)
diag(D) <- 0

print(round(D, 3))


# ============================================================================
#  3. 対称成分 S と歪対称成分 A に分ける
# ============================================================================
#     S = (D + D') / 2   … 二方向の平均。今回は使わないが確認用に作る
#     A = (D - D') / 2   … 二方向の差。これが本研究の対象

S <- (D + t(D)) / 2
A <- (D - t(D)) / 2

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

APM2 <- matrix(0, n, n, dimnames = list(nm, nm))

for (i in 1:n) {
  for (j in 1:n) {
    if (i != j) {
      APM2[i, j] <- sum((A[, i] - A[, j])^2)   # A[, i] は第 i 列
    }
  }
}

APM <- sqrt(APM2)

print(round(APM, 3))


# ============================================================================
#  5. ADD を計算する
# ============================================================================
#     d_ADD^2(i,j) = a_ij^2
#
# 対象ペア自身の非対称量。A の各要素を二乗するだけ。第三者は関与しない。

ADD2 <- A^2
ADD  <- sqrt(ADD2)      # abs(A) と同じ

print(round(ADD, 3))


# ============================================================================
#  6. TPD を計算する
# ============================================================================
#     d_TPD^2(i,j) = sum_{k != i,j} (a_ki - a_kj)^2
#
# APM の和から k = i と k = j の項を除いたもの。
# setdiff(1:n, c(i,j)) が「i でも j でもない k」の並びを作る。

TPD2 <- matrix(0, n, n, dimnames = list(nm, nm))

for (i in 1:n) {
  for (j in 1:n) {
    if (i != j) {
      k <- setdiff(1:n, c(i, j))          # 第三者のインデックス（n-2 個）
      TPD2[i, j] <- sum((A[k, i] - A[k, j])^2)
    }
  }
}

TPD <- sqrt(TPD2)

print(round(TPD, 3))


# ============================================================================
#  7. 分解の検算
# ============================================================================
# 命題：d_APM^2 = 2 * d_ADD^2 + d_TPD^2
# ここがずれたらどこかで間違えている。

cat("分解定理が成立するか :", all(abs(APM2 - (2 * ADD2 + TPD2)) < 1e-10), "\n")

# 各ペアで直接項が占める割合（本文の表1になる）
share <- 2 * ADD2 / APM2
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
# ADD は三角不等式の保証が無い。cmdscale だと負の固有値が出て解が歪むので、
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
# ADD と同じ手順。TPD も三角不等式の保証は無い。

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

# --- 図1：APM-MDS ---
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

# --- 図2：ADD-MDS ---
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

# --- 図3：TPD-MDS ---
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
#  13. 表1のもとになる一覧
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
      share = round(2 * ADD2[i, j] / APM2[i, j] * 100, 1)
    ))
  }
}

tbl <- tbl[order(-tbl$share), ]     # 直接項の寄与が大きい順
print(tbl, row.names = FALSE)

cat("直接項の寄与 : 最小", min(tbl$share), "%  最大", max(tbl$share),
    "%  平均", round(mean(tbl$share), 1), "%  中央値", round(median(tbl$share), 1), "%\n")

cat("相関 ADD-APM :", round(cor(tbl$ADD, tbl$APM), 3), "\n")
cat("相関 TPD-APM :", round(cor(tbl$TPD, tbl$APM), 3), "\n")
cat("相関 ADD-TPD :", round(cor(tbl$ADD, tbl$TPD), 3), "\n")


# ============================================================================
#  14. 布置の解釈に使う量
# ============================================================================
# 第1次元と第2次元が何に対応しているかを確かめるための補助的な計算。
#
#   列和     sum_k a_ki   … 対数ベースの純輸出ポジション。輸出超過なら正
#   絶対値和 sum_k |a_ki| … その国が関与する非対称性の総量

net   <- colSums(A)
total <- colSums(abs(A))

cat("\n国別 純輸出ポジション（列和）\n")
print(round(sort(net, decreasing = TRUE), 3))

cat("\n国別 非対称性の総量（絶対値の列和）\n")
print(round(sort(total, decreasing = TRUE), 3))

cat("\n第1次元との相関\n")
cat("  APM :", round(cor(net, X_apm[, 1]), 3), "\n")
cat("  ADD :", round(cor(net, X_add[, 1]), 3), "\n")
cat("  TPD :", round(cor(net, X_tpd[, 1]), 3), "\n")

cat("\n第2次元と非対称性の総量との相関\n")
cat("  APM :", round(cor(total, X_apm[, 2]), 3), "\n")
cat("  TPD :", round(cor(total, X_tpd[, 2]), 3), "\n")
