# ============================================================================
#  歪対称行列に基づく非対称MDS：APM / ADD / TPD
#
#  論文「歪対称行列に基づく非対称多次元尺度構成法の提案」の分析コード。
#  本文の数値（説明率、ストレス、割合、三角不等式の違反数、固有値比など）は
#  すべてこのスクリプトの出力から取っている。
#
#  前提：data/trade_2023_comtrade.csv（行=輸出国、列=輸出先、単位=10億米ドル）
#        作業ディレクトリをこのファイルのあるフォルダにして実行する。
#  必要：install.packages(c("smacof", "vegan", "ggplot2", "ggrepel"))
#
#  構成：1. データの読み込みと表示
#        2. 関数定義（分解、3つの行列、MDS、作図、確認用）
#        3. 分析（対数変換 → 3つの行列 → MDS → 図 → 一覧 → 各種の確認）
#        効率よりも読みやすさを優先して書いてある。
#
#  記号の対応（論文 → このコード）
#    T_ij        国 i から国 j への輸出額                → exports[i, j]
#    Δ = (δ_ij)  非対称行列（輸出額の対数）              → M
#    A = (a_ij)  歪対称成分 (Δ - Δ') / 2                 → A（または mats$A）
#    d_APM(i,j)  A の第 i 列と第 j 列のユークリッド距離  → APM[i, j]
#    δ_ADD(i,j)  |a_ij|                                  → ADD[i, j]
#    δ_TPD(i,j)  sqrt( Σ_{k≠i,j} (a_ki - a_kj)^2 )       → TPD[i, j]
#    式(3)/(7)   d_APM^2 = 2 δ_ADD^2 + δ_TPD^2          → stopifnot で検算
# ============================================================================

library(smacof)    # SMACOF（ストレス最小化による MDS）。ADD と TPD の布置に使う
library(vegan)     # procrustes()。ADD・TPD の布置を APM の布置に向きだけ揃える
library(ggplot2)   # 図1〜3
library(ggrepel)   # 図の国名ラベルが重ならないようにずらす


# ============================================================================
#  1. データの読み込みと表示（論文3.1節・表1）
# ============================================================================

# CSV の1列目は国名、2列目以降は各輸出先への輸出額。
# 対角（自国への輸出）は空欄で、NA として読み込まれる。
trade <- read.csv("data/trade_2023_comtrade.csv", fileEncoding = "UTF-8-BOM")

# 国名を行名にして、数値部分だけを行列にする。
rownames(trade) <- trade$country
exports <- data.matrix(trade[, -1])   # exports[i, j] = 国 i から国 j への輸出額（10億米ドル）

nm <- rownames(exports)               # 国名コード（USA, CHN, …）
n  <- nrow(exports)                   # 対象数。ここでは 8

cat("\n=== 二国間輸出額（10億米ドル、行=輸出国、列=輸出先）===\n")
print(exports)
cat("対象数 :", n, "か国\n")


# ============================================================================
#  2. 関数定義
# ============================================================================

# ----------------------------------------------------------------------------
#  2.1 対称成分と歪対称成分への分解（論文の式(1)）
# ----------------------------------------------------------------------------
# 非対称行列 M を対称成分 S と歪対称成分 A に分ける。
#   S = (M + M') / 2 … 往復の平均。本研究では使わない
#   A = (M - M') / 2 … 往復の差の半分。以下の計算はすべてこの A だけを使う
# 対角は無視する。先に 0 を入れておくが、A の対角は M の対角に何が入っていても
# (m_ii - m_ii) / 2 = 0 になるので、結果には影響しない。
# 引数 M：n×n の数値行列。戻り値：list(S, A)。
skew_decompose <- function(M) {
  M <- as.matrix(M)
  diag(M) <- 0
  S <- (M + t(M)) / 2
  A <- (M - t(M)) / 2
  list(S = S, A = A)
}

# ----------------------------------------------------------------------------
#  2.2 APM・ADD・TPD の3つの非類似度行列（論文2節）
# ----------------------------------------------------------------------------
# A から3つの行列を作り、あわせて d_APM^2 に占める 2 δ_ADD^2 の割合（論文7節）を返す。
# 引数 M：非対称行列。戻り値：list(A, APM, ADD, TPD, share)。いずれも n×n。
apm_matrices <- function(M) {
  A <- skew_decompose(M)$A
  n <- nrow(A)

  # APM（式(2)）：A の列ベクトルどうしのユークリッド距離。
  # dist() は行どうしの距離を計算するので、転置して列を行にしてから渡す。
  APM <- as.matrix(dist(t(A)))

  # ADD：対象ペア自身に対応する要素の絶対値 |a_ij|。
  # 符号を落とすので、どちら向きの輸出が大きいかは A の符号を見る必要がある（論文5節）。
  ADD <- abs(A)

  # TPD：式(3)右辺第2項を定義どおりに計算する。恒等式 d_APM^2 = 2 δ_ADD^2 + δ_TPD^2
  # から TPD^2 = APM^2 - 2 ADD^2 と逆算すれば1行で済むが、それだと下の stopifnot が
  # 恒等式を恒等式で確かめることになり、検算にならない。
  #
  # 1ペア (i, j) の TPD：i と j 以外の対象 k について、第 i 列と第 j 列の差を取り、
  # その二乗和の平方根。
  tpd_pair <- function(i, j) {
    others <- setdiff(seq_len(n), c(i, j))   # i と j 以外の対象 k（8か国なら6か国）
    diff_k <- A[others, i] - A[others, j]    # k ごとの偏りの差 a_ki - a_kj
    sqrt(sum(diff_k^2))                      # 差の二乗和の平方根 = δ_TPD(i, j)
  }

  # 全ペアについて tpd_pair を呼んで行列に埋める。TPD は対称なので、
  # i < j のペアを1回ずつ計算して (i, j) と (j, i) の両方に入れる。対角は 0 のまま。
  TPD   <- matrix(0, n, n, dimnames = dimnames(A))
  pairs <- combn(n, 2)                       # 全ペア (i, j), i < j を列に並べた 2 x 28 の行列
  for (p in seq_len(ncol(pairs))) {
    i <- pairs[1, p]
    j <- pairs[2, p]
    TPD[i, j] <- tpd_pair(i, j)
    TPD[j, i] <- TPD[i, j]
  }

  # 式(7)の検算：全ペアで d_APM^2 = 2 δ_ADD^2 + δ_TPD^2 が成り立つこと。
  # 左辺と右辺の差は理論上 0 だが、小数計算の丸め誤差がわずかに残るので、
  # 差が 1e-10（0.0000000001）未満なら「等しい」とみなす。
  # 1つでもそれ以上ずれているペアがあれば、エラーを出してここで止まる。
  identity_gap <- abs(APM^2 - (2 * ADD^2 + TPD^2))
  if (any(identity_gap >= 1e-10)) {
    stop("式(7)の恒等式が成り立っていません。最大のずれ: ", max(identity_gap))
  }

  # d_APM^2 のうち a_ij に由来する部分 2 a_ij^2 が占める割合（0〜1）。対角は 0/0 なので NA。
  share <- 2 * A^2 / APM^2
  diag(share) <- NA

  list(A = A, APM = APM, ADD = ADD, TPD = TPD, share = share)
}

# ----------------------------------------------------------------------------
#  2.3 3つの布置（論文3.3節）
# ----------------------------------------------------------------------------
#   APM … ユークリッド距離行列なので古典的 MDS（cmdscale）
#   ADD・TPD … ユークリッド距離行列である保証がないので SMACOF（smacof::mds, ratio）
# ADD・TPD の布置は vegan::procrustes で APM の布置に向きを揃える（回転・鏡映・平行移動
# のみ。scale = FALSE なので拡大縮小はしない）。3枚の図を見比べやすくするためで、
# 布置の形そのものは変わらない。
# 引数 M：非対称行列、ndim：次元数（論文は2）、seed：SMACOF の初期値用の乱数の種。
# 戻り値：list(matrices, conf, eig, gof, stress, fits)
#   matrices … apm_matrices(M) の結果
#   conf     … 3つの布置を縦につないだ data.frame（method, country, Dim1, Dim2）
#   eig      … APM の古典的 MDS の固有値（負がなければユークリッド距離行列）
#   gof      … APM の2次元説明率（上位 ndim 個の固有値の和 / 固有値の絶対値の和）
#   stress   … ADD と TPD の stress-1
#   fits     … 各手法の当てはめオブジェクトそのもの（追加の確認用）
apm_mds <- function(M, ndim = 2, seed = 123) {
  mats <- apm_matrices(M)
  nm   <- rownames(mats$A)

  # APM：古典的 MDS。eig = TRUE で固有値も受け取る。
  fit_apm <- cmdscale(as.dist(mats$APM), k = ndim, eig = TRUE)

  # ADD・TPD：SMACOF。type = "ratio" は非類似度を定数倍するだけで、順位への置き換え
  # などの変換はしない（論文3.3節・5節）。set.seed は結果を再現可能にするため。
  set.seed(seed)
  fit_add <- mds(as.dist(mats$ADD), ndim = ndim, type = "ratio")
  set.seed(seed)
  fit_tpd <- mds(as.dist(mats$TPD), ndim = ndim, type = "ratio")

  # APM の座標は重心が原点になるよう中心化し、これを基準にして他の2つを回転で合わせる。
  X_apm <- scale(fit_apm$points, scale = FALSE)
  X_add <- procrustes(X_apm, fit_add$conf, scale = FALSE)$Yrot
  X_tpd <- procrustes(X_apm, fit_tpd$conf, scale = FALSE)$Yrot

  # 作図用に3つの布置を1つの data.frame にまとめる。
  as_df <- function(X, method) {
    data.frame(method = method, country = nm, Dim1 = X[, 1], Dim2 = X[, 2],
               row.names = NULL)
  }
  conf <- rbind(as_df(X_apm, "APM"), as_df(X_add, "ADD"), as_df(X_tpd, "TPD"))
  conf$method <- factor(conf$method, levels = c("APM", "ADD", "TPD"))  # 図の並び順

  # 2次元説明率：上位 ndim 個の固有値の和を、固有値の絶対値の和で割る。
  # 負の固有値がなければ分母の取り方で値が変わることはない（論文4節）。
  ev  <- fit_apm$eig
  gof <- sum(ev[1:ndim]) / sum(abs(ev))

  list(matrices = mats,
       conf     = conf,
       eig      = ev,
       gof      = gof,
       stress   = c(ADD = fit_add$stress, TPD = fit_tpd$stress),
       fits     = list(apm = fit_apm, add = fit_add, tpd = fit_tpd))
}

# ----------------------------------------------------------------------------
#  2.4 布置の図（論文の図1〜3）
# ----------------------------------------------------------------------------
# methods に1つ渡せば1枚、省略すれば3枚を横に並べる。
# coord_equal() で縦横の縮尺を同じにしている（点間の距離を目で比べるため）。
# seed は ggrepel のラベル配置用。同じ seed なら同じ配置になる。
plot_apm <- function(fit, methods = c("APM", "ADD", "TPD"), seed = 123) {
  df <- subset(fit$conf, method %in% methods)
  ggplot(df, aes(Dim1, Dim2, label = country)) +
    geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +   # 原点を通る薄い補助線
    geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_point(size = 3.2, colour = "steelblue") +
    geom_text_repel(size = 4.2, seed = seed) +                           # 重ならない国名ラベル
    facet_wrap(~ method, nrow = 1) +                                     # 手法ごとに1枚
    coord_equal() +
    labs(x = "Dimension 1", y = "Dimension 2") +
    theme_minimal(base_size = 12)
}

# ----------------------------------------------------------------------------
#  2.5 全ペアの一覧（論文7節で引用している個々のペアの値はここから読む）
# ----------------------------------------------------------------------------
# 全ペア（8か国なら28ペア）の APM・ADD・TPD の値と、d_APM^2 に占める 2 δ_ADD^2 の
# 割合（%）。割合の大きい順に並べる。
#   pair   … "USA-CHN" のようにペアを並べたもの（左が i、右が j）
#   a      … 符号つきの a_ij。正なら m_ij > m_ji、すなわち i から j への輸出のほうが大きい
#   larger … 大きい向きを "X → Y" で示したもの（m_XY > m_YX）。ペアを書く順に依らず向きが読める
#   share  … 100 * 2 a_ij^2 / d_APM^2
# 引数 mats：apm_matrices() の結果。戻り値：data.frame（1ペア1行）。
pair_table <- function(mats) {
  nm <- rownames(mats$A)

  # 上三角の (行, 列) 番号。各ペアを1回ずつ取り出すために使う。
  u <- which(upper.tri(mats$APM), arr.ind = TRUE)
  from <- nm[u[, 1]]    # ペアの左側の国 i
  to   <- nm[u[, 2]]    # ペアの右側の国 j
  a    <- mats$A[u]     # 符号つきの a_ij

  # 大きい向き。a > 0 なら i → j、a < 0 なら j → i、0 なら均衡。
  larger <- ifelse(a > 0, paste(from, "→", to),
            ifelse(a < 0, paste(to, "→", from), "均衡"))

  tbl <- data.frame(pair   = paste(from, to, sep = "-"),
                    a      = a,
                    larger = larger,
                    APM    = mats$APM[u],
                    ADD    = mats$ADD[u],
                    TPD    = mats$TPD[u],
                    share  = 100 * mats$share[u])
  tbl[order(-tbl$share), ]
}

# ----------------------------------------------------------------------------
#  2.6 三角不等式の違反（論文5節・6節）
# ----------------------------------------------------------------------------
# 三角不等式 M[i,j] <= M[i,l] + M[l,j] の違反を数える。
# 順序三つ組 (i, j, l) をすべて調べるので、8か国なら 8*7*6 = 336 通り。
# 違反のうち、超過分 gap が最も大きい三つ組も返す（本文の「中国とインドの 0.992 は …
# の和 0.241 を上回る」はこれ）。
# 引数 M：対称な非類似度行列。戻り値：list(count = 違反数, i, j, l = 最悪の三つ組の番号)。
tri_violations <- function(M) {
  n     <- nrow(M)
  count <- 0        # 違反の数
  worst_gap <- 0    # これまでで最大の超過分
  worst_i <- NA     # 最大の超過分を出した三つ組
  worst_j <- NA
  worst_l <- NA

  for (i in 1:n) {
    for (j in 1:n) {
      for (l in 1:n) {
        if (i == j || j == l || i == l) {
          next                                   # 3つが異なる組だけ調べる
        }
        gap <- M[i, j] - (M[i, l] + M[l, j])     # 正なら三角不等式が破れている
        if (gap > 1e-12) {                       # 丸め誤差ぶんは違反と数えない
          count <- count + 1
          if (gap > worst_gap) {
            worst_gap <- gap
            worst_i <- i
            worst_j <- j
            worst_l <- l
          }
        }
      }
    }
  }

  list(count = count, i = worst_i, j = worst_j, l = worst_l)
}

# ----------------------------------------------------------------------------
#  2.7 二重中心化行列の固有値（論文2節・5節・6節）
# ----------------------------------------------------------------------------
# 二重中心化行列 B = -1/2 J M^2 J の固有値（J = I - 11'/n は中心化行列）。
# M がユークリッド距離行列であることは、B が半正定値（固有値がすべて非負）であることと
# 同値。最小固有値が負なら、何次元に埋め込んでもその距離は再現できない。
# 本文の「最小固有値の絶対値は最大固有値の 48.0%（ADD）、16.5%（TPD）」はここから。
# 引数 M：対称な非類似度行列。戻り値：固有値（降順）。
dc_eigen <- function(M) {
  n <- nrow(M)
  J <- diag(n) - matrix(1 / n, n, n)       # 中心化行列
  B <- -0.5 * J %*% (M^2) %*% J            # 二重中心化行列
  eigen(B, symmetric = TRUE, only.values = TRUE)$values
}


# ============================================================================
#  3. 分析
# ============================================================================

# ----------------------------------------------------------------------------
#  3.1 対数変換（論文3.2節・表2・式(8)）
# ----------------------------------------------------------------------------
# 輸出額の対数をそのまま非対称行列として渡す（非類似度への反転はしない）。
# このとき歪対称成分は a_ij = (1/2) log(T_ij / T_ji) となり（式(9)）、往復の比だけで
# 決まる。単位を変えても、対数に定数を足しても a_ij は変わらない。
# 符号を反転した -log(T) を渡しても A の符号が変わるだけで、APM・ADD・TPD は不変。
# 対角は NA のままだと log で NA になるので 0 にしておく（A の対角には影響しない）。
M <- log(exports)
diag(M) <- 0

cat("\n=== 輸出額の対数（表2）===\n")
print(round(M, 3))

# ----------------------------------------------------------------------------
#  3.2 3つの行列と MDS
# ----------------------------------------------------------------------------
fit  <- apm_mds(M, ndim = 2)
mats <- fit$matrices

cat("\n=== 歪対称成分 a_ij ===\n")
print(round(mats$A, 3))
cat("\n=== APM（列ベクトル間距離）===\n")
print(round(mats$APM, 3))
cat("\n=== ADD（|a_ij|）===\n")
print(round(mats$ADD, 3))
cat("\n=== TPD（第三者項の平方根）===\n")
print(round(mats$TPD, 3))
cat("\n=== d_APM^2 に占める 2 a_ij^2 の割合（%）===\n")
print(round(100 * mats$share, 1))

# ----------------------------------------------------------------------------
#  3.3 適合度（論文4節〜6節）
# ----------------------------------------------------------------------------
# APM は古典的 MDS なので2次元説明率と最小固有値（0 ならユークリッド距離行列）、
# ADD・TPD は SMACOF なので stress-1（0 に近いほど当てはまりがよい）。
cat("\n=== 適合度 ===\n")
cat("APM 2次元説明率 :", round(100 * fit$gof, 1), "%\n")
cat("APM 最小固有値  :", round(min(fit$eig), 6), "\n")
cat("ADD stress-1    :", round(fit$stress["ADD"], 4), "\n")
cat("TPD stress-1    :", round(fit$stress["TPD"], 4), "\n")

# ----------------------------------------------------------------------------
#  3.4 図1〜3
# ----------------------------------------------------------------------------
# 1枚ずつ出したものを論文に使った。最後の1行は3枚を横に並べた確認用。
print(plot_apm(fit, "APM"))
print(plot_apm(fit, "ADD"))
print(plot_apm(fit, "TPD"))
print(plot_apm(fit))                       # 3枚並べ

# ----------------------------------------------------------------------------
#  3.5 全ペアの一覧と要約（論文7節）
# ----------------------------------------------------------------------------
tbl <- pair_table(mats)
cat("\n=== 全ペアの一覧（割合の大きい順）===\n")
print(format(tbl, digits = 3, nsmall = 1), row.names = FALSE)

# 各ペアを1回ずつ数えるため、上三角だけを取り出す。
upper <- upper.tri(mats$APM)
share_pct <- 100 * mats$share[upper]

cat("\n=== 割合 2δ²ADD/d²APM の要約 ===\n")
cat("最小   :", round(min(share_pct), 1), "%\n")
cat("最大   :", round(max(share_pct), 1), "%\n")
cat("平均   :", round(mean(share_pct), 1), "%\n")
cat("中央値 :", round(median(share_pct), 1), "%\n")

cat("\n=== 3つの行列どうしの相関（28ペア）===\n")
cat("ADD と APM :", round(cor(mats$ADD[upper], mats$APM[upper]), 3), "\n")
cat("TPD と APM :", round(cor(mats$TPD[upper], mats$APM[upper]), 3), "\n")
cat("ADD と TPD :", round(cor(mats$ADD[upper], mats$TPD[upper]), 3), "\n")

# ----------------------------------------------------------------------------
#  3.6 三角不等式とユークリッド性（論文5節・6節）
# ----------------------------------------------------------------------------
# 三角不等式を満たすこと（距離であること）と、ユークリッド距離行列であることは別の
# 性質で、前者は後者を含意しない。APM は両方満たし、TPD は三角不等式は満たすが
# ユークリッドではなく、ADD は三角不等式も破る、という3段階を数値で確かめる。
cat("\n=== 三角不等式とユークリッド性 ===\n")
for (label in c("APM", "ADD", "TPD")) {
  D  <- mats[[label]]
  v  <- tri_violations(D)
  ev <- dc_eigen(D)
  ratio_pct <- 100 * abs(min(ev)) / max(ev)   # 最小固有値の絶対値 / 最大固有値（%）

  cat(sprintf("%s : 三角不等式違反 %d 組 / 最小固有値 %.4f（最大固有値比 %.1f%%）\n",
              label, v$count, min(ev), ratio_pct))
  if (v$count > 0) {
    cat(sprintf("     最悪 : %s-%s = %.3f > %s-%s + %s-%s = %.3f\n",
                nm[v$i], nm[v$j], D[v$i, v$j],
                nm[v$i], nm[v$l], nm[v$l], nm[v$j], D[v$i, v$l] + D[v$l, v$j]))
  }
}

# ----------------------------------------------------------------------------
#  3.7 ADD の適合の悪さは次元不足ではない（論文5節）
# ----------------------------------------------------------------------------
# 次元を 2 から n-1 まで増やしても stress-1 が下がらないことを示す。
# 三角不等式を破る値は何次元の布置でも再現できないため。
cat("\n=== ADD の次元別 stress-1 ===\n")
for (k in 2:(n - 1)) {
  set.seed(123)
  fit_k <- mds(as.dist(mats$ADD), ndim = k, type = "ratio")
  cat(sprintf("  %d次元 : %.4f\n", k, fit_k$stress))
}

# ----------------------------------------------------------------------------
#  3.8 APM の古典的 MDS の構造（確認用。論文には数値を載せていない）
# ----------------------------------------------------------------------------
# APM の二重中心化行列 B は J (A'A) J に一致する。つまり APM の古典的 MDS は、
# A の列ベクトルを中心化して主成分分析するのと同じ。歪対称行列の特異値は対で現れる。
# 行平均は、その国が相手全体に対して輸出超過（正）か輸入超過（負）かの目安で、
# 図2の横軸の並び順と一致する（論文5節末尾）。
A <- mats$A
J <- diag(n) - matrix(1 / n, n, n)
B_from_apm <- -0.5 * J %*% (mats$APM^2) %*% J   # APM から作った二重中心化行列
B_from_A   <- J %*% (t(A) %*% A) %*% J           # A の列を中心化した内積行列

cat("\n=== APM の構造 ===\n")
cat("  B と J(A'A)J の最大差 :", format(max(abs(B_from_apm - B_from_A)), digits = 3), "\n")
cat("  A の特異値（対で現れる）:", paste(round(svd(A)$d, 3), collapse = " "), "\n")
cat("  A の行平均（正なら相手全体に対して輸出超過）\n")
print(round(rowMeans(A), 3))

# ----------------------------------------------------------------------------
#  3.9 TPD 布置の重心からの距離（確認用）
# ----------------------------------------------------------------------------
# 中心に近い国ほど、他国に対する偏りの向きと大きさが「平均的」であることを示す。
X_tpd <- as.matrix(subset(fit$conf, method == "TPD")[, c("Dim1", "Dim2")])
X_tpd_centered <- scale(X_tpd, scale = FALSE)                 # 重心を原点に
dist_from_center <- sqrt(rowSums(X_tpd_centered^2))           # 各国の原点からの距離
names(dist_from_center) <- nm

cat("\n=== TPD 布置の重心からの距離 ===\n")
print(round(sort(dist_from_center), 3))
