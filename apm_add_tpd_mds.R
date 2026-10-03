# ============================================================================
#  歪対称行列に基づく非対称MDS：APM / ADD / TPD（本体）
#
#  論文「歪対称行列に基づく非対称多次元尺度構成法の提案」の分析コード。
#  このファイルは論文の流れどおりに進む：
#    1. データを読む（表1）
#    2. 関数を定義する
#    3. 対数をとる（表2）→ 3つの行列を作る → MDS をかける → 図1〜3 → 全ペアの一覧
#  本文の裏付けに使った確認計算（三角不等式、固有値、次元別ストレスなど）は
#  別ファイル apm_add_tpd_checks.R にある。
#
#  前提：data/trade_2023_comtrade.csv（行=輸出国、列=輸出先、単位=10億米ドル）
#        作業ディレクトリをこのファイルのあるフォルダにして実行する。
#  必要：install.packages(c("smacof", "vegan", "ggplot2", "ggrepel"))
#
#  記号の対応（論文 → このコード）
#    T_ij        国 i から国 j への輸出額                → exports[i, j]
#    Δ = (δ_ij)  非対称行列（輸出額の対数）              → Delta
#    A = (a_ij)  歪対称成分 (Δ - Δ') / 2                 → A
#    d_APM(i,j)  A の第 i 列と第 j 列のユークリッド距離  → APM[i, j]
#    δ_ADD(i,j)  |a_ij|                                  → ADD[i, j]
#    δ_TPD(i,j)  sqrt( Σ_{k≠i,j} (a_ki - a_kj)^2 )       → TPD[i, j]
#    式(3)/(7)   d_APM^2 = 2 δ_ADD^2 + δ_TPD^2          → check_identity で検算
#
#  ---------------------------------------------------------------------------
#  R の文法メモ（本文中の初出の箇所にも説明を付けてある）
#  ---------------------------------------------------------------------------
#    x <- 値            代入。「x に値を入れる」
#    #                  行末までコメント
#    f(a, b = 2)        関数呼び出し。b = 2 のように引数名を付けて渡せる
#    function(x) {...}  関数を作る。return(値) が戻り値
#    c(1, 2, 3)         ベクトル（1次元の並び）を作る
#    seq_len(n)         1 から n までの整数の並び
#    M[i, j]            行列の i 行 j 列。M[i, ] は第 i 行全部、M[, j] は第 j 列全部
#    M[, -1]            マイナスは「それ以外」。1列目を除いた全列
#    list(a = 1, b = 2) 名前付きリスト。種類の違うものをまとめて入れる箱
#    lst$a              リストの要素 a を取り出す
#    t(M)               転置。行と列を入れ替える
#    M^2, abs(M), M / N 行列の要素ごとの二乗、絶対値、割り算
#    diag(M) <- 0       対角成分に 0 を入れる
#    for (i in 1:n) {}  繰り返し。i に 1, 2, …, n を順に入れて { } を実行
#    if (条件) {}       条件分岐
#    any(x)             TRUE/FALSE の並び x に、1つでも TRUE があるか
#    cat(...)           文字や数値をそのまま画面に出す。"\n" は改行
#    print(x)           x を R の標準の形で画面に出す
#    round(x, 3)        小数第3位に丸める
#    set.seed(123)      乱数の種を固定する。乱数を使う計算を毎回同じ結果にするため
# ============================================================================

# library(パッケージ名)：追加のパッケージを読み込んで、その関数を使えるようにする。
library(smacof)    # SMACOF（ストレス最小化による MDS）。ADD と TPD の布置に使う
library(vegan)     # procrustes()。ADD・TPD の布置を APM の布置に向きだけ揃える
library(ggplot2)   # 図1〜3
library(ggrepel)   # 図の国名ラベルが重ならないようにずらす


# ============================================================================
#  1. データの読み込みと表示（論文3.1節・表1）
# ============================================================================

# read.csv(ファイル名)：CSV を読み込んで data.frame（表）にする。
# fileEncoding = "UTF-8-BOM" は、先頭に BOM の付いた UTF-8 ファイルを正しく読む指定。
# CSV の1列目は国名、2列目以降は各輸出先への輸出額。
# 対角（自国への輸出）は空欄で、NA（欠損）として読み込まれる。
trade <- read.csv("data/trade_2023_comtrade.csv", fileEncoding = "UTF-8-BOM")

# 国名を行名にして、数値部分だけを行列にする。
# trade$country は列 country。trade[, -1] は1列目（国名）を除いた全列。
rownames(trade) <- trade$country
exports <- data.matrix(trade[, -1])   # exports[i, j] = 国 i から国 j への輸出額（10億米ドル）

nm <- rownames(exports)               # 国名コード（USA, CHN, …）。図や表のラベルに使う
n  <- nrow(exports)                   # 対象数。ここでは 8

cat("\n=== 二国間輸出額（10億米ドル、行=輸出国、列=輸出先）===\n")
print(exports)
cat("対象数 :", n, "か国\n")


# ============================================================================
#  2. 関数定義
#     function(引数) { 処理 } で関数を作り、名前 <- で名前を付ける。
#     定義しただけでは何も起きない。3. の分析部で呼び出したときに動く。
# ============================================================================

# ----------------------------------------------------------------------------
#  2.1 歪対称成分を取り出す（論文の式(1)）
# ----------------------------------------------------------------------------
# 非対称行列 Δ は、対称成分 S = (Δ + Δ') / 2 と歪対称成分 A = (Δ - Δ') / 2 の和に分けられる。
# 本研究で使うのは A だけなので、A だけ返す。
# 対角は先に 0 にしておく。A の対角は Δ の対角が何であっても (δ_ii - δ_ii)/2 = 0 になる。
# 引数 Delta：n×n の数値行列。戻り値：歪対称成分 A（n×n）。
skew_part <- function(Delta) {
  Delta <- as.matrix(Delta)        # data.frame で渡されても行列に直す
  diag(Delta) <- 0                 # 対角成分を 0 にする
  A <- (Delta - t(Delta)) / 2      # t() は転置。引き算と / 2 は要素ごとに計算される
  return(A)
}

# ----------------------------------------------------------------------------
#  2.2 APM・ADD・TPD の3つの非類似度行列（論文2節）
#      1行列1関数。最後の apm_matrices がそれらを順に呼んでまとめる。
# ----------------------------------------------------------------------------

# APM（式(2)）：A の列ベクトルどうしのユークリッド距離。
# dist(X) は X の「行」どうしの距離をまとめて計算するので、t(A) で転置して
# 列を行にしてから渡す。dist の結果は下三角だけの特殊な形なので、
# as.matrix() で普通の n×n 行列に直す。
apm_from <- function(A) {
  APM <- as.matrix(dist(t(A)))
  return(APM)
}

# ADD：対象ペア自身に対応する要素の絶対値 |a_ij|。abs() は要素ごとの絶対値。
# 符号を落とすので、どちら向きの輸出が大きいかは A の符号を見る（論文5節）。
add_from <- function(A) {
  ADD <- abs(A)
  return(ADD)
}

# TPD：式(3)右辺第2項の平方根。「i と j 以外の対象 k について、第 i 列と第 j 列の差を
# とり、その二乗和の平方根」を、定義どおりに計算する。
# 恒等式から TPD^2 = APM^2 - 2 ADD^2 と逆算すれば1行で済むが、それだと check_identity が
# 恒等式を恒等式で確かめることになり、検算にならない。
tpd_from <- function(A) {
  n <- nrow(A)

  # 1ペア (i, j) ぶんの計算。関数の中で定義した関数は、外側の A と n をそのまま使える。
  tpd_pair <- function(i, j) {
    others <- setdiff(seq_len(n), c(i, j))   # setdiff(全体, 除くもの)：i と j 以外の対象 k
    diff_k <- A[others, i] - A[others, j]    # k ごとの偏りの差 a_ki - a_kj（ベクトル）
    value  <- sqrt(sum(diff_k^2))            # 差の二乗和の平方根 = δ_TPD(i, j)
    return(value)
  }

  # 全ペアについて tpd_pair を呼んで行列に埋める。
  # matrix(0, n, n)：全要素 0 の n×n 行列。dimnames = で行名・列名を A と同じにする。
  # combn(n, 2)：1〜n から2つ選ぶ組み合わせを全部列挙し、1組を1列にした 2×28 の行列。
  # TPD は対称なので、i < j のペアを1回ずつ計算して (i, j) と (j, i) の両方に入れる。
  TPD   <- matrix(0, n, n, dimnames = dimnames(A))
  pairs <- combn(n, 2)
  for (p in seq_len(ncol(pairs))) {          # p = 1, 2, …, 28（ペアの番号）
    i <- pairs[1, p]                         # p 番目のペアの1つ目の国
    j <- pairs[2, p]                         # p 番目のペアの2つ目の国
    TPD[i, j] <- tpd_pair(i, j)
    TPD[j, i] <- TPD[i, j]
  }

  return(TPD)
}

# 式(7)の検算：全ペアで d_APM^2 = 2 δ_ADD^2 + δ_TPD^2 が成り立つこと。
# これは A が歪対称でありさえすれば必ず成り立つ恒等式なので、ここで止まるとしたら
# 上の関数のどれかにバグがある。
# 左辺と右辺の差は理論上 0 だが、小数計算の丸め誤差がわずかに残るので、
# 差が 1e-10（0.0000000001）未満なら「等しい」とみなす。
# any(条件)：1つでも条件を満たす要素があれば TRUE。stop(メッセージ)：エラーで止める。
check_identity <- function(APM, ADD, TPD) {
  gap <- abs(APM^2 - (2 * ADD^2 + TPD^2))
  if (any(gap >= 1e-10)) {
    stop("式(7)の恒等式が成り立っていません。最大のずれ: ", max(gap))
  }
}

# d_APM^2 のうち a_ij に由来する部分 2 a_ij^2 が占める割合（0〜1。論文7節）。
# 対角は 0/0 で意味がないので NA にする。
share_from <- function(A, APM) {
  share <- 2 * A^2 / APM^2
  diag(share) <- NA
  return(share)
}

# 上の関数を順に呼んで、3つの行列と割合をひとまとめにする。
# 引数 Delta：非対称行列。戻り値：list(A, APM, ADD, TPD, share)。
# 呼び出し側では mats <- apm_matrices(Delta) と受け取り、mats$APM, mats$TPD のように使う。
apm_matrices <- function(Delta) {
  A   <- skew_part(Delta)
  APM <- apm_from(A)
  ADD <- add_from(A)
  TPD <- tpd_from(A)
  check_identity(APM, ADD, TPD)
  share <- share_from(A, APM)
  return(list(A = A, APM = APM, ADD = ADD, TPD = TPD, share = share))
}

# ----------------------------------------------------------------------------
#  2.3 3つの布置（論文3.3節）
# ----------------------------------------------------------------------------
#   APM … ユークリッド距離行列なので古典的 MDS（cmdscale）
#   ADD・TPD … ユークリッド距離行列である保証がないので SMACOF（smacof の mds, ratio）
# ADD・TPD の布置は procrustes で APM の布置に向きを揃える（回転・鏡映・平行移動のみ。
# scale = FALSE なので拡大縮小はしない）。3枚の図を見比べやすくするためで、
# 布置の形そのものは変わらない。
# 引数 mats：apm_matrices() の結果、ndim：次元数（論文は2）、seed：乱数の種。
#   ndim = 2 のように書いた引数は「省略したときの既定値」。
# 戻り値：list(conf, eig, gof, stress)
#   conf   … 3つの布置を縦につないだ表（method, country, Dim1, Dim2）
#   eig    … APM の古典的 MDS の固有値（負がなければユークリッド距離行列）
#   gof    … APM の2次元説明率
#   stress … ADD と TPD の stress-1
apm_mds <- function(mats, ndim = 2, seed = 123) {
  nm <- rownames(mats$A)

  # APM：古典的 MDS。as.dist() は行列を「距離行列」の型に変える（cmdscale はこの型を受け取る）。
  # eig = TRUE で固有値も一緒に返してもらう。
  fit_apm <- cmdscale(as.dist(mats$APM), k = ndim, eig = TRUE)

  # ADD・TPD：SMACOF。type = "ratio" は非類似度を定数倍するだけで、順位への置き換えなどの
  # 変換はしない（論文3.3節）。乱数で初期値を決めるので、set.seed で結果を再現可能にする。
  set.seed(seed)
  fit_add <- mds(as.dist(mats$ADD), ndim = ndim, type = "ratio")
  set.seed(seed)
  fit_tpd <- mds(as.dist(mats$TPD), ndim = ndim, type = "ratio")

  # APM の座標は scale(scale = FALSE) で重心を原点に移す。
  # procrustes(基準, 合わせたい座標)$Yrot が、基準に向きを合わせたあとの座標。
  X_apm <- scale(fit_apm$points, scale = FALSE)
  X_add <- procrustes(X_apm, fit_add$conf, scale = FALSE)$Yrot
  X_tpd <- procrustes(X_apm, fit_tpd$conf, scale = FALSE)$Yrot

  # 作図用に3つの布置を1つの表にまとめる。
  # data.frame(列名 = 値, …) で表を作る。X[, 1] は座標の1列目（第1次元）。rbind() は表を縦につなぐ。
  as_table <- function(X, method) {
    data.frame(method = method, country = nm, Dim1 = X[, 1], Dim2 = X[, 2], row.names = NULL)
  }
  conf <- rbind(as_table(X_apm, "APM"), as_table(X_add, "ADD"), as_table(X_tpd, "TPD"))

  # factor(…, levels = …)：文字列を「順序の決まったカテゴリ」にする。図を APM, ADD, TPD の順に並べるため。
  conf$method <- factor(conf$method, levels = c("APM", "ADD", "TPD"))

  # 2次元説明率：上位 ndim 個の固有値の和を、固有値の絶対値の和で割る（論文4節）。
  ev  <- fit_apm$eig
  gof <- sum(ev[1:ndim]) / sum(abs(ev))

  # c(ADD = …, TPD = …) は名前付きのベクトル。fit$stress["ADD"] で取り出せる。
  stress <- c(ADD = fit_add$stress, TPD = fit_tpd$stress)

  return(list(conf = conf, eig = ev, gof = gof, stress = stress))
}

# ----------------------------------------------------------------------------
#  2.4 布置の図（論文の図1〜3）
# ----------------------------------------------------------------------------
# methods に1つ渡せば1枚、省略すれば3枚を横に並べる。
# ggplot2 は「ggplot(データ, aes(x, y)) + 層 + 層 + …」と + でつないで図を組み立てる。
# coord_equal() で縦横の縮尺を同じにしている（点間の距離を目で比べるため）。
plot_apm <- function(fit, methods = c("APM", "ADD", "TPD"), seed = 123) {
  # subset(表, 条件)：条件に合う行だけ取り出す。%in% は「左が右の中に含まれるか」。
  df <- subset(fit$conf, method %in% methods)
  ggplot(df, aes(Dim1, Dim2, label = country)) +                        # x, y, ラベルの列を指定
    geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +   # 原点を通る薄い補助線
    geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_point(size = 3.2, colour = "steelblue") +                       # 点
    geom_text_repel(size = 4.2, seed = seed) +                           # 重ならない国名ラベル
    facet_wrap(~ method, nrow = 1) +                                     # method ごとに1枚ずつ横に
    coord_equal() +
    labs(x = "Dimension 1", y = "Dimension 2") +                         # 軸ラベル
    theme_minimal(base_size = 12)                                        # 見た目のテーマ
}

# ----------------------------------------------------------------------------
#  2.5 全ペアの一覧（論文7節で引用している個々のペアの値はここから読む）
# ----------------------------------------------------------------------------
# 全ペア（8か国なら28ペア）の APM・ADD・TPD の値と、d_APM^2 に占める 2 δ_ADD^2 の
# 割合（%）。割合の大きい順に並べる。
#   pair   … "USA-CHN" のようにペアを並べたもの（左が i、右が j）
#   a      … 符号つきの a_ij。正なら i から j への輸出のほうが大きい
#   larger … 大きい向きを "X → Y" で示したもの。ペアを書く順に依らず向きが読める
#   share  … 100 * 2 a_ij^2 / d_APM^2
pair_table <- function(mats) {
  nm <- rownames(mats$A)

  # upper.tri(M)：上三角（行番号 < 列番号）の位置が TRUE の行列。
  # which(…, arr.ind = TRUE)：TRUE の位置を (行, 列) の番号の表にして返す。
  # これで各ペアを1回ずつ取り出せる（28行 × 2列）。
  u    <- which(upper.tri(mats$APM), arr.ind = TRUE)
  from <- nm[u[, 1]]    # ペアの左側の国 i（u の1列目は行番号）
  to   <- nm[u[, 2]]    # ペアの右側の国 j（u の2列目は列番号）
  a    <- mats$A[u]     # 行列を (行, 列) の表で添字付けすると、その位置の値の並びが取れる

  # 大きい向き。ifelse(条件, 真のときの値, 偽のときの値) を要素ごとに評価する。
  larger <- ifelse(a > 0, paste(from, "→", to),
            ifelse(a < 0, paste(to, "→", from), "均衡"))

  tbl <- data.frame(pair   = paste(from, to, sep = "-"),   # "USA-CHN" のような文字列
                    a      = a,
                    larger = larger,
                    APM    = mats$APM[u],
                    ADD    = mats$ADD[u],
                    TPD    = mats$TPD[u],
                    share  = 100 * mats$share[u])

  # order(-x)：x の大きい順に並べたときの行番号。tbl[その順, ] で行を並べ替える。
  tbl <- tbl[order(-tbl$share), ]
  return(tbl)
}


# ============================================================================
#  3. 分析
#     ここから先が実際に動く部分。上で定義した関数を順に呼ぶ。
# ============================================================================

# ----------------------------------------------------------------------------
#  3.1 対数変換（論文3.2節・表2・式(8)）
# ----------------------------------------------------------------------------
# 輸出額の対数をそのまま非対称行列 Δ として使う。
# このとき歪対称成分は a_ij = (1/2) log(T_ij / T_ji) となり（式(9)）、往復の比だけで決まる。
# log() は要素ごとの自然対数。対角は NA のままだと log で NA になるので 0 にしておく
# （A の対角には影響しない）。
Delta <- log(exports)
diag(Delta) <- 0

cat("\n=== 輸出額の対数 Δ（表2）===\n")
print(round(Delta, 3))

# ----------------------------------------------------------------------------
#  3.2 3つの行列
# ----------------------------------------------------------------------------
mats <- apm_matrices(Delta)

cat("\n=== 歪対称成分 A ===\n")
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
#  3.3 MDS と適合度（論文4節〜6節）
# ----------------------------------------------------------------------------
# APM は古典的 MDS なので2次元説明率と最小固有値（0 ならユークリッド距離行列）、
# ADD・TPD は SMACOF なので stress-1（0 に近いほど当てはまりがよい）。
fit <- apm_mds(mats, ndim = 2)

cat("\n=== 適合度 ===\n")
cat("APM 2次元説明率 :", round(100 * fit$gof, 1), "%\n")
cat("APM 最小固有値  :", round(min(fit$eig), 6), "\n")
cat("ADD stress-1    :", round(fit$stress["ADD"], 4), "\n")
cat("TPD stress-1    :", round(fit$stress["TPD"], 4), "\n")

# ----------------------------------------------------------------------------
#  3.4 図1〜3
# ----------------------------------------------------------------------------
# ggplot の図はスクリプト実行時には print() しないと描かれない。
# 1枚ずつ出したものを論文に使った。最後の1行は3枚を横に並べた確認用。
print(plot_apm(fit, "APM"))
print(plot_apm(fit, "ADD"))
print(plot_apm(fit, "TPD"))
print(plot_apm(fit))

# ----------------------------------------------------------------------------
#  3.5 全ペアの一覧と要約（論文7節）
# ----------------------------------------------------------------------------
# format(表, digits = 3, nsmall = 1)：数値を見やすい桁に整える。row.names = FALSE：行番号を出さない。
tbl <- pair_table(mats)
cat("\n=== 全ペアの一覧（割合の大きい順）===\n")
print(format(tbl, digits = 3, nsmall = 1), row.names = FALSE)

cat("\n=== 割合 2δ²ADD/d²APM の要約（28ペア）===\n")
cat("最小   :", round(min(tbl$share), 1), "%\n")
cat("最大   :", round(max(tbl$share), 1), "%\n")
cat("平均   :", round(mean(tbl$share), 1), "%\n")
cat("中央値 :", round(median(tbl$share), 1), "%\n")

# cor(x, y)：ピアソンの積率相関係数。
cat("\n=== 3つの行列どうしの相関（28ペア）===\n")
cat("ADD と APM :", round(cor(tbl$ADD, tbl$APM), 3), "\n")
cat("TPD と APM :", round(cor(tbl$TPD, tbl$APM), 3), "\n")
cat("ADD と TPD :", round(cor(tbl$ADD, tbl$TPD), 3), "\n")
