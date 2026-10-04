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
#  書き方の方針：効率より読みやすさ。1行に1つの処理。入れ子の式は中間変数に分ける。
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
#    x[条件]            条件が TRUE の要素だけ取り出す。x[条件] <- 値 でそこだけ書き換える
#    list(a = 1, b = 2) 名前付きリスト。種類の違うものをまとめて入れる箱
#    lst$a              リストの要素 a を取り出す
#    t(M)               転置。行と列を入れ替える
#    M^2, abs(M), M / N 行列の要素ごとの二乗、絶対値、割り算
#    diag(M) <- 0       対角成分に 0 を入れる
#    for (i in 1:n) {}  繰り返し。i に 1, 2, …, n を順に入れて { } を実行
#    if (条件) {}       条件分岐
#    cat(...)           文字や数値をそのまま画面に出す。"\n" は改行
#    print(x)           x を R の標準の形で画面に出す
#    round(x, 3)        小数第3位に丸める
#    set.seed(123)      乱数の種を固定する。乱数を使う計算を毎回同じ結果にするため
#    pkg::f()           パッケージ pkg の関数 f を、library なしで名前を明示して呼ぶ
# ============================================================================

# library(パッケージ名)：追加のパッケージを読み込んで、その関数を使えるようにする。
library(smacof)    # SMACOF（ストレス最小化による MDS）。ADD と TPD の布置に使う
library(vegan)     # procrustes()。ADD・TPD の布置を APM の布置に向きだけ揃える
library(ggplot2)   # 図1〜3
library(ggrepel)   # 図の国名ラベルが重ならないようにずらす


# ============================================================================
#  1. データの読み込みと表示（論文3.1節・表1）
# ============================================================================

# CSV の場所。"data/…" は相対パスで、「R の作業ディレクトリから見て data フォルダの中」の意味。
# R はこのスクリプトがどこに置いてあるかを知らないので、作業ディレクトリ次第で見つからない。
# よくあるのは、作業ディレクトリが data フォルダそのものになっている場合（CSV を data の中から
# 読み込んだあとなど）。そのときは data/data/… を探してしまうので、いまのフォルダ直下も見る。
# どちらにもなければ、いま R がどこを見ているか（getwd()）と対処を表示して止まる。
#   file.exists(パス)：そのファイルがあるか。getwd()：いまの作業ディレクトリ。
#   stop(…)：エラーを出して止める。paste0(…)：文字列を隙間なくつなぐ。
csv_file <- "data/trade_2023_comtrade.csv"          # 通常：data フォルダの中
if (!file.exists(csv_file)) {
  csv_file <- "trade_2023_comtrade.csv"             # 作業ディレクトリが data のとき
}
if (!file.exists(csv_file)) {
  message_text <- paste0(
    "データファイル trade_2023_comtrade.csv が見つかりません。\n",
    "  いまの R の作業ディレクトリ: ", getwd(), "\n",
    "  このスクリプトと data フォルダがあるフォルダを作業ディレクトリにしてから実行してください。\n",
    "  例: setwd(\"C:/.../APM-ADD-TPD論文\")")
  stop(message_text)
}
cat("読み込むファイル :", normalizePath(csv_file), "\n")   # 実際に読んだ場所を表示しておく

# read.csv(ファイル名)：CSV を読み込んで data.frame（表）にする。
# fileEncoding = "UTF-8-BOM" は、先頭に BOM の付いた UTF-8 ファイルを正しく読む指定。
# CSV の1列目は国名、2列目以降は各輸出先への輸出額。
# 対角（自国への輸出）は空欄で、NA（欠損）として読み込まれる。
trade <- read.csv(csv_file, fileEncoding = "UTF-8-BOM")

# 国名を行名にして、数値部分だけを行列にする。
rownames(trade) <- trade$country        # trade$country は列 country（国名）
values_only <- trade[, -1]              # 1列目（国名）を除いた全列
exports <- data.matrix(values_only)     # 表を数値の行列に。exports[i, j] = 国 i から国 j への輸出額

nm <- rownames(exports)                 # 国名コード（USA, CHN, …）。図や表のラベルに使う
n  <- nrow(exports)                     # 対象数。ここでは 8

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
# 計算には R に同梱の Matrix パッケージの skewpart() を使う（(Δ - Δ') / 2 と同じもの）。
# Matrix::skewpart と書くと、library(Matrix) をしなくてもその関数を呼べる。
# 結果は Matrix パッケージ独自の型で返るので、as.matrix() で普通の行列に戻す。
# 対角は先に 0 にしておく。A の対角は Δ の対角が何であっても (δ_ii - δ_ii)/2 = 0 になる。
# 引数 Delta：n×n の数値行列。戻り値：歪対称成分 A（n×n）。
skew_part <- function(Delta) {
  Delta <- as.matrix(Delta)                   # data.frame で渡されても行列に直す
  diag(Delta) <- 0                            # 対角成分を 0 にする
  skew <- Matrix::skewpart(Delta)             # 歪対称成分 (Δ - Δ') / 2（Matrix 型）
  A <- as.matrix(skew)                        # 普通の行列に戻す
  dimnames(A) <- dimnames(Delta)              # 行名・列名（国名）を付け直す
  return(A)
}

# ----------------------------------------------------------------------------
#  2.2 APM・ADD・TPD の3つの非類似度行列（論文2節）
#      1行列1関数。分析部（3.2）でこれらを順に呼ぶ。
# ----------------------------------------------------------------------------

# APM（式(2)）：A の列ベクトルどうしのユークリッド距離。
# dist(X) は X の「行」どうしの距離をまとめて計算するので、先に転置して列を行にする。
# dist の結果は下三角だけの特殊な形なので、as.matrix() で普通の n×n 行列に直す。
apm_from <- function(A) {
  columns_as_rows <- t(A)                     # 転置。A の列（国）が行になる
  distances <- dist(columns_as_rows)          # 行どうしのユークリッド距離（dist 型）
  APM <- as.matrix(distances)                 # 普通の n×n 行列に
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
    others  <- setdiff(seq_len(n), c(i, j))  # setdiff(全体, 除くもの)：i と j 以外の対象 k
    diff_k  <- A[others, i] - A[others, j]   # k ごとの偏りの差 a_ki - a_kj（ベクトル）
    squared <- diff_k^2                      # 差の二乗
    total   <- sum(squared)                  # 二乗和
    value   <- sqrt(total)                   # 平方根 = δ_TPD(i, j)
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
# 左辺と右辺は理論上一致するが、小数計算の丸め誤差がわずかに残るので、
# R 標準の all.equal() で比べる。all.equal(x, y) は、x と y が丸め誤差の範囲で等しければ
# TRUE を、違えば「どう違うか」の文字列を返す。isTRUE() でその結果が TRUE かどうかを見る。
# stop(メッセージ)：エラーを出して実行を止める。
check_identity <- function(APM, ADD, TPD) {
  lhs <- APM^2                     # 左辺 d_APM^2
  rhs <- 2 * ADD^2 + TPD^2         # 右辺 2 δ_ADD^2 + δ_TPD^2
  comparison <- all.equal(lhs, rhs)
  if (!isTRUE(comparison)) {
    stop("式(7)の恒等式が成り立っていません: ", comparison)
  }
}

# d_APM^2 のうち a_ij に由来する部分 2 a_ij^2 が占める割合（0〜1。論文7節）。
# 対角は 0/0 で意味がないので NA にする。
share_from <- function(A, APM) {
  numerator   <- 2 * A^2           # 分子：2 a_ij^2
  denominator <- APM^2             # 分母：d_APM^2
  share <- numerator / denominator # 要素ごとの割り算
  diag(share) <- NA
  return(share)
}

# ----------------------------------------------------------------------------
#  2.3 MDS（論文3.3節）。手法ごとに1関数。
#      APM は古典的 MDS、ADD と TPD は SMACOF。
# ----------------------------------------------------------------------------

# APM の MDS：古典的多次元尺度構成法（R 標準の cmdscale）。
# APM はユークリッド距離行列なので、固有値分解で座標が直接求まる。
# as.dist() は行列を「距離行列」の型に変える（cmdscale はこの型を受け取る）。
# eig = TRUE で固有値も一緒に返してもらう。
# 引数 APM：APM 行列、ndim：次元数（論文は2）。
# 戻り値：list(points = 座標（n×ndim、重心が原点）, eig = 固有値（n 個、大きい順）,
#              gof = 2次元説明率)
#   gof は cmdscale が GOF として返す2通りのうちの1つ目で、
#   「上位 ndim 個の固有値の和 / 固有値の絶対値の和」。本文で使っている値（論文4節）。
mds_apm <- function(APM, ndim = 2) {
  distances <- as.dist(APM)                            # 距離行列の型に
  result <- cmdscale(distances, k = ndim, eig = TRUE)  # 古典的 MDS
  points <- scale(result$points, scale = FALSE)        # 各列の平均を引き、重心を原点に
  eig    <- result$eig                                 # 固有値
  gof    <- result$GOF[1]                              # 2次元説明率
  return(list(points = points, eig = eig, gof = gof))
}

# SMACOF による MDS（smacof パッケージの mds）。ADD と TPD はこれを使う。
# ユークリッド距離行列である保証がない行列に対して、布置上の距離と非類似度のずれ
# （ストレス）を反復計算で最小にする。
# type = "ratio" は非類似度を定数倍するだけで、順位への置き換えなどの変換はしない。
# 乱数で初期値を決めるので、set.seed で種を固定して結果を再現可能にする。
# 引数 D：非類似度行列、ndim：次元数、seed：乱数の種。
# 戻り値：list(points = 座標（n×ndim）, stress = stress-1（0 に近いほど当てはまりがよい）)
mds_smacof <- function(D, ndim = 2, seed = 123) {
  distances <- as.dist(D)                              # 距離行列の型に
  set.seed(seed)
  result <- mds(distances, ndim = ndim, type = "ratio")
  points <- result$conf                                # 座標
  stress <- result$stress                              # stress-1
  return(list(points = points, stress = stress))
}

# ADD の MDS と TPD の MDS。どちらも SMACOF なので、中身は mds_smacof を呼ぶだけ。
# 名前を分けてあるのは、分析部で「ADD にはこれ、TPD にはこれ」と読めるようにするため。
mds_add <- function(ADD, ndim = 2, seed = 123) {
  result <- mds_smacof(ADD, ndim, seed)
  return(result)
}

mds_tpd <- function(TPD, ndim = 2, seed = 123) {
  result <- mds_smacof(TPD, ndim, seed)
  return(result)
}

# 布置の向きを揃える（プロクラステス回転）。
# MDS の座標は回転・鏡映・平行移動しても距離が変わらないので、向きは任意に決まる。
# 3枚の図を見比べやすくするため、ADD と TPD の布置を APM の布置に向きだけ合わせる。
# vegan の procrustes(基準, 合わせたい座標) を使う。scale = FALSE なので拡大縮小はしない
# （布置の形そのものは変わらない）。$Yrot が、向きを合わせたあとの座標。
# 引数 X_base：基準の座標、X：合わせたい座標。戻り値：向きを合わせた X。
align_to <- function(X_base, X) {
  fit     <- procrustes(X_base, X, scale = FALSE)   # 回転・鏡映・平行移動を求める
  aligned <- fit$Yrot                               # 合わせたあとの座標
  return(aligned)
}

# ----------------------------------------------------------------------------
#  2.4 布置の図（論文の図1〜3）。手法ごとに1関数。
# ----------------------------------------------------------------------------

# 座標 X（n×2、行名が国名）を1枚の散布図にする。plot_apm / plot_add / plot_tpd の共通部分。
# ggplot2 は「ggplot(データ, aes(x, y)) + 層 + 層 + …」と + でつないで図を組み立てる。
# coord_equal() で縦横の縮尺を同じにしている（点間の距離を目で比べるため）。
# seed は ggrepel のラベル配置用。同じ seed なら同じ配置になる。
# 引数 X：座標、title：図の見出し。戻り値：ggplot の図（print() すると描かれる）。
plot_configuration <- function(X, title, seed = 123) {
  # data.frame(列名 = 値, …) で作図用の表を作る。X[, 1] は座標の1列目（第1次元）。
  df <- data.frame(country = rownames(X),
                   Dim1    = X[, 1],
                   Dim2    = X[, 2])

  # 図を層ごとに組み立てる。
  figure <- ggplot(df, aes(Dim1, Dim2, label = country))             # x, y, ラベルの列を指定
  figure <- figure + geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3)  # 原点を通る横の補助線
  figure <- figure + geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3)  # 原点を通る縦の補助線
  figure <- figure + geom_point(size = 3.2, colour = "steelblue")   # 点
  figure <- figure + geom_text_repel(size = 4.2, seed = seed)       # 重ならない国名ラベル
  figure <- figure + coord_equal()                                   # 縦横の縮尺を同じに
  figure <- figure + labs(title = title, x = "Dimension 1", y = "Dimension 2")  # 見出しと軸ラベル
  figure <- figure + theme_minimal(base_size = 12)                   # 見た目のテーマ
  return(figure)
}

# 図1〜3。見出しを変えて plot_configuration を呼ぶだけ。
plot_apm <- function(X_apm) {
  figure <- plot_configuration(X_apm, "APM-MDS")
  return(figure)
}

plot_add <- function(X_add) {
  figure <- plot_configuration(X_add, "ADD-MDS")
  return(figure)
}

plot_tpd <- function(X_tpd) {
  figure <- plot_configuration(X_tpd, "TPD-MDS")
  return(figure)
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
# 引数：A（歪対称成分）、APM、ADD、TPD、share（割合）。いずれも n×n。
pair_table <- function(A, APM, ADD, TPD, share) {
  nm <- rownames(A)

  # 各ペアを1回ずつ取り出すため、上三角（行番号 < 列番号）の位置を使う。
  is_upper   <- upper.tri(APM)                      # 上三角の位置が TRUE の行列
  positions  <- which(is_upper, arr.ind = TRUE)     # TRUE の位置を (行, 列) の番号の表に（28行×2列）
  row_index  <- positions[, 1]                      # 行番号 = ペアの左側の国 i
  col_index  <- positions[, 2]                      # 列番号 = ペアの右側の国 j
  from       <- nm[row_index]                       # 国名に直す
  to         <- nm[col_index]

  # 行列を (行, 列) の表で添字付けすると、その位置の値の並びが取れる。
  a         <- A[positions]                         # 符号つきの a_ij
  apm_value <- APM[positions]
  add_value <- ADD[positions]
  tpd_value <- TPD[positions]
  share_pct <- 100 * share[positions]               # 割合を % に

  # 大きい向き。まず全部「均衡」にしておき、a の符号で書き換える。
  # x[条件] <- 値：条件が TRUE の要素だけを書き換える。
  larger <- rep("均衡", length(a))                              # rep(値, 個数)：同じ値を並べる
  larger[a > 0] <- paste(from[a > 0], "→", to[a > 0])          # a > 0：i から j への輸出が大きい
  larger[a < 0] <- paste(to[a < 0], "→", from[a < 0])          # a < 0：j から i への輸出が大きい

  pair <- paste(from, to, sep = "-")                            # "USA-CHN" のような文字列

  tbl <- data.frame(pair   = pair,
                    a      = a,
                    larger = larger,
                    APM    = apm_value,
                    ADD    = add_value,
                    TPD    = tpd_value,
                    share  = share_pct)

  # 割合の大きい順に並べ替える。order(x, decreasing = TRUE)：大きい順に並べたときの行番号。
  order_by_share <- order(tbl$share, decreasing = TRUE)
  tbl <- tbl[order_by_share, ]
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
A   <- skew_part(Delta)          # 歪対称成分（式(1)）
APM <- apm_from(A)               # 列ベクトル間距離（式(2)）
ADD <- add_from(A)               # |a_ij|
TPD <- tpd_from(A)               # 第三者項の平方根
check_identity(APM, ADD, TPD)    # 式(7)の検算。成り立たなければここで止まる
share <- share_from(A, APM)      # d_APM^2 に占める 2 a_ij^2 の割合

cat("\n=== 歪対称成分 A ===\n")
print(round(A, 3))
cat("\n=== APM（列ベクトル間距離）===\n")
print(round(APM, 3))
cat("\n=== ADD（|a_ij|）===\n")
print(round(ADD, 3))
cat("\n=== TPD（第三者項の平方根）===\n")
print(round(TPD, 3))
cat("\n=== d_APM^2 に占める 2 a_ij^2 の割合（%）===\n")
print(round(100 * share, 1))

# ----------------------------------------------------------------------------
#  3.3 MDS（論文4節〜6節）
# ----------------------------------------------------------------------------
# APM には古典的 MDS、ADD と TPD には SMACOF。
fit_apm <- mds_apm(APM, ndim = 2)
fit_add <- mds_add(ADD, ndim = 2)
fit_tpd <- mds_tpd(TPD, ndim = 2)

# 図を見比べやすいよう、ADD と TPD の布置の向きを APM の布置に揃える。
X_apm <- fit_apm$points
X_add <- align_to(X_apm, fit_add$points)
X_tpd <- align_to(X_apm, fit_tpd$points)

# 適合度。APM は2次元説明率と最小固有値（0 ならユークリッド距離行列）、
# ADD・TPD は stress-1（0 に近いほど当てはまりがよい）。
gof_pct    <- round(100 * fit_apm$gof, 1)     # 説明率を % に
min_eig    <- round(min(fit_apm$eig), 6)      # 最小固有値
stress_add <- round(fit_add$stress, 4)
stress_tpd <- round(fit_tpd$stress, 4)

cat("\n=== 適合度 ===\n")
cat("APM 2次元説明率 :", gof_pct, "%\n")
cat("APM 最小固有値  :", min_eig, "\n")
cat("ADD stress-1    :", stress_add, "\n")
cat("TPD stress-1    :", stress_tpd, "\n")

# ----------------------------------------------------------------------------
#  3.4 図1〜3
# ----------------------------------------------------------------------------
# ggplot の図はスクリプト実行時には print() しないと描かれない。1枚ずつ別の図として出す。
figure_apm <- plot_apm(X_apm)   # 図1
figure_add <- plot_add(X_add)   # 図2
figure_tpd <- plot_tpd(X_tpd)   # 図3
print(figure_apm)
print(figure_add)
print(figure_tpd)

# ----------------------------------------------------------------------------
#  3.5 全ペアの一覧と要約（論文7節）
# ----------------------------------------------------------------------------
tbl <- pair_table(A, APM, ADD, TPD, share)

# format(表, digits = 3, nsmall = 1)：数値を見やすい桁に整える。row.names = FALSE：行番号を出さない。
tbl_shown <- format(tbl, digits = 3, nsmall = 1)
cat("\n=== 全ペアの一覧（割合の大きい順）===\n")
print(tbl_shown, row.names = FALSE)

cat("\n=== 割合 2δ²ADD/d²APM の要約（28ペア）===\n")
cat("最小   :", round(min(tbl$share), 1), "%\n")
cat("最大   :", round(max(tbl$share), 1), "%\n")
cat("平均   :", round(mean(tbl$share), 1), "%\n")
cat("中央値 :", round(median(tbl$share), 1), "%\n")

# cor(x, y)：ピアソンの積率相関係数。
cor_add_apm <- cor(tbl$ADD, tbl$APM)
cor_tpd_apm <- cor(tbl$TPD, tbl$APM)
cor_add_tpd <- cor(tbl$ADD, tbl$TPD)

cat("\n=== 3つの行列どうしの相関（28ペア）===\n")
cat("ADD と APM :", round(cor_add_apm, 3), "\n")
cat("TPD と APM :", round(cor_tpd_apm, 3), "\n")
cat("ADD と TPD :", round(cor_add_tpd, 3), "\n")
