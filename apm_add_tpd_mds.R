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
#        2. 関数定義（分解、3つの行列と検算、MDS、作図、確認用）
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
#    式(3)/(7)   d_APM^2 = 2 δ_ADD^2 + δ_TPD^2          → 検算で確かめる
#
#  ---------------------------------------------------------------------------
#  このスクリプトで使う R の文法メモ（初出の箇所にも同じ説明を付けてある）
#  ---------------------------------------------------------------------------
#    x <- 値            代入。「x に値を入れる」。= でも書けるが R では <- が慣例
#    #                  行末までコメント。R は無視する
#    f(a, b = 2)        関数呼び出し。b = 2 のように引数名を付けて渡せる
#    function(x) {...}  関数を作る。{ } の中の最後の値、または return(値) が戻り値
#    c(1, 2, 3)         ベクトル（1次元の並び）を作る。c は combine の c
#    1:n, seq_len(n)    1 から n までの整数の並び。n が 0 のとき 1:n は 1, 0 になるので
#                       seq_len(n) のほうが安全
#    x[i]               ベクトルの i 番目。行列なら M[i, j] で i 行 j 列
#    M[i, ]  M[, j]     行列の第 i 行、第 j 列（空けた側は「全部」）
#    M[-1]  M[, -1]     マイナスは「それ以外」。M[, -1] は1列目を除いた全列
#    list(a = 1, b = 2) 名前付きリスト。種類の違うものをまとめて入れる箱
#    lst$a              リストの要素 a を取り出す。lst[["a"]] と同じ
#    t(M)               転置。行と列を入れ替える
#    M %*% N            行列の積。* だけだと要素ごとの掛け算になる
#    M^2, abs(M)        行列の要素ごとの二乗、絶対値
#    diag(M) <- 0       対角成分に 0 を入れる。diag(n) は n 次の単位行列
#    nrow(M), ncol(M)   行数、列数
#    for (i in 1:n) {}  繰り返し。i に 1, 2, …, n を順に入れて { } を実行
#    if (条件) {} else {}  条件分岐。next は「この回は飛ばして次へ」
#    &&, ||, !          かつ、または、否定（条件式で使う）
#    any(x), all(x)     TRUE/FALSE の並び x に、1つでも TRUE があるか / 全部 TRUE か
#    cat(...)           文字や数値をそのまま画面に出す。"\n" は改行
#    print(x)           x を R の標準の形で画面に出す
#    sprintf("%.3f", x) 書式付きの文字列を作る。%d 整数、%s 文字列、%.3f 小数3桁
#    paste(a, b)        文字列をつなぐ（間に空白）。sep = "-" で区切りを変えられる
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

# read.csv(ファイル名)：CSV を読み込んで data.frame（表）にする。
# fileEncoding = "UTF-8-BOM" は、先頭に BOM の付いた UTF-8 ファイルを正しく読む指定。
# CSV の1列目は国名、2列目以降は各輸出先への輸出額。
# 対角（自国への輸出）は空欄で、NA（欠損）として読み込まれる。
trade <- read.csv("data/trade_2023_comtrade.csv", fileEncoding = "UTF-8-BOM")

# rownames(x) <- 値：表の行名を付け替える。trade$country は列 country（国名）。
# trade[, -1] は1列目（国名）を除いた全列。data.matrix() はそれを数値の行列にする。
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
#  2.1 対称成分と歪対称成分への分解（論文の式(1)）
# ----------------------------------------------------------------------------
# 非対称行列 M を対称成分 S と歪対称成分 A に分ける。
#   S = (M + M') / 2 … 往復の平均。本研究では使わない
#   A = (M - M') / 2 … 往復の差の半分。以下の計算はすべてこの A だけを使う
# 対角は無視する。先に 0 を入れておくが、A の対角は M の対角に何が入っていても
# (m_ii - m_ii) / 2 = 0 になるので、結果には影響しない。
# 引数 M：n×n の数値行列。戻り値：list(S, A)。
skew_decompose <- function(M) {
  M <- as.matrix(M)        # data.frame で渡されても行列に直す
  diag(M) <- 0             # 対角成分を 0 にする
  S <- (M + t(M)) / 2      # t(M) は転置。行列どうしの + と / 2 は要素ごとに計算される
  A <- (M - t(M)) / 2

  # 2つの行列を名前付きリストにまとめて返す。呼び出し側では $S, $A で取り出す。
  return(list(S = S, A = A))
}

# ----------------------------------------------------------------------------
#  2.2 APM・ADD・TPD の3つの非類似度行列（論文2節）
#      小さな関数に分けてある。apm_from / add_from / tpd_from が1つずつ行列を作り、
#      check_identity が恒等式を検算し、share_from が割合を作る。
#      最後の apm_matrices はそれらを順に呼んでまとめるだけ。
# ----------------------------------------------------------------------------

# APM（式(2)）：A の列ベクトルどうしのユークリッド距離。
# dist(X) は X の「行」どうしの距離をまとめて計算するので、t(A) で転置して
# 列を行にしてから渡す。dist の結果は下三角だけの特殊な形なので、
# as.matrix() で普通の n×n 行列に直す。
# 引数 A：歪対称成分。戻り値：n×n の APM 行列。
apm_from <- function(A) {
  return(as.matrix(dist(t(A))))
}

# ADD：対象ペア自身に対応する要素の絶対値 |a_ij|。abs() は要素ごとの絶対値。
# 符号を落とすので、どちら向きの輸出が大きいかは A の符号を見る必要がある（論文5節）。
# 引数 A：歪対称成分。戻り値：n×n の ADD 行列。
add_from <- function(A) {
  return(abs(A))
}

# TPD：式(3)右辺第2項を定義どおりに計算する。恒等式 d_APM^2 = 2 δ_ADD^2 + δ_TPD^2
# から TPD^2 = APM^2 - 2 ADD^2 と逆算すれば1行で済むが、それだと check_identity が
# 恒等式を恒等式で確かめることになり、検算にならない。
# 引数 A：歪対称成分。戻り値：n×n の TPD 行列。
tpd_from <- function(A) {
  n <- nrow(A)

  # 1ペア (i, j) の TPD を計算する小さな関数。関数の中で関数を定義すると、
  # 外側の変数（ここでは A と n）をそのまま参照できる。
  tpd_pair <- function(i, j) {
    others <- setdiff(seq_len(n), c(i, j))   # setdiff(全体, 除くもの)：i と j 以外の対象 k
    diff_k <- A[others, i] - A[others, j]    # A[行の並び, 列]：k ごとの偏りの差 a_ki - a_kj
    sqrt(sum(diff_k^2))                      # 差の二乗和の平方根 = δ_TPD(i, j)。これが戻り値
  }

  # 全ペアについて tpd_pair を呼んで行列に埋める。TPD は対称なので、
  # i < j のペアを1回ずつ計算して (i, j) と (j, i) の両方に入れる。対角は 0 のまま。
  # matrix(0, n, n)：全要素 0 の n×n 行列。dimnames = で行名・列名を A と同じにする。
  # combn(n, 2)：1〜n から2つ選ぶ組み合わせを全部列挙し、1組を1列にした 2×28 の行列。
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
# 左辺と右辺の差は理論上 0 だが、小数計算の丸め誤差がわずかに残るので、
# 差が 1e-10（0.0000000001）未満なら「等しい」とみなす。
# any(条件)：1つでも条件を満たす要素があれば TRUE。
# stop(メッセージ)：エラーを出して実行を止める。
# 引数：3つの行列。戻り値：なし（成り立たなければエラーで止まる）。
# invisible(NULL) は「返すものはない」を明示する書き方。
check_identity <- function(APM, ADD, TPD) {
  identity_gap <- abs(APM^2 - (2 * ADD^2 + TPD^2))
  if (any(identity_gap >= 1e-10)) {
    stop("式(7)の恒等式が成り立っていません。最大のずれ: ", max(identity_gap))
  }
  invisible(NULL)
}

# d_APM^2 のうち a_ij に由来する部分 2 a_ij^2 が占める割合（0〜1。論文7節）。
# 行列どうしの / は要素ごとの割り算。対角は 0/0 で意味がないので NA にする。
# 引数 A：歪対称成分、APM：APM 行列。戻り値：n×n の割合の行列。
share_from <- function(A, APM) {
  share <- 2 * A^2 / APM^2
  diag(share) <- NA
  return(share)
}

# 上の関数を順に呼んで、3つの行列と割合をひとまとめにする。
# 引数 M：非対称行列。戻り値：list(A, APM, ADD, TPD, share)。いずれも n×n。
apm_matrices <- function(M) {
  A   <- skew_decompose(M)$A   # 歪対称成分だけ取り出す
  APM <- apm_from(A)
  ADD <- add_from(A)
  TPD <- tpd_from(A)
  check_identity(APM, ADD, TPD)
  share <- share_from(A, APM)

  # 5つの行列を名前付きリストにまとめて返す（左が名前、右が中身）。
  # 呼び出し側では mats <- apm_matrices(M) と受け取り、mats$APM, mats$TPD のように使う。
  return(list(A = A, APM = APM, ADD = ADD, TPD = TPD, share = share))
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
#   ndim = 2 のように書いた引数は「省略したときの既定値」。apm_mds(M) と呼べば 2 になる。
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

  # APM：古典的 MDS。as.dist() は行列を dist 型（距離行列であることを示す型）に変える。
  # cmdscale はこの型を受け取る。eig = TRUE で固有値も一緒に返してもらう。
  fit_apm <- cmdscale(as.dist(mats$APM), k = ndim, eig = TRUE)

  # ADD・TPD：SMACOF（smacof パッケージの mds 関数）。type = "ratio" は非類似度を
  # 定数倍するだけで、順位への置き換えなどの変換はしない（論文3.3節・5節）。
  # SMACOF は乱数で初期値を決めるので、set.seed で種を固定して結果を再現可能にする。
  set.seed(seed)
  fit_add <- mds(as.dist(mats$ADD), ndim = ndim, type = "ratio")
  set.seed(seed)
  fit_tpd <- mds(as.dist(mats$TPD), ndim = ndim, type = "ratio")

  # APM の座標 fit_apm$points は、scale(scale = FALSE) で各列の平均を引き、
  # 重心が原点になるよう中心化する。これを基準にして他の2つを回転で合わせる。
  # procrustes(基準, 合わせたい座標)$Yrot が、回転後の座標。
  X_apm <- scale(fit_apm$points, scale = FALSE)
  X_add <- procrustes(X_apm, fit_add$conf, scale = FALSE)$Yrot
  X_tpd <- procrustes(X_apm, fit_tpd$conf, scale = FALSE)$Yrot

  # 作図用に3つの布置を1つの data.frame にまとめる。
  # data.frame(列名 = 値, …) で表を作る。X[, 1] は座標行列の1列目（第1次元）。
  # rbind() は表を縦につなぐ。
  as_df <- function(X, method) {
    data.frame(method = method, country = nm, Dim1 = X[, 1], Dim2 = X[, 2],
               row.names = NULL)
  }
  conf <- rbind(as_df(X_apm, "APM"), as_df(X_add, "ADD"), as_df(X_tpd, "TPD"))

  # factor(…, levels = …)：文字列を「順序の決まったカテゴリ」にする。
  # ここでは図を APM, ADD, TPD の順に並べるため。
  conf$method <- factor(conf$method, levels = c("APM", "ADD", "TPD"))

  # 2次元説明率：上位 ndim 個の固有値の和を、固有値の絶対値の和で割る。
  # ev[1:ndim] は固有値ベクトルの先頭 ndim 個。
  # 負の固有値がなければ分母の取り方で値が変わることはない（論文4節）。
  ev  <- fit_apm$eig
  gof <- sum(ev[1:ndim]) / sum(abs(ev))

  # 結果を名前付きリストにまとめて返す。呼び出し側では fit$conf, fit$stress のように使う。
  # c(ADD = …, TPD = …) は名前付きのベクトル。fit$stress["ADD"] で取り出せる。
  return(list(matrices = mats,
              conf     = conf,
              eig      = ev,
              gof      = gof,
              stress   = c(ADD = fit_add$stress, TPD = fit_tpd$stress),
              fits     = list(apm = fit_apm, add = fit_add, tpd = fit_tpd)))
}

# ----------------------------------------------------------------------------
#  2.4 布置の図（論文の図1〜3）
# ----------------------------------------------------------------------------
# methods に1つ渡せば1枚、省略すれば3枚を横に並べる。
# ggplot2 は「ggplot(データ, aes(x, y)) + 層 + 層 + …」と + でつないで図を組み立てる。
# coord_equal() で縦横の縮尺を同じにしている（点間の距離を目で比べるため）。
# seed は ggrepel のラベル配置用。同じ seed なら同じ配置になる。
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
#   a      … 符号つきの a_ij。正なら m_ij > m_ji、すなわち i から j への輸出のほうが大きい
#   larger … 大きい向きを "X → Y" で示したもの（m_XY > m_YX）。ペアを書く順に依らず向きが読める
#   share  … 100 * 2 a_ij^2 / d_APM^2
# 引数 mats：apm_matrices() の結果。戻り値：data.frame（1ペア1行）。
pair_table <- function(mats) {
  nm <- rownames(mats$A)

  # upper.tri(M)：上三角（行番号 < 列番号）の位置が TRUE の行列。
  # which(…, arr.ind = TRUE)：TRUE の位置を (行, 列) の番号の表にして返す。
  # これで各ペアを1回ずつ取り出せる（28行 × 2列）。
  u <- which(upper.tri(mats$APM), arr.ind = TRUE)
  from <- nm[u[, 1]]    # ペアの左側の国 i（u の1列目は行番号 = i）
  to   <- nm[u[, 2]]    # ペアの右側の国 j（u の2列目は列番号 = j）
  a    <- mats$A[u]     # 行列を (行, 列) の表で添字付けすると、その位置の値の並びが取れる

  # 大きい向き。ifelse(条件, 条件が真のときの値, 偽のときの値) を要素ごとに評価する。
  # a > 0 なら i → j、a < 0 なら j → i、0 なら均衡。
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
  return(tbl[order(-tbl$share), ])
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
  worst_i <- NA     # 最大の超過分を出した三つ組。NA は「まだない」
  worst_j <- NA
  worst_l <- NA

  # for を3つ重ねて、i, j, l のすべての組み合わせを回す。
  for (i in 1:n) {
    for (j in 1:n) {
      for (l in 1:n) {
        if (i == j || j == l || i == l) {       # || は「または」、== は「等しい」
          next                                   # 3つが異なる組だけ調べる。next で次の回へ
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

  # 違反の数と、最悪の三つ組の番号を名前付きリストにまとめて返す。
  return(list(count = count, i = worst_i, j = worst_j, l = worst_l))
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
  J <- diag(n) - matrix(1 / n, n, n)       # diag(n) は単位行列、matrix(1/n, n, n) は全要素 1/n
  B <- -0.5 * J %*% (M^2) %*% J            # %*% は行列の積。M^2 は要素ごとの二乗
  # eigen() は固有値と固有ベクトルを求める。symmetric = TRUE で対称行列向けの計算、
  # only.values = TRUE で固有値だけ求める。$values でその固有値（大きい順）を取り出す。
  return(eigen(B, symmetric = TRUE, only.values = TRUE)$values)
}


# ============================================================================
#  3. 分析
#     ここから先が実際に動く部分。上で定義した関数を順に呼ぶ。
# ============================================================================

# ----------------------------------------------------------------------------
#  3.1 対数変換（論文3.2節・表2・式(8)）
# ----------------------------------------------------------------------------
# 輸出額の対数をそのまま非対称行列として渡す（非類似度への反転はしない）。
# このとき歪対称成分は a_ij = (1/2) log(T_ij / T_ji) となり（式(9)）、往復の比だけで
# 決まる。単位を変えても、対数に定数を足しても a_ij は変わらない。
# 符号を反転した -log(T) を渡しても A の符号が変わるだけで、APM・ADD・TPD は不変。
# log() は要素ごとの自然対数。対角は NA のままだと log で NA になるので 0 にしておく
# （A の対角には影響しない）。
M <- log(exports)
diag(M) <- 0

cat("\n=== 輸出額の対数（表2）===\n")
print(round(M, 3))

# ----------------------------------------------------------------------------
#  3.2 3つの行列と MDS
# ----------------------------------------------------------------------------
# apm_mds() が中で apm_matrices() も呼ぶので、この1行で3つの行列と3つの布置が揃う。
fit  <- apm_mds(M, ndim = 2)
mats <- fit$matrices          # 3つの行列（A, APM, ADD, TPD, share のリスト）を短い名前で持つ

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
# fit$stress["ADD"] は名前付きベクトルから名前で値を取り出す書き方。
cat("\n=== 適合度 ===\n")
cat("APM 2次元説明率 :", round(100 * fit$gof, 1), "%\n")
cat("APM 最小固有値  :", round(min(fit$eig), 6), "\n")
cat("ADD stress-1    :", round(fit$stress["ADD"], 4), "\n")
cat("TPD stress-1    :", round(fit$stress["TPD"], 4), "\n")

# ----------------------------------------------------------------------------
#  3.4 図1〜3
# ----------------------------------------------------------------------------
# ggplot の図は print() しないと（スクリプト実行時には）描かれない。
# 1枚ずつ出したものを論文に使った。最後の1行は3枚を横に並べた確認用。
print(plot_apm(fit, "APM"))
print(plot_apm(fit, "ADD"))
print(plot_apm(fit, "TPD"))
print(plot_apm(fit))                       # 3枚並べ

# ----------------------------------------------------------------------------
#  3.5 全ペアの一覧と要約（論文7節）
# ----------------------------------------------------------------------------
# format(表, digits = 3, nsmall = 1)：数値を有効数字3桁・小数点以下1桁以上に整えて表示用にする。
# row.names = FALSE：行番号を表示しない。
tbl <- pair_table(mats)
cat("\n=== 全ペアの一覧（割合の大きい順）===\n")
print(format(tbl, digits = 3, nsmall = 1), row.names = FALSE)

# 各ペアを1回ずつ数えるため、上三角だけを取り出す。
# upper.tri() で作った TRUE/FALSE の行列で添字付けすると、TRUE の位置の値だけが並びで取れる。
upper <- upper.tri(mats$APM)
share_pct <- 100 * mats$share[upper]

cat("\n=== 割合 2δ²ADD/d²APM の要約 ===\n")
cat("最小   :", round(min(share_pct), 1), "%\n")
cat("最大   :", round(max(share_pct), 1), "%\n")
cat("平均   :", round(mean(share_pct), 1), "%\n")
cat("中央値 :", round(median(share_pct), 1), "%\n")

# cor(x, y)：ピアソンの積率相関係数。
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
# for (label in c(…))：label に "APM", "ADD", "TPD" を順に入れて繰り返す。
# mats[[label]]：リストから、変数 label に入っている名前の要素を取り出す（$ は名前を
# 直接書くとき、[[ ]] は名前が変数に入っているとき）。
cat("\n=== 三角不等式とユークリッド性 ===\n")
for (label in c("APM", "ADD", "TPD")) {
  D  <- mats[[label]]
  v  <- tri_violations(D)
  ev <- dc_eigen(D)
  ratio_pct <- 100 * abs(min(ev)) / max(ev)   # 最小固有値の絶対値 / 最大固有値（%）

  # sprintf：書式の中の %s に文字列、%d に整数、%.4f に小数4桁の数値が順に入る。
  # %% は % の文字そのもの。
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
B_from_A   <- J %*% (t(A) %*% A) %*% J           # A の列を中心化した内積行列（t(A) %*% A = A'A）

# svd(A)$d：特異値分解の特異値。paste(…, collapse = " ")：並びを空白区切りの1つの文字列に。
# rowMeans(A)：各行の平均。
cat("\n=== APM の構造 ===\n")
cat("  B と J(A'A)J の最大差 :", format(max(abs(B_from_apm - B_from_A)), digits = 3), "\n")
cat("  A の特異値（対で現れる）:", paste(round(svd(A)$d, 3), collapse = " "), "\n")
cat("  A の行平均（正なら相手全体に対して輸出超過）\n")
print(round(rowMeans(A), 3))

# ----------------------------------------------------------------------------
#  3.9 TPD 布置の重心からの距離（確認用）
# ----------------------------------------------------------------------------
# 中心に近い国ほど、他国に対する偏りの向きと大きさが「平均的」であることを示す。
# subset(…, method == "TPD")[, c("Dim1", "Dim2")]：TPD の行だけ取り、座標の2列だけ残す。
# rowSums(X^2)：各行の二乗和。names(x) <- nm で並びに国名を付ける。sort() は小さい順。
X_tpd <- as.matrix(subset(fit$conf, method == "TPD")[, c("Dim1", "Dim2")])
X_tpd_centered <- scale(X_tpd, scale = FALSE)                 # 重心を原点に
dist_from_center <- sqrt(rowSums(X_tpd_centered^2))           # 各国の原点からの距離
names(dist_from_center) <- nm

cat("\n=== TPD 布置の重心からの距離 ===\n")
print(round(sort(dist_from_center), 3))
