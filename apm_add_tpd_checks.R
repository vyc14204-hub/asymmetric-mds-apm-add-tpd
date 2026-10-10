# ============================================================================
#  歪対称行列に基づく非対称MDS：APM / ADD / TPD（確認計算）
#
#  本体 apm_add_tpd_mds.R を実行したあとに動かす。本文の次の記述の裏付けを出す。
#    1. 三角不等式の違反数と、二重中心化行列の最小固有値（論文5節・6節）
#         「ADD は28ペア中18ペアで三角不等式が破れ … 最小固有値の絶対値は最大固有値の 48.0%」
#         「TPD は三角不等式はすべて成り立つが … 16.5%」
#    2. ADD の stress-1 は次元を増やしても下がらない（論文5節）
#    3. APM の古典的 MDS は A の列の主成分分析と同じこと、A の行平均の順（論文5節末尾）
#    4. TPD 布置の重心からの距離
#
#  計算はできるだけ R 標準の関数（cmdscale, prcomp, svd）に任せ、自前で書くのは
#  同じことをする標準関数がない三角不等式の数え上げだけにしてある。
#  書き方の方針は本体と同じ：1行に1つの処理。入れ子の式は中間変数に分ける。
#
#  使う変数（本体で作られる）：A, APM, ADD, TPD（行列）、coords_xy_tpd（TPD の布置の座標）、nm（国名）、n（対象数）
#  本体がまだ実行されていなければ、ここで実行する。
#    exists("TPD")：TPD という変数があるか。! は否定。source(ファイル)：そのファイルを実行。
# ============================================================================

if (!exists("TPD")) {
  source("apm_add_tpd_mds.R", encoding = "UTF-8")
}


# ============================================================================
#  関数
# ============================================================================

# ----------------------------------------------------------------------------
#  三角不等式の違反を数える
# ----------------------------------------------------------------------------
# 距離の公理のひとつ、三角不等式 D[i,j] <= D[i,l] + D[l,j] が破れている三つ組を数える。
# 順序三つ組 (i, j, l) をすべて調べるので、8か国なら 8*7*6 = 336 通り。
# 違反のうち超過分 gap が最も大きい三つ組も返す（本文の「中国とインドの 0.992 は
# 中国とフランスの値とフランスとインドの値の和 0.241 を上回る」はこれ）。
# 引数 D：対称な非類似度行列。戻り値：list(count = 違反数, i, j, l = 最悪の三つ組の番号)。
# 使い方：v <- tri_violations(ADD)
#         違反数は v$count、最悪の三つ組の番号は v$i, v$j, v$l で取り出す。
tri_violations <- function(D) {
  n <- nrow(D)
  count     <- 0     # 違反の数
  worst_gap <- 0     # これまでで最大の超過分
  worst_i   <- NA    # 最大の超過分を出した三つ組。NA は「まだない」
  worst_j   <- NA
  worst_l   <- NA

  # for を3つ重ねて、i, j, l のすべての組み合わせを回す。
  for (i in 1:n) {
    for (j in 1:n) {
      for (l in 1:n) {
        if (i == j || j == l || i == l) {       # || は「または」、== は「等しい」
          next                                   # 3つが異なる組だけ調べる。next で次の回へ
        }
        direct  <- D[i, j]                       # i から j へ直接
        via_l   <- D[i, l] + D[l, j]             # l を経由
        gap     <- direct - via_l                # 正なら三角不等式が破れている
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

  return(list(count = count, i = worst_i, j = worst_j, l = worst_l))
}

# ----------------------------------------------------------------------------
#  三角不等式が破れるペアを数える（無順序のペア単位）
# ----------------------------------------------------------------------------
# 本文の「28ペア中18ペア」は、無順序のペア (i, j) のうち、どれかの l で三角不等式が
# 破れているものの数。tri_violations は順序三つ組で数える（ADD では 92）ので、
# ペア単位で数え直す。
# 引数 D：対称な非類似度行列。戻り値：破れているペアの数。
# 使い方：n_violated_pairs <- count_violated_pairs(ADD)
count_violated_pairs <- function(D) {
  n <- nrow(D)
  pairs <- combn(n, 2)                           # 全ペア (i, j), i < j（2×28 の行列）
  violated_count <- 0

  for (p in seq_len(ncol(pairs))) {
    i <- pairs[1, p]
    j <- pairs[2, p]
    others <- setdiff(seq_len(n), c(i, j))       # 経由点 l の候補：i と j 以外
    via_each_l <- D[i, others] + D[others, j]    # l ごとの D[i,l] + D[l,j]（ベクトル）
    shortest_via <- min(via_each_l)              # 最も短い経由
    if (D[i, j] > shortest_via + 1e-12) {        # 直接のほうが長ければ、このペアは破れている
      violated_count <- violated_count + 1
    }
  }

  return(violated_count)
}

# ----------------------------------------------------------------------------
#  二重中心化行列の固有値
# ----------------------------------------------------------------------------
# 二重中心化行列 B = -1/2 J D^2 J の固有値（J = I - 11'/n は中心化行列）。
# D がユークリッド距離行列であることは、B が半正定値（固有値がすべて非負）であることと
# 同値（論文2節）。最小固有値が負なら、何次元に埋め込んでもその距離は再現できない。
# この固有値は古典的 MDS の計算そのものなので、自前で中心化せず、R 標準の cmdscale() に
# eig = TRUE を付けて返してもらう。$eig に n 個の固有値が大きい順に入る（k の値によらない）。
# 引数 D：対称な非類似度行列。戻り値：固有値（大きい順）。
# 使い方：ev <- dc_eigen(TPD)　　最小固有値は min(ev)、最大固有値は max(ev)。
dc_eigen <- function(D) {
  distances <- as.dist(D)                                  # 距離行列の型に
  result    <- cmdscale(distances, k = 2, eig = TRUE)      # 古典的 MDS
  ev        <- result$eig                                  # 固有値
  return(ev)
}


# ============================================================================
#  確認計算
# ============================================================================

# ----------------------------------------------------------------------------
#  1. 三角不等式とユークリッド性（論文5節・6節）
# ----------------------------------------------------------------------------
# 三角不等式を満たすこと（距離であること）と、ユークリッド距離行列であることは別の
# 性質で、前者は後者を含意しない。APM は両方満たし、TPD は三角不等式は満たすが
# ユークリッドではなく、ADD は三角不等式も破る、という3段階を数値で確かめる。
# for (label in c(…))：label に "APM", "ADD", "TPD" を順に入れて繰り返す。
# matrices[[label]]：リストから、変数 label に入っている名前の要素を取り出す
# （$ は名前を直接書くとき、[[ ]] は名前が変数に入っているとき）。
matrices <- list(APM = APM, ADD = ADD, TPD = TPD)   # 3つの行列を名前付きリストに
cat("\n=== 三角不等式とユークリッド性 ===\n")
for (label in c("APM", "ADD", "TPD")) {
  D <- matrices[[label]]

  # 三角不等式
  v <- tri_violations(D)

  # 二重中心化行列の固有値
  ev        <- dc_eigen(D)
  min_eig   <- min(ev)                       # 最小固有値（負ならユークリッドでない）
  max_eig   <- max(ev)                       # 最大固有値
  ratio_pct <- 100 * abs(min_eig) / max_eig  # 最小固有値の絶対値 / 最大固有値（%）

  # sprintf：書式の %s に文字列、%d に整数、%.4f に小数4桁の数値が順に入る。%% は % そのもの。
  line <- sprintf("%s : 三角不等式違反 %d 組 / 最小固有値 %.4f（最大固有値比 %.1f%%）",
                  label, v$count, min_eig, ratio_pct)
  cat(line, "\n")

  # 最悪の三つ組（違反があるときだけ）
  if (v$count > 0) {
    name_i <- nm[v$i]
    name_j <- nm[v$j]
    name_l <- nm[v$l]
    direct <- D[v$i, v$j]
    via_l  <- D[v$i, v$l] + D[v$l, v$j]
    worst  <- sprintf("     最悪 : %s-%s = %.3f > %s-%s + %s-%s = %.3f",
                      name_i, name_j, direct, name_i, name_l, name_l, name_j, via_l)
    cat(worst, "\n")
  }
}

n_violated_pairs <- count_violated_pairs(ADD)
cat("ADD で三角不等式が破れるペアの数（無順序）:", n_violated_pairs, "\n")

# ----------------------------------------------------------------------------
#  2. ADD の適合の悪さは次元不足ではない（論文5節）
# ----------------------------------------------------------------------------
# 次元を 2 から n-1 まで増やしても stress-1 が下がらないことを示す。
# 三角不等式を破る値は何次元の布置でも再現できないため。
cat("\n=== ADD の次元別 stress-1 ===\n")
for (k in 2:(n - 1)) {
  fit_k <- mds(as.dist(ADD), ndim = k, type = "ratio")
  line  <- sprintf("  %d次元 : %.4f", k, fit_k$stress)
  cat(line, "\n")
}

# ----------------------------------------------------------------------------
#  3. APM の古典的 MDS の構造と A の行平均（論文5節末尾）
# ----------------------------------------------------------------------------
# APM は A の列ベクトル間のユークリッド距離なので、APM の古典的 MDS は、A の列ベクトル
# （各国を1つの点とみなす）を中心化して主成分分析するのと同じになる。
# これを、R 標準の prcomp()（主成分分析）と cmdscale()（古典的 MDS）の結果を突き合わせて
# 確かめる。prcomp(t(A)) は A の列（国）を観測、行（相手国）を変数とした主成分分析。
# 主成分の分散 sdev^2 に (n - 1) を掛けたものが、古典的 MDS の固有値に一致するはず。
# 主成分分析側
columns_as_rows <- t(A)                                        # A の列（国）を行に
pca <- prcomp(columns_as_rows, center = TRUE, scale. = FALSE)  # 主成分分析
pc_variances <- pca$sdev^2                                     # 各主成分の分散
eig_from_pca <- pc_variances * (n - 1)                         # 分散 × (n-1) = 固有値

# 古典的 MDS 側（同じ個数だけ取り出す）
eig_all      <- dc_eigen(APM)                                  # n 個の固有値
how_many     <- length(eig_from_pca)
eig_from_mds <- eig_all[seq_len(how_many)]

# 突き合わせ
max_gap <- max(abs(eig_from_pca - eig_from_mds))

# paste(…, collapse = " ")：並びを空白区切りの1つの文字列に。
cat("\n=== APM の構造 ===\n")
cat("  主成分分析の固有値 :", paste(round(eig_from_pca, 4), collapse = " "), "\n")
cat("  古典的 MDS の固有値:", paste(round(eig_from_mds, 4), collapse = " "), "\n")
cat("  両者の最大差       :", format(max_gap, digits = 3), "\n")

# 歪対称行列の特異値は対で現れる。svd(A)$d：特異値分解の特異値。
singular_values <- svd(A)$d
cat("  A の特異値（対で現れる）:", paste(round(singular_values, 3), collapse = " "), "\n")

# 行平均 rowMeans(A)：その国が相手全体に対して輸出超過（正）か輸入超過（負）かの目安。
# 図2の横軸の並び順と一致する。
row_means <- rowMeans(A)
cat("  A の行平均（正なら相手全体に対して輸出超過）\n")
print(round(row_means, 3))

# ----------------------------------------------------------------------------
#  4. TPD 布置の重心からの距離
# ----------------------------------------------------------------------------
# 中心に近い国ほど、他国に対する偏りの向きと大きさが「平均的」であることを示す。
# coords_xy_tpd は本体で作った TPD の布置の座標（n×2）。
coords_xy_tpd_centered <- scale(coords_xy_tpd, scale = FALSE)   # 重心を原点に
squared          <- coords_xy_tpd_centered^2                    # 座標の二乗
sum_of_squares   <- rowSums(squared)                            # 各国の二乗和
dist_from_center <- sqrt(sum_of_squares)                        # 平方根 = 原点からの距離
names(dist_from_center) <- nm                                   # 並びに国名を付ける
dist_sorted <- sort(dist_from_center)                           # 小さい順

cat("\n=== TPD 布置の重心からの距離 ===\n")
print(round(dist_sorted, 3))
