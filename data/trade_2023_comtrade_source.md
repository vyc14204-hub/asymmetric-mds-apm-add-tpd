# trade_2023_comtrade.csv の出典

## 書誌情報（論文に書く形）

United Nations Statistics Division. *UN Comtrade Database*.
https://comtrade.un.org/ (accessed 18 August 2026).

## 取得条件

| 項目 | 値 |
|---|---|
| エンドポイント | `https://comtradeapi.un.org/public/v1/preview/C/A/HS` |
| 取得日時 (UTC) | 2026-08-18（`comtrade_2023_raw.json` の `pulled_utc` に秒まで記録） |
| 年 | 2023 |
| 貿易フロー | `flowCode=X`（**輸出**） |
| 品目 | `cmdCode=TOTAL`（全品目合計） |
| 分類 | HS |
| 頻度 | A（年次） |
| 単位 | 10億米ドル（小数第1位に丸め）。論文の表1と同一の値であり、分析はこの値に対して行う。取得時の生値は名目米ドル |
| 対象 | 財のみ（サービスを含まない） |

## 報告国の統一

**行＝輸出国、列＝輸出先。全セルが「輸出国自身の報告による輸出額」。**
Comtrade の reporter を輸出国、flow を X に固定して取得したので、
全方向が **FOB 建て**で揃っている。輸入側報告（CIF、運賃・保険料込み）は一切混ざっていない。

これは本研究にとって死活的に重要である。CIF/FOB を混ぜると評価基準の差が
そのまま見かけの非対称性 $a_{ij}$ に化けるため、非対称性そのものを分析対象とする
本研究では結果が汚染される。

## 行の絞り込み

Comtrade は同一の reporter–partner–year に対し、輸送手段（`motCode`）、
第二パートナー（`partner2Code`）、税関手続（`customsCode`）別の内訳行を返す。
合計行を一意に取るため、以下で絞り込んだ：

```
customsCode == "C00"  かつ  motCode == 0  かつ  partner2Code == 0
```

## 国コード（Comtrade M49）

USA 842 / CHN 156 / DEU 276 / JPN 392 / KOR 410 / GBR 826 / FRA 251 / IND 699

注：France は 250 ではなく **251**、India は 356 ではなく **699**。

## 検証

USA→CHN = 147,805,519,996 USD。米国の対中財輸出2023年の公表値（約1,478億ドル）と一致。
CHN→DEU = 100,570,000,000 USD 規模。UN Comtrade の公表値（約1,005.7億ドル）と一致。

## 旧ファイル `trade_2023_8economies.csv` について

**出典不明のため使用しない。** 2026-08-17 に生成されたが、CSV・R Canvas ワークブック
(`非対称.rfa1`) の内部 sqlite・Downloads のいずれにも出所の記録がなく、
LLM の記憶由来の可能性が高い。値は概ね正しいが整数に丸められており、
実データとの乖離が大きいセルがある（下表、10億ドル）：

| 方向 | 旧CSV | Comtrade実データ | 乖離 |
|---|---|---|---|
| DEU→FRA | 94 | 130.3 | −28% |
| USA→DEU | 89 | 76.5 | +16% |
| GBR→USA | 60 | 71.9 | −17% |
| DEU→USA | 161 | 173.1 | −7% |
| CHN→DEU | 107 | 100.6 | +6% |

とくに旧CSVは DEU→CHN と CHN→DEU を**ともに 107** としており、
$a_{ij}=0$ という完全な対称を人工的に作っていた。実データでは 106.9 対 100.6 で差がある。
$a_{ij}=0$ は本研究が主題とする現象そのものであるため、この丸めは致命的だった。
旧ファイルは経緯の記録として残すが、分析には用いない。
