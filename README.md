# mult_SVD — SVD + lifting による多変数近似 GCD の実験

Maple code for the paper "An SVD-and-Lifting Based Framework for Approximate Multivariate GCD Computation"
(M. Sanuki). It contains the implementation (`svdlift.mpl`), the examples and experiment scripts that
produce every table of the paper (`exp_*.mpl`, output in `tables/`), and the scripts that check the
statements of the paper against the computation (`verify_claims.mpl`, `verify_section4.py`, with logs).

以下は日本語の説明。論文（SVD + lifting による多変数近似 GCD）の実験コード一式（2026-11-23 作成）。

## 使用例

```
read "common.mpl":    # svdlift.mpl, appgcd.mpl も読み込まれる
tv := [t1, t2]:
F := expand(C26*A26): G := expand(C26*B26):        # 論文の Example 1
Cc, Ft, Gt, sig := SVDLift:-Solve(F, G, tv, 3, 4, x, jstar = 1):   # k = 3, W = 4
SVDLift:-RelErr(Cc, C26, [x, t1, t2]);              # GCD の相対誤差
F2, G2, At, Bt := Ex2(1.0003, 1.0005):              # Example 2（初期因子が近似共通因子をもつ）
SVDLift:-Solve(F2, G2, tv, 3, 6, x, mode = "ill", jstar = 1)[1];   # 2 次元自由空間リフティング
```

## 実行

```
"/Applications/Maple 2025/maple" -q run_all.mpl
```

表が `tables/*.tex` に書き出される。環境変数 `MULTSVD_PAPERTAB` に論文の `experiments/tables/` を指定するとそこにも複製される
（論文は `\exptable{...}` でそれを `\input` する）。個別に回すときは
`read "exp_examples.mpl":` などをワークシートで読み込めばよい。

## ファイル

| ファイル | 内容 |
|---|---|
| `svdlift.mpl` | 本体（module `SVDLift`）。`Lift`（標準）, `LiftIll`（悪条件：自由方向を span{v_{K-1}, v_K} に広げ 1 段遅れて決める）, `LiftDense`（密な場合：乗数 c(t) を SVD で復元）, `GcdFromCof`, `DirectNull`（比較用の直接法）, `ChooseShift`（展開点の選択）など |
| `appgcd.mpl` | 神戸 2026-09 の道具箱のコピー。Bezout lifting [SNC09] との比較に使う |
| `common.mpl` | 例 1, 2 の定義、表の入出力 |
| `exp_examples.mpl` | 表 1：Example 1, 2 |
| `exp_illcond.mpl` | 表 2, 3：初期因子の近似共通因子の距離 delta を振る／補正パスの収束 |
| `exp_dense.mpl` | 表 4：密な場合（乗数の復元） |
| `exp_random.mpl` | 表 5：ランダム例、直接 SVD との比較、展開点の規則 |
| `verify_claims.mpl` | 論文 2, 3 節の記述を検証（27 項目を PASS/FAIL で出力、結果は `verify_claims.log`） |
| `verify_section4.py` | 論文 4 節の本文の数値を `tables/*.tex` と突き合わせる（29 項目、結果は `verify_section4.log`） |
| `exp_vectors.mpl` | 本文の例 1, 2 に載せる途中経過（v_K, δz, δq, δr, 残差）を出力（表は作らない） |
| `prototype_python/` | アルゴリズム検証用の NumPy 試作（`python3 runexp.py <出力先>` で同じ形式の暫定表） |

## 状況（2026-11-23）

- Maple 2025 で `run_all.mpl` 相当（exp_examples / exp_illcond / exp_dense / exp_random）を実行済み。
  論文の表 1〜5 はすべて Maple の出力（表の下に "Maple 2025, Digits=15"）。所要時間は全体で数分。
- すべての実験で入力に相対ノイズ eps を加えている（`SVDLift:-AddNoise`、同じ台の上の乱数摂動）。
- 密な場合の乗数次数 D の判定は許容誤差 tol = max(1e-8, 1e4*eps)。ノイズ下ではここが最も不安定（16/20 正解）。
- `prototype_python/` は Maple が動かなかった時期の NumPy 試作。Maple 版と同じアルゴリズムで、数値もほぼ一致した。
- Maple の注意：スクリプトのトップレベルで `D` に代入しないこと（微分演算子で保護されている）。
