# exp_examples.mpl -- Example 1, 2（入力に相対ノイズ eps を加える）（表 tab_examples.tex）
read "common.mpl":
randomize(11):
tv := [t1, t2]: vars := [x, t1, t2]:
NAME := "tab_examples.tex": fd := OpenTab(NAME,
  "\\begin{tabular}{lllccc}\\hline\nexample & $\\varepsilon$ & method & cofactor error & GCD error & backward error\\\\ \\hline\n"):
for eps in [0, 1e-10, 1e-8] do
  F := SVDLift:-AddNoise(expand(C26*A26), x, tv, eps):
  G := SVDLift:-AddNoise(expand(C26*B26), x, tv, eps):
  W := max(SVDLift:-TDeg(F, tv), SVDLift:-TDeg(G, tv)):
  Cc, Ft, Gt, sig := SVDLift:-Solve(F, G, tv, 3, W, x, jstar = 1):
  if eps = 0 then printf("Ex.1 singular values: %a\n", convert(sig, list)) end if;
  fprintf(fd, "Ex.~1 & %s & standard lifting & %s & %s & %s\\\\\n", SVDLift:-Fmt(eps),
    SVDLift:-Fmt(SVDLift:-RelErr(Ft, A26, vars)), SVDLift:-Fmt(SVDLift:-RelErr(Cc, C26, vars)),
    SVDLift:-Fmt(SVDLift:-BackErr(F, G, Cc, Ft, Gt, x, tv)));
  Cb := AppGCD:-BezoutLift(F, G, tv, 3, W, x)[1]:
  fprintf(fd, "Ex.~1 & %s & B\\'ezout lifting~\\cite{sanuki09} & -- & %s & --\\\\\n", SVDLift:-Fmt(eps),
    SVDLift:-Fmt(SVDLift:-RelErr(Cb, C26, vars)));
end do:
fprintf(fd, "\\hline\n"):
for eps in [0, 1e-10, 1e-8] do
  F, G, At, Bt := Ex2(1.0003, 1.0005):
  F := SVDLift:-AddNoise(F, x, tv, eps): G := SVDLift:-AddNoise(G, x, tv, eps):
  W := max(SVDLift:-TDeg(F, tv), SVDLift:-TDeg(G, tv)):
  for md in ["std", "ill"] do
    Cc, Ft, Gt, sig := SVDLift:-Solve(F, G, tv, 3, W, x, mode = md, jstar = 1);
    fprintf(fd, "Ex.~2 & %s & %s & %s & %s & %s\\\\\n", SVDLift:-Fmt(eps),
      `if`(md = "std", "standard lifting", "two-dim.\\ free space"),
      SVDLift:-Fmt(SVDLift:-RelErr(Ft, At, vars)), SVDLift:-Fmt(SVDLift:-RelErr(Cc, C26, vars)),
      SVDLift:-Fmt(SVDLift:-BackErr(F, G, Cc, Ft, Gt, x, tv)));
  end do;
  if eps = 0 then printf("Ex.2 singular values: %a\n", convert(sig, list)) end if;
  Cb := AppGCD:-BezoutLift(F, G, tv, 3, W, x)[1]:
  fprintf(fd, "Ex.~2 & %s & B\\'ezout lifting~\\cite{sanuki09} & -- & %s & --\\\\\n", SVDLift:-Fmt(eps),
    SVDLift:-Fmt(SVDLift:-RelErr(Cb, C26, vars)));
end do:
CloseTab(fd, NAME):
