# exp_illcond.mpl -- 初期因子が近似共通因子を持つ場合（入力に相対ノイズ）
#   表 tab_illcond.tex, tab_passes.tex
read "common.mpl":
randomize(13):
tv := [t1, t2]: vars := [x, t1, t2]:
NAME := "tab_illcond.tex": fd := OpenTab(NAME,
  "\\begin{tabular}{cccccccc}\\hline\n$\\varepsilon$ & $\\delta$ & $\\sigma_{K-1}/\\sigma_1$ & $\\sigma_K/\\sigma_1$ & standard & two-dim.\\ free space & adaptive & B\\'ezout~\\cite{sanuki09}\\\\ \\hline\n"):
for eps in [0, 1e-10] do
  for dl in [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6] do
    F, G, At, Bt := Ex2(1 + dl, 1 + 2*dl);
    F := SVDLift:-AddNoise(F, x, tv, eps); G := SVDLift:-AddNoise(G, x, tv, eps);
    W := max(SVDLift:-TDeg(F, tv), SVDLift:-TDeg(G, tv));
    C1, F1, G1, sig := SVDLift:-Solve(F, G, tv, 3, W, x, mode = "std", jstar = 1);
    C2, F2, G2, sig := SVDLift:-Solve(F, G, tv, 3, W, x, mode = "ill", jstar = 1);
    Ca, Fa, Ga, sig, which := SVDLift:-SolveAdaptive(F, G, tv, 3, W, x, jstar = 1);
    Cb := AppGCD:-BezoutLift(F, G, tv, 3, W, x)[1];
    K := LinearAlgebra:-Dimension(sig);
    fprintf(fd, "%s & %s & %s & %s & %s & %s & %s (%s) & %s\\\\\n",
      SVDLift:-Fmt(eps), SVDLift:-Fmt(dl), SVDLift:-Fmt(sig[K-1]/sig[1]), SVDLift:-Fmt(sig[K]/sig[1]),
      SVDLift:-Fmt(SVDLift:-RelErr(C1, C26, vars)), SVDLift:-Fmt(SVDLift:-RelErr(C2, C26, vars)),
      SVDLift:-Fmt(SVDLift:-RelErr(Ca, C26, vars)), which,
      SVDLift:-Fmt(SVDLift:-RelErr(Cb, C26, vars)));
  end do;
  if eps = 0 then fprintf(fd, "\\hline\n") end if;
end do:
CloseTab(fd, NAME):

NAME := "tab_passes.tex": fd := OpenTab(NAME,
  "\\begin{tabular}{ccccc}\\hline\n& \\multicolumn{2}{c}{$\\varepsilon=0$} & \\multicolumn{2}{c}{$\\varepsilon=10^{-10}$}\\\\\npasses & GCD error & backward error & GCD error & backward error\\\\ \\hline\n"):
F0, G0, At, Bt := Ex2(1.0003, 1.0005):
F1 := SVDLift:-AddNoise(F0, x, tv, 1e-10): G1 := SVDLift:-AddNoise(G0, x, tv, 1e-10):
W := max(SVDLift:-TDeg(F0, tv), SVDLift:-TDeg(G0, tv)):
for np from 1 to 6 do
  Cc, Ft, Gt, sig := SVDLift:-Solve(F0, G0, tv, 3, W, x, mode = "ill", jstar = 1, npass = np);
  Cn, Fn, Gn, sig := SVDLift:-Solve(F1, G1, tv, 3, W, x, mode = "ill", jstar = 1, npass = np);
  fprintf(fd, "%d & %s & %s & %s & %s\\\\\n", np,
    SVDLift:-Fmt(SVDLift:-RelErr(Cc, C26, vars)), SVDLift:-Fmt(SVDLift:-BackErr(F0, G0, Cc, Ft, Gt, x, tv)),
    SVDLift:-Fmt(SVDLift:-RelErr(Cn, C26, vars)), SVDLift:-Fmt(SVDLift:-BackErr(F1, G1, Cn, Fn, Gn, x, tv)));
end do:
CloseTab(fd, NAME):
