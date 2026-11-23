# exp_dense.mpl -- 密な場合：乗数の復元（表 tab_dense.tex）
read "common.mpl":
randomize(2026):
NTRIAL := 5:
NAME := "tab_dense.tex": fd := OpenTab(NAME,
  "\\begin{tabular}{ccccccc}\\hline\n$\\ell$ & $d_{\\vect{t}}$ & $\\varepsilon$ & $D$ correct & fixed normalization & multiplier reconstruction & backward error\\\\ \\hline\n"):
for ll in [2, 3] do
  tv := [seq(t||i, i = 1 .. ll)]: vars := [x, op(tv)]:
  for dt in [1, 2] do for eps in [0, 1e-8] do
    en := []: er := []: bw := []: ok := 0:
    for trial to NTRIAL do
      Cp := SVDLift:-RandPoly(x, tv, 3, dt, dense = true, monic = false);
      Ap := SVDLift:-RandPoly(x, tv, 3, dt, dense = true, monic = false);
      Bp := SVDLift:-RandPoly(x, tv, 3, dt, dense = true, monic = false);
      F := SVDLift:-AddNoise(expand(Cp*Ap), x, tv, eps): G := SVDLift:-AddNoise(expand(Cp*Bp), x, tv, eps):
      Wb := max(SVDLift:-TDeg(F, tv), SVDLift:-TDeg(G, tv));
      # 定数成分を仮定した正規化（誤り）
      R, sig, js, nu, nv := SVDLift:-Lift(F, G, tv, 3, Wb, x);
      Ft2, Gt2 := SVDLift:-VecToCof(R, Wb, nu, nv, x);
      en := [op(en), SVDLift:-RelErr(Ft2, Ap, vars)];
      # 乗数の復元
      R, sig, js, nu, nv, Dd, ratio := SVDLift:-LiftDense(F, G, tv, 3, Wb, x, tol = max(1e-8, 1e4*eps));
      Ft, Gt := SVDLift:-VecToCof(R, Dd, nu, nv, x);
      Cc := SVDLift:-GcdFromCof(F, G, Ft, Gt, tv, 3, Dd, x);
      er := [op(er), SVDLift:-RelErr(Ft, Ap, vars)];
      bw := [op(bw), SVDLift:-BackErr(F, G, Cc, Ft, Gt, x, tv)];
      if Dd = max(SVDLift:-TDeg(Ap, tv), SVDLift:-TDeg(Bp, tv)) then ok := ok + 1 end if;
    end do;
    fprintf(fd, "%d & %d & %s & %d/%d & %s & %s & %s\\\\\n", ll, dt, SVDLift:-Fmt(eps), ok, NTRIAL,
      SVDLift:-Fmt(Median(en)), SVDLift:-Fmt(Median(er)), SVDLift:-Fmt(Median(bw)));
  end do; end do;
end do:
CloseTab(fd, NAME):
