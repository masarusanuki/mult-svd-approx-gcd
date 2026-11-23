# exp_random.mpl -- ランダム例：直接法との比較、展開点の規則（表 tab_random.tex）
read "common.mpl":
randomize(7):
NTRIAL := 10: TAU := 1e-2: k := 3:
NAME := "tab_random.tex": fd := OpenTab(NAME, cat(
  "\\begin{tabular}{ccc|ccc|ccc|cc}\\hline\n",
  "&&& \\multicolumn{3}{c|}{direct SVD} & \\multicolumn{3}{c|}{proposed, $\\vect{s}=\\vect{0}$} & \\multicolumn{2}{c}{proposed, shift rule}\\\\\n",
  "$\\ell$ & $d_{\\vect{t}}$ & $\\varepsilon$ & size & time (s) & error (max) & time (s) & error (med.) & error (max) & error (med.) & error (max)\\\\ \\hline\n")):
for cfg in [[2,1],[2,2],[2,3],[3,1],[3,2],[4,1],[4,2]] do
  ll, dt := op(cfg);
  tv := [seq(t||i, i = 1 .. ll)]: vars := [x, op(tv)]:
  for eps in [0, 1e-8, 1e-6] do
    ed := []: td := []: e0 := []: t0 := []: es := []:
    for trial to NTRIAL do
      Cp := SVDLift:-RandPoly(x, tv, k, dt);
      Ap := SVDLift:-RandPoly(x, tv, 3, dt);
      Bp := SVDLift:-RandPoly(x, tv, 3, dt);
      F := SVDLift:-AddNoise(expand(Cp*Ap), x, tv, eps);
      G := SVDLift:-AddNoise(expand(Cp*Bp), x, tv, eps);
      Wb := max(SVDLift:-TDeg(F, tv), SVDLift:-TDeg(G, tv));
      # 直接法（余因子の真の t-次数を与える）
      st := time();
      Ftd, Gtd, shp, tsvd := SVDLift:-DirectNull(F, G, tv, k,
                                max(SVDLift:-TDeg(Ap, tv), SVDLift:-TDeg(Bp, tv)), x);
      td := [op(td), time() - st];
      ed := [op(ed), SVDLift:-RelErr(Ftd, Ap, vars)];
      # 提案法（s = 0）
      st := time();
      Cc, Ft, Gt, sig := SVDLift:-Solve(F, G, tv, k, Wb, x, jstar = 1);
      t0 := [op(t0), time() - st];
      e0 := [op(e0), SVDLift:-RelErr(Cc, Cp, vars)];
      # 展開点の規則：sigma_{K-1}/sigma_1 < TAU のときだけ動かす
      if SVDLift:-SigmaRatio(F, G, tv, k, x) < TAU then
        s, rs := SVDLift:-ChooseShift(F, G, tv, k, x);
        Cs, Fs, Gs, sig := SVDLift:-Solve(SVDLift:-Shift(F, tv, s), SVDLift:-Shift(G, tv, s),
                                          tv, k, Wb, x, jstar = 1);
        es := [op(es), SVDLift:-RelErr(SVDLift:-Shift(Cs, tv, -s), Cp, vars)];
      else
        es := [op(es), e0[-1]];
      end if;
    end do;
    fprintf(fd, "%d & %d & %s & $%d\\times %d$ & %.3f & %s & %.3f & %s & %s & %s & %s\\\\\n",
      ll, dt, SVDLift:-Fmt(eps), shp[1], shp[2], Median(td), SVDLift:-Fmt(max(op(ed))),
      Median(t0), SVDLift:-Fmt(Median(e0)), SVDLift:-Fmt(max(op(e0))),
      SVDLift:-Fmt(Median(es)), SVDLift:-Fmt(max(op(es))));
    printf("l=%d dt=%d eps=%g done\n", ll, dt, eps);
  end do;
end do:
CloseTab(fd, NAME):
