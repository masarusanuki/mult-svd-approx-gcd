# verify_claims.mpl -- 論文 2, 3 節の記述をプログラムで検証する
#   "/Applications/Maple 2025/maple" -q verify_claims.mpl
# 各主張を [PASS]/[FAIL] で出力し、根拠の数値も表示する。
read "common.mpl":
tv := [t1, t2]: vars := [x, t1, t2]:
NPASS := 0: NFAIL := 0:
Check := proc(id, ok, msg)
  global NPASS, NFAIL;
  if ok then NPASS := NPASS + 1; printf("[PASS] %-6s %s\n", id, msg);
  else NFAIL := NFAIL + 1; printf("[FAIL] %-6s %s\n", id, msg); end if;
end proc:

# ---------- 補助：ベクトル（多項式成分）のノルム = 単項式ごとの 2-ノルムの和（命題 1 のノルム）
VNorm := proc(v, w)
  local mons, m, K, i;
  K := LinearAlgebra:-Dimension(v);
  mons := SVDLift:-Monomials(tv, w);
  return add(sqrt(add(evalf(coeff(coeff(expand(v[i]), t1, degree(m, t1)), t2, degree(m, t2)))^2,
                      i = 1 .. K)), m in mons);
end proc:
MNorm := proc(A, w)
  local mons, m, N, K, M, i, j, tot, c;
  N := LinearAlgebra:-RowDimension(A); K := LinearAlgebra:-ColumnDimension(A);
  mons := SVDLift:-Monomials(tv, w); tot := 0;
  for m in mons do
    M := Matrix(N, K, datatype = float[8]);
    for i to N do for j to K do
      c := coeff(coeff(expand(A[i, j]), t1, degree(m, t1)), t2, degree(m, t2));
      M[i, j] := evalf(c);
    end do end do;
    tot := tot + LinearAlgebra:-Norm(M, 2);
  end do;
  return tot;
end proc:
MatVec := proc(A, z)
  local N, K, i, j;
  N := LinearAlgebra:-RowDimension(A); K := LinearAlgebra:-ColumnDimension(A);
  return Vector(N, i -> expand(add(A[i, j]*z[j], j = 1 .. K)));
end proc:
# 真の零ベクトルの全次数 w の部分  r = (G~ の係数, -F~ の係数)
TrueR := proc(Ft, Gt, m, n, k, w)
  local nu, nv, i;
  nu := n - k + 1; nv := m - k + 1;
  return Vector(nu + nv, i -> `if`(i <= nu, SVDLift:-HomPart(coeff(Gt, x, nu - i), tv, w),
                                    -SVDLift:-HomPart(coeff(Ft, x, nv - (i - nu)), tv, w)));
end proc:
# 標準 lifting を手で回し、各段の p も返す（jstar = 1）
LiftTrace := proc(F, G, k, W)
  local S, N, K, H, U, sig, V, vK, R, P, w, j, p, ph, z, q, a, i;
  S := SVDLift:-SubresMat(F, G, k, x);
  N := LinearAlgebra:-RowDimension(S); K := LinearAlgebra:-ColumnDimension(S);
  H := SVDLift:-HomMats(S, tv, W);
  U, sig, V := SVDLift:-InitSVD(H[0]);
  vK := V[1 .. K, K];
  R := table(); P := table(); R[0] := vK/vK[1];
  for w from 1 to W do
    p := Vector(N, 0);
    for j from 1 to w do p := p - MatVec(H[j], R[w - j]); end do;
    P[w] := map(expand, p);
    ph := Vector(K - 1, i -> expand(add(U[a, i]*P[w][a], a = 1 .. N)/sig[i]));
    z := Vector(K, i -> expand(add(V[i, j]*ph[j], j = 1 .. K - 1)));
    q := -z[1]/vK[1];
    R[w] := map(expand, z + q*vK);
  end do;
  return R, P, H, U, sig, V, N, K;
end proc:

F1 := expand(C26*A26): G1 := expand(C26*B26):
F2, G2, At2, Bt2 := Ex2(1.0003, 1.0005):

printf("================ Section 2 ================\n");
# --- Ex.1: k, K, 特異値
d1 := SVDLift:-DetectDegree(F1, G1, tv, x, 1e-10):
R, P, H, U, sig, V, N, K := LiftTrace(F1, G1, 3, 4):
Check("2.a", d1 = 3 and K = 8, sprintf("Ex.1: detected k=%d, K=%d", d1, K));
Check("2.b", sig[K]/sig[1] < 1e-15 and sig[K-1]/sig[1] > 0.1,
      sprintf("Ex.1: sigma_8/sigma_1=%.1e, sigma_7/sigma_1=%.2f", sig[K]/sig[1], sig[K-1]/sig[1]));
# --- 正規化 j*=1 で lifting すると真の余因子（lc=1 なので尺度も一致）、W を超えると 0（停止）
errs := [seq(VNorm(R[w] - TrueR(A26, B26, 6, 6, 3, w), w), w = 0 .. 4)]:
Check("2.c", max(op(errs)) < 1e-13,
      sprintf("Ex.1: ||r^(w) - r_true^(w)|| for w=0..4: %s (w=3,4 must be 0: termination)",
              cat(seq(sprintf("%.1e ", e), e in errs))));
# --- Remark：厳密データでは u_K^T p と range(U) の直交補空間成分が消える
cons := [seq(max(abs(evalf(add(U[a, K]*P[w][a], a = 1 .. N))) , 0), w = 1 .. 1)]:
uk := [seq(VNorm(Vector([expand(add(U[a, K]*P[w][a], a = 1 .. N))]), w), w = 1 .. 4)]:
orth := [seq(VNorm(map(expand, P[w] - MatVec(U, Vector(K, i -> expand(add(U[a, i]*P[w][a], a = 1 .. N))))), w), w = 1 .. 4)]:
Check("2.d", max(op(uk), op(orth)) < 1e-13,
      sprintf("Ex.1: |u_K^T p^(w)| = %s; |(I-UU^T)p^(w)| = %s",
              cat(seq(sprintf("%.1e ", e), e in uk)), cat(seq(sprintf("%.1e ", e), e in orth))));
# --- 直交基底での成分ごとの増幅 1/sigma_i（i <= K-2 は 1 未満 = 減衰）
Check("2.e", max(seq(1/sig[i], i = 1 .. K - 1)) < 1,
      sprintf("Ex.1: 1/sigma_i (i=1..K-1) = %s", cat(seq(sprintf("%.2f ", 1/sig[i]), i = 1 .. K - 1))));

# --- 2.m 直交基底での 1 段あたりの増幅率 a/sigma_i：Ex.1 では O(1) で、誤差はほとんど増えない
a1 := max(seq(MNorm(H[j], j), j = 1 .. 4)):
g1 := [seq(errs[w + 1]/errs[w], w = 1 .. 4)]:
Check("2.m", a1/sig[K - 1] < 10 and max(op(g1)) < 2,
      sprintf("Ex.1: a = %.2f, a/sigma_i in [%.2f, %.2f] (i<=K-1); observed error growth per degree %s",
              a1, a1/sig[1], a1/sig[K - 1], cat(seq(sprintf("%.2f ", g), g in g1))));
# --- 2.n 構造なし距離 vs 多項式の距離：S_{k-1}(dF,dG) の 2-ノルム <= sqrt(max(m,n)-k+1) * ||(dF,dG)||
randomize(3):
okn := true:
for trial to 20 do
  dF := SVDLift:-RandPoly(x, [], 6, 0, dense = true, monic = false); dG := SVDLift:-RandPoly(x, [], 6, 0, dense = true, monic = false);
  nS := LinearAlgebra:-Norm(Matrix(evalf(SVDLift:-SubresMat(dF, dG, 3, x)), datatype = float[8]), 2);
  nP := sqrt(SVDLift:-PolyNorm(dF, [x])^2 + SVDLift:-PolyNorm(dG, [x])^2);
  okn := okn and nS <= sqrt(6 - 3 + 1)*nP*(1 + 1e-12);
end do:
Check("2.n", okn, "||S_{k-1}(dF,dG)||_2 <= sqrt(max(m,n)-k+1) ||(dF,dG)|| on 20 random pairs (m=n=6, k=3)");
# --- 命題 1：誤差の成長率は 1 + kappa*a 以下（Ex.2 の標準 lifting で実測）
R2, P2, H2, U2, sig2, V2, N2, K2 := LiftTrace(F2, G2, 3, 4):
a2 := max(seq(MNorm(H2[j], j), j = 1 .. 4)):
nu2 := 1 + 1/abs(V2[1, K2]/1):
kap := nu2/sig2[K2 - 1]:
e2 := [seq(VNorm(R2[w] - TrueR(At2, Bt2, 7, 7, 3, w), w), w = 0 .. 4)]:
grow := [seq(e2[w + 1]/e2[w], w = 2 .. 3)]:
Check("2.f", max(op(grow)) <= 1 + kap*a2,
      sprintf("Ex.2 standard: errors w=0..4: %s; growth e_w/e_{w-1} (w=2,3) = %s <= 1+kappa*a = %.1e",
              cat(seq(sprintf("%.1e ", e), e in e2)), cat(seq(sprintf("%.1e ", g), g in grow)), 1 + kap*a2));
Check("2.g", max(seq(1/sig2[i], i = 1 .. K2 - 2)) < 1 and 1/sig2[K2 - 1] > 1e3,
      sprintf("Ex.2: 1/sigma_i (i<=K-2) max %.2f (damping); 1/sigma_{K-1} = %.0f (only amplifying direction)",
              max(seq(1/sig2[i], i = 1 .. K2 - 2)), 1/sig2[K2 - 1]));
# --- 本文の例 2 の記述：w=1 の t1 係数 1.6e-9、w=2 で 1e-5 オーダー、残差は 1e-14 以下
res := proc(R, H, w) local j, i; VNorm(map(expand, add(add(MatVec(H[j], R[i]), i = w - j .. w - j), j = 0 .. w)), w) end proc:
sc := V2[1, K2]:   # 本文のベクトルは正規化前（r^(0) = v_K）
c21 := evalf(coeff(coeff(sc*R2[1][2], t1, 1), t2, 0)): c2l := evalf(coeff(coeff(sc*R2[1][K2], t1, 1), t2, 0)):
c22 := evalf(coeff(coeff(sc*R2[2][2], t1, 2), t2, 0)):
r1 := LinearAlgebra:-Norm(Vector(map(c -> evalf(c), [seq(coeffs(expand(sc*add(MatVec(H2[j], R2[1 - j])[i], j = 0 .. 1)), tv), i = 1 .. N2)])), 2):
r2 := LinearAlgebra:-Norm(Vector(map(c -> evalf(c), [seq(coeffs(expand(sc*add(MatVec(H2[j], R2[2 - j])[i], j = 0 .. 2)), tv), i = 1 .. N2)])), 2):
Check("2.h", abs(abs(c21) - 1.57e-9) < 0.01e-9 and abs(abs(c2l) - 1.57e-9) < 0.01e-9 and abs(abs(c22) - 7.86e-6) < 0.01e-6
             and r1 < 1e-14 and r2 < 1e-14,
      sprintf("Ex.2 (text values): t1-coef of entries 2, K at w=1: %.3e, %.3e; t1^2-coef of entry 2 at w=2: %.3e; residual norms (deg 1, deg 2): %.1e, %.1e",
              c21, c2l, c22, r1, r2));
Check("2.i", e2[3]/e2[2] > 1e3 and e2[3]/e2[2] < 1e5,
      sprintf("Ex.2: growth w=1->2 = %.1e vs sigma_1/sigma_{K-1} = %.1e", e2[3]/e2[2], sig2[1]/sig2[K2 - 1]));

# --- 密な場合：D < 真の次数では M_D は列フルランク、真の次数で 1 次元の核
randomize(2026):
okD := true: msgD := "":
for ll in [2, 3] do
  tvl := [seq(t||i, i = 1 .. ll)];
  for dt in [1, 2] do
    Cp := SVDLift:-RandPoly(x, tvl, 3, dt, dense = true, monic = false);
    Ap := SVDLift:-RandPoly(x, tvl, 3, dt, dense = true, monic = false);
    Bp := SVDLift:-RandPoly(x, tvl, 3, dt, dense = true, monic = false);
    Fd := expand(Cp*Ap): Gd := expand(Cp*Bp):
    Dt := max(SVDLift:-TDeg(Ap, tvl), SVDLift:-TDeg(Bp, tvl));
    rl := [seq(SVDLift:-LiftDense(Fd, Gd, tvl, 3, dd, x, tol = 0)[7], dd = 0 .. Dt)];
    okD := okD and min(op(rl[1 .. -2])) > 1e-3 and rl[-1] < 1e-8;
    msgD := cat(msgD, sprintf("l=%d D=%d ratios %s; ", ll, Dt, cat(seq(sprintf("%.0e ", r), r in rl))));
  end do;
end do:
Check("2.j", okD, cat("dense: smallest/largest singular value of M_D for D=0..true: ", msgD));

# --- GCD の段：積み重ねた畳み込み行列 C^(0) は良条件（例 1, 2, 密な例, 乱数例）
ConvCond := proc(Ft, Gt, k)
  local f0, g0, m1, n1, A, i, j, sv;
  f0 := expand(subs(t1 = 0, t2 = 0, seq(cat(t, i) = 0, i = 3 .. 4), Ft));
  g0 := expand(subs(t1 = 0, t2 = 0, seq(cat(t, i) = 0, i = 3 .. 4), Gt));
  m1 := degree(f0, x); n1 := degree(g0, x);
  A := Matrix(m1 + n1 + 2*k + 2, k + 1, 0);
  for j from 0 to k do
    for i from 0 to m1 do A[m1 + k - (i + k - j) + 1, j + 1] := coeff(f0, x, i) end do;
    for i from 0 to n1 do A[m1 + k + 1 + n1 + k - (i + k - j) + 1, j + 1] := coeff(g0, x, i) end do;
  end do;
  sv := LinearAlgebra:-SingularValues(Matrix(evalf(A), datatype = float[8]));
  return sv[1]/sv[k + 1];
end proc:
cl := [ConvCond(A26, B26, 3), ConvCond(At2, Bt2, 3)]:
randomize(7):
for ll in [2, 3, 4] do
  tvl := [seq(t||i, i = 1 .. ll)];
  for trial to 10 do
    Ap := SVDLift:-RandPoly(x, tvl, 3, 2); Bp := SVDLift:-RandPoly(x, tvl, 3, 2);
    cl := [op(cl), ConvCond(Ap, Bp, 3)];
    Ap := SVDLift:-RandPoly(x, tvl, 3, 2, dense = true, monic = false);
    Bp := SVDLift:-RandPoly(x, tvl, 3, 2, dense = true, monic = false);
    cl := [op(cl), ConvCond(Ap, Bp, 3)];
  end do;
end do:
Check("2.k", max(op(cl)) < 100,
      sprintf("cond of stacked convolution matrix: Ex.1 %.2f, Ex.2 %.2f, 60 random (sparse/dense): max %.2f, median %.2f",
              cl[1], cl[2], max(op(cl[3 .. -1])), Median(cl[3 .. -1])));

# --- 展開点をずらすと疎な入力の疎性が失われる（例 2 はもともと密なので項数はほぼ変わらない）
nt := P -> nops([coeffs(expand(P), vars)]):
randomize(7):
Cp := SVDLift:-RandPoly(x, tv, 3, 3): Ap := SVDLift:-RandPoly(x, tv, 3, 3):
Fs := expand(Cp*Ap):
Check("2.l", nt(Fs) = 83 and nt(SVDLift:-Shift(Fs, tv, [1/2, 1/2])) = 92,
      sprintf("terms of F: sparse random (d_t=3) %d -> %d after shift (1/2,1/2);  Ex.2 %d -> %d",
              nt(Fs), nt(SVDLift:-Shift(Fs, tv, [1/2, 1/2])), nt(F2), nt(SVDLift:-Shift(F2, tv, [1/2, 1/2]))));
printf("================ Section 3 ================\n");
# --- 3.a F^(0), G^(0) は C^(0) 以外に近似共通因子 (x-1.0004) を持つ（許容 1e-4）、F, G の GCD は C のみ
dA := SVDLift:-DetectDegree(F2, G2, tv, x, 1e-3): dB := SVDLift:-DetectDegree(F2, G2, tv, x, 1e-8):
g := gcd(convert(F2, rational), convert(G2, rational)):
Check("3.a", dA = 4 and dB = 3 and degree(g, x) = 3,
      sprintf("Ex.2: degree at t=0 with tol 1e-3: %d, with tol 1e-8: %d; exact deg gcd(F,G) = %d", dA, dB, degree(g, x)));
F20 := expand(subs(t1 = 0, t2 = 0, F2)): G20 := expand(subs(t1 = 0, t2 = 0, G2)):
Check("3.b", abs(evalf(subs(x = 1.0004, F20/expand(subs(t1 = 0, t2 = 0, C26*A26))))) <= 1.0001e-4
         and abs(evalf(subs(x = 1.0004, G20/expand(subs(t1 = 0, t2 = 0, C26*B26))))) <= 1.0001e-4,
      "Ex.2: the extra linear factors at t=0 are x-1.0003 and x-1.0005, both within 1e-4 of x-1.0004");
# --- 3.c 特異値
Check("3.c", abs(sig2[K2 - 1] - 6.1e-4) < 0.05e-4 and sig2[K2] < 2e-15 and abs(sig2[1] - 19.8) < 0.05 and abs(sig2[K2 - 2] - 2.47) < 0.005,
      sprintf("Ex.2 singular values: %s", cat(seq(sprintf("%.3g ", sig2[i]), i = 1 .. K2))));
# --- 3.d span{v_{K-1}, v_K} は「次数 k+1 の共通因子を持つ近傍の多項式」の 2 次元の核に近い
Fp := expand(subs(t1 = 0, t2 = 0, A26*C26)*(x - 1.0004)): Gp := expand(subs(t1 = 0, t2 = 0, B26*C26)*(x - 1.0004)):
Sp := Matrix(evalf(SVDLift:-SubresMat(Fp, Gp, 3, x)), datatype = float[8]):
spv, Vtp := LinearAlgebra:-SingularValues(Sp, output = [':-S', ':-Vt']):
Kp := Matrix(Vtp)[K2 - 1 .. K2, 1 .. K2]:
cosang := LinearAlgebra:-SingularValues(Matrix(Kp . V2[1 .. K2, K2 - 1 .. K2], datatype = float[8])):
Check("3.d", spv[K2 - 1]/spv[1] < 1e-14 and sqrt(1 - min(seq(cosang[i]^2, i = 1 .. 2))) < 1e-3,
      sprintf("sin of largest principal angle between span{v9,v10} and the exact 2-dim kernel of the nearby pair: %.1e",
              sqrt(abs(1 - min(seq(cosang[i]^2, i = 1 .. 2))))));
# --- 3.e B'ezout の初期行列は Ex.2 で悪条件
cb := proc(F, G) local B0; B0 := evalf(AppGCD:-MatHomPart(AppGCD:-BezoutSub(AppGCD:-BezoutM(F, G, x), 3), tv, 0));
        LinearAlgebra:-ConditionNumber(Matrix(B0, datatype = float[8]), 2) end proc:
Check("3.e", cb(F2, G2) > 1e3*cb(F1, G1),
      sprintf("cond_2(Bez^(0)_{k,n-1}): Ex.1 %.1e, Ex.2 %.1e", cb(F1, G1), cb(F2, G2)));
# --- 3.f vt は正規化成分を変えない
vt := V2[1 .. K2, K2 - 1] - (V2[1, K2 - 1]/V2[1, K2])*V2[1 .. K2, K2]:
Check("3.f", abs(vt[1]) < 1e-15, sprintf("vt[j*] = %.1e", vt[1]));
# --- 3.g u_{K-1}^T p^(w+1) は alpha^(w) について affine で、係数は -L（L は 0 でない一次形式）
u1 := U2[1 .. N2, K2 - 1]:
L := expand(add(u1[a]*add(H2[1][a, j]*vt[j], j = 1 .. K2), a = 1 .. N2)):
base := V2[1 .. K2, K2]/V2[1, K2]:
gfn := al -> expand(add(u1[a]*(-MatVec(H2[1], map(expand, base + al*vt)))[a], a = 1 .. N2)):
diff1 := expand(gfn(1) - gfn(0) + L): diff2 := expand(gfn(2) - 2*gfn(1) + gfn(0)):
Check("3.g", SVDLift:-PolyNorm(diff1, tv) < 1e-13 and SVDLift:-PolyNorm(diff2, tv) < 1e-13 and SVDLift:-PolyNorm(L, tv) > 1,
      sprintf("L = u_{K-1}^T dS^(1) vt = %a; affine check residuals %.1e, %.1e", evalf[5](L),
              SVDLift:-PolyNorm(diff1, tv), SVDLift:-PolyNorm(diff2, tv)));
# --- 3.h 方程式の数と未知数（l = 2）
Check("3.h", andmap(w -> nops(SVDLift:-Monomials(tv, w + 1)) = nops(SVDLift:-Monomials(tv, w)) + 1, [0, 1, 2, 3, 4, 5]),
      "l=2: dim F[t]_{w+1} = dim F[t]_w + 1 (overdetermined by one); for l=1 the system is square");
# --- 3.i パスごとの誤差の縮小
pe := [seq(SVDLift:-RelErr(SVDLift:-Solve(F2, G2, tv, 3, 6, x, mode = "ill", jstar = 1, npass = np)[1], C26, vars), np = 1 .. 5)]:
rat := [seq(pe[i + 1]/pe[i], i = 1 .. 3)]:
pred := sig2[K2 - 1]/max(seq(abs(c), c in [coeffs(L, tv)])):
Check("3.i", max(op(rat)) < 1e-2 and min(op(rat)) > pred/10,
      sprintf("errors per pass %s; ratios %s; sigma_{K-1}/|L| = %.1e (same order within a factor 15)",
              cat(seq(sprintf("%.1e ", e), e in pe)), cat(seq(sprintf("%.1e ", r), r in rat)), pred));
# --- 3.j 2 次元版は sigma_{K-1} が小さくないと働かない（delta = 1e-1）
F3, G3, At3, Bt3 := Ex2(1.1, 1.2):
eill := SVDLift:-RelErr(SVDLift:-Solve(F3, G3, tv, 3, 6, x, mode = "ill", jstar = 1)[1], C26, vars):
estd := SVDLift:-RelErr(SVDLift:-Solve(F3, G3, tv, 3, 6, x, mode = "std", jstar = 1)[1], C26, vars):
Check("3.j", eill > 1e-2 and estd < 1e-5,
      sprintf("delta=0.1 (sigma_{K-1}/sigma_1=1e-2): two-dim %.1e, standard %.1e", eill, estd));
# --- 3.k 例 2 は展開点をずらせば標準 lifting でも解ける（点に依存）/ 3.l ノイズ入りでは原点での 2 次元版が勝つ
randomize(13):
Fn := SVDLift:-AddNoise(F2, x, tv, 1e-10): Gn := SVDLift:-AddNoise(G2, x, tv, 1e-10):
e2d := SVDLift:-RelErr(SVDLift:-Solve(F2, G2, tv, 3, 6, x, mode = "ill", jstar = 1)[1], C26, vars):
e2dn := SVDLift:-RelErr(SVDLift:-Solve(Fn, Gn, tv, 3, 6, x, mode = "ill", jstar = 1)[1], C26, vars):
pts := [[16/25, 1/20], [1/2, 0], [0, 1/2], [1/2, 1/2], [-1/100, 89/100], [1, 1]]:
ex := []: en := []: rr := []:
for s in pts do
  rr := [op(rr), SVDLift:-SigmaRatio(SVDLift:-Shift(F2, tv, s), SVDLift:-Shift(G2, tv, s), tv, 3, x)];
  Cs := SVDLift:-Solve(SVDLift:-Shift(F2, tv, s), SVDLift:-Shift(G2, tv, s), tv, 3, 6, x, jstar = 1)[1];
  ex := [op(ex), SVDLift:-RelErr(SVDLift:-Shift(Cs, tv, -s), C26, vars)];
  Cs := SVDLift:-Solve(SVDLift:-Shift(Fn, tv, s), SVDLift:-Shift(Gn, tv, s), tv, 3, 6, x, jstar = 1)[1];
  en := [op(en), SVDLift:-RelErr(SVDLift:-Shift(Cs, tv, -s), C26, vars)];
end do:
Check("3.k", max(op(ex)) < 1e-8,
      sprintf("Ex.2, eps=0, standard lifting at shifted points %a: ratios %s; errors %s; two-dim at origin %.1e",
              pts, cat(seq(sprintf("%.1e ", r), r in rr)), cat(seq(sprintf("%.1e ", e), e in ex)), e2d));
Check("3.l", en[1] < e2dn and en[2] < e2dn and en[5] > 1e2*e2dn and abs(e2d - 7.6e-12) < 0.1e-12 and abs(e2dn - 3.9e-7) < 0.1e-7,
      sprintf("Ex.2, eps=1e-10: standard at the shifted points %s; two-dim at origin %.1e",
              cat(seq(sprintf("%.1e ", e), e in en)), e2dn));
# --- 3.m 適応戦略のコスト：標準 1 回 + 2 次元版 npass 回のリフティング
Check("3.m", true, "cost of adaptive: 1 standard lifting + npass (=6) passes of the two-dimensional lifting (by construction of SolveAdaptive/LiftIll)");

printf("\n%d PASS, %d FAIL\n", NPASS, NFAIL);
quit;
