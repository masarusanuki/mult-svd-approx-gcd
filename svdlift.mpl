# svdlift.mpl -- SVD + lifting による多変数近似 GCD（論文 paper-full-multsvd.tex の実装）
#
#   read "svdlift.mpl":
#
# 記法は論文に合わせる。
#   F, G   : x と t = [t1,...,tl] の多項式（浮動小数係数）
#   k      : 近似 GCD の x についての次数
#   S      : 部分終結式行列 S_{k-1}(F,G)  (u,v) -> u*F + v*G,
#            deg u <= n-k, deg v <= m-k。(m+n-k+1) x K,  K = m+n-2k+2。
#            零ベクトルは r = (G~_{n-k},...,G~_0, -F~_{m-k},...,-F~_0)。
#   S^(w)  : S の t についての全次数 w の斉次部分
#
# 手続き
#   Lift        : 標準の lifting（自由方向は v_K のみ。成分 jstar を固定して正規化）
#   LiftIll     : 初期因子が近似共通因子を持つ場合（sigma_{K-1} が小さい）。
#                 自由方向を span{v_{K-1}, v_K} に広げ、v_{K-1} の係数を 1 段遅れて決める。
#   LiftDense   : どの成分も t について定数でない場合。冪級数として持ち上げてから
#                 多項式の乗数 c(t) を（次数 D を上げながら）SVD で求める。
#   GcdFromCof  : 余因子から GCD を求める（積み重ねた畳み込み行列の lifting）
#   DirectNull  : 比較用。多変数部分終結式を数値行列として作り、SVD で零ベクトルを求める
#
# Maple 2025 / Digits = 15（ハードウェア浮動小数）を想定。

SVDLift := module()
  option package;
  export HomPart, TDeg, Truncate, Monomials, SubresMat, HomMats, InitSVD,
         Lift, LiftIll, LiftDense, VecToCof, GcdFromCof, Solve, SolveAdaptive,
         DirectNull, DetectDegree, SigmaRatio, ChooseShift, Shift,
         PolyNorm, RelErr, BackErr, RandPoly, AddNoise, Fmt;
  local MatVec, Rhs, CoeffEqs;

  # ------------------------------------------------------------------
  # 斉次部分・全次数・打ち切り
  # ------------------------------------------------------------------
  HomPart := proc(p, tv::list, w::nonnegint)
    local L, q, v;
    if nops(tv) = 0 then
      if w = 0 then return p else return 0 end if
    end if;
    q := expand(subs({seq(v = L*v, v in tv)}, expand(p)));
    return expand(coeff(collect(q, L), L, w));
  end proc;

  TDeg := proc(p, tv::list)
    local q;
    q := expand(p);
    if nops(tv) = 0 or q = 0 then return 0 end if;
    return degree(q, {op(tv)});
  end proc;

  Truncate := proc(p, tv::list, T::nonnegint)
    local w;
    return expand(add(HomPart(p, tv, w), w = 0 .. T));
  end proc;

  # 全次数ちょうど w の単項式のリスト
  Monomials := proc(tv::list, w::nonnegint)
    local l, L;
    l := nops(tv);
    if l = 0 then return [1] end if;
    if w = 0 then return [1] end if;
    if l = 1 then return [tv[1]^w] end if;
    L := [seq(seq(tv[1]^i * mm, mm in Monomials(tv[2 .. l], w - i)), i = 0 .. w)];
    return L;
  end proc;

  # ------------------------------------------------------------------
  # 部分終結式行列 S_{k-1}(F,G)：行は x の高次から、列は u ブロック（高次から）→ v ブロック
  # ------------------------------------------------------------------
  SubresMat := proc(F, G, k::nonnegint, x)
    local m, n, f, g, nu, nv, rows, S, i, j;
    m := degree(F, x); n := degree(G, x);
    f := [seq(coeff(expand(F), x, i), i = 0 .. m)];
    g := [seq(coeff(expand(G), x, i), i = 0 .. n)];
    nu := n - k + 1; nv := m - k + 1;
    rows := m + n - k + 1;
    S := Matrix(rows, nu + nv, 0);
    for j from 1 to nu do          # 列 j <-> x^(n-k-j+1) * F
      for i from 0 to m do
        S[(m + n - k) - (i + n - k - j + 1) + 1, j] := f[i+1];
      end do;
    end do;
    for j from 1 to nv do          # 列 nu+j <-> x^(m-k-j+1) * G
      for i from 0 to n do
        S[(m + n - k) - (i + m - k - j + 1) + 1, nu + j] := g[i+1];
      end do;
    end do;
    return S;
  end proc;

  # 行列の斉次部分 S^(0), ..., S^(wmax)
  HomMats := proc(S, tv::list, wmax::nonnegint)
    local H, w;
    H := table();
    for w from 0 to wmax do
      H[w] := map(HomPart, S, tv, w);
    end do;
    return H;
  end proc;

  # S^(0) の薄い SVD：U (N x K), sig (K), V (K x K)
  InitSVD := proc(S0)
    local N, K, U, sv, Vt;
    N := LinearAlgebra:-RowDimension(S0);
    K := LinearAlgebra:-ColumnDimension(S0);
    U, sv, Vt := LinearAlgebra:-SingularValues(
                   Matrix(evalf(S0), datatype = float[8]), output = [':-U', ':-S', ':-Vt']);
    return Matrix(U[1 .. N, 1 .. K]), Vector([seq(sv[i], i = 1 .. K)]),
           LinearAlgebra:-Transpose(Matrix(Vt));
  end proc;

  # t の単項式ごとの係数（0 なら空リスト）
  CoeffEqs := proc(g, tv::list)
    local q;
    q := expand(g);
    if q = 0 then return [] end if;
    if nops(tv) = 0 then return [q] end if;
    return [coeffs(collect(q, tv, 'distributed'), tv)];
  end proc;

  # 行列（多項式成分）×ベクトル（多項式成分）
  MatVec := proc(A, z)
    local N, K, i;
    N := LinearAlgebra:-RowDimension(A); K := LinearAlgebra:-ColumnDimension(A);
    return Vector(N, i -> expand(add(A[i, j]*z[j], j = 1 .. K)));
  end proc;

  # p^(w) = - sum_{j=1}^{min(w,wS)} S^(j) r^(w-j)
  Rhs := proc(H, wS, R, w, N)
    local p, j;
    p := Vector(N, 0);
    for j from 1 to min(w, wS) do
      p := p - MatVec(H[j], R[w - j]);
    end do;
    return map(expand, p);
  end proc;

  # ------------------------------------------------------------------
  # 標準の lifting。jstar = 0 なら |v_K| 最大の成分で正規化。
  # 戻り値: R（table, R[w] は全次数 w の斉次部分のベクトル）, sig, jstar, nu, nv
  # ------------------------------------------------------------------
  Lift := proc(F, G, tv::list, k::posint, W::nonnegint, x, {jstar::nonnegint := 0})
    local S, N, K, wS, H, U, sig, V, vK, js, R, w, p, ph, z, beta, i, j, a, m, n;
    m := degree(F, x); n := degree(G, x);
    S := SubresMat(F, G, k, x);
    N := LinearAlgebra:-RowDimension(S); K := LinearAlgebra:-ColumnDimension(S);
    wS := max(seq(seq(TDeg(S[i, j], tv), j = 1 .. K), i = 1 .. N));
    H := HomMats(S, tv, min(W, wS));
    U, sig, V := InitSVD(H[0]);
    vK := V[1 .. K, K];
    js := jstar;
    if js = 0 then
      js := 1;
      for i from 2 to K do if abs(vK[i]) > abs(vK[js]) then js := i end if end do;
    end if;
    R := table();
    R[0] := Vector(K, i -> vK[i]/vK[js]);
    for w from 1 to W do
      p := Rhs(H, min(W, wS), R, w, N);
      ph := Vector(K - 1, i -> expand(add(U[a, i]*p[a], a = 1 .. N)/sig[i]));
      z := Vector(K, i -> expand(add(V[i, j]*ph[j], j = 1 .. K - 1)));
      beta := -z[js]/vK[js];
      R[w] := map(expand, z + beta*vK);
    end do;
    return R, sig, js, n - k + 1, m - k + 1;
  end proc;

  # ------------------------------------------------------------------
  # 悪条件の場合：自由方向を span{v_{K-1}, v_K} に広げる。
  #   r^(w) = y^(w) + alpha^(w) vt + (正規化のための v_K 成分),
  #   vt = v_{K-1} - (v_{K-1}[js]/v_K[js]) v_K   （r[js] を動かさない方向）
  # alpha^(w-1) は次の段の u_{K-1} 方程式
  #   u_{K-1}^T p^(w) = sigma_{K-1} alpha^(w)
  # から最小二乗で決める。1 パス目は右辺を 0、以後は前のパスの alpha^(w) を使う。
  # ------------------------------------------------------------------
  LiftIll := proc(F, G, tv::list, k::posint, W::nonnegint, x,
                  {jstar::nonnegint := 0, npass::posint := 6})
    local S, N, K, wS, H, U, sig, V, vK, v1, u1, s1, vt, js, R, w, p, ph, z, beta,
          i, j, a, m, n, ps, alphas, aprev, D, mons, av, alpha, g, eqs, A, b, sol, Wh;
    m := degree(F, x); n := degree(G, x);
    S := SubresMat(F, G, k, x);
    N := LinearAlgebra:-RowDimension(S); K := LinearAlgebra:-ColumnDimension(S);
    wS := max(seq(seq(TDeg(S[i, j], tv), j = 1 .. K), i = 1 .. N));
    Wh := min(W + 1, wS);
    H := HomMats(S, tv, Wh);
    U, sig, V := InitSVD(H[0]);
    vK := V[1 .. K, K]; v1 := V[1 .. K, K - 1]; u1 := U[1 .. N, K - 1]; s1 := sig[K - 1];
    js := jstar;
    if js = 0 then
      js := 1;
      for i from 2 to K do if abs(vK[i]) > abs(vK[js]) then js := i end if end do;
    end if;
    vt := v1 - (v1[js]/vK[js])*vK;
    aprev := table();
    for ps from 1 to npass do
      R := table(); alphas := table();
      D := Vector(K, i -> vK[i]/vK[js]);
      for w from 1 to W + 1 do
        mons := Monomials(tv, w - 1);
        av := [seq(cat(`_a`, i), i = 1 .. nops(mons))];
        alpha := add(av[i]*mons[i], i = 1 .. nops(mons));
        R[w - 1] := map(expand, D + alpha*vt);
        p := Rhs(H, Wh, R, w, N);
        g := expand(add(u1[a]*p[a], a = 1 .. N));
        if ps > 1 and w <= W then
          g := expand(g - s1*aprev[w]);
        end if;
        eqs := CoeffEqs(g, tv);
        if nops(eqs) = 0 then
          alpha := 0;
        else
          A, b := LinearAlgebra:-GenerateMatrix(eqs, av);
          sol := LinearAlgebra:-LeastSquares(Matrix(evalf(A), datatype = float[8]),
                                             Vector(evalf(b), datatype = float[8]));
          alpha := add(sol[i]*mons[i], i = 1 .. nops(mons));
        end if;
        alphas[w - 1] := alpha;
        R[w - 1] := map(expand, D + alpha*vt);
        if w <= W then
          p := Rhs(H, Wh, R, w, N);
          ph := Vector(K - 2, i -> expand(add(U[a, i]*p[a], a = 1 .. N)/sig[i]));
          z := Vector(K, i -> expand(add(V[i, j]*ph[j], j = 1 .. K - 2)));
          beta := -z[js]/vK[js];
          D := map(expand, z + beta*vK);
        end if;
      end do;
      aprev := copy(alphas);
    end do;
    return R, sig, js, n - k + 1, m - k + 1;
  end proc;

  # ------------------------------------------------------------------
  # 密な場合：|v_K| 最大成分で正規化して 2*Dmax+1 次まで持ち上げ（冪級数）、
  # tdeg c <= D の乗数 c(t) で c*r^ の D+1..2D+1 次の部分が消えるものを SVD で探す。
  # D = 0,1,... と上げ、最小特異値比が tol を切った最初の D を採る。
  # ------------------------------------------------------------------
  LiftDense := proc(F, G, tv::list, k::posint, Dmax::nonnegint, x, {tol := 1e-8})
    local R, sig, js, nu, nv, K, D, cm, cv, cpoly, Rs, eqs, w, comp, Q, A, b, sv, Vt, c,
          rr, i, ratio, out;
    R, sig, js, nu, nv := Lift(F, G, tv, k, 2*Dmax + 1, x);
    K := nu + nv;
    Rs := Vector(K, comp -> add(R[w][comp], w = 0 .. 2*Dmax + 1));
    for D from 0 to Dmax do
      cm := [seq(op(Monomials(tv, w)), w = 0 .. D)];
      cv := [seq(cat(`_c`, i), i = 1 .. nops(cm))];
      cpoly := add(cv[i]*cm[i], i = 1 .. nops(cm));
      eqs := [];
      for comp from 1 to K do
        Q := expand(cpoly*Truncate(Rs[comp], tv, 2*D + 1));
        for w from D + 1 to 2*D + 1 do
          eqs := [op(eqs), op(CoeffEqs(HomPart(Q, tv, w), tv))];
        end do;
      end do;
      A, b := LinearAlgebra:-GenerateMatrix(eqs, cv);
      sv, Vt := LinearAlgebra:-SingularValues(Matrix(evalf(A), datatype = float[8]),
                                               output = [':-S', ':-Vt']);
      ratio := sv[nops(cv)]/sv[1];
      if ratio < tol or D = Dmax then
        c := add(Vt[nops(cv), i]*cm[i], i = 1 .. nops(cm));
        break;
      end if;
    end do;
    out := table();
    for w from 0 to D do
      out[w] := Vector(K, comp -> HomPart(expand(c*Rs[comp]), tv, w));
    end do;
    return out, sig, js, nu, nv, D, ratio;
  end proc;

  # r（table of 斉次部分）-> [F~, G~]
  VecToCof := proc(R, W, nu, nv, x)
    local Gt, Ft, j, w;
    Gt := expand(add(add(R[w][j]*x^(nu - j), j = 1 .. nu), w = 0 .. W));
    Ft := expand(-add(add(R[w][nu + j]*x^(nv - j), j = 1 .. nv), w = 0 .. W));
    return Ft, Gt;
  end proc;

  # ------------------------------------------------------------------
  # 余因子から GCD：[Conv(F~); Conv(G~)] c = [F; G] を lifting で解く。
  # ------------------------------------------------------------------
  GcdFromCof := proc(F, G, Ft, Gt, tv::list, k::nonnegint, W::nonnegint, x)
    local m, n, A, rows, i, j, ft, gt, wA, H, A0p, Cw, w, b, fF, fG, cc, a;
    m := degree(F, x); n := degree(G, x);
    rows := m + n + 2;
    ft := [seq(coeff(Ft, x, i), i = 0 .. m - k)];
    gt := [seq(coeff(Gt, x, i), i = 0 .. n - k)];
    A := Matrix(rows, k + 1, 0);
    for j from 0 to k do               # 列 j+1 <-> x^(k-j)
      for i from 0 to m - k do A[m - (i + k - j) + 1, j + 1] := ft[i + 1]; end do;
      for i from 0 to n - k do A[m + 1 + n - (i + k - j) + 1, j + 1] := gt[i + 1]; end do;
    end do;
    wA := max(seq(seq(TDeg(A[i, j], tv), j = 1 .. k + 1), i = 1 .. rows));
    H := HomMats(A, tv, min(W, wA));
    A0p := LinearAlgebra:-MatrixInverse(Matrix(evalf(H[0]), datatype = float[8]),
                                        method = pseudo);
    fF := expand(F); fG := expand(G);
    Cw := table();
    for w from 0 to W do
      b := Vector(rows, i -> `if`(i <= m + 1, HomPart(coeff(fF, x, m - i + 1), tv, w),
                                   HomPart(coeff(fG, x, n - (i - m - 1) + 1), tv, w)));
      for j from 1 to min(w, wA) do
        b := b - MatVec(H[j], Cw[w - j]);
      end do;
      Cw[w] := Vector(k + 1, i -> expand(add(A0p[i, a]*b[a], a = 1 .. rows)));
    end do;
    cc := expand(add(add(Cw[w][j + 1]*x^(k - j), j = 0 .. k), w = 0 .. W));
    return cc;
  end proc;

  # 一式：lifting → 余因子 → GCD
  #   mode = "std" | "ill"
  Solve := proc(F, G, tv::list, k::posint, W::nonnegint, x,
                {mode::string := "std", jstar::nonnegint := 0, npass::posint := 6})
    local R, sig, js, nu, nv, Ft, Gt, C;
    if mode = "ill" then
      R, sig, js, nu, nv := LiftIll(F, G, tv, k, W, x, ':-jstar' = jstar, ':-npass' = npass);
    else
      R, sig, js, nu, nv := Lift(F, G, tv, k, W, x, ':-jstar' = jstar);
    end if;
    Ft, Gt := VecToCof(R, W, nu, nv, x);
    C := Truncate(GcdFromCof(F, G, Ft, Gt, tv, k, W, x), tv, W);
    return C, Truncate(Ft, tv, W), Truncate(Gt, tv, W), sig;
  end proc;

  # 標準と悪条件用の両方を実行し、後退誤差の小さい方を採る
  SolveAdaptive := proc(F, G, tv::list, k::posint, W::nonnegint, x, {jstar::nonnegint := 0})
    local r1, r2, b1, b2;
    r1 := [Solve(F, G, tv, k, W, x, ':-mode' = "std", ':-jstar' = jstar)];
    r2 := [Solve(F, G, tv, k, W, x, ':-mode' = "ill", ':-jstar' = jstar)];
    b1 := BackErr(F, G, r1[1], r1[2], r1[3], x, tv);
    b2 := BackErr(F, G, r2[1], r2[2], r2[3], x, tv);
    if b1 <= b2 then return op(r1), "std" else return op(r2), "ill" end if;
  end proc;

  # ------------------------------------------------------------------
  # 比較用：多変数部分終結式を数値行列として作り、SVD で零ベクトル。
  # 余因子の t についての全次数は Wc 以下と仮定（真の値を与えれば直接法に有利な条件）。
  # 戻り値: Ft, Gt, [行数, 列数], SVD の CPU 時間
  # ------------------------------------------------------------------
  DirectNull := proc(F, G, tv::list, k::posint, Wc::nonnegint, x)
    local m, n, cm, gv, fv, Gt, Ft, res, eqs, A, b, sv, Vt, t0, tsvd, nc, i, j, w, vals, allv;
    m := degree(F, x); n := degree(G, x);
    cm := [seq(op(Monomials(tv, w)), w = 0 .. Wc)];
    nc := nops(cm);
    gv := [seq(seq(cat(`_g`, i, `_`, j), j = 1 .. nc), i = 0 .. n - k)];
    fv := [seq(seq(cat(`_f`, i, `_`, j), j = 1 .. nc), i = 0 .. m - k)];
    Gt := add(add(cat(`_g`, i, `_`, j)*cm[j], j = 1 .. nc)*x^(n - k - i), i = 0 .. n - k);
    Ft := add(add(cat(`_f`, i, `_`, j)*cm[j], j = 1 .. nc)*x^(m - k - i), i = 0 .. m - k);
    res := expand(F*Gt - G*Ft);
    eqs := [coeffs(collect(res, [x, op(tv)], 'distributed'), [x, op(tv)])];
    allv := [op(gv), op(fv)];
    A, b := LinearAlgebra:-GenerateMatrix(eqs, allv);
    A := Matrix(evalf(A), datatype = float[8]);
    t0 := time();
    sv, Vt := LinearAlgebra:-SingularValues(A, output = [':-S', ':-Vt']);
    tsvd := time() - t0;
    vals := {seq(allv[i] = Vt[nops(allv), i], i = 1 .. nops(allv))};
    return expand(subs(vals, Ft)), expand(subs(vals, Gt)),
           [LinearAlgebra:-RowDimension(A), LinearAlgebra:-ColumnDimension(A)], tsvd;
  end proc;

  # ------------------------------------------------------------------
  # 次数判定・展開点
  # ------------------------------------------------------------------
  DetectDegree := proc(F, G, tv::list, x, tol)
    local F0, G0, d, S0, sv, K, v;
    F0 := expand(subs({seq(v = 0, v in tv)}, F));
    G0 := expand(subs({seq(v = 0, v in tv)}, G));
    for d from min(degree(F0, x), degree(G0, x)) by -1 to 1 do
      S0 := Matrix(evalf(SubresMat(F0, G0, d, x)), datatype = float[8]);
      sv := LinearAlgebra:-SingularValues(S0);
      K := LinearAlgebra:-ColumnDimension(S0);
      if sv[K]/sv[1] < tol then return d end if;
    end do;
    return 0;
  end proc;

  SigmaRatio := proc(F, G, tv::list, k::posint, x)
    local S0, sv, K, v;
    S0 := SubresMat(expand(subs({seq(v = 0, v in tv)}, F)),
                    expand(subs({seq(v = 0, v in tv)}, G)), k, x);
    S0 := Matrix(evalf(S0), datatype = float[8]);
    sv := LinearAlgebra:-SingularValues(S0);
    K := LinearAlgebra:-ColumnDimension(S0);
    return sv[K - 1]/sv[1];
  end proc;

  Shift := proc(P, tv::list, s::list)
    local i;
    return expand(subs({seq(tv[i] = tv[i] + s[i], i = 1 .. nops(tv))}, P));
  end proc;

  # 候補点 ncand 個から sigma_{K-1}/sigma_1 最大の点を選ぶ（単変数の SVD だけで判定）
  ChooseShift := proc(F, G, tv::list, k::posint, x, {ncand::posint := 10})
    local best, bs, i, j, s, r, rr;
    rr := rand(-100 .. 100);
    best := SigmaRatio(F, G, tv, k, x); bs := [seq(0, j = 1 .. nops(tv))];
    for i from 1 to ncand do
      s := [seq(rr()/100, j = 1 .. nops(tv))];
      r := SigmaRatio(Shift(F, tv, s), Shift(G, tv, s), tv, k, x);
      if r > best then best := r; bs := s end if;
    end do;
    return bs, best;
  end proc;

  # ------------------------------------------------------------------
  # 誤差の尺度
  # ------------------------------------------------------------------
  PolyNorm := proc(P, vars::list)
    local q, c;
    q := expand(P);
    if q = 0 then return 0.0 end if;
    return sqrt(add(evalf(c)^2, c in [coeffs(q, vars)]));
  end proc;

  # min_s || s*Pc - Pt || / || Pt ||
  RelErr := proc(Pc, Pt, vars::list)
    local a, b, ab, aa, s;
    a := expand(evalf(Pc)); b := expand(evalf(Pt));
    aa := PolyNorm(a, vars)^2;
    if aa = 0 then return 1.0 end if;
    # <a,b> = (||a+b||^2 - ||a-b||^2)/4
    ab := (PolyNorm(a + b, vars)^2 - PolyNorm(a - b, vars)^2)/4;
    s := ab/aa;
    return PolyNorm(s*a - b, vars)/PolyNorm(b, vars);
  end proc;

  BackErr := proc(F, G, C, Ft, Gt, x, tv::list)
    local vars;
    vars := [x, op(tv)];
    return max(PolyNorm(F - C*Ft, vars)/PolyNorm(F, vars),
               PolyNorm(G - C*Gt, vars)/PolyNorm(G, vars));
  end proc;

  # ------------------------------------------------------------------
  # ランダムな多項式
  #   dense = true : x の各係数は全次数 dt 以下の全単項式を持つ
  #   dense = false: 定数項 + (nterms-1) 個のランダムな単項式
  #   monic = true : x の最高次係数は 1
  # ------------------------------------------------------------------
  RandPoly := proc(x, tv::list, dx::nonnegint, dt::nonnegint,
                   {dense::truefalse := false, monic::truefalse := true, nterms::posint := 4})
    local r, ms, P, e, c, i, w, pick, sg;
    r := rand(-10^6 .. 10^6);
    pick := rand(1 .. 10^9);
    ms := [seq(op(Monomials(tv, w)), w = 0 .. dt)];
    P := 0;
    for e from 0 to dx do
      if e = dx and monic then
        c := 1;
      elif dense then
        c := add(evalf(r()/10^6)*ms[i], i = 1 .. nops(ms));
      else
        sg := `if`(r() >= 0, 1, -1);
        c := sg*(0.5 + 0.5*abs(evalf(r()/10^6)));
        for i from 1 to nterms - 1 do
          c := c + evalf(r()/10^6)*ms[(pick() mod nops(ms)) + 1];
        end do;
      end if;
      P := P + c*x^e;
    end do;
    return expand(P);
  end proc;

  # 同じ台の上に相対ノルム eps の摂動を加える
  AddNoise := proc(P, x, tv::list, eps)
    local r, vars, cl, ml, E, i;
    if eps = 0 then return P end if;
    r := rand(-10^6 .. 10^6);
    vars := [x, op(tv)];
    cl := [coeffs(expand(P), vars, 'ml')];
    ml := [ml];
    E := add(evalf(r()/10^6)*ml[i], i = 1 .. nops(ml));
    return expand(P + eps*PolyNorm(P, vars)/PolyNorm(E, vars)*E);
  end proc;

  # LaTeX 用の数値表記  $1.2\times 10^{-3}$
  Fmt := proc(v)
    local e, m;
    if v = 0 then return "$0$" end if;
    e := floor(log10(abs(evalf(v))));
    m := evalf(v/10^e);
    if abs(m) >= 9.95 then m := m/10; e := e + 1 end if;
    return sprintf("$%.1f\\times 10^{%d}$", m, e);
  end proc;

end module:
