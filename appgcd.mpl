# appgcd.mpl -- 近似 GCD のための Sylvester / Bezout 行列とリフティング（Maple 版）
#
# Python 版 appgcd.py と同じ内容。もとになる論文は ../PDF/ :
#   [SNC09] M. Sanuki, Computing Multivariate Approximate GCD Based on
#           Barnett's Theorem, Proc. SNC 2009, 149-157.
#   [SVD26] M. Sanuki, An SVD-and-Lifting Based Framework for Approximate
#           Multivariate GCD Computation, 2026.
#
# 使い方:   read "appgcd.mpl":
#           AppGCD:-BezoutLift(F, G, [u], 2, 8, x);

AppGCD := module()
  option package;
  export
    HomPart, TDeg, Truncate, MatHomPart, Chop,
    BezoutM, BezoutCheck, BezoutSub, SylvesterSub, SplitNullVec,
    BezoutLift, BezoutLiftAdjoint, SylvesterLift,
    RuizSym, SymSolve,
    GammaOf, GcdRelErr, CoeffMax,
    Examples;
  local `&lam`;

  # ------------------------------------------------------------------
  # 全次数 w の斉次成分を取り出す。
  # t -> L*t と置いて L の w 次の係数を見る、という古典的な手。
  # ------------------------------------------------------------------
  HomPart := proc(p, tv::list, w::nonnegint)
    local L, q, v;
    if nops(tv) = 0 then
      if w = 0 then return p else return 0 end if
    end if;
    q := expand(subs({seq(v = L*v, v in tv)}, expand(p)));
    return expand(coeff(collect(q, L), L, w));
  end proc;

  # 行列・ベクトルの各成分に HomPart を適用
  MatHomPart := proc(M, tv::list, w::nonnegint)
    return map(q -> HomPart(q, tv, w), M);
  end proc;

  # t についての全次数
  TDeg := proc(p, tv::list)
    if nops(tv) = 0 or p = 0 then return 0 end if;
    return degree(expand(p), {op(tv)});
  end proc;

  # 全次数 T で打ち切る
  Truncate := proc(p, tv::list, T::nonnegint)
    local w;
    return expand(add(HomPart(p, tv, w), w = 0 .. T));
  end proc;

  # 係数の最大絶対値
  CoeffMax := proc(p, x, tv::list)
    local q, cl, c;
    q := expand(p);
    if q = 0 then return 0.0 end if;
    cl := [coeffs(q, [x, op(tv)])];
    return max(seq(abs(evalf(c)), c in cl));
  end proc;

  # 微小な係数を落とす。浮動小数の最小二乗から出る 1e-14 程度のゴミを消す。
  Chop := proc(p, x, tv::list, tol := 1e-9, dig := 12)
    local q, mx, cl, ml, i, res;
    q := expand(evalf(p));
    if q = 0 then return 0 end if;
    mx := CoeffMax(q, x, tv);
    if mx = 0 then return 0 end if;
    cl := [coeffs(q, [x, op(tv)], 'ml')];
    ml := [ml];
    res := add(`if`(abs(evalf(cl[i])) > tol*mx,
                    round(evalf(cl[i])*10^dig)/10^dig*ml[i], 0), i = 1 .. nops(cl));
    return expand(res);
  end proc;

  # ------------------------------------------------------------------
  # Bezout 行列   B(F,G) = (F(x)G(y)-F(y)G(x))/(x-y) = sum b_ij x^i y^j
  # 閉形式: b_ij = sum_{l=0}^{min(i,j)} ( f_{i+j+1-l} g_l - f_l g_{i+j+1-l} )
  # 対称。rank Bez = n - deg gcd。
  # ------------------------------------------------------------------
  BezoutM := proc(F, G, x, nn := -1)
    local n, f, g, fc, gc, B, i, j, l, s;
    n := `if`(nn > 0, nn, max(degree(F, x), degree(G, x)));
    f := [seq(coeff(expand(F), x, i), i = 0 .. n)];
    g := [seq(coeff(expand(G), x, i), i = 0 .. n)];
    fc := r -> `if`(r >= 0 and r <= n, f[r+1], 0);
    gc := r -> `if`(r >= 0 and r <= n, g[r+1], 0);
    B := Matrix(n, n, 0);
    for i from 0 to n-1 do
      for j from i to n-1 do
        s := expand(add(fc(i+j+1-l)*gc(l) - fc(l)*gc(i+j+1-l), l = 0 .. min(i,j)));
        B[i+1, j+1] := s;
        B[j+1, i+1] := s;
      end do;
    end do;
    return B;
  end proc;

  # 定義どおり割って作ったものと突き合わせる（検算用）。
  # Maple 組み込みの LinearAlgebra:-BezoutMatrix とは全体の符号が異なることがある
  # （B(F,G) = -B(G,F) なので階数も零空間も GCD も変わらない）。
  BezoutCheck := proc(F, G, x)
    local y, n, Bdef, Bcf, i, j, P;
    n := max(degree(F, x), degree(G, x));
    P := expand(normal((F*subs(x = y, G) - subs(x = y, F)*G)/(x - y)));
    Bdef := Matrix(n, n, (i, j) -> coeff(coeff(P, x, i-1), y, j-1));
    Bcf := BezoutM(F, G, x, n);
    return simplify(expand(Bdef - Bcf));
  end proc;

  # Bez_{k,n-1} : 最後の n-k 行・n-k 列（[SNC09] 3 節）
  BezoutSub := proc(B, k::nonnegint)
    local n;
    n := LinearAlgebra:-RowDimension(B);
    return B[k+1 .. n, k+1 .. n];
  end proc;

  # ------------------------------------------------------------------
  # 部分終結式行列  (u,v) -> u*F + v*G,  deg u <= n-d, deg v <= m-d
  # 形は (m+n-d+1) x (m+n-2d+2)。deg gcd = d なら零空間はちょうど 1 次元で、
  # 零ベクトルは (G~_{n-d},...,G~_0, -F~_{m-d},...,-F~_0)。
  # ------------------------------------------------------------------
  SylvesterSub := proc(F, G, d::nonnegint, x)
    local m, n, f, g, nu, nv, rows, S, i, j;
    m := degree(F, x); n := degree(G, x);
    f := [seq(coeff(expand(F), x, i), i = 0 .. m)];
    g := [seq(coeff(expand(G), x, i), i = 0 .. n)];
    nu := n - d + 1; nv := m - d + 1;
    rows := m + n - d + 1;
    S := Matrix(rows, nu + nv, 0);
    for j from 0 to nu-1 do
      for i from 0 to m do S[i + (nu-1-j) + 1, j+1] := f[i+1]; end do;
    end do;
    for j from 0 to nv-1 do
      for i from 0 to n do S[i + (nv-1-j) + 1, nu+j+1] := g[i+1]; end do;
    end do;
    # 高次の行を上に（Python 版と同じ並び）
    return Matrix(rows, nu+nv, (i, j) -> S[rows+1-i, j]);
  end proc;

  # 零ベクトル -> 余因子 (F~, G~)
  SplitNullVec := proc(r::list, m, n, d, x)
    local nu, Gt, Ft, j;
    nu := n - d + 1;
    Gt := expand(add(r[j+1]*x^(nu-1-j), j = 0 .. nu-1));
    Ft := expand(-add(r[nu+j+1]*x^(m-d-j), j = 0 .. m-d));
    return [Ft, Gt];
  end proc;

  # ------------------------------------------------------------------
  # 対称平衡化（Ruiz）。Bez は対称なので D M D は対称性を保つ。
  # [SNC09] (5.5) の「Bez^{-1} の (i,j) 成分が O(1/gamma^{i+j-2})」という評価は
  # 対角スケーリングで消せる形をしている。
  # ------------------------------------------------------------------
  RuizSym := proc(M, it::posint := 30)
    local n, D, A, i, j, r, k, kk;
    n := LinearAlgebra:-RowDimension(M);
    D := Vector(n, 1.0);
    A := Matrix(evalf(M));
    for kk to it do
      r := Vector(n, i -> sqrt(max(seq(abs(A[i, j]), j = 1 .. n), 1e-300)));
      for i to n do D[i] := D[i]/r[i]; end do;
      A := Matrix(n, n, (i, j) -> A[i, j]/(r[i]*r[j]));
    end do;
    return D;
  end proc;

  # 対称行列 M0 の解法。scale="ruiz" で対称平衡化、solver="eigh" で対称固有分解。
  # 戻り値は rhs を受け取る手続き。
  SymSolve := proc(M0, {scale::string := "none", solver::string := "lu"})
    local n, D, Ms, i, j, ev, Q, lam, solve1;
    local zz;
    n := LinearAlgebra:-RowDimension(M0);
    if scale = "ruiz" then D := RuizSym(M0) else D := Vector(n, 1.0) end if;
    Ms := Matrix(n, n, (i, j) -> evalf(M0[i, j])*D[i]*D[j]);
    if solver = "eigh" then
      ev := LinearAlgebra:-Eigenvectors(Matrix(Ms, shape = symmetric));
      lam := ev[1]; Q := ev[2];
      solve1 := proc(rhs)
        local z, w2, y2;
        z := Vector(n, i -> D[i]*rhs[i]);
        w2 := LinearAlgebra:-Transpose(Q) . z;
        y2 := Vector(n, i -> w2[i]/lam[i]);
        z := Q . y2;
        return Vector(n, i -> D[i]*z[i]);
      end proc;
    else
      solve1 := proc(rhs)
        local z, y2;
        z := Vector(n, i -> D[i]*rhs[i]);
        y2 := LinearAlgebra:-LinearSolve(Ms, z);
        return Vector(n, i -> D[i]*y2[i]);
      end proc;
    end if;
    return solve1, Ms;
  end proc;

  # ------------------------------------------------------------------
  # Bezout lifting（[SNC09] Algorithm 1）
  #   Bez^(0)_{k,n-1} dc^(w) = db~^(w) - sum_{j>=1} dBez^(j)_{k,n-1} dc^(w-j)
  # 解の第 1 成分が c_i/c_k（Barnett の定理）。
  # ------------------------------------------------------------------
  BezoutLift := proc(F, G, tv::list, k::posint, T::nonnegint, x,
                     {scale::string := "none", solver::string := "lu"})
    local n, B, Bsub, Brhs, M0, wmax, C, w, j, R, dM, gcdpoly, i, co, sol, Ms, a, b, c;
    n := max(degree(F, x), degree(G, x));
    B := BezoutM(F, G, x, n);
    Bsub := BezoutSub(B, k);
    Brhs := B[k+1 .. n, 1 .. k];
    M0 := MatHomPart(Bsub, tv, 0);
    wmax := max(seq(seq(TDeg(Bsub[i, j], tv), j = 1 .. n-k), i = 1 .. n-k));
    sol, Ms := SymSolve(M0, ':-scale' = scale, ':-solver' = solver);
    C := table();
    for w from 0 to T do
      R := Matrix(MatHomPart(Brhs, tv, w));
      for j to min(w, wmax) do
        dM := MatHomPart(Bsub, tv, j);
        R := map(expand, R - LinearAlgebra:-Multiply(dM, C[w-j]));
      end do;
      C[w] := Matrix(n-k, k, 0);
      for i to k do
        C[w][1 .. n-k, i] := sol(Vector(n-k, a -> R[a, i]));
      end do;
    end do;
    gcdpoly := x^k;
    for i to k do
      co := add(C[w][1, i], w = 0 .. T);
      gcdpoly := gcdpoly + expand(co)*x^(i-1);
    end do;
    return expand(gcdpoly), C, M0;
  end proc;

  # ------------------------------------------------------------------
  # 随伴リフティング：第 1 成分しか要らないので B y = e_1 を 1 本だけ持ち上げる。
  #   e_1^T x_i = (B^{-T} e_1)^T b~_i = y^T b~_i
  # Bez は対称なので B^T = B、随伴系は元と同じ行列（分解を使い回せる）。
  # 持ち上げるベクトルが k 本から 1 本になる。
  # ------------------------------------------------------------------
  BezoutLiftAdjoint := proc(F, G, tv::list, k::posint, T::nonnegint, x,
                            {scale::string := "none", solver::string := "lu"})
    local n, B, Bsub, Brhs, M0, wmax, Y, w, j, R, dM, e1, gcdpoly, i, co, sol, Ms, N, a, c;
    n := max(degree(F, x), degree(G, x));
    B := BezoutM(F, G, x, n);
    Bsub := BezoutSub(B, k);
    Brhs := B[k+1 .. n, 1 .. k];
    N := n - k;
    M0 := MatHomPart(Bsub, tv, 0);
    wmax := max(seq(seq(TDeg(Bsub[i, j], tv), j = 1 .. N), i = 1 .. N));
    sol, Ms := SymSolve(M0, ':-scale' = scale, ':-solver' = solver);
    e1 := Vector(N, 0); e1[1] := 1;
    Y := table();
    for w from 0 to T do
      if w = 0 then R := Vector(e1) else R := Vector(N, 0) end if;
      for j to min(w, wmax) do
        dM := MatHomPart(Bsub, tv, j);
        R := map(expand, R - LinearAlgebra:-Multiply(dM, Y[w-j]));
      end do;
      Y[w] := sol(R);
    end do;
    gcdpoly := x^k;
    for i to k do
      co := Truncate(expand(add(add(Y[w][a]*Brhs[a, i], a = 1 .. N), w = 0 .. T)), tv, T);
      gcdpoly := gcdpoly + co*x^(i-1);
    end do;
    return expand(gcdpoly), Y;
  end proc;

  # ------------------------------------------------------------------
  # 部分終結式の零空間を lifting で求める（r[1] = 1 と正規化して 1 列目を右辺へ）。
  # 縦長なので最小二乗で解く。
  # ------------------------------------------------------------------
  SylvesterLift := proc(F, G, tv::list, d::nonnegint, T::nonnegint, x)
    local m, n, S, rows, cols, A, b, A0, wmax, Z, w, j, R, dA, db, r, i, ncol, a, c;
    local ii;
    m := degree(F, x); n := degree(G, x);
    S := SylvesterSub(F, G, d, x);
    rows := LinearAlgebra:-RowDimension(S); cols := LinearAlgebra:-ColumnDimension(S);
    ncol := cols - 1;
    A := S[1 .. rows, 2 .. cols];
    b := Vector(rows, i -> -S[i, 1]);
    A0 := MatHomPart(A, tv, 0);
    wmax := max(seq(seq(TDeg(A[i, j], tv), j = 1 .. ncol), i = 1 .. rows),
                seq(TDeg(b[i], tv), i = 1 .. rows));
    Z := table();
    for w from 0 to T do
      R := Vector(rows, i -> HomPart(b[i], tv, w));
      for j to min(w, wmax) do
        dA := MatHomPart(A, tv, j);
        R := map(expand, R - LinearAlgebra:-Multiply(dA, Z[w-j]));
      end do;
      Z[w] := LinearAlgebra:-LeastSquares(evalf(A0), evalf(R));
    end do;
    r := [1, seq(Chop(add(Z[w][i], w = 0 .. T), x, tv), i = 1 .. ncol)];
    return r, S;
  end proc;

  # ------------------------------------------------------------------
  # gamma = |lc(C(x,0))| / ||C(x,0)||_inf   （[SNC09] 5 節）
  # ------------------------------------------------------------------
  GammaOf := proc(C, tv::list, x)
    local C0, k, cl, v, i;
    C0 := expand(subs({seq(v = 0, v in tv)}, C));
    k := degree(C0, x);
    cl := [seq(abs(evalf(coeff(C0, x, i))), i = 0 .. k)];
    return abs(evalf(coeff(C0, x, k)))/max(op(cl));
  end proc;

  # モニック正規化したうえでの相対誤差
  GcdRelErr := proc(Cc, Ctrue, tv::list, x)
    local a, b, d;
    a := expand(Cc/lcoeff(expand(Cc), x));
    b := expand(Ctrue/lcoeff(expand(Ctrue), x));
    d := expand(a - b);
    return evalf(CoeffMax(d, x, tv)/CoeffMax(b, x, tv));
  end proc;

  # ------------------------------------------------------------------
  # 例の一覧（Python 版と同じ）
  # ------------------------------------------------------------------
  Examples := proc(x, u, t1, t2, t3)
    local E, A1, B1, C26;
    A1 := x^3 + (u^3+u-1)*x + 1;
    B1 := x^3 - (u^3+u+2)*x + 1;
    C26 := x^3 + (1+t2-2*t1+t1^2)*x + 3;
    E := table();
    E["snc-ex1"] := [x^2-u^2*x+u+1, A1, B1, [u], 2, "SNC09 Ex.1  gamma=1"];
    E["snc-ex3"] := [x^2-100*u^2*x+u+1, A1, B1, [u], 2, "SNC09 Ex.3  gamma=1"];
    E["snc-ex4"] := [x^2-100*(u^2+1)*x+u+1, A1, B1, [u], 2, "SNC09 Ex.4  gamma=0.01"];
    E["snc-ex5"] := [x^2-100*(u^2+1)*x+u+1, A1, x^3-(u^3+u+2)*x+100, [u], 2,
                     "SNC09 Ex.5  gamma=0.01"];
    E["svd-ex1"] := [C26, x^3+(t2^2+t1+t2-2)*x^2-1, x^3+(2*t2^2-t1+3)*x^2-1,
                     [t1, t2], 3, "SVD26 Ex.1  l=2"];
    E["new-l2"]  := [x^2+(1+t1-t2)*x+(2+t1*t2),
                     x^3+(t1+2*t2)*x^2-x+(1-t1^2),
                     x^3-(2*t1-t2)*x^2+(3+t2^2)*x+1, [t1, t2], 2, "new  l=2"];
    E["new-l3"]  := [x^2+(1+t1)*x+(2+t2-t3),
                     x^3+(t1+t2)*x^2-x+(1+t3),
                     x^3-(t2-t3)*x^2+(3+t1)*x+2, [t1, t2, t3], 2, "new  l=3"];
    return E;
  end proc;

end module:
