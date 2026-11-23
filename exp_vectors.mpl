# exp_vectors.mpl -- 論文本文の例 1, 2 に載せる途中経過（v_K, delta z, delta r）を出力する
read "common.mpl":
tv := [t1, t2]:
Show := proc(name, v)
  local i;
  printf("%s:\n", name);
  for i to LinearAlgebra:-Dimension(v) do
    printf("  [%2d] %a\n", i, evalf[10](expand(v[i])));
  end do;
end proc:
Steps := proc(F, G, k, js, wmax)
  local S, N, K, H, U, sig, V, vK, R, w, p, j, ph, z, q, a, i, res;
  S := SVDLift:-SubresMat(F, G, k, x);
  N := LinearAlgebra:-RowDimension(S); K := LinearAlgebra:-ColumnDimension(S);
  H := SVDLift:-HomMats(S, tv, wmax);
  U, sig, V := SVDLift:-InitSVD(H[0]);
  printf("singular values: %a\n", convert(map(evalf[4], sig), list));
  vK := V[1 .. K, K];
  Show("v_K", vK);
  R := table(); R[0] := vK;
  for w from 1 to wmax do
    p := Vector(N, 0);
    for j from 1 to w do p := p - LinearAlgebra:-Multiply(H[j], R[w - j]); end do;
    p := map(expand, p);
    ph := Vector(K - 1, i -> expand(add(U[a, i]*p[a], a = 1 .. N)/sig[i]));
    z := Vector(K, i -> expand(add(V[i, j]*ph[j], j = 1 .. K - 1)));
    Show(sprintf("delta z^(%d)  (delta q = 0)", w), z);
    q := -z[js]/vK[js];
    printf("delta q^(%d) = %a\n", w, evalf[10](expand(q)));
    R[w] := map(expand, z + q*vK);
    Show(sprintf("delta r^(%d)", w), R[w]);
    # 残差 S r^(w) mod <t>^(w+1)
    res := map(expand, add(add(LinearAlgebra:-Multiply(H[j], R[i]), i = 0 .. w - j), j = 0 .. w));
    printf("residual (degree <= %d) = %.2e\n", w, SVDLift:-PolyNorm(add(res[i]^2, i = 1 .. N), tv)^(1/2));
  end do;
end proc:
printf("===== Example 1 =====\n");
Steps(expand(C26*A26), expand(C26*B26), 3, 1, 1);
printf("===== Example 2 =====\n");
F2, G2, At, Bt := Ex2(1.0003, 1.0005):
Steps(F2, G2, 3, 1, 2);
printf("true cofactor coefficients (Ex.2, x-degree 4 cofactor of G, deg-1 part):\n");
print(SVDLift:-HomPart(Bt, tv, 1), SVDLift:-HomPart(At, tv, 1));
quit;
