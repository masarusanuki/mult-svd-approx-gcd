"""Prototype of the SVD + lifting approximate GCD (numpy only).

Polynomials in x, t_1..t_l are dicts {(ex, (e1,..,el)): coeff}.
Vectors of homogeneous t-polynomials of degree w are arrays (len, M_w),
columns indexed by the monomials of total degree w.
"""
import itertools, math, time
import numpy as np

# ---------------------------------------------------------------- monomials
class Mon:
    def __init__(self, l, wmax):
        self.l = l
        self.mons = []
        self.idx = []
        for w in range(wmax + 1):
            ms = [e for e in itertools.product(range(w + 1), repeat=l) if sum(e) == w]
            ms.sort(reverse=True)
            self.mons.append(ms)
            self.idx.append({e: i for i, e in enumerate(ms)})
        self._mul = {}

    def M(self, w):
        return len(self.mons[w])

    def multab(self, a, b):
        key = (a, b)
        if key not in self._mul:
            T = np.empty((self.M(a), self.M(b)), dtype=int)
            for i, ea in enumerate(self.mons[a]):
                for j, eb in enumerate(self.mons[b]):
                    T[i, j] = self.idx[a + b][tuple(p + q for p, q in zip(ea, eb))]
            self._mul[key] = T
        return self._mul[key]


def hmul_mat_vec(Mo, A, a, z, b):
    """A: (N,K,M_a) homogeneous deg a matrix, z: (K,M_b) -> (N, M_{a+b})."""
    T = Mo.multab(a, b)
    N = A.shape[0]
    out = np.zeros((N, Mo.M(a + b)))
    # out[:, T[i,j]] += A[:,:,i] @ z[:,j]
    P = np.einsum('nki,kj->nij', A, z)
    np.add.at(out.T, T.ravel(), P.reshape(N, -1).T)
    return out


def hmul_scal_vec(Mo, c, a, z, b):
    """c: (M_a,) homogeneous scalar of deg a, z: (K,M_b) -> (K, M_{a+b})."""
    T = Mo.multab(a, b)
    K = z.shape[0]
    out = np.zeros((K, Mo.M(a + b)))
    P = c[:, None, None] * z.T[None, :, :]  # (M_a, M_b, K)
    np.add.at(out.T, T.ravel(), P.reshape(-1, K))
    return out

# ---------------------------------------------------------------- polynomials
def pmul(P, Q):
    R = {}
    for (ea, ta), ca in P.items():
        for (eb, tb), cb in Q.items():
            k = (ea + eb, tuple(p + q for p, q in zip(ta, tb)))
            R[k] = R.get(k, 0.0) + ca * cb
    return R

def padd(P, Q, s=1.0):
    R = dict(P)
    for k, v in Q.items():
        R[k] = R.get(k, 0.0) + s * v
    return R

def pnorm(P):
    return math.sqrt(sum(v * v for v in P.values()))

def xdeg(P):
    return max(e for (e, _), v in P.items() if v != 0)

def tdeg(P):
    return max(sum(t) for (_, t), v in P.items() if v != 0)

def coeff_x(P, i, l):
    return {t: v for (e, t), v in P.items() if e == i}

# ---------------------------------------------------------------- subresultant
def sylv_hom(F, G, k, Mo, W):
    """Homogeneous parts dS^(w) (N,K,M_w), w=0..W, of S_k(F,G).
    Columns: u (deg n-k) high->low, then v (deg m-k). r=(Gt,-Ft) is null."""
    l = Mo.l
    m, n = xdeg(F), xdeg(G)
    nu, nv = n - k + 1, m - k + 1
    N, K = m + n - k + 1, nu + nv
    S = [np.zeros((N, K, Mo.M(w))) for w in range(W + 1)]
    # column j of u-block: x^(n-k-j) * F ; row index 0 = highest power m+n-k
    for (e, t), c in F.items():
        w = sum(t)
        if w > W: continue
        mi = Mo.idx[w][t]
        for j in range(nu):
            p = e + (n - k - j)
            S[w][(m + n - k) - p, j, mi] += c
    for (e, t), c in G.items():
        w = sum(t)
        if w > W: continue
        mi = Mo.idx[w][t]
        for j in range(nv):
            p = e + (m - k - j)
            S[w][(m + n - k) - p, nu + j, mi] += c
    return S, N, K, nu, nv


def vec_to_cof(R, Mo, nu, nv, n, m, k):
    """R: list of (K,M_w) -> (Ft, Gt) dicts."""
    Ft, Gt = {}, {}
    for w, Z in enumerate(R):
        for mi, t in enumerate(Mo.mons[w]):
            for j in range(nu):
                if Z[j, mi] != 0: Gt[(n - k - j, t)] = Gt.get((n - k - j, t), 0) + Z[j, mi]
            for j in range(nv):
                if Z[nu + j, mi] != 0: Ft[(m - k - j, t)] = Ft.get((m - k - j, t), 0) - Z[nu + j, mi]
    return Ft, Gt

# ---------------------------------------------------------------- degree detection
def detect_degree(F, G, l, tol):
    """largest d with sigma_min(S_d(F0,G0)) / ||S_d|| < tol."""
    Mo = Mon(l, 0)
    F0 = {k: v for k, v in F.items() if sum(k[1]) == 0}
    G0 = {k: v for k, v in G.items() if sum(k[1]) == 0}
    m, n = xdeg(F0), xdeg(G0)
    for d in range(min(m, n), 0, -1):
        S = sylv_hom(F0, G0, d, Mo, 0)[0][0][:, :, 0]
        s = np.linalg.svd(S, compute_uv=False)
        if s[-1] / s[0] < tol:
            return d, s
    return 0, None

# ---------------------------------------------------------------- lifting
def svd_lift(F, G, k, l, W, jstar=None, mode='standard', npass=3, verbose=False):
    """Lift null vector of S_k up to degree W.
    mode='standard': free direction v_K only, normalisation r[jstar] const.
    mode='ill':      free directions v_{K-1}, v_K; alpha determined one step late.
    returns list R[w] (K,M_w), info."""
    Mo = Mon(l, W)
    S, N, K, nu, nv = sylv_hom(F, G, k, Mo, W)
    U, sig, Vt = np.linalg.svd(S[0][:, :, 0], full_matrices=False)
    V = Vt.T
    vK = V[:, -1]
    if jstar is None:
        jstar = int(np.argmax(np.abs(vK)))
    wS = max(w for w in range(W + 1) if np.any(S[w] != 0))

    def rhs(R, w):
        p = np.zeros((N, Mo.M(w)))
        for j in range(1, min(w, wS) + 1):
            p -= hmul_mat_vec(Mo, S[j], j, R[w - j], w - j)
        return p

    if mode == 'standard':
        R = [vK[:, None] / vK[jstar]]
        for w in range(1, W + 1):
            p = rhs(R, w)
            ph = U.T @ p                       # (K, M_w)
            z = V[:, :-1] @ (ph[:-1] / sig[:-1, None])
            beta = -z[jstar] / vK[jstar]
            R.append(z + vK[:, None] * beta[None, :])
        return R, dict(sigma=sig, jstar=jstar, nu=nu, nv=nv, N=N, K=K)

    # ---- ill-conditioned: two-dimensional free space span{v_{K-1}, v_K}
    # alpha^(w-1) (coefficient of v_{K-1}) is fixed one step late from the
    # u_{K-1}-equation  sigma_{K-1} alpha^(w) = u_{K-1}^T p^(w),
    # first with sigma_{K-1} alpha^(w) := 0, then corrected in later passes.
    Mo = Mon(l, W + 1)
    S = sylv_hom(F, G, k, Mo, W + 1)[0]
    wS = max(w for w in range(W + 2) if np.any(S[w] != 0))
    def rhs(R, w):
        p = np.zeros((N, Mo.M(w)))
        for j in range(1, min(w, wS) + 1):
            p -= hmul_mat_vec(Mo, S[j], j, R[w - j], w - j)
        return p
    v1, u1, s1 = V[:, -2], U[:, -2], sig[-2]
    vt = v1 - (v1[jstar] / vK[jstar]) * vK      # keeps r[jstar] unchanged
    alpha_prev = None
    for ps in range(npass):
        R = [None] * (W + 1)
        alphas = [None] * (W + 1)
        D = (vK / vK[jstar])[:, None]           # particular part of degree w-1
        for w in range(1, W + 2):
            Mw1 = Mo.M(w - 1)
            R[w - 1] = D
            g0 = u1 @ rhs(R, w)
            Gm = np.zeros((Mo.M(w), Mw1))
            for mi in range(Mw1):
                E = np.zeros((K, Mw1)); E[:, mi] = vt
                Gm[:, mi] = -u1 @ hmul_mat_vec(Mo, S[1], 1, E, w - 1)
            target = -g0
            if alpha_prev is not None and w <= W:
                target = target + s1 * alpha_prev[w]
            a, *_ = np.linalg.lstsq(Gm, target, rcond=None)
            alphas[w - 1] = a
            R[w - 1] = D + vt[:, None] * a[None, :]
            if w <= W:
                ph = U.T @ rhs(R, w)
                z = V[:, :-2] @ (ph[:-2] / sig[:-2, None])
                beta = -z[jstar] / vK[jstar]
                D = z + vK[:, None] * beta[None, :]
        alpha_prev = alphas
    return R, dict(sigma=sig, jstar=jstar, nu=nu, nv=nv, N=N, K=K)


def naive_ill(F, G, k, l, W, jstar=None):
    return svd_lift(F, G, k, l, W, jstar, mode='standard')

# ---------------------------------------------------------------- GCD from cofactors
def conv_mat(P, d, l, Mo, W):
    """homogeneous parts of the convolution matrix of P (x-deg dp) acting on deg-d poly."""
    dp = xdeg(P)
    rows, cols = dp + d + 1, d + 1
    A = [np.zeros((rows, cols, Mo.M(w))) for w in range(W + 1)]
    for (e, t), c in P.items():
        w = sum(t)
        if w > W: continue
        for j in range(cols):
            A[w][(dp + d) - (e + d - j), j, Mo.idx[w][t]] += c
    return A

def gcd_from_cof(F, G, Ft, Gt, k, l, W):
    """solve [Conv(Ft);Conv(Gt)] c = [F;G] by lifting with SVD of the deg-0 part."""
    Mo = Mon(l, W + max(tdeg(F), tdeg(G)))
    WW = W
    A1 = conv_mat(Ft, k, l, Mo, WW); A2 = conv_mat(Gt, k, l, Mo, WW)
    A = [np.vstack([a1, a2]) for a1, a2 in zip(A1, A2)]
    m, n = xdeg(F), xdeg(G)
    rows = A[0].shape[0]
    def bvec(P, deg, w):
        b = np.zeros((deg + 1, Mo.M(w)))
        for (e, t), c in P.items():
            if sum(t) == w: b[deg - e, Mo.idx[w][t]] += c
        return b
    wA = max(w for w in range(WW + 1) if np.any(A[w] != 0))
    Ap = np.linalg.pinv(A[0][:, :, 0])
    Cw = []
    for w in range(WW + 1):
        b = np.vstack([bvec(F, m, w), bvec(G, n, w)])
        for j in range(1, min(w, wA) + 1):
            b -= hmul_mat_vec(Mo, A[j], j, Cw[w - j], w - j)
        Cw.append(Ap @ b)
    C = {}
    for w, Z in enumerate(Cw):
        for mi, t in enumerate(Mo.mons[w]):
            for j in range(k + 1):
                if Z[j, mi] != 0: C[(k - j, t)] = Z[j, mi]
    return C

# ---------------------------------------------------------------- error measures
def vec_err(Pc, Pt):
    """distance after optimal scalar: min_s ||s Pc - Pt|| / ||Pt||."""
    keys = set(Pc) | set(Pt)
    a = np.array([Pc.get(k, 0.0) for k in keys]); b = np.array([Pt.get(k, 0.0) for k in keys])
    s = (a @ b) / (a @ a)
    return np.linalg.norm(s * a - b) / np.linalg.norm(b)

def chop(P, tol=1e-300):
    return {k: v for k, v in P.items() if abs(v) > tol}

def trunc(P, W):
    return {k: v for k, v in P.items() if sum(k[1]) <= W}

def backward(F, G, C, Ft, Gt):
    return max(pnorm(padd(F, pmul(C, Ft), -1)) / pnorm(F), pnorm(padd(G, pmul(C, Gt), -1)) / pnorm(G))

# ---------------------------------------------------------------- direct (big matrix) method
def direct_nullvec(F, G, k, l, Wc):
    """null vector of the full numerical multivariate subresultant:
    unknown Gt (x-deg n-k), Ft (x-deg m-k) with t-tdeg <= Wc each."""
    m, n = xdeg(F), xdeg(G)
    WF = max(tdeg(F), tdeg(G))
    Mo = Mon(l, Wc + WF)
    cmon = [(w, t) for w in range(Wc + 1) for t in Mo.mons[w]]
    rmon = {}
    for w in range(Wc + WF + 1):
        for t in Mo.mons[w]:
            for e in range(m + n - k + 1):
                rmon[(e, t)] = len(rmon)
    cols = []
    for deg, P, sgn in ((n - k, F, 1.0), (m - k, G, 1.0)):
        for e in range(deg, -1, -1):
            for (w, t) in cmon:
                col = {}
                for (ep, tp), c in P.items():
                    key = (ep + e, tuple(a + b for a, b in zip(tp, t)))
                    col[rmon[key]] = col.get(rmon[key], 0) + sgn * c
                cols.append(col)
    A = np.zeros((len(rmon), len(cols)))
    for j, col in enumerate(cols):
        for i, c in col.items(): A[i, j] = c
    t0 = time.time()
    _, s, Vt = np.linalg.svd(A, full_matrices=False)
    tsvd = time.time() - t0
    v = Vt[-1]
    Gt, Ft = {}, {}
    j = 0
    nc = len(cmon)
    for e in range(n - k, -1, -1):
        for (w, t) in cmon:
            Gt[(e, t)] = v[j]; j += 1
    for e in range(m - k, -1, -1):
        for (w, t) in cmon:
            Ft[(e, t)] = -v[j]; j += 1
    return Ft, Gt, A.shape, tsvd, s

# ---------------------------------------------------------------- random generation
def rand_tpoly(rng, l, d, dense=True, nterms=3, const=None):
    ms = [t for w in range(d + 1) for t in itertools.product(range(d + 1), repeat=l) if sum(t) == w]
    ms = sorted(set(ms))
    if dense:
        P = {t: rng.uniform(-1, 1) for t in ms}
    else:
        ch = rng.choice(len(ms), size=min(nterms, len(ms)), replace=False)
        P = {ms[i]: rng.uniform(-1, 1) for i in ch}
        P[(0,) * l] = rng.uniform(0.5, 1) * rng.choice([-1, 1])
    if const is not None:
        P = {(0,) * l: const}
    return P

def rand_poly(rng, l, dx, dt, dense=True, monic=True, nterms=3):
    P = {}
    for e in range(dx + 1):
        if e == dx and monic:
            c = {(0,) * l: 1.0}
        else:
            c = rand_tpoly(rng, l, dt, dense, nterms)
        for t, v in c.items():
            P[(e, t)] = v
    return P

def add_noise(rng, P, eps):
    if eps == 0: return dict(P)
    E = {k: rng.uniform(-1, 1) for k in P}
    s = eps * pnorm(P) / pnorm(E)
    return {k: P[k] + s * E[k] for k in P}


# ---------------------------------------------------------------- dense case
def multiplier(R, Mo, D):
    """homogeneous system for c (tdeg<=D): parts of degree D+1..2D of c*R vanish.
    returns (coeff matrix, monomial list)."""
    K = R[0].shape[0]
    cmons = [(d, mi) for d in range(D + 1) for mi in range(Mo.M(d))]
    blocks = []
    for w in range(D + 1, 2 * D + 2):
        Bw = np.zeros((K * Mo.M(w), len(cmons)))
        for j, (d, mi) in enumerate(cmons):
            if w - d > len(R) - 1: continue
            e = np.zeros(Mo.M(d)); e[mi] = 1.0
            Bw[:, j] = hmul_scal_vec(Mo, e, d, R[w - d], w - d).ravel()
        blocks.append(Bw)
    return (np.vstack(blocks) if blocks else np.zeros((0, len(cmons)))), cmons

def dense_lift(F, G, k, l, Dmax, tol=1e-8):
    """normalise at the largest entry of v_K, lift as power series, then find the
    polynomial multiplier c(t) of minimal degree D (detected by SVD)."""
    Mo = Mon(l, 2 * Dmax + 1)
    R, info = svd_lift(F, G, k, l, 2 * Dmax + 1)
    for D in range(0, Dmax + 1):
        A, cmons = multiplier(R, Mo, D)
        if A.shape[0] == 0:
            c = np.array([1.0]); s = [1.0, 0.0]
        else:
            _, s, Vt = np.linalg.svd(A, full_matrices=False)
            c = Vt[-1]
        if A.shape[0] == 0 or s[-1] / s[0] < tol:
            break
    # r = trunc(c * R, D)
    out = []
    for w in range(D + 1):
        z = np.zeros((R[0].shape[0], Mo.M(w)))
        for j, (d, mi) in enumerate(cmons):
            if d > w: continue
            e = np.zeros(Mo.M(d)); e[mi] = c[j]
            z += hmul_scal_vec(Mo, e, d, R[w - d], w - d)
        out.append(z)
    info['D'] = D; info['smin'] = s[-1] / s[0]
    return out, info


# ---------------------------------------------------------------- expansion point
def shift(P, s):
    """P(x, t + s)."""
    R = {}
    for (e, t), c in P.items():
        # expand prod (t_i + s_i)^{a_i}
        terms = [((), c)]
        for a, si in zip(t, s):
            new = []
            for (pre, cc) in terms:
                for b in range(a + 1):
                    new.append((pre + (b,), cc * math.comb(a, b) * si ** (a - b)))
            terms = new
        for (tt, cc) in terms:
            R[(e, tt)] = R.get((e, tt), 0.0) + cc
    return R

def sigma_ratio(F, G, k, l):
    Mo = Mon(l, 0)
    S = sylv_hom(F, G, k, Mo, 0)[0][0][:, :, 0]
    s = np.linalg.svd(S, compute_uv=False)
    return s[-2] / s[0]

def choose_shift(F, G, k, l, rng, ncand=10, radius=1.0):
    best = (sigma_ratio(F, G, k, l), (0.0,) * l)
    for _ in range(ncand):
        s = tuple(np.round(rng.uniform(-radius, radius, l), 2))
        r = sigma_ratio(shift(F, s), shift(G, s), k, l)
        if r > best[0]: best = (r, s)
    return best
