"""Numerical evidence for the convexity claims in notes/tau_sign_structure.md §3:

  (i)  s |-> ||diag(s) O diag(s)^-1||_2 is NOT convex in s (midpoint tests fail), so the
       "not convex" half of the question is correct; what survives is quasi-convexity of
       the LMI formulation in x = s^2 (Proposition 3.2(a));
  (ii) the diagonal-scaling optimum for the 2-norm of a SYMMETRIC O is ||O||_2 (no gain,
       Proposition 2.1), while for the infinity norm of a NONNEGATIVE A it is rho(A)
       (Theorem 3.4).
"""
import math, random

def jacobi(A, iters=400):
    n = len(A); a = [r[:] for r in A]
    for _ in range(iters):
        if sum(a[i][j] ** 2 for i in range(n) for j in range(n) if i != j) < 1e-30:
            break
        for p in range(n - 1):
            for q in range(p + 1, n):
                if abs(a[p][q]) < 1e-18:
                    continue
                th = (a[q][q] - a[p][p]) / (2 * a[p][q])
                t = (1 if th >= 0 else -1) / (abs(th) + math.sqrt(th * th + 1))
                c = 1 / math.sqrt(t * t + 1); s = t * c
                for k in range(n):
                    akp, akq = a[k][p], a[k][q]
                    a[k][p] = c * akp - s * akq; a[k][q] = s * akp + c * akq
                for k in range(n):
                    apk, aqk = a[p][k], a[q][k]
                    a[p][k] = c * apk - s * aqk; a[q][k] = s * apk + c * aqk
    return sorted(a[i][i] for i in range(n))

def specnorm(A):
    n = len(A)
    G = [[sum(A[k][i] * A[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
    return math.sqrt(max(jacobi(G)))

def f(s, O):
    n = len(O)
    return specnorm([[s[i] * O[i][j] / s[j] for j in range(n)] for i in range(n)])

def rand_sym(n):
    O = [[random.gauss(0, 1) for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(i):
            O[i][j] = O[j][i]
    return O

random.seed(4)
viol, worst = 0, 0.0
for _ in range(60):
    O = rand_sym(3)
    for _ in range(300):
        s = [math.exp(random.gauss(0, 1)) for _ in range(3)]
        t = [math.exp(random.gauss(0, 1)) for _ in range(3)]
        mid = [(s[i] + t[i]) / 2 for i in range(3)]
        gap = f(mid, O) - (f(s, O) + f(t, O)) / 2
        if gap > 1e-9:
            viol += 1; worst = max(worst, gap)
print(f"(i)  midpoint-convexity violations in s: {viol}/18000, worst excess {worst:.3f}")

O = [[2, 3, -1], [3, 1, 4], [-1, 4, 0.5]]
base = specnorm(O); lny = [0.0, 0.0, 0.0]; best = base
for it in range(3000):
    step = 0.5 / (1 + it / 200)
    for k in range(3):
        for d in (step, -step):
            cand = lny[:]; cand[k] += d
            v = f([math.exp(z) for z in cand], O)
            if v < best - 1e-15:
                best, lny = v, cand
print(f"(ii) symmetric O: ||O||_2 = {base:.6f}, best after descent = {best:.6f}, gain = {base-best:.2e}")

A = [[0, 1000.0, 0], [0, 0, 1000.0], [0.001, 0, 0]]
def perron(B, it=20000):
    n = len(B); x = [1.0] * n
    for _ in range(it):
        y = [sum(B[i][j] * x[j] for j in range(n)) + x[i] for i in range(n)]
        m = max(y); x = [v / m for v in y]
    return x, sum(sum(B[i][j] * x[j] for j in range(n)) / x[i] for i in range(n)) / n
v, r = perron(A)
def infn(M): return max(sum(abs(z) for z in row) for row in M)
S = [1.0 / z for z in v]
print(f"      nonnegative A: rho = {r:.6f}, ||A||_inf = {infn(A):.3f}, "
      f"||diag(1/v) A diag(1/v)^-1||_inf = {infn([[S[i]*A[i][j]/S[j] for j in range(3)] for i in range(3)]):.6f}")
