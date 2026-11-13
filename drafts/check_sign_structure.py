import math, random
random.seed(7)

def mm(A,B):
    n=len(A); m=len(B[0]); k=len(B)
    return [[sum(A[i][t]*B[t][j] for t in range(k)) for j in range(m)] for i in range(n)]
def T(A): return [list(r) for r in zip(*A)]
def eye(n): return [[1.0 if i==j else 0.0 for j in range(n)] for i in range(n)]

def jacobi(A, iters=200):
    """eigenvalues of real symmetric matrix (list of lists)."""
    n=len(A); a=[r[:] for r in A]
    for _ in range(iters):
        off=sum(a[i][j]**2 for i in range(n) for j in range(n) if i!=j)
        if off < 1e-26: break
        for p in range(n-1):
            for q in range(p+1,n):
                if abs(a[p][q])<1e-18: continue
                th=(a[q][q]-a[p][p])/(2*a[p][q])
                t=(1 if th>=0 else -1)/(abs(th)+math.sqrt(th*th+1))
                c=1/math.sqrt(t*t+1); s=t*c
                for k in range(n):
                    akp=a[k][p]; akq=a[k][q]
                    a[k][p]=c*akp-s*akq; a[k][q]=s*akp+c*akq
                for k in range(n):
                    apk=a[p][k]; aqk=a[q][k]
                    a[p][k]=c*apk-s*aqk; a[q][k]=s*apk+c*aqk
    return sorted(a[i][i] for i in range(n))

def specnorm(A):
    return math.sqrt(max(jacobi(mm(T(A),A))))

def rho_nonneg(A):
    """Perron root of entrywise nonnegative matrix via power iteration."""
    n=len(A); x=[1.0]*n
    r=0.0
    for _ in range(2000):
        y=[sum(A[i][j]*x[j] for j in range(n)) for i in range(n)]
        nrm=sum(y); x=[v/nrm for v in y]; r=nrm
    return r

def infnorm(A): return max(sum(abs(v) for v in row) for row in A)
def diagsim2(s,A):
    n=len(A); return [[s[i]*A[i][j]/s[j] for j in range(n)] for i in range(n)]

def rand_sym(n, lo=-1, hi=1):
    A=[[random.uniform(lo,hi) for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(i): A[i][j]=A[j][i]
    return A

def rand_nonneg(n):
    return [[abs(random.gauss(0,1)) for _ in range(n)] for _ in range(n)]

print("=== (1) symmetric O: random search over positive diagonal S cannot beat ||O||_2")
O=[[2,3,-1],[3,1,4],[-1,4,0.5]]
nO=specnorm(O)
best=1e9
for _ in range(4000):
    s=[math.exp(random.gauss(0,2)) for _ in range(3)]
    best=min(best,specnorm(diagsim2(s,O)))
print("  eig:",[round(v,4) for v in jacobi(O)],"||O||_2=",round(nO,6),"min_searched=",round(best,6),
      ">= ||O||_2:", best>=nO-1e-9)

print("=== (2) non-symmetric: diagonal similarity collapses the 2-norm")
A=[[1.0,1000.0],[0.0,1.0]]
print("  ||A||_2=",round(specnorm(A),3),"   best of ||S A S^-1||_2 over s=t grid:",
      round(min(specnorm(diagsim2([1.0,t],A)) for t in [10**(-6+12*k/400) for k in range(401)]),4),
      " rho(A)=",1.0)
B=[[0,1000.0,0],[0,0,1000.0],[0.001,0,0]]
print("  cyclic 3x3: ||B||_2=",round(specnorm(B),3)," rho=",round(rho_nonneg([[abs(v) for v in r] for r in B]),3),
      " min over diag S of ||.||_inf:",round(min(infnorm(diagsim2([1.0,t,t*t],B)) for t in [10**(-3+6*k/600) for k in range(601)]),3))

print("=== (3) A>=0: min over positive diagonal S of ||S A S^-1||_inf == rho(A)")
for trial in range(3):
    A=rand_nonneg(5); r=rho_nonneg(A); best=1e9
    for _ in range(1500):
        s=[math.exp(random.gauss(0,2.5)) for _ in range(5)]
        best=min(best,infnorm(diagsim2(s,A)))
    print("   rho=",round(r,6),"  min_inf_searched=",round(best,6),"  rowsum(I)=",round(infnorm(A),4))

print("=== (4) ||O||_2 <= rho(|O|) <= max row sum; strictness vs balance")
O2=[[0,1,-1],[1,0,1],[-1,1,0]]   # unbalanced triangle
print("  unbalanced:  ||O||_2=",round(specnorm(O2),6)," rho(|O|)=",round(rho_nonneg([[abs(v) for v in r] for r in O2]),6),
      " maxrowsum=",infnorm(O2))
Dg=[[1,0,0],[0,1,0],[0,0,-1]]
Ob=mm(mm(Dg,O2),Dg)
print("  switched(balanced): ||Ob||_2=",round(specnorm(Ob),6)," rho(|O|)=",round(rho_nonneg([[abs(v) for v in r] for r in Ob]),6))
Z=[[-1.0,-2.0,-0.5],[-2.0,-1.5,-3.0],[-0.5,-3.0,-2.2]]  # Z-matrix
print("  Z-matrix:    ||Z||_2=",round(specnorm(Z),6)," rho(|Z|)=",round(rho_nonneg([[abs(v) for v in r] for r in Z]),6),
      " maxrowsum=",infnorm(Z))

print("=== (5) pointwise |x^T O x| <= |x|^T |O| |x|")
Ok=rand_sym(6); worst=-1e9
for _ in range(3000):
    x=[random.gauss(0,1) for _ in range(6)]
    lhs=abs(sum(Ok[i][j]*x[i]*x[j] for i in range(6) for j in range(6)))
    rhs=sum(abs(Ok[i][j])*abs(x[i])*abs(x[j]) for i in range(6) for j in range(6))
    worst=max(worst,lhs-rhs)
print("  max(lhs-rhs) =",worst)

print("=== (6) weighted Collatz-Wielandt: tau_v = max_i (|O|v)_i/v_i >= ||O||_2")
v=[math.exp(random.gauss(0,1)) for _ in range(6)]
AbsO=[[abs(w) for w in r] for r in Ok]
tv=max(sum(AbsO[i][j]*v[j] for j in range(6))/v[i] for i in range(6))
print("  tau_v=",round(tv,6),"  ||O||_2=",round(specnorm(Ok),6),"  rho(|O|)=",round(rho_nonneg(AbsO),6))
