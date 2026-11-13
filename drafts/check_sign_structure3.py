import math, random
random.seed(11)
def mm(A,B): return [[sum(A[i][t]*B[t][j] for t in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]
def T(A): return [list(r) for r in zip(*A)]
def jacobi(A,iters=300):
    n=len(A); a=[r[:] for r in A]
    for _ in range(iters):
        if sum(a[i][j]**2 for i in range(n) for j in range(n) if i!=j)<1e-28: break
        for p in range(n-1):
            for q in range(p+1,n):
                if abs(a[p][q])<1e-18: continue
                th=(a[q][q]-a[p][p])/(2*a[p][q])
                t=(1 if th>=0 else -1)/(abs(th)+math.sqrt(th*th+1))
                c=1/math.sqrt(t*t+1); s=t*c
                for k in range(n):
                    akp=a[k][p]; akq=a[k][q]; a[k][p]=c*akp-s*akq; a[k][q]=s*akp+c*akq
                for k in range(n):
                    apk=a[p][k]; aqk=a[q][k]; a[p][k]=c*apk-s*aqk; a[q][k]=s*apk+c*aqk
    return sorted(a[i][i] for i in range(n))
def specnorm(A): return math.sqrt(max(jacobi(mm(T(A),A))))
def rho_sym(A): return max(abs(v) for v in jacobi(A))
def infnorm(A): return max(sum(abs(v) for v in r) for r in A)
def sim2(s,A): return [[s[i]*A[i][j]/s[j] for j in range(len(A))] for i in range(len(A))]      # S A S^-1
def cong(s,A): return [[s[i]*A[i][j]*s[j] for j in range(len(A))] for i in range(len(A))]      # S A S

def perron(A,it=20000):
    n=len(A); x=[1.0]*n
    for _ in range(it):
        y=[sum(A[i][j]*x[j] for j in range(n))+x[i] for i in range(n)]
        m=max(y); x=[v/m for v in y]
    r=sum(sum(A[i][j]*x[j] for j in range(n))/x[i] for i in range(n))/n
    return x,r

print("=== D'. A>=0 : min_s ||diag(s) A diag(s)^-1||_inf = rho(A), witness s = 1/Perron")
for tag,A in [("rand5a",[[abs(random.gauss(1,0.6)) for _ in range(5)] for _ in range(5)]),
              ("cyclic3",[[0,1000.,0],[0,0,1000.],[.001,0,0]])]:
    v,r=perron(A); s=[1.0/t for t in v]
    print(f"   {tag}: rho={r:.6f}  ||.||_inf(s=1/v)={infnorm(sim2(s,A)):.6f}  ||.||_inf(s=1)={infnorm(A):.4f}"
          f"  ||.||_2={specnorm(A):.4f}")

print("=== G'. congruence can shrink (S O S), similarity cannot")
O=[[0,10.],[10.,0.]]
print("   ||O||_2=",specnorm(O)," ||SOS||_2 with s=(1,0.1):",specnorm(cong([1,0.1],O)),
      " ||SOS^-1||_2:",specnorm(sim2([1,0.1],O)))

print("=== B'. equality criterion iff: zero gap <=> some eigenvector at +-rho has no row cancellation")
def perron_vec_of(A,lam):
    # solve (A-lam I) w = 0 via inverse iteration-ish: use (A - lam I + tiny shift) power iteration on inverse
    n=len(A); M=[[A[i][j]-(lam if i==j else 0) + (1e-9 if i==j else 0) for j in range(n)] for i in range(n)]
    # gaussian invert
    aug=[M[i][:]+[1.0 if i==k else 0.0 for k in range(n)] for i in range(n)]
    for c in range(n):
        p=max(range(c,n),key=lambda r:abs(aug[r][c])); aug[c],aug[p]=aug[p],aug[c]
        pv=aug[c][c]
        if abs(pv)<1e-14: return None
        aug[c]=[v/pv for v in aug[c]]
        for r in range(n):
            if r!=c and aug[r][c]!=0:
                f=aug[r][c]; aug[r]=[aug[r][k]-f*aug[c][k] for k in range(2*n)]
    Minv=[row[n:] for row in aug]
    x=[1.0]*n
    for _ in range(300):
        y=[sum(Minv[i][j]*x[j] for j in range(n)) for i in range(n)]
        nm=math.sqrt(sum(v*v for v in y)); x=[v/nm for v in y]
    return x
def is_balanced(O):
    n=len(O); sig=[None]*n
    for st in range(n):
        if sig[st] is not None: continue
        sig[st]=1; stk=[st]
        while stk:
            i=stk.pop()
            for j in range(n):
                if i!=j and O[i][j]!=0:
                    want=sig[i]*(1 if O[i][j]>0 else -1)
                    if sig[j] is None: sig[j]=want; stk.append(j)
                    elif sig[j]!=want: return False
    return True
bad=0; tested=0
for t in range(120):
    n=4; A=[[0.0]*n for _ in range(n)]; Mag=[[0.0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1,n):
            m=abs(random.gauss(1,0.3)); sg=random.choice([1,-1])
            A[i][j]=A[j][i]=m*sg; Mag[i][j]=Mag[j][i]=m
    if min(sum(1 for j in range(n) if Mag[i][j]>0) for i in range(n))<2: continue
    tested+=1
    gap=rho_sym(Mag)-rho_sym(A)
    found=False
    for lam in (rho_sym(A),-rho_sym(A)):
        x=perron_vec_of(A,lam)
        if x is None: continue
        if all(abs(abs(sum(A[i][j]*x[j] for j in range(n)))-sum(Mag[i][j]*abs(x[j]) for j in range(n)))<1e-6 for i in range(n)):
            found=True
    if (gap<1e-7)!=found: bad+=1
print(f"   tested {tested} signed 4x4: criterion mismatches = {bad}")
print("   balanced=>zero gap check:", all((rho_sym(Mag)-rho_sym(A))<1e-9 for _ in [] ) )

print("=== F'. clean examples")
C4=[[0,1,0,-1],[1,0,1,0],[0,1,0,1],[-1,0,1,0]]
print("   signed C4 (one negative edge): spectrum",[round(v,4) for v in jacobi(C4)],
      " rho=",round(rho_sym(C4),6)," rho(|O|)=",round(rho_sym([[abs(w) for w in r] for r in C4]),6))
print("=== E'. weighted CW: tau_v >= rho(|O|) with equality iff v is the Perron vector")
Ok=[[2,3,-1],[3,1,4],[-1,4,0.5]]; Mag=[[abs(w) for w in r] for r in Ok]; v,r=perron(Mag)
print("   rho(|O|)=",round(r,6)," tau_Perron=",round(max(sum(Mag[i][j]*v[j] for j in range(3))/v[i] for i in range(3)),6),
      " tau_allones=",round(infnorm(Mag),4)," ||O||_2=",round(specnorm(Ok),6))
