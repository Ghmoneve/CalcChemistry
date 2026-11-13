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
def ds2(s,A): return [[s[i]*A[i][j]/s[j] for j in range(len(A))] for i in range(len(A))]
def ds_inf(s,A): return [[s[i]*A[i][j]/s[j] for j in range(len(A))] for i in range(len(A))]

def perron(A,it=5000):           # A >= 0 irreducible: Perron vector via (A+I) shift
    n=len(A); x=[1.0]*n; it_use=it
    for _ in range(it_use):
        y=[sum(A[i][j]*x[j] for j in range(n))+x[i] for i in range(n)]
        m=max(y); x=[v/m for v in y]
    r=sum(sum(A[i][j]*x[j] for j in range(n))/x[i] for i in range(n))/n
    return x,r

def is_balanced(O):
    n=len(O); sgn=[[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i!=j and O[i][j]!=0: sgn[i][j]=1 if O[i][j]>0 else -1
    sig=[None]*n
    for st in range(n):
        if sig[st] is not None: continue
        sig[st]=1; stk=[st]
        while stk:
            i=stk.pop()
            for j in range(n):
                if sgn[i][j]!=0:
                    want=sig[i]*sgn[i][j]
                    if sig[j] is None: sig[j]=want; stk.append(j)
                    elif sig[j]!=want: return False,(i,j)
    return True,None

print("=== A. signature balance vs spectral gap  rho(|O|)-rho(O)")
print("   (unbalanced triangle, spectrum {1,1,-2} as computed by hand)")
O2=[[0,1,-1],[1,0,1],[-1,1,0]]
print("   eigen:",[round(v,6) for v in jacobi(O2)],"rho=",rho_sym(O2),"rho|O|=",rho_sym([[abs(v) for v in r] for r in O2]),
      "balanced:",is_balanced(O2)[0])
print("   -> strict gap NOT implied by imbalance; test random signed patterns:")
found_gap=0; found_nogap_unbal=0
for t in range(400):
    n=4
    A=[[0.0]*n for _ in range(n)]
    mag=[[0.0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1,n):
            m=abs(random.gauss(1,0.3)); s=random.choice([1,-1])
            A[i][j]=A[j][i]=m*s; mag[i][j]=mag[j][i]=m
    if min(sum(1 for j in range(n) if mag[i][j]>0) for i in range(n))==0: continue
    bal,_=is_balanced(A)
    g=rho_sym(mag)-rho_sym(A)
    if g>1e-9 and bal: print("   COUNTEREXAMPLE to 'balanced => no gap'!")
    if g>1e-9: found_gap+=1
    if g<=1e-9 and not bal: found_nogap_unbal+=1
    if found_gap==1 and g>1e-9:
        print("   example strict gap: rho(O)=",round(rho_sym(A),5),"rho(|O|)=",round(rho_sym(mag),5),"balanced=",bal)
        found_gap=2
print("   #strict-gap samples:",found_gap,"  #imbalanced-but-zero-gap samples:",found_nogap_unbal)

print("=== B. equality criterion: no row cancellation at a Perron vector")
A=O2; lam=-2.0; x=[1.0,-1.0,1.0]
Ax=[sum(A[i][j]*x[j] for j in range(3)) for i in range(3)]
print("   Ax =",Ax," (= lam*x?)",all(abs(Ax[i]-lam*x[i])<1e-12 for i in range(3)))
print("   |Ax| == |A||x| (no cancellation):",all(abs(abs(Ax[i])-sum(abs(A[i][j])*abs(x[j]) for j in range(3)))<1e-12 for i in range(3)))

print("=== C. Prop2 by descent: min_s ||diag(s) O diag(s)^-1||_2 for symmetric O stays at ||O||_2")
O=[[2,3,-1],[3,1,4],[-1,4,0.5]]; base=specnorm(O)
lny=[0.0,0.0,0.0]; best=base
for it in range(3000):
    step=0.5/(1+it/200)
    for k in range(3):
        for d in (step,-step):
            cand=lny[:]; cand[k]+=d
            v=specnorm(ds2([math.exp(z) for z in cand],O))
            if v<best-1e-15: best=v; lny=cand
print("   ||O||_2=",round(base,6),"  descend result=",round(best,6),"  gain=",round(base-best,12))

print("=== D. A>=0: ||diag(v) A diag(v)^-1||_inf at Perron vector == rho(A)")
for _ in range(3):
    A=[[abs(random.gauss(1,0.6)) for _ in range(5)] for _ in range(5)]
    v,r=perron(A)
    print("   rho=",round(r,6),"  ||.||_inf at Perron v=",round(infnorm(ds_inf(v,A)),6),"  ||.||_inf at v=1:",round(infnorm(A),4))
for B,name in [([[0,1000.,0],[0,0,1000.],[.001,0,0]],"cyclic")]:
    v,r=perron(B)
    print("   ",name,"rho=",round(r,6)," inf at Perron v=",round(infnorm(ds_inf(v,B)),6)," spectral norm:",round(specnorm(B),2))

print("=== E. weighted CW bound tau_v with v = Perron vector of |O| equals rho(|O|)")
Ok=[[2,3,-1],[3,1,4],[-1,4,0.5]]; Mag=[[abs(w) for w in r] for r in Ok]
v,r=perron(Mag)
tav=max(sum(Mag[i][j]*v[j] for j in range(3))/v[i] for i in range(3))
print("   rho(|O|)=",round(r,6)," tau_v=",round(tav,6)," ||O||_2=",round(specnorm(Ok),6),
      " rowsum=",infnorm(Mag))
# C4 with one negative edge: clean strict-gap exhibit
C4=[[0,1,0,-1],[1,0,1,0],[0,1,0,1],[-1,0,1,0]]
print("=== F. signed C4 with one negative edge: eigen",[round(v,4) for v in jacobi(C4)],
      "rho=",round(rho_sym(C4),6),"rho|O|=",round(rho_sym([[abs(w) for w in r] for r in C4]),6),
      "balanced?",is_balanced(C4)[0])
print("=== G. congruence (not similarity) can shrink: S O S with O=[[0,10],[10,0]], S=diag(1,0.1)")
print("   ||SOS||_2 =",round(specnorm(ds2([1,0.1],[[0,10.],[10.,0.]])),4))
print("=== H. naive 'drop negative entries' bound is FALSE: O=[[0,1],[1,0]], rho(O^-)=0 but lam_min=-1")
print("   lam_min=",jacobi([[0,1],[1,0]])[0])
