"""Independent Hensel-series replay of all 29 genus-46 RR jet matrices."""
from sage.all import *
from pathlib import Path
import ast, json, time, hashlib

ROOT = Path('work/deg10/continuation_20261030b/C81')
START = time.monotonic()
def read(name):
    return ast.literal_eval((ROOT/name).read_text().replace('[*','[').replace('*]',']'))
def freeze(a):
    return tuple(map(freeze,a)) if isinstance(a,list) else a
MODEL=read('y46_kummer_model.m')
F2=GF(2)
R2=PolynomialRing(F2,'T')
F4=GF(4,'h',modulus=R2([1,1,1])); h=F4.gen()
PREC=64

class LocalArc:
    def __init__(self,meta):
        self.meta=meta
        self.mult=int(meta[1])
        self.K=GF(2**(len(meta[2])-1),'a',modulus=R2(meta[2]))
        K=self.K; a=K.gen()
        def absolute(v):
            return sum((K(c)*a**i for i,c in enumerate(v)),K(0))
        self.absolute=absolute
        self.eta=absolute(meta[3]); self.z=absolute(meta[4])
        assert self.eta**2+self.eta+1==0
        self.e=(len(meta[2])-1)//2
        basis=[v for i in range(self.e) for v in [self.z**i,self.eta*self.z**i]]
        self.binarybasis=matrix(F2,[list(v.vector()) for v in basis])
        assert self.binarybasis.is_invertible()
        self.binverse=self.binarybasis.inverse()
        self.S=LaurentSeriesRing(K,'s',default_prec=PREC)
        s=self.S.gen(); self.s=s
        t0,y0=map(absolute,meta[5])
        self.t0=t0;self.y0=y0
        # The curve polynomial has polynomial coefficients at these affine points.
        def cv(bits):return K(bits[0])+K(bits[1])*self.eta
        self.cv=cv
        self.pcache={};self.ecache={}
        def rat(c,t):
            n,d=c
            return self.poly(n,t)/self.poly(d,t)
        def eq(t,y):return sum((rat(c,t)*y**i for i,c in enumerate(MODEL[0])),self.S(0))
        def dy(t,y):return sum((rat(c,t)*y**(i-1) for i,c in enumerate(MODEL[0]) if i%2),self.S(0))
        def dp(cs,t):return sum((cv(c)*t**(i-1) for i,c in enumerate(cs) if i%2),self.S(0))
        def dt(t,y):
            ans=self.S(0)
            for i,(n,d) in enumerate(MODEL[0]):
                nn=self.poly(n,t);dd=self.poly(d,t)
                ans+=(dp(n,t)*dd-nn*dp(d,t))/dd**2*y**i
            return ans
        assert eq(self.S(t0),self.S(y0))==0
        if dy(self.S(t0),self.S(y0))!=0:
            t=self.S(t0)+s+O(s**PREC);y=self.S(y0)+O(s**PREC)
            for _ in range(8):y-=eq(t,y)/dy(t,y)
            self.parameter='t'
        else:
            assert dt(self.S(t0),self.S(y0))!=0
            y=self.S(y0)+s+O(s**PREC);t=self.S(t0)+O(s**PREC)
            for _ in range(8):t-=eq(t,y)/dt(t,y)
            self.parameter='y'
        assert eq(t,y).valuation()>=PREC-2
        self.t=t;self.y=y
        self.tau=self.ef(meta[6]);assert self.tau.valuation()==1
        self.f=self.ef(MODEL[2]);assert self.f.valuation()==0
        self.roots=list(map(absolute,meta[7]));self.actual=absolute(meta[8])
        assert len(set(self.roots))==3 and self.actual in self.roots
        self.wseries={}
        for w0 in self.roots:
            assert w0**3==self.f[0]
            w=self.S(w0)+O(s**PREC)
            for _ in range(8):w-=(w**3-self.f)/w**2
            assert (w**3-self.f).valuation()>=PREC-2
            self.wseries[w0]=w
        # Independently construct the exact elliptic curve and Velu x-image.
        d=len(SOURCEKEY)//2
        b=sum((K(v)*self.z**i for i,v in enumerate(SOURCEKEY[:d])),K(0))
        c=sum((K(v)*self.z**i for i,v in enumerate(SOURCEKEY[d:])),K(0))
        E=EllipticCurve(K,[1-c,-b,-b,0,0]);P=E(0,0)
        assert 81*P==E(0) and 27*P!=E(0)
        U=P if self.mult==1 else 2*P
        T=3*U;x0,yy=T.xy();ai=E.a_invariants()
        a3=2*yy+ai[0]*x0+ai[2];a2=3*x0+ai[1]
        a4=3*x0*x0+2*ai[1]*x0+ai[3]-ai[0]*yy
        tangent=a4/a3;a2-=tangent*tangent+ai[0]*tangent
        kap=a3/a2;bb=-a2**3/a3**2;cc=1-(ai[0]+2*tangent)*a2/a3
        ET=EllipticCurve(K,[1-cc,-bb,-bb,0,0])
        V=ET((U[0]-x0)/kap**2,(U[1]-yy-tangent*(U[0]-x0))/kap**3)
        assert 3*V==ET(0,0)
        assert [self.ef(v)[0] for v in MODEL[9]]==list(ET.a_invariants())
        Q=9*ET(0,0)
        vx=V[0]+(V+Q)[0]+(V-Q)[0]-2*Q[0]
        assert vx==absolute(meta[9])
        aval=self.ef(MODEL[7])[0];aa=self.ef(MODEL[5])[0];bb=self.ef(MODEL[6])[0]
        assert [w for w in self.roots if aval+aa*w+bb*w*w==vx]==[self.actual]
        nx=[self.ef(c)[0] for c in MODEL[3]]
        assert sum((v*vx**i for i,v in enumerate(nx)),K(0))==0
        # Check that native local and relative-field encodings agree exactly.
        for v in self.roots+[self.actual,self.eta,self.z]:
            co=self.relative(v)
            assert sum((self.cv(q)*self.z**i for i,q in enumerate(co)),K(0))==v

    def poly(self,cs,t):
        ans=t.parent()(0)
        for c in reversed(cs):ans=ans*t+self.cv(c)
        return ans

    def rat(self,c):
        key=freeze(c)
        if key not in self.pcache:
            n,d=c;self.pcache[key]=self.poly(n,self.t)/self.poly(d,self.t)
        return self.pcache[key]

    def ef(self,v):
        key=freeze(v)
        if key not in self.ecache:
            ans=self.S(0)
            for c in reversed(v):ans=ans*self.y+self.rat(c)
            self.ecache[key]=ans
        return self.ecache[key]

    def relative(self,v):
        row=vector(F2,list(self.K(v).vector()))*self.binverse
        return [[int(row[2*i]),int(row[2*i+1])] for i in range(self.e)]

alltests=[];arcs=[]
for caseid in [1,2,3,7,10]:
    meta=read('rr_case%d_localmeta.m'%caseid)
    SOURCEKEY=meta[2]
    locals_=[LocalArc(v) for v in meta[5]]
    for n,L in enumerate(locals_):
        arcs.append({'case':caseid,'block':n,'base_parameter':L.parameter,'relative_degree':L.e,'multiplicity':L.mult,'tau_valuation':int(L.tau.valuation()),'hensel_precision':PREC-2,'exact_velu_branch':True})
    cusps=range(1,10) if caseid in [1,2,3] else [0]
    for ci in cusps:
        proof=read('rr_case%d_cusp%d.m'%(caseid,ci))
        entries=[(j,v) for j,bs in proof[0] for v in bs]
        rows=[[] for _ in entries]
        for L in locals_:
            for w0 in L.roots:
                if w0==L.actual:continue
                w=L.wseries[w0]
                for order in range(L.mult):
                    for i,(j,v) in enumerate(entries):
                        sec=L.tau**L.mult*L.ef(v)*w**j
                        assert sec.valuation()>=0 and sec.precision_absolute()>=4
                        val=sec[0] if order==0 else sec[1]/L.tau[1]
                        rows[i].extend(L.relative(val))
        assert rows==proof[1],('entrywise matrix mismatch',caseid,ci)
        M=matrix(F4,[[F4(v[0])+h*F4(v[1]) for v in row] for row in rows])
        rank=int(M.rank());assert M.nrows()-rank==1
        record={'case':caseid,'cusp':ci,'rows':int(M.nrows()),'columns':int(M.ncols()),'rank':rank,'h0':1,'matrix_entrywise_equal':True}
        alltests.append(record)
        print('SAGE_C81_INDEPENDENT_LOCAL_ARCS',caseid,ci,M.nrows(),M.ncols(),rank,flush=True)
out={'status':'PASS','method':'independent Sage Hensel arcs, explicit Velu branches, entrywise full jet matrix comparison','series_precision':PREC,'local_arcs':arcs,'tests':alltests,'elapsed_seconds':time.monotonic()-START}
assert len(alltests)==29
(ROOT/'independent_sage_local_arcs.json').write_text(json.dumps(out,indent=2)+'\n')
print('C81_SAGE_LOCAL_ARCS_ALL_29_PASS',time.monotonic()-START,flush=True)
