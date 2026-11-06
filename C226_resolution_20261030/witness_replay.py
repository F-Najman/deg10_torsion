from pathlib import Path
import ast,json,math,time
from sage.all import QQ,ZZ,PolynomialRing,NumberField,EllipticCurve,GF,prime_range
P=Path('work/deg10/C226_resolution_20261030');tt=time.monotonic()
def read(p):return ast.literal_eval(p.read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
def q(v):return QQ(v[0])/QQ(v[1])
def check_point(idx):
 rows=read(P/f'exact_points_{idx}.m');assert len(rows)==1;fc,bc,cc,mult=rows[0];assert mult==1;R=PolynomialRing(QQ,'x');f=R([q(v) for v in fc]);assert f.degree()==10 and f.is_irreducible();K=NumberField(f,'a');a=K.gen();b=K([q(v) for v in bc]);c=K([q(v) for v in cc]);E=EllipticCurve(K,[0,c,0,b*(1-b-c),0]);T=E(0,0);Q=E(b,b);assert E.discriminant()!=0 and 2*T==E(0) and T!=E(0) and 26*Q==E(0) and 13*Q!=E(0) and 2*Q!=E(0) and T!=13*Q
 g=0;recs=[]
 for p in prime_range(17,2000):
  p=int(p)
  try:
   kp=GF(p);S=PolynomialRing(kp,'x');fp=S(f);bp=S([kp(q(v)) for v in bc]);cp=S([kp(q(v)) for v in cc])
  except (ZeroDivisionError,ValueError):continue
  if fp.gcd(fp.derivative()).degree()!=0:continue
  for root,_ in fp.roots():
   bb=bp(root);ccc=cp(root);ee=EllipticCurve(kp,[0,ccc,0,bb*(1-bb-ccc),0])
   if ee.discriminant()==0:continue
   aa=bb*(1-bb-ccc);count=1+sum(1 if (z:=xx**3+ccc*xx**2+aa*xx)==0 else (2 if z.is_square() else 0) for xx in kp);assert count==ee.cardinality() and count%52==0;g=math.gcd(g,int(count));recs.append([p,int(root),int(bb),int(ccc),int(count)])
   if g==52:break
  if g==52:break
 assert g==52
 print('C226_FULL_TORSION_INDEPENDENTLY_PROVED',idx,recs,time.monotonic()-tt,flush=True)
 r=(13*Q)[0];roots=[K(0),r,-c-r];assert len(set(roots))==3 and all(x**3+c*x*x+b*(1-b-c)*x==0 for x in roots);squares=[];isogenous=[]
 for j,rr in enumerate(roots):
  oth=[s for s in roots if s!=rr]
  for v,other in [(oth[0]-rr,oth[1]-rr),(oth[1]-rr,oth[0]-rr)]:
   yes=bool(v.is_square());squares.append([j,[str(t) for t in v.list()],yes])
   if yes:
    sq=v.sqrt();A=c+3*rr;B=3*rr*rr+2*c*rr+b*(1-b-c);EP=EllipticCurve(K,[0,-2*A,0,A*A-4*B,0]);xx=v-other;R4=EP(xx,2*xx*sq);S13=2*Q;u=S13[0]-rr;W=EP(u+A+B/u,S13[1]*(1-B/u**2));P52=R4+W;assert 4*R4==EP(0) and 2*R4!=EP(0) and 13*W==EP(0) and W!=EP(0) and 52*P52==EP(0) and 26*P52!=EP(0) and 4*P52!=EP(0)
    en=lambda x:[[int(t.numerator()),int(t.denominator())] for t in K(x).list()]
    isogenous.append({'kernel_root':en(rr),'a_invariants':[en(v) for v in EP.a_invariants()],'point52':[en(P52[0]),en(P52[1])]})
 report={'status':'PASS','class':idx,'field_polynomial':[str(t) for t in f.list()],'full_torsion':[2,26],'irreducible_degree':10,'reductions':recs,'torsion_upper_bound':g,'j_degree':int(E.j_invariant().minpoly().degree()),'two_isogeny_square_tests':squares,'cyclic52_witnesses':isogenous,'elapsed_seconds':time.monotonic()-tt}
 (P/f'witness_replay_{idx}.json').write_text(json.dumps(report,indent=2)+'\n');print('C226_EXACT_WITNESS_REPLAY_COMPLETE',idx,'C52 isogenous witnesses',len(isogenous),'j degree',report['j_degree'],flush=True)
 return report
