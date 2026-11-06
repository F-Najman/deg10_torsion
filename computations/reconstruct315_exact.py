from pathlib import Path
import json,ast,re
from sage.all import QQ,ZZ,PolynomialRing,NumberField,EllipticCurve
H=Path(__file__).resolve().parents[1];s=(H/'exclusions/315d5/reconstruction_residues.txt').read_text();rows=[ast.literal_eval(x.replace('[*','[').replace('*]',']')) for x in re.findall(r'\[\*.*?\*\]',s,re.S)];pairs={}
for row in rows:pairs.setdefault(row[0],[]).append(row)
R=PolynomialRing(QQ,'z');F=NumberField(R.gen()**2+R.gen()+1,'z');z=F.gen();T=PolynomialRing(F,'t');acc=[ZZ(0)]*20;mod=ZZ(1);used=[];found=False
for p,rr in sorted(pairs.items()):
 if len(rr)!=2:continue
 _,z1,f1,q1=rr[0];_,z2,f2,q2=rr[1];vals=[]
 for a,b in [(f1[:5],f2[:5]),(q1,q2)]:
  for x,y in zip(a,b):
   v=(x-y)*pow(z1-z2,-1,p)%p;u=(x-v*z1)%p;vals.extend([u,v])
 assert len(vals)==20
 for j in range(20):acc[j]+=mod*(((vals[j]-int(acc[j]))*pow(int(mod),-1,p))%p)
 mod*=p;used.append(p)
 if len(used)<4:continue
 try:cs=[x.rational_reconstruction(mod) for x in acc]
 except (ValueError,ArithmeticError):continue
 fc=[cs[2*i]+cs[2*i+1]*z for i in range(5)]+[F(1)];qc=[cs[10+2*i]+cs[11+2*i]*z for i in range(5)];f=T(fc)
 if not f.is_irreducible():continue
 K=F.extension(f,'w');w=K.gen();q=sum((K(c)*w**i for i,c in enumerate(qc)),K(0));E=EllipticCurve(K,[(1+q)*w+(2-q),0,q*w*w+(1-q)*w,0,0]);P=E([0,0]);Q=E([-w,w*w]);
 if E.discriminant()==0:continue
 if 15*Q!=E(0) or 3*Q==E(0) or 5*Q==E(0):
  print('DECIC315_RECONSTRUCTION_NOT_YET_EXACT primes',used,flush=True);continue
 assert 3*P==E(0) and P!=E(0) and 5*Q not in [P,-P]
 enc=lambda c:[str(a) for a in F(c).list()]+['0']*(2-len(F(c).list()))
 d={'status':'EXACT POINTS VERIFIED; FULL TORSION BOUND PENDING','base_polynomial':[1,1,1],'relative_polynomial':[enc(c) for c in fc],'q_coefficients':[enc(c) for c in qc],'reconstruction_primes':used,'reconstruction_modulus':str(mod),'contains':[3,15]}
 (H/'occurrences/new315_exact_candidate.json').write_text(json.dumps(d,indent=2)+'\n');print('DECIC315_EXACT_CONTAINS_3x15 field_degree=10 f=',f,'q=',T(qc),'j_degree',E.j_invariant().absolute_minpoly().degree(),flush=True);found=True;break
print('DECIC315_EXACT_RECONSTRUCTION_COMPLETE found='+str(found),flush=True)
