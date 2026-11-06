from sage.all import *
from pathlib import Path
from collections import defaultdict,Counter
import ast,json,itertools,time
P=Path('work/deg10/continuation_20261031b/M3N18');T=time.monotonic()
def rd(n):return ast.literal_eval((P/n).read_text().replace('[*','[').replace('*]',']'))
g7=rd('paired_p7_group_data.m');g13=rd('paired_p13_group_data.m');assert g7[1:3]==g13[1:3]
assert g7==rd('inputs/M3N18_9multiple_fresh_p7_group_data.m') and g13==rd('inputs/M3N18_9multiple_aut27_p13_group_data.m')
r7=[v[:-1] for v in g7[4]];r13=[v[:-1] for v in g13[4]];c7=[v[:-1] for v in g7[5]];c13=[v[:-1] for v in g13[5]];n7=len(c7[0]);n13=len(c13[0]);ds=g7[2];base=ds.index(1);assert c7[base]==[0]*n7

def hnf(vs):return matrix(ZZ,vs).hermite_form(transformation=True,include_zero_rows=False)
R7,_=hnf(r7);R13,_=hnf(r13);HC,UHC=hnf(r7+c7);assert R7.det()/HC.det()==20253807
DG,_=hnf([v+[0]*n13 for v in r7]+[[0]*n7+v for v in r13]+[v+w for v,w in zip(c7,c13)])
assert R7.det()*R13.det()/DG.det()==20253807

def keyer(H):
 V=H.inverse();den=V.denominator();A=matrix(ZZ,den*V)
 def fun(v):return tuple(int(z%den) for z in vector(ZZ,v)*A)
 return fun
kr7=keyer(R7);key=keyer(DG);key7=lambda v:key(v+[0]*n13);key13=lambda v:key([0]*n7+v)
a=rd('split_hecke_endomorphisms_arithmetic.txt');assert a[2]==g7[3] and a[3]==g7[4] and a[4]==g7[5]
ex=[v[:-1] for v in a[6]];HH,UU=hnf(r7+c7+ex);assert HH==identity_matrix(ZZ,n7)
Amats=[]
for q,imgs,out,idx in a[8]:
 assert idx==27 and len(imgs)==27
 Y=matrix(ZZ,[[0]*n7 for v in r7]+[v[:-1] for v in imgs]);AA=UU*Y
 assert all(kr7(list(vector(ZZ,v)*AA))==kr7([0]*n7) for v in r7+c7)
 Amats.append((q,AA))
def allowed(v):return all(kr7(list(vector(ZZ,v)*A))==kr7([0]*n7) for q,A in Amats)

def canon(qt,p,pol):
 F=GF(p**(len(pol)-1),'z',modulus=PolynomialRing(GF(p),'u')(pol));z=F.gen();v=[sum(F(a)*z**i for i,a in enumerate(w)) for w in qt];return min(tuple(int((a**(p**j)).polynomial()[i]) for a in v for i in range(F.degree())) for j in range(F.degree()))
def own(p,d):
 a=rd(f'paired_p{p}_degree{d}_data.m');assert a[0:2]==[p,d]
 if d<=3:
  name=(f'inputs/M3N18_9multiple_aut27_p13_degree{d}_data.m' if p==13 else f'inputs/M3N18_9multiple_fresh_p7_degree{d}_data.m');old=rd(name);assert a[2]==old[2];lookup={canon(v,p,old[2]):cl for v,cl in zip(old[4],old[3])};repl=[];kr=kr7 if p==7 else keyer(R13)
  for qt,cl in zip(a[4],a[3]):
   known=lookup[canon(qt,p,a[2])];assert kr([9*(x-y) for x,y in zip(cl[:-1],known[:-1])])==kr([0]*(len(cl)-1));repl.append(known)
  a[3]=repl
 return a

report=[]
for p in [7,13]:
 G=g7 if p==7 else g13;R=R7 if p==7 else R13;kr=keyer(R);vc=[v[:-1] for v in G[5]];bind=[]
 for d in[2,3]:
  a=rd(f'paired_p{p}_degree{d}_data.m');old=rd(f'inputs/M3N18_9multiple_fresh_p7_degree{d}_data.m'if p==7 else f'inputs/M3N18_9multiple_aut27_p13_degree{d}_data.m');lookup={canon(qt,p,old[2]):cl for qt,cl in zip(old[4],old[3])}
  for qt,cl in zip(a[4],a[3]):
   known=lookup[canon(qt,p,a[2])];assert kr([x-y for x,y in zip(cl[:-1],known[:-1])])==kr([0]*(len(cl)-1));vc.append(cl[:-1]);bind.append([d,canon(qt,p,a[2])])
 H,_=hnf(([v[:-1]for v in G[4]])+vc);report.append({'p':p,'exact_point_bindings':len(bind),'generated_index':int(abs(H.det())),'all_class_coordinates_exact':True})
print(report);(P/'class_frame_audit.json').write_text(json.dumps(report,indent=2)+'\n');print('D31B_M318_CLASS_FRAME_AUDIT_COMPLETE')
