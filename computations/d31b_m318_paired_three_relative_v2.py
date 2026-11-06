from sage.all import *
from pathlib import Path
from collections import defaultdict,Counter
import ast,json,itertools,time,re
P=Path('work/deg10/continuation_20261031b/M3N18');T=time.monotonic()
def rd(n):
 s=(P/n).read_text();s=re.sub(r'\[\s*(-?\d+)\s*\.\.\s*(-?\d+)\s*\]',lambda m:str(list(range(int(m[1]),int(m[2])+1))),s);s=re.sub(r'\[\s*\*','[',s);s=re.sub(r'\*\s*\]',']',s);return ast.literal_eval(s)
assert rd('kernel_structure.m')[1:]==[[3,3,3],3]
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
a=rd('split_hecke_framed_arithmetic.txt');assert a[2]==g7[3] and a[3]==g7[4] and a[4]==g7[5]
gp=rd('split_hecke_generator_points.m');old2=rd('inputs/M3N18_9multiple_fresh_p7_degree2_data.m')
def simplecanon(qt,p,pol):
 F=GF(p**(len(pol)-1),'z',modulus=PolynomialRing(GF(p),'u')(pol));z=F.gen();v=[sum(F(a)*z**i for i,a in enumerate(w)) for w in qt];return min(tuple(int((a**(p**j)).polynomial()[i]) for a in v for i in range(F.degree())) for j in range(F.degree()))
lookup={simplecanon(qt,7,old2[2]):cl for qt,cl in zip(old2[4],old2[3])}
for pr,cl in zip(gp,a[6]):
 known=lookup[simplecanon(pr[2],7,pr[1])];assert kr7([x-y for x,y in zip(cl[:-1],known[:-1])])==kr7([0]*n7)
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
own7={d:own(7,d) for d in range(2,6)};own13={d:own(13,d) for d in range(2,5)}
all7={d:rd(f'inputs/M3N18_9multiple_fresh_p7_degree{d}_data.m')[3] for d in[2,3]}
all13={d:rd(f'inputs/M3N18_9multiple_aut27_p13_degree{d}_data.m')[3] for d in[2,3]}
own13[4]=rd('inputs/M3N18_9multiple_aut27_p13_degree4_data.m')
rat=[i for i,d in enumerate(ds) if d==1];cub=[i for i,d in enumerate(ds) if d==3]

def add(vs,n):return [sum(v[i] for v in vs) for i in range(n)]
def comp(c,d,n):
 if d==0:return [([], [0]*n)]
 ans=[]
 for ids in itertools.combinations_with_replacement(rat,d):ans.append(([i+1 for i in ids],add([c[i] for i in ids],n)))
 if d==3:ans +=[([i+1],c[i]) for i in cub]
 return ans

def patterns(p,ownrows,allrows,c,n):
 out=[]
 for d in sorted(ownrows):
  for i,v in enumerate(ownrows[d][3],1):
   for ids,cv in comp(c,5-d,n):out.append(([d,i,ids],add([v[:-1],cv],n)))
 for i,v in enumerate(ownrows[2][3],1):
  for j,w in enumerate(allrows[2],1):
   for ids,cv in comp(c,1,n):out.append(([22,i,j,ids],add([v[:-1],w[:-1],cv],n)))
  for j,w in enumerate(allrows[3],1):out.append(([23,i,j],add([v[:-1],w[:-1]],n)))
 return out
source=patterns(7,own7,all7,c7,n7);assert len(source)==3622
sourcegood=[(i,desc,v) for i,(desc,v) in enumerate(source,1) if allowed(v)]
r7geom=rd('relative_cubic_rank.m');r13geom=rd('relative_cubic_rank_p13.m');assert r7geom[-1]==r13geom[-1]==[]
clear7={b[0] for b in r7geom[4]};clear13={b[0] for b in r13geom[4]};assert clear7=={4,5,7,8} and all(len(b[5])==171 and all(r[2]==4 for r in b[5]) for b in r7geom[4]+r13geom[4])
sourcegood=[r for r in sourcegood if not(r[1][0]==3 and r[1][1]in clear7)];assert len(sourcegood)==242
per=rd('cusp_automorphisms.m')[3];assert len(per)==27
HCinv=HC.inverse();sourcekeys=defaultdict(list)
for ix,desc,cl in sourcegood:
 coef=vector(ZZ,[3*v for v in cl])*HCinv;assert all(z in ZZ for z in coef);coef=vector(ZZ,coef)*UHC;vc=list(coef[len(r7):]);assert kr7(add([[vc[i]*z for z in c7[i]] for i in range(24)],n7))==kr7([3*z for z in cl])
 for h,perm in enumerate(per,1):
  v=add([[vc[i]*z for z in c7[perm[i]-1]] for i in range(24)],n7);adjust=15-sum(ds[i]*vc[i] for i in range(24));v=add([v,[adjust*z for z in c7[perm[base]-1]]],n7)
  sourcekeys[key7(v)].append([ix,h])
print('M318_SOURCE_COMPLETE',len(source),len(sourcegood),len(sourcekeys),flush=True)
# If all four streamed quintic shards exist, merge their disjoint orbit coverage.
shards=[]
for h in range(4):
 f=P/f'paired_lean_s{h}_p13_degree5_data.m'
 if f.exists():shards.append(rd(f.name))
full=len(shards)==4
if full:
 vals=[];coords=[];keys=[];labels=[];images=[];offset=0
 for sh in shards:
  assert sh[:3]==[13,5,shards[0][2]] and len(sh[6])==27*len(sh[3]);vals+=sh[3];coords+=sh[4];keys+=sh[6];labels +=[v+offset for v in sh[7]];images+=sh[8];offset+=len(sh[3])
 assert len(vals)==2792 and len(set(map(tuple,keys)))==len(keys)==75384
 own13[5]=[13,5,shards[0][2],vals,coords,list(range(1,len(vals)+1)),keys,labels,images]
 raw=json.dumps(own13[5]);(P/'paired_merged_p13_degree5_data.m').write_text('[*'+raw[1:-1]+'*]')
target=patterns(13,own13,all13,c13,n13);matches=[]
for i,(desc,cl) in enumerate(target,1):
 if desc[0]==3 and desc[1]in clear13:continue
 ss=sourcekeys.get(key13([-3*z for z in cl]))
 if ss:matches.append({'target_id':i,'target_descriptor':desc,'target_class':cl,'source_matches':ss})
report={'status':'PASS','full_target_quintics_included':full,'source_patterns':len(source),'source_after_split_Hecke':len(sourcegood),'split_Hecke_filter_applied':True,'source_Aut_stable_three_multiple_keys':len(sourcekeys),'target_patterns':len(target),'matched_target_patterns':len(matches),'matched_pattern_types':dict(Counter(str(m['target_descriptor'][0]) for m in matches)),'source_genus':10,'relative_degree':5,'cuspidal_saturation_assumed':False,'elapsed_seconds':time.monotonic()-T}
(P/('paired_three_relative_arithmetic_complete.json' if full else 'paired_three_relative_arithmetic_lower.json')).write_text(json.dumps(report,indent=2)+'\n');(P/('paired_three_relative_matches_complete.json' if full else 'paired_three_relative_matches_lower.json')).write_text(json.dumps(matches,separators=(',',':'))+'\n')
print(json.dumps(report,indent=2));print('D31B_M318_PAIRED_ARITHMETIC_COMPLETE',full,flush=True)
