from pathlib import Path
from math import gcd
from collections import Counter
import json,hashlib
BASE=Path(__file__).resolve().parent

def rd(p):return json.loads(p.read_text().replace('[*','[').replace('*]',']'))
def mag(v):
 if isinstance(v,int):return str(v)
 if all(isinstance(a,int)for a in v):return '['+','.join(map(str,v))+']'
 return '[*'+','.join(mag(a)for a in v)+'*]'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
for N in[48,54,66]:
 root=BASE/f'C{N}';fp=rd(root/'frobenius_p5.m');g=(len(fp)-1)//2;tr=[];closed=[]
 for d in range(1,11):
  s=-d*fp[2*g-d]-sum(fp[2*g-i]*tr[d-i-1]for i in range(1,d));tr.append(s);nn=5**d+1-s-sum(i*closed[i-1]for i in range(1,d)if d%i==0);assert nn%d==0;closed.append(nn//d)
 labels={(a,b)for a in range(N)for b in range(N)if gcd(gcd(a,b),N)==1};unseen=set(labels);cs=[]
 while unseen:
  a,b=min(unseen);o={(s*(a+k*b)%N,s*b%N)for k in range(N)for s in[-1,1]};unseen-=o;cs.append(sorted(o))
 ix={x:i for i,C in enumerate(cs)for x in C}
 def perm(u):return[ix[(u*C[0][0]%N,C[0][1])]for C in cs]
 def orbits(perms):
  unseen=set(range(len(cs)));out=[]
  while unseen:
   x=min(unseen);todo=[x];o={x}
   while todo:
    a=todo.pop()
    for pp in perms:
     b=pp[a]
     if b not in o:o.add(b);todo.append(b)
   unseen-=o;out.append(tuple(sorted(o)))
  return set(out)
 glob=orbits([perm(u)for u in range(1,N)if gcd(u,N)==1]);fin=orbits([perm(5)]);sig=Counter(map(len,fin));nc=[closed[d-1]-sig[d]for d in range(1,11)];fields=[];records=[];counts=[];orbc=[];inputs=[]
 for d in range(1,11):
  ns=1 if d<=7 else 6;keys=set();cnt=0;pol=None
  for s in range(ns):
   p=root/f'source_torsion_p5_d{d}_s{s}.m';a=rd(p);assert a[0]==d and a[2]is True and a[3]==s and a[4]==ns
   assert pol is None or pol==a[1];pol=a[1];cnt+=a[8][d-1]
   for e,key in a[9]:
    if e!=d:continue
    assert len(key)==2*d and all(0<=v<5 for v in key)and tuple(key)not in keys;keys.add(tuple(key))
   inputs.append({'path':str(p.relative_to(BASE.parent.parent.parent)),'sha256':sha(p),'closed_full_degree':a[8][d-1],'diamond_orbits_full_degree':sum(e==d for e,key in a[9])})
  assert cnt==nc[d-1],(N,d,cnt,nc[d-1]);counts.append(cnt);orbc.append(len(keys));fields.append([d,pol]);records.extend([d,list(v)]for v in sorted(keys))
 out=[fields,nc,records];(root/f'source{N}_p5.m').write_text(mag(out)+'\n');(root/f'source{N}_p5.json').write_text(json.dumps(out,separators=(',',':'))+'\n')
 cusp={'level':N,'geometric_cusps':len(cs),'primitive_label_orbits':cs,'global_orbits':[list(v)for v in sorted(glob)],'Frob5_orbits':[list(v)for v in sorted(fin)],'global_orbits_intact_at5':glob==fin,'global_degree_signature':dict(Counter(map(len,glob))),'finite_degree_signature':dict(sig)};(root/'cusp_actions_p5.json').write_text(json.dumps(cusp,indent=2)+'\n')
 summary={'status':'PASS','degree':10,'level':N,'genus':g,'noncuspidal_closed_counts':nc,'source_diamond_orbits_by_degree':orbc,'source_diamond_orbits_total':len(records),'fresh_Frobenius_Newton_counts_agree':True,'census_algorithm':'Every finite j modulo Frobenius, all twists; exact halving/triplication or division-polynomial roots, no elliptic point counts','inputs':inputs};(root/'source_census_verified.json').write_text(json.dumps(summary,indent=2)+'\n');print('D31A_F5_TORSION_CENSUS_VERIFIED',N,nc,orbc,len(records))
