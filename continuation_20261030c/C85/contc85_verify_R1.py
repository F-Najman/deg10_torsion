#!/usr/bin/env python3
"""Complete exact R1 multiplication sieve for actual Hecke target diamond orbits."""
from pathlib import Path
import json,gzip,time,random,sys
from arithmetic10 import *
H=Path(__file__).resolve().parent;t0=time.monotonic();REPLAY='--readonly'in sys.argv;rng=random.Random(20261026)
a=(H/'inputs/canonical_unit_stage1_data.m').read_text();b=(H/'jets10.m').read_text();D=matbits(a,'Diamond2Matrix');D3=matbits(a,'Diamond3Matrix');labels=[min((pow(2,-aa,85)*pow(3,-bb,85))%85,85-(pow(2,-aa,85)*pow(3,-bb,85))%85)for bb in range(16)for aa in range(2)];assert len(set(labels))==32;index={u:i for i,u in enumerate(labels)};rep=lambda u:min(u%85,85-u%85);sigma=lambda i:index[rep(labels[i]*43)];tau=index[3];st=sigma(tau);si=sigma(0);sr=lit(b,'SourceRecords');tr=lit(b,'SourceHeckeRecords');assert len(sr)==len(tr)==4;T=json.loads((H/'inputs/canonical_ring_R1_tensor.json').read_text());assert T['dimensions']==[225,672];TT=[[int(v,16)for v in row]for row in T['product_1_1_into_2']]
class Span:
 def __init__(self,rows=()):
  self.p={}
  for v in rows:self.add(v)
 def add(self,v):
  while v:
   k=(v&-v).bit_length()-1
   if k in self.p:v^=self.p[k]
   else:self.p[k]=v;return True
  return False
 def quotient(self,v):
  for k in sorted(self.p):
   if v>>k&1:v^=self.p[k]
  return v
 def __len__(self):return len(self.p)
def kernel(rows):
 p={};out=[]
 for i,v in enumerate(rows):
  c=1<<i
  while v:
   k=(v&-v).bit_length()-1
   if k in p:a,b=p[k];v^=a;c^=b
   else:p[k]=(v,c);break
  else:out.append(c)
 return out
families={};keys={};SOURCE_ORBITS=[];TARGET_ORBITS=[]
def orbit(rec,degree):
 v=decode_record(rec);F=v['F'];z=v['jets'];assert F.normalize(z)==z;original=z;oo=[]
 for b in range(16):
  for aa in range(2):
   za=z if aa==0 else transform(F,z,D);key=F.key(za,degree)
   if key not in keys:
    fid=len(keys)+1;keys[key]=fid;families[fid]={'F':F,'jets':za,'degree':degree,'evals':{}}
   fid=keys[key];assert families[fid]['degree']==degree;oo.append(fid)
  z=transform(F,z,D3)
 assert z==original and len(set(oo))==32;return oo
for sid,(s,targets)in enumerate(zip(sr,tr),1):
 SOURCE_ORBITS.append(orbit(s,10));TARGET_ORBITS.append([(mult,orbit(r,deg))for e,deg,mult,r in targets]);assert sum(m*families[o[0]]['degree']for m,o in TARGET_ORBITS[-1])==40
 print('DEG10_C85_DIAMOND_JET_ORBITS_READY source',sid,'families_so_far',len(families),'elapsed',time.monotonic()-t0,flush=True)
assert len(set(sum(SOURCE_ORBITS,[])))==128

def evaluation(fid,m):
 v=families[fid];assert 1<=m<=5
 if m not in v['evals']:
  F=v['F'];v['evals'][m]=[sum(a[n]<<(F.d*n)for n in range(m))for a in v['jets']]
 return v['evals'][m],v['F'].d*m

def conditions(dd,positive):
 rows=[0]*225;width=0;degree=0;support=[]
 for fid,c in sorted(dd.items()):
  c=c if positive else-c
  if c<=0:continue
  ev,w=evaluation(fid,c);rows=[x|(y<<width)for x,y in zip(rows,ev)];width+=w;degree+=families[fid]['degree']*c;support.append([fid,c])
 return rows,width,degree,support

from sage.all import GF,Matrix,vector
F2=GF(2);saved=json.loads(gzip.open(H/'canonical_R1_matrices.json.gz','rt').read());results=json.loads((H/'canonical_R1_sieve_verified.json').read_text())['results'];assert len(saved)==len(results)==4
reports=[]
def matrix_bits(rows,width):return Matrix(F2,[[v>>j&1 for j in range(width)]for v in rows])
def productmap(bits):return matrix_bits([xorrows(bits,[row[j]for row in TT])for j in range(225)],672)
for sid,(so,to,old,entry)in enumerate(zip(SOURCE_ORBITS,TARGET_ORBITS,saved,results),1):
 dd={}
 def add(fid,n):dd[fid]=dd.get(fid,0)+n
 for mult,o in to:add(o[sigma(0)],mult);add(o[0],-mult)
 add(so[si],-1);add(so[0],1);add(so[st],-3);add(so[tau],3)
 bp,pw,pd,ps=conditions(dd,True);bm,nw,nd,ns=conditions(dd,False);assert pd==nd==80
 CP=matrix_bits(bp,pw);CM=matrix_bits(bm,nw);assert CP.rank()==CM.rank()==80
 ff=int(old['F'],16);fv=vector(F2,[ff>>i&1 for i in range(225)]);assert fv*CP==0
 MF=productmap(ff);assert MF.rank()==225;N=MF.right_kernel().basis_matrix().transpose();assert N.ncols()==447 and MF*N==0
 A=CM
 for value in old['H']:
  hh=int(value,16);hv=vector(F2,[hh>>i&1 for i in range(225)]);assert hv*CP==0;A=A.augment(productmap(hh)*N)
 assert A.rank()==225
 stored=matrix_bits([int(v,16)for v in old['rows']],old['columns']);assert stored.rank()==225
 reports.append(dict(source=sid,positive_rank=80,negative_rank=80,multiplication_rank=225,independent_joint_rank=225))
out=dict(status='PASS',cases=4,necessary_convention=2,simultaneous_target_shift=0,domain_dimension=225,product_dimension=672,reports=reports,seconds=time.monotonic()-t0);(H/'R1_independent_sage.json').write_text(json.dumps(out,indent=2)+'\n');print('C85_DEG10_CANONICAL_MULTIPLICATION_INDEPENDENT_SAGE_PASS',time.monotonic()-t0,flush=True)
