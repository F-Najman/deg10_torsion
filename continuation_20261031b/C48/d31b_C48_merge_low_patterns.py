from pathlib import Path
from itertools import product,combinations_with_replacement
from collections import Counter,defaultdict
import json,re,time
T=time.monotonic();P=Path('work/deg10/continuation_20261031b/C48')
def rd(n):return json.loads((P/n).read_text().replace('[*','[').replace('*]',']'))
old={'__builtins__':{}};exec((P/'inputs/old_target_actions.py').read_text(),old);rels=old['quotient_relations'];mods=[rels[i][i]for i in range(4)];assert rels==[[mods[i]if i==j else 0 for j in range(4)]for i in range(4)];assert mods==[2,2,24,48]
cl={}
for sh in range(3):
 a=rd(f'low_class_s{sh}.m');assert a[2]==rels and a[3]==old['quotient_images_of_G_generators']
 for row in a[0]:assert row[0]not in cl;cl[row[0]]=tuple(row[5])
assert len(cl)==3550
A={tuple(label):im for label,im in zip(old['labels'],old['action_images'])};inv={}
for label,M in A.items():
 back={}
 for v in product(*(range(n)for n in mods)):
  w=tuple(sum(v[i]*M[i][j]for i in range(4))%mods[j]for j in range(4));assert w not in back;back[w]=v
 assert len(back)==4608;inv[label]=back
labels=rd('geometric_action_labels.m');geom=rd('low_geometry.m');records=geom[1];byd=defaultdict(list);anchors=defaultdict(list);actual={}
for row in records:
 rid,d,key,rep,ti,h=row;actual[rid]=inv[tuple(labels[h-1])][cl[ti]];byd[d].append(rid)
 if rep:anchors[d].append(rid)
assert[len(byd[d])for d in range(3,9)]==[16,64,656,2712,10784,6172]
assert[len(anchors[d])for d in range(3,9)]==[2,8,82,342,1348,6172]
s=(P/'inputs/portable_inputs.m').read_text()
def portable(n):return json.loads(re.search(n+r'\s*:=\s*(.*?);',s,re.S)[1].replace('[*','[').replace('*]',']'))
qimg=old['quotient_images_of_G_generators']
def q(v):return tuple(sum(v[i]*qimg[i][j]for i in range(9))%mods[j]for j in range(4))
allowed=[set()for _ in range(8)];allowed[0].add((0,0,0,0))
for d,v in zip(portable('PortableCuspDegrees'),portable('PortableCuspClasses')):
 w=q(v)
 for j in range(d,8):allowed[j]|={tuple((a+b)%m for a,b,m in zip(w,x,mods))for x in allowed[j-d]}
assert list(map(len,allowed))==[1,5,23,75,178,378,661,1065]
tests=[];left=[]
def test(ds,ids):
 cd=10-sum(ds);v=tuple(sum(actual[i][j]for i in ids)%mods[j]for j in range(4));yes=tuple(-z%m for z,m in zip(v,mods))in allowed[cd];row=[len(tests)+1,ds,list(ids),cd,list(v),yes];tests.append(row)
 if yes:left.append(row)
for d in range(3,9):
 for i in anchors[d]:test([d],[i])
for d,e in[(3,3),(3,4),(3,5),(3,6),(3,7),(4,4),(4,5),(4,6),(5,5)]:
 for i in anchors[d]:
  for j in byd[e]:test([d,e],[i,j])
for i in anchors[3]:
 for j,k in combinations_with_replacement(byd[3],2):test([3,3,3],[i,j,k])
 for j in byd[3]:
  for k in byd[4]:test([3,3,4],[i,j,k])
summary={'status':'PASS','support_cases':len(tests),'actual_source_points':len(records),'canonical_target_queries':3550,'all_mixed_actions_undone_before_summing':True,'group_action_frame_exactly_matched':True,'allowed_cusp_sizes':list(map(len,allowed)),'remaining':left,'remaining_count':len(left),'remaining_types':[[list(k),n]for k,n in Counter(tuple(v[1])for v in left).items()],'geometric_labels':labels,'elapsed_seconds':time.monotonic()-T}
(P/'low_patterns_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
def mag(v):
 if isinstance(v,bool):return'true'if v else'false'
 if isinstance(v,int):return str(v)
 flat=all(isinstance(z,int)and not isinstance(z,bool)for z in v);return('['if flat else'[*')+','.join(mag(z)for z in v)+(']'if flat else'*]')
(P/'low_patterns_remaining.m').write_text(mag([geom[0],records,left])+'\n');(P/'low_pattern_tests.json').write_text(json.dumps(tests,separators=(',',':'))+'\n');print('D31B_C48_LOW_PATTERN_SIEVE_COMPLETE',len(tests),len(left),summary['remaining_types'],time.monotonic()-T,flush=True)
