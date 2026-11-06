from sage.all import *
from pathlib import Path
import json,re,ast
D=Path(__file__).resolve().parents[1];C=D/'continuation_20261031b/C66';lit=lambda v:ast.literal_eval(v.replace('[*','[').replace('*]',']').replace('<','[').replace('>',']'))
def ass(s,k):return lit(re.search(r'\b'+k+r':=(.*?);',s,re.S)[1])
s=(C/'F5_common_basis58.m').read_text();ac=ass(s,'ACoefficients');bf=ass(s,'BasisFunctions');inp=(C/'inputs/F5_rr/input33.m').read_text();blocks=ass(inp,'GlobalCuspBlocks');old=ass(inp,'CuspRepresentatives');rats=[i for i,b in enumerate(blocks)if b[1]==1];targets=[]
for v in old:
 for j in rats:w=v[:];w[j]+=1;targets.append(w)
assert len(targets)==100 and ac==[max([0]+[v[j]for v in targets])for j in range(17)];assert sum(a*b[1]for a,b in zip(ac,blocks))==58
k=GF(5);R=PolynomialRing(k,'t');K=R.fraction_field();ff=[[K(R(z[0]))/K(R(z[1]))for z in row]for row in bf];assert len(ff)==38 and all(len(z)==12 for z in ff);dn=lcm([c.denominator()for row in ff for c in row]);md=max((dn*c).numerator().degree()for row in ff for c in row);IM=matrix(k,[[R(dn*c)[i]for c in row for i in range(md+1)]for row in ff]);assert IM.rank()==38
raw=(C/'F5_common_basis58_jets.jsonl').read_text();data=[];dec=json.JSONDecoder()
while raw.strip():raw=raw.lstrip();v,i=dec.raw_decode(raw);data.append(v);raw=raw[i:]
used=set()
for z in data:
 j=z['block_index']-1;deg,e,finite,ib=blocks[j][2][z['place_index']-1];assert deg==z['degree'];assert z['A_order']==ac[j]*e;assert z['orders']==[(ac[j]-v[j])*e for v in targets];assert len(z['rows'])==38 and all(len(r)==max(z['orders'])*deg for r in z['rows']);used.add(j)
assert used=={j for j in range(17)if any(ac[j]-v[j]>0 for v in targets)}
ranks=[]
for i in range(100):
 M=matrix(k,[sum([z['rows'][j][:z['orders'][i]*z['degree']]for z in data],[])for j in range(38)]);assert M.ncols()==48 and M.rank()==38;ranks.append(38)
rows=[]
for z in data:rows.append('<'+','.join(map(str,[z['block_index'],z['place_index'],z['degree'],z['A_order']]))+','+repr(z['orders'])+',Matrix(GF(5),38,'+str(len(z['rows'][0]))+','+repr(sum(z['rows'],[]))+')>')
(C/'F5_common_basis58_jet_tables.m').write_text('JetTables:=[*'+',\n'.join(rows)+'*];\n')
out={'status':'PASS','A_degree':58,'basis_dimension_RR':38,'independent_coefficient_rank':38,'coefficient_columns':int(IM.ncols()),'maximum_polynomial_degree':int(md),'jet_tables':len(data),'test_count':100,'target_degrees':[10]*100,'restriction_degrees':[48]*100,'independent_ranks':ranks,'all_h0_zero':True};(C/'F5_common_basis58_independent.json').write_text(json.dumps(out,indent=2)+'\n');print('D31B_C66_COMMON_BASIS58_INDEPENDENT_COMPLETE 100 38',flush=True)
