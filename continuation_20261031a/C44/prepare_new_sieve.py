from pathlib import Path
import json,re
H=Path(__file__).resolve().parent
s=(H/'new_finite_jets.raw.jsonseq').read_text();s=re.sub(r'\[\s*\*','[',s);s=re.sub(r'\*\s*\]',']',s)
dec=json.JSONDecoder();i=0;data=[]
while i<len(s):
 while i<len(s)and s[i].isspace():i+=1
 if i==len(s):break
 a,i=dec.raw_decode(s,i);data.append(a)
assert data[-1]['type']=='complete'and data[-1]['cases']==data[0]['requested_cases']==956
assert data[-1]['distinct_source10']==data[0]['source_places10']==5854 and data[-1]['distinct_source5']==60
(H/'new_finite_jets.normalized.jsonl').write_text(''.join(json.dumps(a,separators=(',',':'))+'\n'for a in data))
J={i:{s:[]for s in ['T','U']}for i in range(1,957)};degrees={}
for a in data:
 if a['type']=='case':degrees[a['case']]=a['positive_degree'];assert a['positive_degree']==a['negative_degree']
 if a['type']=='jet':J[a['case']][a['side']].append(a)
assert len(degrees)==956
path=H/'new_sieve_inputs.m'
with path.open('w')as f:
 f.write('SieveOutput:="'+str(H/'new_matrix_results.raw.jsonseq')+'";\nJetCases:=[* *];\n')
 for i in range(1,957):
  vals=[]
  for side in ['T','U']:
   rr=J[i][side];assert sum(a['degree']*a['multiplicity']for a in rr)==degrees[i]
   vals.append('[*'+','.join('[*'+','.join(map(str,[a['degree'],a['multiplicity'],a['field_polynomial'],a['canonical_jets']]))+'*]'for a in rr)+'*]')
  f.write('Append(~JetCases,[*'+str(i)+','+','.join(vals)+'*]);\n')
print('D31A_C44_NEW_SIEVE_INPUTS_READY 956 cases')
