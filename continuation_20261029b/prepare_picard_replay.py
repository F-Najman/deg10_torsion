from pathlib import Path
import ast
H=Path(__file__).resolve().parents[1];P=Path(__file__).resolve().parent;C=P/'C240';D=C/'portable'
def m(x):
 if isinstance(x,int):return str(x)
 assert isinstance(x,list)
 if all(isinstance(y,int) for y in x):return '['+','.join(map(str,x))+']'
 return '[*'+','.join(m(y) for y in x)+'*]'
for s in range(8):
 w=ast.literal_eval(ast.parse((D/f'shard{s}_RR_identity.log.witnesses.py').read_text()).body[0].value)
 (D/f'replay{s}_data.m').write_text('ReplayWitnesses:='+m(w)+';\n')
 (H/'computations'/f'late240_replay_s{s}.m').write_text(f'SHARD:={s};\nload "work/deg10/continuation_20261029b/C240/replay_ideal_identities.m";\n')
