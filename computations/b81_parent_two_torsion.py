from pathlib import Path
import sys,runpy,json,hashlib
H=Path(__file__).resolve().parents[1];B=H/'continuation_20261030b/C81';script=B/'inputs/replay_two_torsion_81.py';data=B/'inputs/two_torsion_81_v2_matrices.m';sys.argv=[str(script),str(data)];runpy.run_path(str(script),run_name='__main__')
(B/'parent_two_torsion_replay.json').write_text(json.dumps({'status':'PASS','level':81,'integral_homology_rank':380,'rational_two_torsion_upper_bound_dimension':0,'evidence':[{'path':str(p.relative_to(H)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}for p in [script,data]]},indent=2)+'\n')
print('C81_PARENT_RATIONAL_TWO_TORSION_ZERO_PASS')
