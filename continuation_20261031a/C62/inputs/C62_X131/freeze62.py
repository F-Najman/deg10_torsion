#!/usr/bin/env python3
"""Freeze the reviewed proof exactly once."""
from pathlib import Path
import json,hashlib
P=Path(__file__).resolve().parent;H=P.parents[2];B=P.parent;out=P/'manifest62.json';assert not out.exists()
paths={p for p in P.rglob('*') if p.is_file() and p.suffix not in ('.pyc',) and '__pycache__' not in p.parts and p.name not in ('manifest62.json','exclusion62_pending.tex')}
paths.add(B/'run_native_job.py')
# Bind the complete local software snapshot, including all loaded specifications.
S=json.loads((H/'shared/snapshot_manifest.json').read_text());paths.add(H/'shared/snapshot_manifest.json')
for e in S['files']:
 q=H/e['path'];assert hashlib.sha256(q.read_bytes()).hexdigest()==e['sha256'];paths.add(q)
for n in ['gonality','DEHMZ.tex','citation_correction_20261026_v2/exclusion35.tex']:paths.add(H/'realizations/continuation/C35'/n)
jobs=[]
for tag,marker in [('C62_N31_presentation','NONIC62_N31_PRESENTATION_COMPLETE N=31 fb_size=3481 relation_rows=4391'),('C62_N31_classify','NONIC62_N31_CLASSIFICATION_COMPLETE tests=2190 survivors=5 complete=true')]:
 status=B/(tag+'_status.json');log=B/(tag+'.log');s=json.loads(status.read_text());assert s['native_exit_code']==0 and s['status']=='EXITED' and marker in log.read_text();paths.update([status,log]);jobs.append({'status':status.relative_to(H).as_posix(),'log':log.relative_to(H).as_posix(),'completion_marker':marker})
M={'degree':9,'excluded_group':[1,62],'complete':True,'method':'Projection to rank-zero J1(31), complete primary class obstruction at3 and positive nonic31 divisor rigidity','source_genus':91,'target_genus':26,'prime':3,'FB_places':3481,'relations':4391,'tests':2190,'positive_survivors':5,'possible62_lifts':0,'native_jobs':jobs,'files':[{'path':q.relative_to(H).as_posix(),'sha256':hashlib.sha256(q.read_bytes()).hexdigest()} for q in sorted(paths)]}
out.write_text(json.dumps(M,indent=2)+'\n');print('C62_FROZEN files='+str(len(M['files']))+' sha256='+hashlib.sha256(out.read_bytes()).hexdigest())
