import os,sys,runpy
os.environ['PYTHONDONTWRITEBYTECODE']='1'
p='work/deg10/continuation_20261030c/C65/inherited_degree9/cyclic/continuation_20261026/C65_completion_20261027/certify.py'
sys.argv=[p,'--deep']
runpy.run_path(p,run_name='__main__')
print('C65_INHERITED_DEGREE_INDEPENDENT_EVIDENCE_ALL_PASS',flush=True)
