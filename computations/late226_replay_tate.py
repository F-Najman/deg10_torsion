from pathlib import Path
import runpy,sys
P=Path("work/deg10/continuation_20261029b/C226/tate_alignment_inputs").resolve()
sys.argv=[str(P/"tate_cusp_replay.py"),str(P)]
runpy.run_path(str(P/"tate_cusp_replay.py"),run_name="__main__")
