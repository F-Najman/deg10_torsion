# Written for this repository: LMFDB labels of the fields of van Hoeij's decic points of Proposition 3.4
# (polredabs taken from the census records in inputs/), using the public LMFDB mirror.
import json, glob, sys
from pathlib import Path
import psycopg2
conn = psycopg2.connect(host='devmirror.lmfdb.xyz', dbname='lmfdb', user='lmfdb', password='lmfdb', connect_timeout=20)
cur = conn.cursor()
from sympy import Poly, symbols
x = symbols('x')
for f in sorted(glob.glob(str(Path(__file__).resolve().parents[1] / 'inputs' / '1.*.10.*.json'))):
    d = json.load(open(f))
    p = Poly(d['field']['polredabs'].replace('^', '**'), x)
    coeffs = [int(c) for c in reversed(p.all_coeffs())]
    cur.execute("select label, disc_abs, r2, class_number from nf_fields where coeffs = %s::numeric[]", (coeffs,))
    print(d['id'], d['field']['disc'], d['field']['signature'], cur.fetchall())
