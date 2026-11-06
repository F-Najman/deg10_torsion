from pathlib import Path
import hashlib,json
H=Path('work/deg10');D=H/'continuation_20261029';D.mkdir(exist_ok=True);(D/'deps').mkdir(exist_ok=True)
pairs=[('shared/fast_hecke.m','fast_hecke.m'),('shared/fast_hecke_11.m','fast_hecke_11.m'),
 ('noncyclic/deps/fast_full_aut_hash_m3_15.m','fast_full_aut_hash_m3_15.m'),
 ('noncyclic/deps/fast_full_aut_m3_21.m','fast_full_aut_m3_21.m'),
 ('noncyclic/deps/fast_full_aut_hash_m3_18_root.m','fast_full_aut_hash_m3_18_root.m'),
 ('noncyclic/deps/fast_full_aut_hash_m4_16_relative.m','fast_full_aut_hash_m4_16_relative.m'),
 ('noncyclic/deps/class_group_direct_RR_m3_18_root.m','class_group_direct_RR_m3_18_root.m'),
 ('noncyclic/deps/fast_diamond_hash_m2.m','fast_diamond_hash_m2.m'),
 ('noncyclic/model_6_12_root.m','model_6_12.m')]
manifest={}
for old,new in pairs:
 src=Path('work/degree8')/old;out=D/'deps'/new;data=src.read_bytes()
 data=data.replace(b'work/degree8/shared/mdmagma',b'work/deg10/shared/mdmagma')
 assert not out.exists();out.write_bytes(data)
 manifest[str(out.relative_to(H))]={'source':str(src),'source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'sha256':hashlib.sha256(data).hexdigest()}
(D/'dependency_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
header='''SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(5*10^9);t0:=Cputime();
AttachSpec("work/deg10/shared/mdmagma/v2/mdmagma.spec");
Attach("work/deg10/continuation_20261029/deps/fast_hecke.m");
Attach("work/deg10/continuation_20261029/deps/fast_hecke_11.m");
'''
for m,n,q,g in [(3,18,7,10),(4,16,5,13),(6,12,7,9),(3,21,4,25),(3,24,7,25),(4,24,5,33),(3,27,4,46),(4,20,9,25),(4,28,9,61)]:
 stem=f'cont_local_{m}_{n}_q{q}'
 text=header
 if m==6:text+='Attach("work/deg10/continuation_20261029/deps/model_6_12.m");\nX:=MDX612(GF(7));\n'
 else:text+=f'X:=MDX11({m},{n},GF({q}));\n'
 text+=f'''C:=Curve(X);assert Genus(C) eq {g};
cs:=Cusps(X);ds:=[Degree(c):c in cs];
print "DECIC_LOCAL_CUSPS",{m},{n},{q},ds,Cputime(t0);
ns:=[];
for d in [1..5] do xs:=NoncuspidalPlaces(X,d);Append(~ns,#xs);print "DECIC_LOCAL_COUNT",{m},{n},{q},d,#xs,Cputime(t0);end for;
J:=JH({m*n},[{n+1}]);assert Dimension(J) eq {g};
p:=Characteristic(GF({q}));e:=Degree(GF({q}));f:=FrobeniusPolynomial(J,p);assert Degree(f) eq {2*g};tr:=[];
for i in [1..5*e] do
 a:=-i*Coefficient(f,{2*g}-i);
 for j in [1..i-1] do a-:=Coefficient(f,{2*g}-j)*tr[i-j];end for;
 Append(~tr,Integers()!a);
end for;
cl:=[];
for d in [1..5] do
 a:={q}^d+1-tr[e*d];for j in Divisors(d) do if j lt d then a-:=j*cl[j];end if;end for;assert a mod d eq 0;Append(~cl,a div d);
end for;
assert cl eq [ns[d]+#[x:x in ds|x eq d]:d in [1..5]];
Write("work/deg10/continuation_20261029/{stem}_census.m",Sprint([*{m},{n},{q},ds,ns,Coefficients(f),cl*]));
print "DECIC_LOCAL_CENSUS_COMPLETE",{m},{n},{q},ns,Cputime(t0);
'''
 if g<=13:
  text+='''AF:=AlgorithmicFunctionField(FunctionField(C));c0:=[c:c in cs|Degree(c) eq 1][1];
G,_,phi:=ClassGroup(AF:ReductionDivisor:=FunctionFieldDivisor(Divisor(c0)),Proof:=true);
print "DECIC_LOCAL_CLASSGROUP_PROVEN",Invariants(G),Cputime(t0);
'''
 text+='quit;\n'
 (H/'computations'/(stem+'.m')).write_text(text)
source=(H/'computations/winding111.m').read_text()
for n in [66,68,69,70,72,74,76,78,80,81,84,87,90,93,99,100,105,135]:
 text=source.replace('N:=111;',f'N:={n};').replace('work/deg10/exclusions/krylov_', 'work/deg10/continuation_20261029/krylov_')
 (H/'computations'/f'cont_winding{n}.m').write_text(text)
print('PREPARED continuation dependencies, 9 local census jobs and 18 exact winding jobs')
