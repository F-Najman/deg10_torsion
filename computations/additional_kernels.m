SetColumns(0);SetQuitOnError(true);SetNthreads(1);SetMemoryLimit(4*10^9);
N:=175;a:=18;bad:={1,5};
checked:=0;badchecked:=0;for chi in GaloisConjugacyRepresentatives(FullDirichletGroup(N)) do
 if IsOdd(chi) or Order(chi) notin bad then continue;end if;
 S:=CuspidalSubspace(ModularSymbols(MinimalBaseRingCharacter(chi),2,1));w:=Dimension(WindingSubmodule(S));checked+:=1;
 if w ne Dimension(S) then assert Evaluate(chi,a) eq 1;badchecked+:=1;end if;end for;
assert badchecked eq 2;
printf "DECIC_ADDITIONAL_KERNEL N=%o a=%o diamond_order=12 potentially_bad_orbits=2 native_replay=true checked=%o\n",N,a,checked;
N:=185;a:=2;bad:={1,4};
checked:=0;badchecked:=0;for chi in GaloisConjugacyRepresentatives(FullDirichletGroup(N)) do
 if IsOdd(chi) or Order(chi) notin bad then continue;end if;
 S:=CuspidalSubspace(ModularSymbols(MinimalBaseRingCharacter(chi),2,1));w:=Dimension(WindingSubmodule(S));checked+:=1;
 if w ne Dimension(S) then assert Evaluate(chi,a) eq 1;badchecked+:=1;end if;end for;
assert badchecked eq 2;
printf "DECIC_ADDITIONAL_KERNEL N=%o a=%o diamond_order=18 potentially_bad_orbits=2 native_replay=true checked=%o\n",N,a,checked;
N:=259;a:=3;bad:={1};
checked:=0;badchecked:=0;for chi in GaloisConjugacyRepresentatives(FullDirichletGroup(N)) do
 if IsOdd(chi) or Order(chi) notin bad then continue;end if;
 S:=CuspidalSubspace(ModularSymbols(MinimalBaseRingCharacter(chi),2,1));w:=Dimension(WindingSubmodule(S));checked+:=1;
 if w ne Dimension(S) then assert Evaluate(chi,a) eq 1;badchecked+:=1;end if;end for;
assert badchecked eq 1;
printf "DECIC_ADDITIONAL_KERNEL N=%o a=%o diamond_order=9 potentially_bad_orbits=1 native_replay=true checked=%o\n",N,a,checked;
N:=407;a:=10;bad:={1,6};
checked:=0;badchecked:=0;for chi in GaloisConjugacyRepresentatives(FullDirichletGroup(N)) do
 if IsOdd(chi) or Order(chi) notin bad then continue;end if;
 S:=CuspidalSubspace(ModularSymbols(MinimalBaseRingCharacter(chi),2,1));w:=Dimension(WindingSubmodule(S));checked+:=1;
 if w ne Dimension(S) then assert Evaluate(chi,a) eq 1;badchecked+:=1;end if;end for;
assert badchecked eq 2;
printf "DECIC_ADDITIONAL_KERNEL N=%o a=%o diamond_order=6 potentially_bad_orbits=2 native_replay=true checked=%o\n",N,a,checked;
N:=481;a:=3;bad:={1};
checked:=0;badchecked:=0;for chi in GaloisConjugacyRepresentatives(FullDirichletGroup(N)) do
 if IsOdd(chi) or Order(chi) notin bad then continue;end if;
 S:=CuspidalSubspace(ModularSymbols(MinimalBaseRingCharacter(chi),2,1));w:=Dimension(WindingSubmodule(S));checked+:=1;
 if w ne Dimension(S) then assert Evaluate(chi,a) eq 1;badchecked+:=1;end if;end for;
assert badchecked eq 1;
printf "DECIC_ADDITIONAL_KERNEL N=%o a=%o diamond_order=18 potentially_bad_orbits=1 native_replay=true checked=%o\n",N,a,checked;
N:=629;a:=6;bad:={1,9};
checked:=0;badchecked:=0;for chi in GaloisConjugacyRepresentatives(FullDirichletGroup(N)) do
 if IsOdd(chi) or Order(chi) notin bad then continue;end if;
 S:=CuspidalSubspace(ModularSymbols(MinimalBaseRingCharacter(chi),2,1));w:=Dimension(WindingSubmodule(S));checked+:=1;
 if w ne Dimension(S) then assert Evaluate(chi,a) eq 1;badchecked+:=1;end if;end for;
assert badchecked eq 2;
printf "DECIC_ADDITIONAL_KERNEL N=%o a=%o diamond_order=16 potentially_bad_orbits=2 native_replay=true checked=%o\n",N,a,checked;
N:=703;a:=3;bad:={1};
checked:=0;badchecked:=0;for chi in GaloisConjugacyRepresentatives(FullDirichletGroup(N)) do
 if IsOdd(chi) or Order(chi) notin bad then continue;end if;
 S:=CuspidalSubspace(ModularSymbols(MinimalBaseRingCharacter(chi),2,1));w:=Dimension(WindingSubmodule(S));checked+:=1;
 if w ne Dimension(S) then assert Evaluate(chi,a) eq 1;badchecked+:=1;end if;end for;
assert badchecked eq 1;
printf "DECIC_ADDITIONAL_KERNEL N=%o a=%o diamond_order=9 potentially_bad_orbits=1 native_replay=true checked=%o\n",N,a,checked;
N:=851;a:=3;bad:={1};
checked:=0;badchecked:=0;for chi in GaloisConjugacyRepresentatives(FullDirichletGroup(N)) do
 if IsOdd(chi) or Order(chi) notin bad then continue;end if;
 S:=CuspidalSubspace(ModularSymbols(MinimalBaseRingCharacter(chi),2,1));w:=Dimension(WindingSubmodule(S));checked+:=1;
 if w ne Dimension(S) then assert Evaluate(chi,a) eq 1;badchecked+:=1;end if;end for;
assert badchecked eq 1;
printf "DECIC_ADDITIONAL_KERNEL N=%o a=%o diamond_order=198 potentially_bad_orbits=1 native_replay=true checked=%o\n",N,a,checked;
N:=1073;a:=3;bad:={1};
checked:=0;badchecked:=0;for chi in GaloisConjugacyRepresentatives(FullDirichletGroup(N)) do
 if IsOdd(chi) or Order(chi) notin bad then continue;end if;
 S:=CuspidalSubspace(ModularSymbols(MinimalBaseRingCharacter(chi),2,1));w:=Dimension(WindingSubmodule(S));checked+:=1;
 if w ne Dimension(S) then assert Evaluate(chi,a) eq 1;badchecked+:=1;end if;end for;
assert badchecked eq 1;
printf "DECIC_ADDITIONAL_KERNEL N=%o a=%o diamond_order=252 potentially_bad_orbits=1 native_replay=true checked=%o\n",N,a,checked;
N:=1147;a:=3;bad:={1};
checked:=0;badchecked:=0;for chi in GaloisConjugacyRepresentatives(FullDirichletGroup(N)) do
 if IsOdd(chi) or Order(chi) notin bad then continue;end if;
 S:=CuspidalSubspace(ModularSymbols(MinimalBaseRingCharacter(chi),2,1));w:=Dimension(WindingSubmodule(S));checked+:=1;
 if w ne Dimension(S) then assert Evaluate(chi,a) eq 1;badchecked+:=1;end if;end for;
assert badchecked eq 1;
printf "DECIC_ADDITIONAL_KERNEL N=%o a=%o diamond_order=45 potentially_bad_orbits=1 native_replay=true checked=%o\n",N,a,checked;
N:=1369;a:=6;bad:={1,3,9};
checked:=0;badchecked:=0;for chi in GaloisConjugacyRepresentatives(FullDirichletGroup(N)) do
 if IsOdd(chi) or Order(chi) notin bad then continue;end if;
 S:=CuspidalSubspace(ModularSymbols(MinimalBaseRingCharacter(chi),2,1));w:=Dimension(WindingSubmodule(S));checked+:=1;
 if w ne Dimension(S) then assert Evaluate(chi,a) eq 1;badchecked+:=1;end if;end for;
assert badchecked eq 3;
printf "DECIC_ADDITIONAL_KERNEL N=%o a=%o diamond_order=74 potentially_bad_orbits=3 native_replay=true checked=%o\n",N,a,checked;
print "DECIC_ADDITIONAL_KERNELS_COMPLETE";quit;
