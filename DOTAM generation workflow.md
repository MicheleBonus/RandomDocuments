## Creating the DOX residue

The **DOX** residue was taken from the file `Mol-ia4_m1-c1.mol2` from PyRED Job `P2260`. However, the `@<TRIPOS>HEADTAIL` and `@<TRIPOS>RESIDUECONNECT` sections in the file were changed from:

```
@<TRIPOS>HEADTAIL
C1 1
C3 1
@<TRIPOS>RESIDUECONNECT
1 C1 C3 C5 0 0 0
```

to

```
@<TRIPOS>HEADTAIL
C1 1
C5 1
@<TRIPOS>RESIDUECONNECT
1 C1 C5 C3 0 0 0
```

and the file was saved as `DOX.mol2`. This ensures that **C5** is the tail atom that connects to the head atom of the next residue. This was changed primarily for aesthetic reasons, so that the carbon with the highest count is the connecting atom.

Another key feature of the **DOX** residue is its 'connect2' entry for its **C3** atom, allowing for additional connectivity beyond the standard head and tail atoms.

In tLEaP, we can then just **add** the residue to the existing `DC_SIGN_Linkers.off` library. Since the _gaff2_ parameters and the additional _frcmod_ for the linkers should cover the atom types, bonds, angles, and dihedrals in **DOX**, we won't need an additional frcmod file. However, since we are naturally cautious, we will perform _parmchk_ analysis on the entire molecule once more at the very end as an extra precaution:

```
source leaprc.gaff2

loadOff DC_SIGN_Linkers.off
loadAmberParams DC_SIGN_Linkers.frcmod

DOX = loadMol3 DOX.mol2

check DOX

saveOff DOX DC_SIGN_Linkers.off

quit
```

## Generating the residues for the DOTAM structure

The molecular structure we are first going to assemble consists of a single chain that can be viewed in two distinct segments: a forward part and a backward part.

![[DOTAM_Plan.svg]]

The forward part of the chain is composed of a sequence of parameterized residues arranged in the following order: **MAN-MMA-MEA-MEO-HEA-DMS-EDA-DOX**. In this arrangement, the head and tail atoms of each residue are set to ensure proper connectivity. Notably, the first residue in this sequence, **MAN**, lacks a head atom, as it represents the starting point of the chain and thus only contains a tail atom. This tail atom connects to the head atom of the subsequent residue, **MMA**. Moving along the chain, each residue, from **MMA** through to **DOX**, is structured with defined head and tail atoms, facilitating correct attachment to adjacent residues in the sequence.

However, constructing the backward part of the chain poses a unique challenge. A straightforward inversion of the residue sequence from the forward chain would not work due to the specific head and tail atom configurations. This issue becomes evident with the **DOX** residue. In **DOX**, the **C1** atom is designated as the head, and the **C5** atom as the tail. If another **DOX** were to be placed in sequence, as per the inverted order, the **C5** of the first **DOX** would incorrectly connect to the **C1** of the second **DOX**, whereas the proper connection should be between two **C5** atoms. To resolve this, the backward part of the chain requires a new set of residues that have the same topology as the original residues but contain swapped head and tail atoms. These new residues will be named **XOD, ADE, SMD, AEH, OEM, AEM, AMM**, and **NAM**.

To generate these residues, we will identify the head and tail atoms of the "foward"  residues in the linkers. To accomplish this, we will use the `desc` command for each residue after loading the library into _tLEaP_.

```
source leaprc.gaff2

loadOff DC_SIGN_Linkers.off
loadAmberParams DC_SIGN_Linkers.frcmod

desc MAN
desc MMA
desc MEA
desc MEO
desc HEA
desc DMS
desc EDA
desc ECA
desc DOX

quit
```

The output of the `desc` commands should look somewhat like this:

```
UNIT name: MAN
Head atom: null
Tail atom: .R<MAN 1>.A<C1 1>

UNIT name: MMA
Head atom: .R<MMA 1>.A<O1 1>
Tail atom: .R<MMA 1>.A<C3 5>

UNIT name: MEA
Head atom: .R<MEA 1>.A<N1 1>
Tail atom: .R<MEA 1>.A<O1 9>

UNIT name: MEO
Head atom: .R<MEO 1>.A<C1 1>
Tail atom: .R<MEO 1>.A<O2 7>

UNIT name: HEA
Head atom: .R<HEA 1>.A<C1 1>
Tail atom: .R<HEA 1>.A<N1 7>

UNIT name: DMS
Head atom: .R<DMS 1>.A<C1 1>
Tail atom: .R<DMS 1>.A<C4 9>

UNIT name: EDA
Head atom: .R<EDA 1>.A<N1 1>
Tail atom: .R<EDA 1>.A<N2 9>

UNIT name: ECA
Head atom: .R<ECA 1>.A<C 5>
Tail atom: null

UNIT name: DOX
Head atom: .R<DOX 1>.A<C1 1>
Tail atom: .R<DOX 1>.A<C5 10>
```

We will use this information as a guide to generate the new set of residues. In general, the generation of these residues is entirely possible in tLEaP and consists of four parts:

1. Generating copies of the original residue Units in Variables with the respective name of the corresponding "backward" residue.
2. Setting the **unit names** inside the variables to the respective name of the corresponding "backward" residue.
3. Setting the **residue names** inside the respective residues in the units to the corresponding "backward" residue
4. Swapping the head and tail atoms in the units and the connect0 and connect1 atoms in the residues.

```
source leaprc.gaff2

loadOff DC_SIGN_Linkers.off
loadAmberParams DC_SIGN_Linkers.frcmod

NAM = copy MAN
AMM = copy MMA
AEM = copy MEA
OEM = copy MEO
AEH = copy HEA
SMD = copy DMS
ADE = copy EDA
XOD = copy DOX

# Correct the names of the copied units
set NAM name "NAM"
set AMM name "AMM"
set AEM name "AEM"
set OEM name "OEM"
set AEH name "AEH"
set SMD name "SMD"
set ADE name "ADE"
set XOD name "XOD"

# Also correct the residue names inside the renamed units
set NAM.1 name "NAM"
set AMM.1 name "AMM"
set AEM.1 name "AEM"
set OEM.1 name "OEM"
set AEH.1 name "AEH"
set SMD.1 name "SMD"
set ADE.1 name "ADE"
set XOD.1 name "XOD"

# Now swap the head and tail atoms in the units and the residues

# Swap head and tail for NAM unit and connect0 and connect1 for NAM residue
# Note: head was null in MAN, so tail is null in NAM
set NAM head NAM.1.C1
set NAM tail null
set NAM.1 connect0 NAM.1.C1
set NAM.1 connect1 null

# Swap head and tail for AMM unit and connect0 and connect1 for AMM residue
set AMM head AMM.1.C3
set AMM tail AMM.1.O1
set AMM.1 connect0 AMM.1.C3
set AMM.1 connect1 AMM.1.C3

# Swap head and tail for AEM unit and connect0 and connect1 for AEM residue
set AEM head AEM.1.O1
set AEM tail AEM.1.N1
set AEM.1 connect0 AEM.1.O1
set AEM.1 connect1 AEM.1.N1

# Swap head and tail for OEM unit and connect0 and connect1 for OEM residue
set OEM head OEM.1.O2
set OEM tail OEM.1.C1
set OEM.1 connect0 OEM.1.O2
set OEM.1 connect1 OEM.1.C1

# Swap head and tail for AEH unit and connect0 and connect1 for AEH residue
set AEH head AEH.1.N1
set AEH tail AEH.1.C1
set AEH.1 connect0 AEH.1.N1
set AEH.1 connect1 AEH.1.C1

# Swap head and tail for SMD unit and connect0 and connect1 for SMD residue
set SMD head SMD.1.C4
set SMD tail SMD.1.C1
set SMD.1 connect0 SMD.1.C4
set SMD.1 connect1 SMD.1.C1

# Swap head and tail for ADE unit and connect0 and connect1 for ADE residue
set ADE head ADE.1.N2
set ADE tail ADE.1.N1
set ADE.1 connect0 ADE.1.N2
set ADE.1 connect1 ADE.1.N1

# Swap head and tail for XOD unit and connect0 and connect1 for XOD residue
# connect2 remains C3, as set in the Mol3 file
set XOD head XOD.1.C5
set XOD tail XOD.1.C1
set XOD.1 connect0 XOD.1.C5
set XOD.1 connect1 XOD.1.C1

saveOff NAM DC_SIGN_Linkers.off
saveOff AMM DC_SIGN_Linkers.off
saveOff AEM DC_SIGN_Linkers.off
saveOff OEM DC_SIGN_Linkers.off
saveOff AEH DC_SIGN_Linkers.off
saveOff SMD DC_SIGN_Linkers.off
saveOff ADE DC_SIGN_Linkers.off
saveOff XOD DC_SIGN_Linkers.off

quit
```

## Constructing a fully symmetric DOTAM structure

Even with all the residues available, constructing the complete structure in tLEaP is difficult: some dihedral angles will be incorrect, particularly those within the macrocycle, even when constructing only half of the molecule.

Do generate an initial full DOTAM structure, the `Mol-ia5_m1-c1.mol2` file from PyRED Job `P2241` was used as a guideline. Linear linker structures were constructed separately in `tLEaP` and then superimposed onto the corresponding parts in `Mol-ia5_m1-c1.mol2`.

If one is lucky, and the constructed symmetry is "close enough", one can then symmetrize the structure to C4 in Maestro or GaussView.

After manually renaming all the residues to their corresponding force field equivalents (those we constructed above), renumbering the residues, and correcting for atom name mismatches, one may get a structure like `DOT.pdb`.

## Reintroducing some convenience

Obviously, for each possible linker, this tedious process would have to be repeated. To circumvent this issue, we will set the dihedral angles _while_ constructing the half-molecule, using the manually constructed structure as a guideline.

```
source leaprc.gaff2
loadOff DC_SIGN_Linkers.off
loadAmberParams DC_SIGN_Linkers.frcmod

chain1 = sequence { MAN MMA MEA MEO HEA DMS EDA DOX XOD ADE SMD AEH OEM AEM AMM NAM }

impose chain1 { 1 2 } { { "H2" "C1" "O1" "C2" -60.0 } }
impose chain1 { 2 3 } { { "C2" "C3" "N1" "C3" 180.00 } }
impose chain1 { 3 4 } { { "C4" "O1" "C1" "C2" 180.00 } }
impose chain1 { 4 5 } { { "C2" "O2" "C1" "C2" 180.00 } }
impose chain1 { 5 6 } { { "C2" "N1" "C1" "C2" 180.00 } }
impose chain1 { 6 7 } { { "C3" "C4" "N1" "C3" 180.00 } }

impose chain1 { 7 8 } { { "N2" "C1" "C2" "N2" -169.83066 } }
impose chain1 { 7 8 } { { "C4" "N2" "C1" "C2" 180.00000 } }
impose chain1 { 8 } { { "C1" "C2" "N2" "C3" -146.59512 } }
impose chain1 { 8 } { { "C1" "C2" "N2" "C5" 89.82430 } }
impose chain1 { 8 9 } { { "C2" "N2" "C5" "C5" -76.21161 } }
impose chain1 { 8 9 } { { "N2" "C5" "C5" "N2" -61.19138 } }
impose chain1 { 8 9 } { { "C5" "C5" "N2" "C2" 159.97023 } }
impose chain1 { 8 9 } { { "C5" "C5" "N2" "C3" -76.0422 } }
impose chain1 { 8 9 } { { "C3", "N2", "C5", "C5" 160.28267 } }
impose chain1 { 9 } { { "C3" "N2" "C2" "C1" 89.82430 } }
impose chain1 { 9 } { { "C5" "N2" "C2" "C1" -146.59512 } }
impose chain1 { 9 10 } { { "N2" "C2" "C1" "N2" -169.83066 } }
impose chain1 { 9 10 } { { "C2" "C1" "N2" "C4" 180.00000 } }

impose chain1 { 10 11 } { { "C3" "N1" "C4" "C3" 180.00000 } }
impose chain1 { 11 12 } { { "C2" "C1" "N1" "C2" 180.00000 } }
impose chain1 { 12 13 } { { "C2" "C1" "O2" "C2" 180.00000 } }
impose chain1 { 13 14 } { { "C2" "C1" "O1" "C4" 180.00000 } }
impose chain1 { 14 15 } { { "C3" "N1" "C3" "C2" 180.00000 } }
impose chain1 { 15 16 } { { "C2" "O1" "C1" "H2" -60.0 } }
```

And now I don't know what to do :D