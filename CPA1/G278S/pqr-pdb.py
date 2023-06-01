import MDAnalysis as mda

propka_output="./c9xs0s622d.pqr"
pdb_out="G278S_protonated.pdb"

system = mda.Universe(propka_output, format="PQR", in_memory=True)

## Remove the segment IDs (w/o this, `SYST` gets annoyingly printed by element)
for atom in system.atoms:
    atom.segment.segid = ''

system.atoms.write(pdb_out)
