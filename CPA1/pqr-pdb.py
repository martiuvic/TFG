import MDAnalysis as mda

propka_output="./gvp6g8gm12.pqr"
pdb_out="5OM9_R122H_protonated.pdb"

system = mda.Universe(propka_output, format="PQR", in_memory=True)

## Remove the segment IDs (w/o this, `SYST` gets annoyingly printed by element)
for atom in system.atoms:
    atom.segment.segid = ''

system.atoms.write(pdb_out)
