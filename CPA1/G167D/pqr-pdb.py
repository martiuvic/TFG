import MDAnalysis as mda

propka_output="./nmes9byl54.pqr"
pdb_out="G167D_protonated.pdb"

system = mda.Universe(propka_output, format="PQR", in_memory=True)

## Remove the segment IDs (w/o this, `SYST` gets annoyingly printed by element)
for atom in system.atoms:
    atom.segment.segid = ''

system.atoms.write(pdb_out)
