#1 download the pdb
wget https://files.rcsb.org/view/2F8A.pdb -O raw.pdb

#2 running PISA https://www.ebi.ac.uk/msd-srv/prot_int/cgi-bin/piserver (to get the proper biological unit of the structures)
oligomer.pdb 

#3 Clean PDB on the oligomer.pdb (it is done to erase all data that is not important to run simulation and can confuse it.)
sh clean-pdb.sh 
#TER lines should be between the protein and any metals, ions, or waters

#3 check pKa of the amino acids 
propka3 oligomer.pdb

#4  Predict the proper protonation states for each residue. The pqr file contains the adjusted residue names.  https://server.poissonboltzmann.org/
60y0z5k439.pqr  
#convert the pqr file to pdb file by the script file name - pqr-pdb.py

#5 Molprobity -  http://molprobity.biochem.duke.edu/  

#6 LEAP

#Remove the flip and new marks that MolProbity prints
#Remove USER and REMARK from PROPKA and MolProbity
vi <pdb file coming from molprobity>
#:%s/  flip//g
#:%s/  new//g
#:%s/ OW  WAT/ O   WAT/g

#7 tLEAP
tleap -f tleap.in

Input file

source leaprc.gaff - general amber force field
source leaprc.protein.ff14SB - specific protein parameter
source leaprc.water.tip3p
c = loadpdb 7fc2_ph7_protonatedFH.pdb
savepdb c 7fc2_vac.pdb
saveamberparm c 7fc2_vac.prmtop  - vacuum parameter and topology  7fc2_vac.inpcrd - vacuum coordinate file
addions c CL 0.0
addions c NA 0.0
solvatebox c TIP3PBOX 12.0
savepdb c 7fc2_wat.pdb
saveamberparm c 7fc2_wat.prmtop  - water parameter and topology  7fc2_wat_init0.rst - water coordinate file
quit
