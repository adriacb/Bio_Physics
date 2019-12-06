# BioPhysics
Energy analysis exercise 2019-20 - Biophysics 	           
Protein-protein interface analysis            	                          
Adrià Cabello Blanque                                                         

## Requirements 
- Python >= 3.6
- Biopython module (>= 1.72)
- Molecular viewer (pymol, chimera)
- biobb_structure_checking python module
- Anaconda
	- Create a new environment: conda create -n name_env
    - Install Biopython: conda install -c anaconda biopython
	- Install Biobb_Structure_checking: conda install -c bioconda biobb_structure_checking
	- Install Numpy: conda install -c anaconda numpy
	- conda config --add channels bioconda
	- Activate environment: source activate name_env

• Clone this git-hub link in your work Folder: $ git clone https://github.com/jlgelpi/Biophysics.git

See also:
- https://github.com/jlgelpi/Biophysics/tree/master/Examples
- http://biopython.org/DIST/docs/tutorial/Tutorial.html
- https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ
- http://biopython.org/DIST/docs/api/

## Biobb_structure_checking

	$ check_structure fixed
	$ chec_structure add_hydrogen

## Alanine_scanning.py


	$ python3 alanine_scanning.py pdb.pdb1 pdb.pdb1 --c1 A --c2 B --maxdist 8.0

Where pdb.pdb1 and pdb.pdb1 is the same biological unit, without heteroatoms, and fixed with biobb_structure_checking
and with hydrogen added with add_hydrogen

External dependencies

    Bio.PDB.NeighborSearch (BioPython)
    Bio.PDB.PDBParser (Biopython)

## Output_files

	2ki5log.gdoc	The general output of the script for 2ki5
	6axg_log.gdoc	The general output of the script for 6axg

	2ki5.pse	Pymol
	6axg.pse	Pymol
	chain_A_2ki5_opt.pdb.pdb	PDB Chain A 2ki5
	chain_A_2ki5_opt.pdb.txt	Output residues (Contacts 1)
	chain_A_6axg_opt.pdb.pdb	PDB Chain A 6axg
	chain_A_6axg_opt.pdb.txt	Output residues (Contacts 2)
	chain_B_2ki5_opt.pdb.pdb	PDB Chain B 2ki5
	chain_B_6axg_opt.pdb.pdb	PDB Chain B 6axg
