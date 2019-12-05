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

## forcefield.py
Module for forcefield parameters management

## residue_library.py
Module for Residue Library management

## data/vdwprm
Simple vdw parameters (based on AMBER ff)

## data/aaLib.lib
Library for obtaining amino acid atom types and partial charges (based on AMBER ff)

## Alanine_scanning.py


	$ python3 alanine_scanning.py pdb.pdb1 pdb.pdb1 --c1 A --c2 B --maxdist 8.0

Where pdb.pdb1 and pdb.pdb1 is the same biological unit, without heteroatoms, and fixed with biobb_structure_checking
and with hydrogen added with add_hydrogen

External dependencies

    Bio.PDB.NeighborSearch (BioPython)
    Bio.PDB.PDBParser (Biopython)

## 1. Objective:

To evaluate the relative contribution of interface residues to the interaction energy in a protein-protein complex.

## 2. Strategy:

- Determine amino acid residues that form the interface between the complex components
- Determine the contribution to the stability of the complex by mimicking an Ala-scanning experiment, i.e. replacing each residue in turn by Ala and evaluating the changes in the inter-complex interaction energy.
- Compare the results obtained between two types of complexes, a permanent interface (2ki5), and a temporary one (6axg)

## 3. Preparation:

- Obtain the required structures from the PDB.
- Check at PDB which is the composition of a “Biological unit”. Remove all chains but those involved in the biological unit, if necessary
- Remove all heteroatoms
- Perform a quality checking on the structures, and add missing side-chains and hydrogen atoms.


