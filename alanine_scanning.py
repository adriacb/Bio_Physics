#! /usr/bin/python3

"""
First Step: Checking the contact residues for each chain and generate a pdb file for each Chain and
    a TXT for each chain containing all the residues that belong to the interface.
Second Step: ΔGA-B = ΔGelectA-B + ΔGvdwA-B + ΔGSolvA-B - ΔGSolvA - ΔGSolvB
Third Step: Check the change of energies for each residue with ALANINE muthagenesis.

Use a pdb File with no HETAM, use Biobb_structure_checking for fixside and add_hydrogen functions.
For more information look at read_me file

USAGE: python3 res_dist.py 6axg_hyd_noHETAM.pdb1 6axg_hyd_noHETAM.pdb1 --c1 A --c2 B --maxdist 10.0
       python3 res_dist.py 2ki5_hyd_noHETAM.pdb1 2ki5_hyd_noHETAM.pdb1 --c1 A --c2 B --maxdist 10.0
"""
import argparse
import os
import sys
import Bio.PDB
import random
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import Selection
from Bio.PDB import PDBIO
import numpy as np
from Bio.PDB.NACCESS import NACCESS_atomic
from forcefield import VdwParamset
from residue_library import ResiduesDataLib
from Bio.PDB.Residue import Residue, DisorderedResidue
''' MAIN - Parsing the PDB FILE '''
parser = argparse.ArgumentParser (
    prog='res_dist',
    description='Getting CA pairs within max dist',
    usage='res_dist.py [options] pdb_file [> output_file]'
)

parser.add_argument(
    '--maxdist',
    action='store',
    dest='max_dist',
    default=20,
    type=float,
    help='Max contact distance (A)'
)

parser.add_argument('pdb_file', help='Input PDB 1', type=open)
parser.add_argument('pdb_file2', help='Input PDB 1', type=open)
parser.add_argument("--c1", help="First molecule chain", dest="c1")
parser.add_argument("--c2", help="Second molecule chain", dest="c2")
parser.add_argument(
    '--rlib',
    action='store',
    dest='reslib_file',
    default='data/aaLib.lib',
    help='Residue Library'
)
parser.add_argument(
    '--vdw',
    action='store',
    dest='vdwprm_file',
    default='data/vdwprm',
    help='Vdw parameters'
)

args = parser.parse_args()

for k, v in vars(args).items():
    print ('{:10}:'.format(k), v)

print("PDB.filename:", args.pdb_file.name)

parser = PDBParser(PERMISSIVE=1)

print ('Parsing', args.pdb_file)

out_folder = ""


print('Setting atom_type, charge and vdw for all the structure')
print()

def parsing(pdb_path):
    # Loading Libraries
    # loading residue library from data/aaLib.lib
    residue_library = ResiduesDataLib(args.reslib_file)

    # loading VdW parameters
    ff_params = VdwParamset(args.vdwprm_file)
    # load structure from PDB file
    parser = PDBParser(PERMISSIVE=1)
    st = parser.get_structure('STR', pdb_path)
    for at in st.get_atoms():
        resname = at.get_parent().get_resname()
        params = residue_library.get_params(resname, at.id)
        try:
            at.xtra['atom_type'] = params.at_type
            at.xtra['charge'] = params.charge
            at.xtra['vdw'] = ff_params.at_types[at.xtra['atom_type']]
        except:
            print("ERROR: residue/atom pair not in library (" + resname + ' ' + at.id + ')')
            #sys.exit(2)
        #if not params:
            #sys.exit("ERROR: residue/atom pair not in library (" + resname + ' ' + at.id + ')')

    srf = NACCESS_atomic(st[0],naccess_binary ='/Users/Adria/Desktop/BioPhy/soft/NACCESS/naccess' )
    return st


st = parsing(str(args.pdb_file.name))

''' Some code for getting contacts and generate a list of residues and atoms, also a pdb files for each chain '''

def is_contact(res_1, other_atoms, cutoff):

	for atom in res_1:
		ns = NeighborSearch(other_atoms)
		center = atom.get_coord()
		neighbors = ns.search(center, cutoff)
		residue_list = Selection.unfold_entities(neighbors, 'R') # R for residues
		if len(residue_list)>0:
			return True
	return False

def get_contacts(st, all_atoms, verbose, cutoff):
    progress = 0
    contacts = []

    for residue in st:
        progress+=1
        if len(verbose)>0:
            print(verbose,progress,"out of",len(st), residue)
        atom_list = st[residue]
        outcome = is_contact(atom_list, all_atoms, cutoff)
        if outcome:
            for at in atom_list:
                res = at.get_parent()
                #chain = res.get_parent()
                #at_id =  "{} {} {}".format (res, res.get_resname(), res.id[1], chain.id)
                if res not in contacts:
                    contacts.append(res)
    return contacts

def get_all_atoms(residue_map):
	all_atoms_out = []
	for residue in residue_map:
		for atom in residue_map[residue]:
			all_atoms_out.append(atom)
	return all_atoms_out

def get_atom_list(st, chains):

	output = dict()
	for chain in st:
		if chain.id in chains:
			for residue in chain.get_residues():
				hetflag, resseq, icode = residue.get_id()
				the_id = (chain.id+"_"+str(resseq)+"_"+icode).strip()
				for atom in residue.get_unpacked_list():
					if hetflag==' ':
						if the_id in output:
							output[the_id].append(atom)
						else:
							output[the_id] = [atom]
	return output

def save_residues(filename, contacts):
    f = open(filename,'w')
    for ele in contacts:
        f.write(str(ele)+"\n")
    f.close()

def save_contacts(structure, chains,out_file):
    '''Separated pdb files by chain'''
    Select = Bio.PDB.Select
    class ConstrSelect(Select):
        def accept_chain(self, chain):
            if chain.id in chains:
                return 1
            else:
                return 0
    w = PDBIO()
    w.set_structure(structure)
    randint = random.randint(0,9)
    w.save("TMP"+str(randint)+".pdb",ConstrSelect())
    #Remove the HETATM
    f_tmp = open("TMP"+str(randint)+".pdb", 'r')
    f_out = open(out_file, 'w')
    for line in f_tmp.readlines():
        if line[0:3]!="TER" and line[0:6]!="HETATM":
            f_out.write(line)
    f_tmp.close()
    f_out.close()
	#Delete the created temporary files
    os.remove("TMP"+str(randint)+".pdb")


cwd = os.getcwd()
os.system("mkdir "+cwd+"/"+out_folder)

str_1 = PDBParser().get_structure('first_st', args.pdb_file)
str_2 = PDBParser().get_structure('second_st', args.pdb_file2)

chains_1 = args.c1
chains_2 = args.c2

# Loading the structures
atoms_1 = Selection.unfold_entities(str_1, 'C')
atoms_2 = Selection.unfold_entities(str_2, 'C')
print('Chain {} : {}'.format(chains_1, len(chains_1)), 'Chain {} : {}'.format(chains_2, len(chains_2)))

# Returns the mapping from chain, residue id to the atoms list
input_1 = get_atom_list(atoms_1, chains_1)
input_2 = get_atom_list(atoms_2, chains_2)

all_atoms_1 = get_all_atoms(input_1)
all_atoms_2 = get_all_atoms(input_2)

contacts_1 = get_contacts(input_1, all_atoms_2, "First chain, res. ", args.max_dist)
contacts_2 = get_contacts(input_2, all_atoms_1, "Second chain, res. ", args.max_dist)

print("# of residues of chain {}: {}".format(chains_1,len(contacts_1)), "# of residues of chain {}: {}".format(chains_2,len(contacts_2)))

#save the interface residues
save_residues(cwd+"/"+out_folder+"/chain_{}_{}.txt".format(chains_1, args.pdb_file.name), contacts_1)
save_residues(cwd+"/"+out_folder+"/chain_{}_{}.txt".format(chains_1, args.pdb_file.name), contacts_2)

#Save the chains
save_contacts(str_1,chains_1,cwd+"/"+out_folder+"/chain_{}_{}.pdb".format(chains_1, args.pdb_file.name))
save_contacts(str_2,chains_2,cwd+"/"+out_folder+"/chain_{}_{}.pdb".format(chains_2, args.pdb_file.name))

'''                                 COMPUTING ENERGIES                                          '''
''' DEFINING FUNCTIONS '''

def solvation(at):
    ''' this function computes the solvation energy for one atom '''
    try:
        if at.element == 'H':# or 'EXP_NACCESS' not in at.xtra or 'vdw' not in at.xtra:
            return 0.
        sigma = at.xtra['vdw'].fsrf
        surface = float(at.xtra['EXP_NACCESS'])
    except:
        return 0
    return sigma * surface

def mehler_solmajer(dist): # This is for computing the dielectric constant using the Mehler-Solmajer approximation
    return ((86.9525)/(1-7.7839 * np.exp(-0.3153 * dist))) - 8.5525


def calc_vdw2(at, at2):
    '''returns Van der Waals interactions of one atom'''
    total = 0.
    try:
        e1 = at.xtra['vdw'].eps
        s1 = at.xtra['vdw'].sig
        d = at2 - at
        if d > 2:
            e2 = at2.xtra['vdw'].eps
            s2 = at2.xtra['vdw'].sig
            epsilon = np.sqrt(e1 * e2)
            sigma = np.sqrt(s1 * s2)
            total += 4*epsilon*((sigma/d)**12 - (sigma/d)**6)
    except:
        total += 0
    return total


def calc_elec_int2(at,at2):
    ''' returns electrostatic interactions of one atom '''
    total = 0.
    try:
        dist = at2-at
        e_r = mehler_solmajer(dist)
        total += 332.16*((at.xtra['charge'] * at2.xtra['charge']) / (e_r * dist))
    except:
        total += 1
    return total


''' PARSING STRUCTURES WITH ELEMENS OF CHAIN A AND B '''

print('Parsing chain A')

st_chA = parsing("chain_{}_{}.pdb".format(chains_1, args.pdb_file.name))

print('Parsing chain B')
st_chB = parsing("chain_{}_{}.pdb".format(chains_2, args.pdb_file.name))

atoms_alanine = ['N', 'H', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'HB3', 'C', 'O']


def delta_G_AB(v_d_w, electrostatic, solvationAB, solvation_A, solvation_B):
    return v_d_w + electrostatic + solvationAB - solvation_A - solvation_B


residues_AB = list()

for res1 in contacts_1:
    residues_AB.append(res1)
for res2 in contacts_2:
    residues_AB.append(res2)


print()
def alanine_scanning(st, sta, stb):
    print('------------ Computing Energies - Alanine Scanning ------------')
    print('ΔGA-B = ΔGelectA-B + ΔGvdwA-B + ΔGSolvA-B - ΔGSolvA - ΔGSolvB')
    energies = dict()

    #Energies ALANINE
    res_vdw_ala1 = 0.0
    res_elect_ala1 = 0.0
    res_sol_ala1 = 0.0
    res_vdw_ala2 = 0.0
    res_elect_ala2 = 0.0
    res_sol_ala2 = 0.0


    vdw_T = 0.0
    elect_T = 0.0
    solv_T = 0.0
    sol_AT = 0.0
    sol_BT = 0.0
    dic_chains = dict()
    first = True
    first2 = True
    print('Computing v_d_w and electrostatic energies of the residues from Chain A that belongs to the interface')
    for r1 in sta.get_residues():
        res_vdw = 0.0
        res_ele = 0.0
        if r1 in contacts_1:
            for r2 in stb.get_residues():
                if r1 != r2 and r2 in contacts_2:
                    for a1 in r1.get_atoms():
                        if a1.xtra:
                            for a2 in r2.get_atoms():
                                if a2.xtra:
                                    res_vdw += calc_vdw2(a1, a2)
                                    res_ele += calc_elec_int2(a1, a2)
                                    if first:
                                        if r1.resname == 'ALA' and a1 in atoms_alanine:
                                            res_vdw_ala1 += calc_vdw2(a2, a1)
                                            res_elect_ala1 += calc_elec_int2(a2, a1)
                                            first = False
            if r1.resname != 'ALA':
                dic_chains[r1] = [[res_vdw, res_ele, 0, 0], [res_vdw_ala1,res_elect_ala1,res_sol_ala1,res_sol_ala1]]
            elif r1.resname == 'ALA':
                dic_chains[r1] = [[res_vdw, res_ele, 0, 0], [0,0,0,0]]
        vdw_T += res_vdw
        elect_T += res_ele

    print('Computing v_d_w and electrostatic energies of the residues from Chain B that belongs to the interface')
    for r2 in st_chB.get_residues():
        res_vdw = 0.0
        res_ele = 0.0
        if r2 in contacts_2:
            for r1 in sta.get_residues():
                if r1 != r2 and r1 in contacts_1:
                    for a2 in r2.get_atoms():
                        if a2.xtra:
                            for a1 in r1.get_atoms():
                                if a1.xtra:
                                    res_vdw += calc_vdw2(a2, a1)
                                    res_ele += calc_elec_int2(a2, a1)
                                    if first2:
                                        if r2.resname == 'ALA' and a2 in atoms_alanine:
                                            res_vdw_ala2 += calc_vdw2(a1, a2)
                                            res_elect_ala2 += calc_elec_int2(a1, a2)
                                            first2 = False
            if r2.resname != 'ALA':
                dic_chains[r2] = [[res_vdw, res_ele, 0, 0], [res_vdw_ala2,res_elect_ala2,res_sol_ala2,res_sol_ala2]]
            elif r2.resname == 'ALA':
                dic_chains[r2] = [[res_vdw, res_ele, 0, 0], [0,0,0,0]]
        vdw_T += res_vdw
        elect_T += res_ele
    f3 = True
    f4 = True
    for res in st.get_residues(): #In this for loop I computed the solvation energies
        if res in residues_AB: # if this residue belongs to the interface
            res_sol = 0.0
            res_sol_A = 0.0
            res_sol_B = 0.0
            for at in res.get_atoms():
                res_sol += solvation(at)
                if res in contacts_1:
                    res_sol_A += solvation(at)
                    if f3:
                        res_sol_ala1 += solvation(at)
                        f3 = False
                if res in contacts_2:
                    res_sol_B += solvation(at)
                    if f4:
                        res_sol_ala2 += solvation(at)
                        f4 = False
            solv_T += res_sol
            sol_AT += res_sol_A
            sol_BT += res_sol_B
            if res in dic_chains:
                dic_chains[res][0][2] = res_sol
                if res in contacts_1:
                    dic_chains[res][0][3] = res_sol_A
                    if res.resname != 'ALA':
                        dic_chains[res][1][2] = res_sol_ala1
                        dic_chains[res][1][3] = res_sol_ala1
                elif res in contacts_2:
                    dic_chains[res][0][3] = res_sol_B
                    if res.resname != 'ALA':
                        dic_chains[res][1][2] = res_sol_ala1
                        dic_chains[res][1][3] = res_sol_ala1

    detla_g = delta_G_AB(vdw_T, elect_T, solv_T, sol_AT, sol_BT)
    print("ΔGVan der Waals A-B: {} \nΔGElect A-B: {} \nΔGSolvA-B: {}\nΔGSolvA: {} \nΔGSolvB: {} \nΔG A-B: {} ".format(vdw_T, elect_T, solv_T, sol_AT, sol_BT, detla_g))

    return dic_chains


residues_alanine_scan = alanine_scanning(st, st_chA, st_chB)

def compute_differences(residues_alanine_scan):
    dicto = dict()
    for res in residues_alanine_scan:
        whole = 0.0
        differences = 0.0
        differences += (residues_alanine_scan[res][0][0] - residues_alanine_scan[res][1][0]) + (residues_alanine_scan[res][0][1] - residues_alanine_scan[res][1][1]) + (residues_alanine_scan[res][0][2] - residues_alanine_scan[res][1][2]) + (residues_alanine_scan[res][0][3] - residues_alanine_scan[res][1][3])
        whole += (residues_alanine_scan[res][0][0]) + (residues_alanine_scan[res][0][1]) + (residues_alanine_scan[res][0][2]) + (residues_alanine_scan[res][0][3])
        dicto[res] = [whole, differences]
    return dicto

diff = compute_differences(residues_alanine_scan)

total_dif = list()
total_diff = list()
range_seq = list()
for residue in diff:
    total_dif.append(diff[residue][0])
    total_diff.append(diff[residue][1])
    range_seq.append(residue.id[1])
    print(residue, diff[residue][1])
################
# Make the plots
import numpy
import matplotlib.pyplot as plt


seq = numpy.arange(0, len(total_diff))

plt.plot(seq,total_diff)
plt.title('All energies Folded State')
plt.ylabel('Delta Delta G')
plt.xlabel('Sequence position')
plt.show()
