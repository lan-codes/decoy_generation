# -*- mode: python; tab-width: 4 -*- 

## 
# gen_decoys.py 
#
# Copyright (C) 2024 Thomas A. Seidel <thomas.seidel@univie.ac.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
##

import sys
import io
import pathlib
import argparse

import CDPL.Chem as Chem
import CDPL.MolProp as MolProp
import CDPL.Descr as Descr
import CDPL.Util as Util
import CDPL.Base as Base


MAX_POOL_SIZE           = 10000
ECFP4_LENGTH            = 6007
MAX_DECOY_INPUT_MOL_SIM = 0.35
MAX_DECOY_DECOY_SIM     = 0.80


class MolPropData():

    def __init__(self):
        pass

    def matches(self, props: object, strict_match: bool, params: argparse.Namespace) -> bool:
        if abs(self.weight - props.weight) > params.mw_tol:
            return False
 
        if abs(self.logP - props.logP) > params.logp_tol:
            return False

        if abs(self.netCharge - props.netCharge) > params.charge_tol:
            return False

        if not strict_match:
            if abs(self.numRotBonds - props.numRotBonds) > (params.rbc_tol + 1):
                return False

            if abs(self.numHBA - props.numHBA) > (params.hba_tol + 1):
                return False

            return (abs(self.numHBD - props.numHBD) <= (params.hbd_tol + 1))
        
        if abs(self.numRotBonds - props.numRotBonds) > params.rbc_tol:
            return False

        if abs(self.numHBA - props.numHBA) > params.hba_tol:
            return False

        return (abs(self.numHBD - props.numHBD) <= params.hbd_tol)

    def calculate(self, mol: Chem.Molecule, mol_idx: int) -> object:
        Chem.initSubstructureSearchTarget(mol, False)

        self.index = mol_idx
        self.name = Chem.getName(mol)
        self.weight = MolProp.calcMass(mol)
        self.logP = MolProp.calcXLogP(mol)
        self.netCharge = MolProp.getNetFormalCharge(mol)
        self.numRotBonds = MolProp.getRotatableBondCount(mol)
        self.numHBA = MolProp.getHBondAcceptorAtomCount(mol)
        self.numHBD = MolProp.getHBondDonorAtomCount(mol)
        self.hashCode = Chem.calcHashCode(mol, Chem.AtomPropertyFlag.TYPE,
                                          Chem.BondPropertyFlag.ORDER | Chem.BondPropertyFlag.AROMATICITY)
        self.molecule = mol

    def read(self, f: io.TextIOBase, skip_header: bool) -> list:
        if skip_header:
             f.readline()

        line = f.readline()

        if not line:
             return None
             
        fields = line.strip().split(',')

        self.index = int(fields[0])
        self.weight = float(fields[1])
        self.logP = float(fields[2])
        self.netCharge = int(fields[3])
        self.numRotBonds = int(fields[4])
        self.numHBA = int(fields[5])
        self.numHBD = int(fields[6])
        self.hashCode = int(fields[7])

        return self

    def write(self, f: io.TextIOBase, write_header: bool) -> None:
        if write_header:
            f.write('Mol. Index,MW,XLogP,Net Charge,#Rot. Bonds,#HBA,#HBD,Hash Code\n')

        f.write(f'{self.index},{self.weight},{self.logP},{self.netCharge},{self.numRotBonds},{self.numHBA},{self.numHBD},{self.hashCode}\n')

    
def parseArguments():
    parser = argparse.ArgumentParser(description='Retrieves decoys for a given set of input molecules.')
    
    parser.add_argument('-i',
                        dest='input',
                        required=False,
                        help='[Required/Optional] File containing the input molecules (e.g. actives).',
                        default=None,
                        metavar='<path>')
    parser.add_argument('-o',
                        dest='output',
                        required=False,
                        help='[Required/Optional] Output file for the retrieved decoy molecules.',
                        default=None,
                        metavar='<path>')
    parser.add_argument('-d',
                        dest='decoy_db',
                        required=True,
                        help='[Required] File providing the decoy source molecules.',
                        metavar='<path>')
    parser.add_argument('-p',
                        dest='decoy_db_props',
                        required=True,
                        help='[Required] CSV-file providing the precalculated decoy database molecule properties.',
                        metavar='<path>')
    parser.add_argument('-n',
                        dest='max_num_decoys',
                        required=False,
                        help='[Optional] Maximum number of decoys per input molecule (default: 50).',
                        metavar='<int>',
                        type=int,
                        default=50)
    parser.add_argument('-m',
                        dest='min_num_decoys',
                        required=False,
                        help='[Optional] Minimum number of decoys per input molecule (default: 20).',
                        metavar='<int>',
                        type=int,
                        default=20)
    parser.add_argument('-x',
                        dest='excl_mols',
                        required=False,
                        help='[Optional] File with molecules that shall not be contained in the output decoy molecule list (e.g. all known actives).',
                        default=None,
                        metavar='<path>')
    parser.add_argument('--mw-tol',
                        dest='mw_tol',
                        required=False,
                        help='[Optional] Moleular weight matching tolerance (default: 20.0).',
                        metavar='<float>',
                        type=float,
                        default=20.0)
    parser.add_argument('--charge-tol',
                        dest='charge_tol',
                        required=False,
                        help='[Optional] Net charge matching tolerance (default: 0).',
                        metavar='<int>',
                        type=int,
                        default=0)
    parser.add_argument('--logp-tol',
                        dest='logp_tol',
                        required=False,
                        help='[Optional] LogP matching tolerance (default: 0.75).',
                        metavar='<float>',
                        type=float,
                        default=0.75)
    parser.add_argument('--rbc-tol',
                        dest='rbc_tol',
                        required=False,
                        help='[Optional] Rotatable bond count matching tolerance (default: 0).',
                        metavar='<int>',
                        type=int,
                        default=0)
    parser.add_argument('--hba-tol',
                        dest='hba_tol',
                        required=False,
                        help='[Optional] H-bond acceptor atom count matching tolerance (default: 0).',
                        metavar='<int>',
                        type=int,
                        default=0)
    parser.add_argument('--hbd-tol',
                        dest='hbd_tol',
                        required=False,
                        help='[Optional] H-bond donor atom count matching tolerance (default: 0).',
                        metavar='<int>',
                        type=int,
                        default=0)
    parser.add_argument('--calc-props',
                        dest='calc_props',
                        required=False,
                        help='[Optional] Calculate and store decoy database molecule properties.',
                        action='store_true',
                        default=False)
 
    return parser.parse_args()

def findPropMatchingDecoyCandidates(input_mols: list, excl_mol_hashes: set, args: argparse.Namespace) -> list:
    print('Searching for property matching molecules in decoy database \'{args.decoy_db_props}\'...', file=sys.stderr)

    matched_decoys = []
    decoy_props = MolPropData()
    i = 1
    
    for input_mol in input_mols:
        num_sel = 0
        skip_header = True

        print(f' Trying to find matching decoys for input molecule {i}...', file=sys.stderr, end='\r')
        
        with open(args.decoy_db_props, 'r') as f:
            while decoy_props.read(f, skip_header):
                skip_header = False

                if decoy_props.hashCode in excl_mol_hashes:
                    continue
            
                if input_mol.matches(decoy_props, True, args):
                    matched_decoys.append(decoy_props)
                    excl_mol_hashes.add(decoy_props.hashCode)

                    decoy_props = MolPropData()
                    num_sel += 1

                    if num_sel >= MAX_POOL_SIZE:
                        break

        i += 1
        
    print(f' -> Found {len(matched_decoys)} decoy candidates                             ', file=sys.stderr)

    return matched_decoys

def filterDecoysByInputMolSim(decoy_mol_reader: Chem.MoleculeReader, input_mols: list, decoy_mols: list) -> list:
    print(f'Removing decoy candidates with input molecule Tanimoto similarity > {MAX_DECOY_INPUT_MOL_SIM}...', file=sys.stderr)

    ecfp_gen = Descr.CircularFingerprintGenerator()
        
    for input_mol in input_mols:
        input_mol.ecfp = Util.BitSet(ECFP4_LENGTH)

        ecfp_gen.generate(input_mol.molecule) 
        ecfp_gen.setFeatureBits(input_mol.ecfp)

        del input_mol.molecule
   
    res_decoys = []
    mol = Chem.BasicMolecule()
    ecfp = Util.BitSet(ECFP4_LENGTH)
    num_proc = 0
    
    for decoy_mol in decoy_mols:
        if not decoy_mol_reader.read(decoy_mol.index, mol):
            print(f' Error: could not read decoy database molecule at index {decoy_mol.index}', file=sys.stderr)
            continue

        if (num_proc % 1000) == 0:
            print(f' -> Processed {num_proc} candidates', file=sys.stderr, end='\r')
            
        Chem.calcBasicProperties(mol, False)

        ecfp_gen.generate(mol) 
        ecfp_gen.setFeatureBits(ecfp)

        max_sim = -1.0

        for input_mol in input_mols:
            max_sim = max(max_sim, Descr.calcTanimotoSimilarity(ecfp, input_mol.ecfp))
            
            if max_sim > MAX_DECOY_INPUT_MOL_SIM:
                break

        num_proc += 1
            
        if max_sim > MAX_DECOY_INPUT_MOL_SIM:
            continue

        decoy_mol.maxInputSim = max_sim
        decoy_mol.ecfp = ecfp

        ecfp = Util.BitSet(ECFP4_LENGTH)

        res_decoys.append(decoy_mol)

    print(f' -> Discarded {len(decoy_mols) - len(res_decoys)} candidates, {len(res_decoys)} remaining', file=sys.stderr)
    
    return res_decoys

def filterDecoysByDiversity(decoy_mols: list) -> list:
    print('Performing diversity picking on decoy candidates...', file=sys.stderr)

    decoy_mols = sorted(decoy_mols, key=lambda entry: entry.weight)
    res_decoys = []
    num_proc = 0
    
    for decoy_mol in decoy_mols:
        if (num_proc % 1000) == 0:
            print(f' -> Clustered {num_proc} candidates', file=sys.stderr, end='\r')
            
        sim_entry = -1
        
        for i in range(len(res_decoys) - 1, -1, -1):
            sim = Descr.calcTanimotoSimilarity(res_decoys[i].ecfp, decoy_mol.ecfp)

            if sim > MAX_DECOY_DECOY_SIM:
                sim_entry = i
                break
            
        num_proc += 1
        
        if sim_entry < 0:
            res_decoys.append(decoy_mol)
            continue
        
        if res_decoys[sim_entry].maxInputSim > decoy_mol.maxInputSim:
            res_decoys[sim_entry] = decoy_mol
            
    print(f' -> Discared {len(decoy_mols) - len(res_decoys)} candidates, {len(res_decoys)} final candidates remaining', file=sys.stderr)

    return res_decoys

def assignDecoys(input_mols: list, decoy_mols: list, first_pass: bool, params: argparse.Namespace) -> tuple:
    unass_decoy_mols = []
    
    for decoy_mol in decoy_mols:
        sel_ipt_mol = None
        
        for input_mol in input_mols:
            if first_pass:
                if len(input_mol.decoys) >= params.max_num_decoys:
                    continue
            else:
                if len(input_mol.decoys) >= params.min_num_decoys:
                    continue
                
            if not input_mol.matches(decoy_mol, first_pass, params):
                continue
            
            if not sel_ipt_mol or len(input_mol.decoys) < len(sel_ipt_mol.decoys):
                sel_ipt_mol = input_mol

        if not sel_ipt_mol:
            unass_decoy_mols.append(decoy_mol)
            continue

        sel_ipt_mol.decoys.append(decoy_mol)
       
    for input_mol in input_mols:
        if len(input_mol.decoys) < params.min_num_decoys:
            return False, unass_decoy_mols
            
    return True, unass_decoy_mols
        
def selectOutputDecoys(input_mols: list, decoy_mols: list, params: argparse.Namespace) -> None:
    print('Selecting output decoys...', file=sys.stderr)

    decoy_mols = sorted(decoy_mols, key=lambda entry: entry.maxInputSim)
 
    for input_mol in input_mols:
        input_mol.decoys = []

    complete, decoy_mols = assignDecoys(input_mols, decoy_mols, True, params)
    
    if not complete:
        assignDecoys(input_mols, decoy_mols, False, params)

    num_selected = 0
            
    for input_mol in input_mols:
        num_selected += len(input_mol.decoys)
        
        if len(input_mol.decoys) < params.min_num_decoys:
            print(f' Warning: insufficient number of decoys for input molecule \'{input_mol.name}\' ({len(input_mol.decoys)})', file=sys.stderr)
        
    print(f' -> Selected {num_selected} final decoys', file=sys.stderr)

def outputDecoys(input_mols: list, decoy_mol_reader: Chem.MoleculeReader, args: argparse.Namespace) -> None:
    print(f'Writing decoy molecules to file \'{args.output}\'...', file=sys.stderr)

    mol_writer = Chem.MolecularGraphWriter(args.output)
    mol = Chem.BasicMolecule()

    Chem.setMultiConfExportParameter(mol_writer, False)

    for input_mol in input_mols:
        for decoy_mol in input_mol.decoys:
            if not decoy_mol_reader.read(decoy_mol.index, mol):
                print(f' Error: could not read decoy database molecule at index {decoy_mol.index}', file=sys.stderr)
                continue

            Chem.calcBasicProperties(mol, False)
        
            if not mol_writer.write(mol):
                print(f' Error: could not output decoy database molecule at index {decoy_mol.index}', file=sys.stderr)
        
    mol_writer.close()
    decoy_mol_reader.close()
    
    print('Done!')   

def calcDecoyDBMolProperties(args: argparse.Namespace) -> None:
    print(f'Calculating properties for molecules in \'{args.decoy_db}\'...', file=sys.stderr)
    print(f' Writing properties to file \'{args.decoy_db_props}\'', file=sys.stderr)
    
    mol_reader = Chem.MoleculeReader(args.decoy_db)
    mol = Chem.BasicMolecule()
    decoy_props = MolPropData()
    write_header = True

    Chem.setMultiConfImportParameter(mol_reader, False)
    
    with open(args.decoy_db_props, 'w') as f:
        while True:
            try:
                while mol_reader.read(mol):
                    if (mol_reader.getRecordIndex() % 1000) == 0:
                        print(f' -> Processed {mol_reader.getRecordIndex()} molecules', file=sys.stderr, end='\r')

                    decoy_props.calculate(mol, mol_reader.getRecordIndex() - 1)
                    decoy_props.write(f, write_header)
                    
                    write_header = False

                break
                    
            except Base.IOError as e:
                print(f' Error: reading molecule at index {mol_reader.getRecordIndex()} failed', file=sys.stderr)

                reader.setRecordIndex(mol_reader.getRecordIndex() + 1)
             
    mol_reader.close()
    
    print(f' -> Calculated properties for {mol_reader.getRecordIndex()} molecules', file=sys.stderr)
    print('Done!')   

def loadInputMolecules(args: argparse.Namespace) -> tuple:
    print(f'Loading input molecules from file \'{args.input}\'...', file=sys.stderr)

    mol_reader = Chem.MoleculeReader(args.input)
    input_mols = []
    mol = Chem.BasicMolecule()
    
    Chem.setMultiConfImportParameter(mol_reader, False)

    while True:
        try:
            while mol_reader.read(mol):
                if (mol_reader.getRecordIndex() % 1000) == 0:
                    print(f' -> Processed {mol_reader.getRecordIndex()} molecules', file=sys.stderr, end='\r')

                props = MolPropData()

                props.calculate(mol, mol_reader.getRecordIndex() - 1)

                input_mols.append(props)
                
                mol = Chem.BasicMolecule()

            break
        
        except Base.IOError as e:
             print(f' Error: reading molecule at index {mol_reader.getRecordIndex()} failed', file=sys.stderr)

             mol_reader.setRecordIndex(mol_reader.getRecordIndex() + 1)

    print(f' -> Read {len(input_mols)} molecules', file=sys.stderr)

    excl_mol_hashes = set()
    
    for entry in input_mols:
        excl_mol_hashes.add(entry.hashCode)
    
    return (input_mols, excl_mol_hashes)

def loadExcludeMolecules(args: argparse.Namespace, excl_mol_hashes: set) -> None:
    print(f'Loading exclude molecule list from file \'{args.excl_mols}\'...', file=sys.stderr)

    mol_reader = Chem.MoleculeReader(args.excl_mols)
    mol = Chem.BasicMolecule()
    
    Chem.setMultiConfImportParameter(mol_reader, False)

    while True:
        try:
            while mol_reader.read(mol):
                Chem.calcBasicProperties(mol, False)
        
                excl_mol_hashes.add(Chem.calcHashCode(mol, Chem.AtomPropertyFlag.TYPE, Chem.BondPropertyFlag.ORDER | Chem.BondPropertyFlag.AROMATICITY))

            break
        
        except Base.IOError as e:
             print(f' Error: reading molecule at index {mol_reader.getRecordIndex()} failed', file=sys.stderr)

             reader.setRecordIndex(mol_reader.getRecordIndex() + 1)
        
    print(f' -> Read {mol_reader.getRecordIndex()} molecules', file=sys.stderr)
    
def main(args):
    if args.calc_props:
        calcDecoyDBMolProperties(args)
        return

    if not args.input:
        sys.exit('Error: missing input file argument')
   
    if not args.output:
        sys.exit('Error: missing output file argument')

    decoy_mol_reader = Chem.MoleculeReader(args.decoy_db)
    input_mols, excl_mol_hashes = loadInputMolecules(args)
    
    Chem.setMultiConfImportParameter(decoy_mol_reader, False)

    if args.excl_mols:
        loadExcludeMolecules(args, excl_mol_hashes)
    
    decoy_mols = findPropMatchingDecoyCandidates(input_mols, excl_mol_hashes, args)
    decoy_mols = filterDecoysByInputMolSim(decoy_mol_reader, input_mols, decoy_mols)
    decoy_mols = filterDecoysByDiversity(decoy_mols)

    selectOutputDecoys(input_mols, decoy_mols, args)
    outputDecoys(input_mols, decoy_mol_reader, args)

if __name__ == '__main__':
    main(parseArguments())
