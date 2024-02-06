# -*- mode: python; tab-width: 4 -*- 

## 
# gen_decoys.py 
#
# Copyright (C) 2024 Thomas A. Seidel <thomas.seidel@univie.ac.at>
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 3 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program; see the file LICENSE. If not, write to
# the Free Software Foundation, Inc., 59 Temple Place - Suite 330,
# Boston, MA 02111-1307, USA.
##

import sys
import io
import pathlib
import argparse

import CDPL.Chem as Chem
import CDPL.MolProp as MolProp
import CDPL.Descr as Descr
import CDPL.Util as Util


MAX_POOL_SIZE            = 10000
ECFP4_LENGTH             = 6007
MAX_INPUT_MOL_SIMILARITY = 0.35
MAX_DECOY_MOL_SIMILARITY = 0.80


class MolPropData():

    def __init__(self):
        pass

    def matches(self, props: object, params: argparse.Namespace) -> bool:
        if abs(self.weight - props.weight) > params.mw_tol:
            return False
 
        if abs(self.logP - props.logP) > params.logp_tol:
            return False

        if abs(self.netCharge - props.netCharge) > params.charge_tol:
            return False

        if abs(self.numRotBonds - props.numRotBonds) > params.rbc_tol:
            return False

        if abs(self.numHBA - props.numHBA) > params.hba_tol:
            return False

        return (abs(self.numHBD - props.numHBD) <= params.hbd_tol)

    @staticmethod
    def calculate(mol: Chem.Molecule, mol_idx: int) -> object:
        Chem.initSubstructureSearchTarget(mol, False)

        props = MolPropData()

        props.index = mol_idx
        props.name = Chem.getName(mol)
        props.weight = MolProp.calcMass(mol)
        props.logP = MolProp.calcXLogP(mol)
        props.numRotBonds = MolProp.getRotatableBondCount(mol)
        props.numHBA = MolProp.getHBondAcceptorAtomCount(mol)
        props.numHBD = MolProp.getHBondDonorAtomCount(mol)
        props.hashCode = Chem.calcHashCode(mol, Chem.AtomPropertyFlag.TYPE, Chem.BondPropertyFlag.ORDER | Chem.BondPropertyFlag.AROMATICITY)
        props.netCharge = 0

        for a in mol.atoms:
            props.netCharge += Chem.getFormalCharge(a)

        return props

    @staticmethod
    def read(f: io.TextIOBase) -> list:
        prop_table = []

        f.readline()

        for line in f:
            fields = line.strip().split(',')

            props = MolPropData()
            props.index = int(fields[0])
            props.weight = float(fields[1])
            props.logP = float(fields[2])
            props.netCharge = int(fields[3])
            props.numRotBonds = int(fields[4])
            props.numHBA = int(fields[5])
            props.numHBD = int(fields[6])
            props.hashCode = int(fields[7])

            prop_table.append(props)

        return prop_table

    @staticmethod
    def write(table: list, f: io.TextIOBase) -> None:
        f.write('Mol. Index,MW,XLogP,Net Charge,#Rot. Bonds,#HBA,#HBD,Hash Code\n')

        for entry in table:
            f.write(f'{entry.index},{entry.weight},{entry.logP},{entry.netCharge},{entry.numRotBonds},{entry.numHBA},{entry.numHBD},{entry.hashCode}\n')

    
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
                        help='[Optional] LogP matching tolerance (default: 1.5).',
                        metavar='<float>',
                        type=float,
                        default=1.5)
    parser.add_argument('--rbc-tol',
                        dest='rbc_tol',
                        required=False,
                        help='[Optional] Rotatable bond count matching tolerance (default: 1).',
                        metavar='<int>',
                        type=int,
                        default=1)
    parser.add_argument('--hba-tol',
                        dest='hba_tol',
                        required=False,
                        help='[Optional] H-bond acceptor atom count matching tolerance (default: 1).',
                        metavar='<int>',
                        type=int,
                        default=1)
    parser.add_argument('--hbd-tol',
                        dest='hbd_tol',
                        required=False,
                        help='[Optional] H-bond donor atom count matching tolerance (default: 1).',
                        metavar='<int>',
                        type=int,
                        default=1)
    parser.add_argument('--calc-props',
                        dest='calc_props',
                        required=False,
                        help='[Optional] Calculate and store decoy database molecule properties.',
                        action='store_true',
                        default=False)
 
    return parser.parse_args()

def findPropMatchingDecoyCandidates(input_mol_props: list, decoy_mol_props: list, excl_mol_set: set, params: argparse.Namespace) -> list:
    print('Searching for property matching molecules in decoy database...', file=sys.stderr)

    sel_mask = Util.BitSet(len(decoy_mol_props))
    matching = []

    for entry in input_mol_props:
        num_sel = 0

        for i in range(len(decoy_mol_props)):
            if sel_mask.test(i):
                continue

            if decoy_mol_props[i].hashCode in excl_mol_set:
                continue
            
            if entry.matches(decoy_mol_props[i], params):
                matching.append(decoy_mol_props[i])
                sel_mask.set(i)

                num_sel += 1

                if num_sel >= MAX_POOL_SIZE:
                    break

    print(f' -> Found {len(matching)} decoy candidates', file=sys.stderr)

    return matching

def filterDecoysByInputMolSim(decoy_mol_reader: Chem.MoleculeReader, input_mol_props: list, decoy_mol_props: list) -> list:
    ecfp_gen = Descr.CircularFingerprintGenerator()
    
    print(f'Removing decoy candidates with input molecule Tanimoto similarity > {MAX_INPUT_MOL_SIMILARITY}...', file=sys.stderr)

    for entry in input_mol_props:
        entry.ecfp = Util.BitSet(ECFP4_LENGTH)

        ecfp_gen.generate(entry.molecule) 
        ecfp_gen.setFeatureBits(entry.ecfp)

        del entry.molecule
   
    remaining = []
    mol = Chem.BasicMolecule()
    ecfp = Util.BitSet(ECFP4_LENGTH)
    
    for decoy_entry in decoy_mol_props:
        if not decoy_mol_reader.read(decoy_entry.index, mol):
            print(f' -> Error: could not read decoy database molecule at index {decoy_entry.index}', file=sys.stderr)
            continue

        Chem.calcBasicProperties(mol, False)

        ecfp_gen.generate(mol) 
        ecfp_gen.setFeatureBits(ecfp)

        max_sim = -1.0

        for ipt_mol_entry in input_mol_props:
            max_sim = max(max_sim, Descr.calcTanimotoSimilarity(ecfp, ipt_mol_entry.ecfp))
            
            if max_sim > MAX_INPUT_MOL_SIMILARITY:
                break

        if max_sim > MAX_INPUT_MOL_SIMILARITY:
            continue

        decoy_entry.maxInputSim = max_sim
        decoy_entry.ecfp = ecfp

        ecfp = Util.BitSet(ECFP4_LENGTH)

        remaining.append(decoy_entry)

    print(f' -> Discarded {len(decoy_mol_props) - len(remaining)} decoy candidates, {len(remaining)} remaining', file=sys.stderr)
    
    return remaining

def filterDecoysByDiversity(decoy_mol_props: list) -> list:
    print('Performing diversity picking on decoy candidates...', file=sys.stderr)

    decoy_mol_props = sorted(decoy_mol_props, key=lambda entry: entry.weight)
    selected = []

    for entry in decoy_mol_props:
        sim_entry_idx = -1
        discard = False
        
        for i in range(len(selected)):
            sim = Descr.calcTanimotoSimilarity(selected[i].ecfp, entry.ecfp)

            if sim > MAX_DECOY_MOL_SIMILARITY:
                if sim_entry_idx >= 0:
                    discard = True
                    break
                
                sim_entry_idx = i

        if discard:
            continue
        
        if sim_entry_idx < 0:
            selected.append(entry)
            continue
        
        if selected[sim_entry_idx].maxInputSim > entry.maxInputSim:
            selected[sim_entry_idx] = entry
            
    print(f' -> Discared {len(decoy_mol_props) - len(selected)} decoy candidates, {len(selected)} final candidates remaining', file=sys.stderr)

    return selected

def selectOutputDecoys(input_mol_props: list, decoy_mol_props: list, params: argparse.Namespace) -> list:
    print('Selecting output decoys...', file=sys.stderr)

    selected = []

    for entry in input_mol_props:
        entry.numDecoys = 0
    
    for entry in decoy_mol_props:
        sel_ipt_mol_idx = -1
        exit_loop = True
        
        for i in range(len(input_mol_props)):
            if input_mol_props[i].numDecoys >= params.max_num_decoys:
                continue

            exit_loop = False

            if not input_mol_props[i].matches(entry, params):
                continue
            
            if sel_ipt_mol_idx < 0 or input_mol_props[i].numDecoys < input_mol_props[sel_ipt_mol_idx].numDecoys:
                sel_ipt_mol_idx = i

        if sel_ipt_mol_idx < 0:
            if exit_loop:
                break

            continue

        selected.append(entry)
        
        input_mol_props[sel_ipt_mol_idx].numDecoys += 1

    for entry in input_mol_props:
        if entry.numDecoys < params.min_num_decoys:
            print(f' Warning: insufficient number of decoys for input molecule \'{entry.name}\' ({entry.numDecoys})', file=sys.stderr)
        
    print(f' -> Selected {len(selected)} final decoys', file=sys.stderr)
    
    return selected

def outputDecoys(decoy_mol_props: list, decoy_mol_reader: Chem.MoleculeReader, args: argparse.Namespace) -> None:
    print(f'Writing decoy molecules to file \'{args.output}\'...', file=sys.stderr)

    mol_writer = Chem.MolecularGraphWriter(args.output)
    mol = Chem.BasicMolecule()
    
    for entry in decoy_mol_props:
        if not decoy_mol_reader.read(entry.index, mol):
            print(f' -> Error: could not read decoy database molecule at index {entry.index}', file=sys.stderr)
            continue

        Chem.calcBasicProperties(mol, False)
        
        if not mol_writer.write(mol):
            print(f' -> Error: could not write decoy database molecule at index {entry.index}', file=sys.stderr)
        
    mol_writer.close()
    decoy_mol_reader.close()
    
    print('Done!')   

def readMolsAndCalcProperties(reader: Chem.MoleculeReader, store_mol: bool) -> list:
    prop_table = []
    mol = Chem.BasicMolecule()

    while reader.read(mol):
        if (reader.getRecordIndex() % 1000) == 0:
            print(f' -> Passed {reader.getRecordIndex()}th molecule', file=sys.stderr, end='\r')

        props = MolPropData.calculate(mol, reader.getRecordIndex() - 1)

        if (store_mol):
            props.molecule = mol
            mol = Chem.BasicMolecule()

        prop_table.append(props)

    return prop_table

def calcDecoyDBMolProperties(args: argparse.Namespace) -> None:
    print(f'Calculating properties for molecules in \'{args.decoy_db}\'...', file=sys.stderr)

    mol_reader = Chem.MoleculeReader(args.decoy_db)
    prop_table = readMolsAndCalcProperties(mol_reader, False)

    print(f' -> Calculated properties for {len(prop_table)} molecules', file=sys.stderr)
    print(f'Writing properties to file \'{args.decoy_db_props}\'...', file=sys.stderr)
        
    with open(args.decoy_db_props, 'w') as f:
        MolPropData.write(prop_table, f)

    mol_reader.close()
    
    print('Done!')   

def loadDecoyDBMolProperties(args: argparse.Namespace) -> list:
    print(f'Loading precalculated decoy database molecule properties from file \'{args.decoy_db_props}\'...', file=sys.stderr)

    with open(args.decoy_db_props, 'r') as f:
        decoy_db_props = MolPropData.read(f)

        print(f' -> Read {len(decoy_db_props)} records', file=sys.stderr)

        return decoy_db_props

def loadInputMolecules(args: argparse.Namespace) -> tuple:
    print(f'Loading input molecules from file \'{args.input}\'...', file=sys.stderr)

    mol_reader = Chem.MoleculeReader(args.input)
    input_mol_props = readMolsAndCalcProperties(mol_reader, True)

    print(f' -> Read {len(input_mol_props)} molecules', file=sys.stderr)

    excl_mol_set = set()

    for entry in input_mol_props:
        excl_mol_set.add(entry.hashCode)
    
    return (input_mol_props, excl_mol_set)

def loadExcludeMolecules(args: argparse.Namespace, excl_mol_set: set) -> None:
    print(f'Loading exclude molecule list from file \'{args.excl_mols}\'...', file=sys.stderr)

    mol_reader = Chem.MoleculeReader(args.excl_mols)
    mol = Chem.BasicMolecule()
    
    while mol_reader.read(mol):
        Chem.calcBasicProperties(mol, False)
        
        excl_mol_set.add(Chem.calcHashCode(mol, Chem.AtomPropertyFlag.TYPE, Chem.BondPropertyFlag.ORDER | Chem.BondPropertyFlag.AROMATICITY))
        
    print(f' -> Read {mol_reader.getNumRecords()} molecules', file=sys.stderr)
    
def process(args):
    if args.calc_props:
        calcDecoyDBMolProperties(args)
        return

    if not args.input:
        sys.exit('Error: missing input file argument')
   
    if not args.output:
        sys.exit('Error: missing output file argument')

    decoy_mol_reader = Chem.MoleculeReader(args.decoy_db)
    input_mol_props, excl_mol_set = loadInputMolecules(args)

    if args.excl_mols:
        loadExcludeMolecules(args, excl_mol_set)
    
    decoy_mol_props = loadDecoyDBMolProperties(args)
    decoy_mol_props = findPropMatchingDecoyCandidates(input_mol_props, decoy_mol_props, excl_mol_set, args)
    decoy_mol_props = filterDecoysByInputMolSim(decoy_mol_reader, input_mol_props, decoy_mol_props)
    decoy_mol_props = filterDecoysByDiversity(decoy_mol_props)
    decoy_mol_props = selectOutputDecoys(input_mol_props, decoy_mol_props, args)

    outputDecoys(decoy_mol_props, decoy_mol_reader, args)

if __name__ == '__main__':
    process(parseArguments())
