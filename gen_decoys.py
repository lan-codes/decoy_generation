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


MAX_POOL_SIZE = 10000
ECFP4_LENGTH  = 8179


class MolPropData():

    def __init__(self):
        pass

    def matches(self, props: object, params: argparse.Namespace) -> bool:
        if abs(self.molWeight - props.molWeight) > params.mw_tol:
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
        props.molIndex = mol_idx
        props.molWeight = MolProp.calcMass(mol)
        props.logP = MolProp.calcXLogP(mol)
        props.numRotBonds = MolProp.getRotatableBondCount(mol)
        props.numHBA = MolProp.getHBondAcceptorAtomCount(mol)
        props.numHBD = MolProp.getHBondDonorAtomCount(mol)
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
            props.molIndex = int(fields[0])
            props.molWeight = float(fields[1])
            props.logP = float(fields[2])
            props.netCharge = int(fields[3])
            props.numRotBonds = int(fields[4])
            props.numHBA = int(fields[5])
            props.numHBD = int(fields[6])

            prop_table.append(props)

        return prop_table

    @staticmethod
    def write(table: list, f: io.TextIOBase) -> None:
        f.write('Mol. Index,MW,XLogP,Net Charge,#Rot. Bonds,#HBA,#HBD\n')

        for entry in table:
            f.write(f'{entry.molIndex},{entry.molWeight},{entry.logP},{entry.netCharge},{entry.numRotBonds},{entry.numHBA},{entry.numHBD}\n')


def readMolsAndCalcProperties(reader: Chem.MoleculeReader, store_mol: bool) -> list:
    prop_table = []
    mol = Chem.BasicMolecule()
    
    while reader.read(mol):
        print(f' - At molecule {reader.getRecordIndex()}', file=sys.stderr, end='\r')

        props = MolPropData.calculate(mol, reader.getRecordIndex() - 1)

        if (store_mol):
            props.molecule = mol
            mol = Chem.BasicMolecule()

        prop_table.append(props)

    return prop_table

def calcInputMolECFPs(input_mol_props: list) -> None:
    print('Calculating input molecule fingerprints...', file=sys.stderr)

    ecfp_gen = Descr.CircularFingerprintGenerator()

    ecfp_gen.setNumIterations(2)

    for ipt_prop_set in input_mol_props:
        ipt_prop_set.ecfp = Util.BitSet(ECFP4_LENGTH)

        ecfp_gen.generate(ipt_prop_set.molecule) 
        ecfp_gen.setFeatureBits(ipt_prop_set.ecfp)

def selectPropMatchingDecoyCandidates(input_mol_props: list, decoy_db_props: list, params: argparse.Namespace) -> list:
    print('Selecting property matched decoy candidates...', file=sys.stderr)

    selected = Util.BitSet(len(decoy_db_props))
    candidates = []

    for ipt_prop_set in input_mol_props:
        num_sel = 0

        for i in range(len(decoy_db_props)):
            if num_sel >= MAX_POOL_SIZE:
                break

            if selected.test(i):
                continue

            if ipt_prop_set.matches(decoy_db_props[i], params):
                candidates.append(decoy_db_props[i])
                selected.set(i)

                num_sel += 1

    print(f' - Selected {len(candidates)} candidates', file=sys.stderr)

    return candidates

def parseArguments():
    parser = argparse.ArgumentParser(description='Retrieves decoys for a given set of input molecules.')
    
    parser.add_argument('-i',
                        dest='input',
                        required=False,
                        help='[Required/Optional] File containing the input molecules (e.g. actives).',
                        default=None,
                        metavar='<path>',
                        nargs=1)
    parser.add_argument('-o',
                        dest='output',
                        required=False,
                        help='[Required/Optional] Output file for the retrieved decoy molecules.',
                        default=None,
                        metavar='<path>',
                        nargs=1)
    parser.add_argument('-d',
                        dest='decoy_db',
                        required=True,
                        help='[Required] File providing the decoy source molecules.',
                        metavar='<path>',
                        nargs=1)
    parser.add_argument('-p',
                        dest='decoy_db_props',
                        required=True,
                        help='[Required] CSV-file providing the precalculated decoy database molecule properties.',
                        metavar='<path>',
                        nargs=1)
    parser.add_argument('--mw-tol',
                        dest='mw_tol',
                        required=False,
                        help='[Optional] Moleular weight matching tolerance (default: 20.0).',
                        metavar='<float>',
                        type=float,
                        default=20.0,
                        nargs=1)
    parser.add_argument('--charge-tol',
                        dest='charge_tol',
                        required=False,
                        help='[Optional] Net charge matching tolerance (default: 0).',
                        metavar='<int>',
                        type=int,
                        default=0,
                        nargs=1)
    parser.add_argument('--logp-tol',
                        dest='logp_tol',
                        required=False,
                        help='[Optional] LogP matching tolerance (default: 1.5).',
                        metavar='<float>',
                        type=float,
                        default=1.5,
                        nargs=1)
    parser.add_argument('--rbc-tol',
                        dest='rbc_tol',
                        required=False,
                        help='[Optional] Rotatable bond count matching tolerance (default: 1).',
                        metavar='<int>',
                        type=int,
                        default=1,
                        nargs=1)
    parser.add_argument('--hba-tol',
                        dest='hba_tol',
                        required=False,
                        help='[Optional] H-bond acceptor atom count matching tolerance (default: 1).',
                        metavar='<int>',
                        type=int,
                        default=1,
                        nargs=1)
    parser.add_argument('--hbd-tol',
                        dest='hbd_tol',
                        required=False,
                        help='[Optional] H-bond donor atom count matching tolerance (default: 1).',
                        metavar='<int>',
                        type=int,
                        default=1,
                        nargs=1)
    parser.add_argument('--calc-props',
                        dest='calc_props',
                        required=False,
                        help='[Optional] Calculate and store decoy database molecule properties.',
                        action='store_true',
                        default=False)
 
    return parser.parse_args()

def process(args):
    if args.calc_props:
        mol_reader = Chem.MoleculeReader(args.decoy_db[0])

        print(f'Calculating properties for molecules in \'{args.decoy_db[0]}\'...', file=sys.stderr)

        prop_table = readMolsAndCalcProperties(mol_reader, False)

        print(f' - Calculated properties for {len(prop_table)} molecules', file=sys.stderr)
        print(f'Writing properties to file \'{args.decoy_db_props[0]}\'...', file=sys.stderr)
        
        with open(args.decoy_db_props[0], 'w') as f:
            MolPropData.write(prop_table, f)

        mol_reader.close()


        return

    if not args.input:
        sys.exit('Error: missing input file argument!')
   
    if not args.output:
        sys.exit('Error: missing output file argument!')

    print(f'Loading precalculated decoy database molecule properties from \'{args.decoy_db_props[0]}\'...', file=sys.stderr)

    decoy_db_props = None

    with open(args.decoy_db_props[0], 'r') as f:
        decoy_db_props = MolPropData.read(f)

    print(f' - Loaded {len(decoy_db_props)} property records', file=sys.stderr)

    mol_reader = Chem.MoleculeReader(args.input[0])
    
    print(f'Loading/preprocessing input molecules from file \'{args.input[0]}\'...', file=sys.stderr)

    input_mol_props = readMolsAndCalcProperties(mol_reader, True)

    print(f' - Loaded {len(input_mol_props)} molecules', file=sys.stderr)
    
    selectPropMatchingDecoyCandidates(input_mol_props, decoy_db_props, args)
    calcInputMolECFPs(input_mol_props)

    print('Done!')

if __name__ == '__main__':
    process(parseArguments())
