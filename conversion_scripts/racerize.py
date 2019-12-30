#!/usr/bin/env python3
'''
Convert atomic PDB file to coarse-grained PDB or TXYZ file
'''
import os
import sys
from collections import defaultdict
from typing import List
import numpy as np
import argparse

BOND_CUTOFF = 5.0
RADIUS = 2.5
RESNAME_RNA = ("A", "C", "G", "U")
RESNAME_OTHER = ('MG','CA','NA','CL','K','TL','BA','NH4','HOH')
# [resname,] name : List[CG atom types], List[CG atom names]
# new atom types, e.g. heteroatoms, should be added here
CG_ATOMS = {
    ('A', 'N6'): ((4,), ('N6',)), 
    ('A', 'C8'): ((3,), ('CG',)), 
    ('A', 'C2'): ((9,), ('CA',)), 
    ('G', 'O6'): ((6,), ('O6',)), 
    ('G', 'N2'): ((5,), ('N2',)), 
    ('G', 'C8'): ((3,), ('CG',)), 
    ('C', 'O2'): ((7,), ('O2',)), 
    ('C', 'C6'): ((8,), ('CU',)), 
    ('C', 'N4'): ((4,), ('N6',)), 
    ('U', 'O2'): ((7,), ('O2',)), 
    ('U', 'C6'): ((8,), ('CU',)), 
    ('U', 'O4'): ((6,), ('O6',)), 
    "C4'": ((2, 11, 12, 13), ('S1', 'S2', 'G1', 'G2')),
    "P": ((1, 10, ), ('P1', 'P2')),
    }
_type_to_name = {}
for _, atypes in CG_ATOMS.items():
    for itype, name in zip(atypes[0], atypes[1]):
        _type_to_name[itype] = name
# used to determine bond cutoff
# this is copied from tinker prm file
_bond_param_str = '''\
bond         1    2         11.12440191     3.8900
bond         2    3         9.784688995     3.7400  
bond         3    4         71.77033493     4.2864
bond         4    9         71.98564593     3.5261
bond         3    9         71.79425837     4.3282
bond         2    8         10.88516746     3.6146
bond         4    8         71.81818182     3.5859
bond         4    7         71.96172249     4.5499
bond         7    8         71.81818182     3.5189
bond         3    6         71.84210526     4.2829
bond         5    6         71.81818182     4.5712
bond         3    5         71.96172249     5.6590 
bond         6    8         71.84210526     3.5544
bond         6    7         71.88995215     4.5255
bond        11    1         11.12440191     3.8900
bond        10   11         11.12440191     3.8900
bond        10    2         11.12440191     3.8900
bond        11    3         9.784688995     3.7400 
bond        11    8         10.88516746     3.6146 

bond         1   12         11.12440191     3.8900
bond        12    3         9.784688995     3.7400 
bond        12    8         10.88516746     3.6146
bond        10   12         11.12440191     3.8900
bond        13    1         11.12440191     3.8900
bond        10   13         11.12440191     3.8900
bond        13    3         9.784688995     3.7400 
bond        13    8         10.88516746     3.6146
'''
def get_bond_param(lines):
    res = []
    for line in lines:
        w = line.split('#')[0].split()
        if len(w) == 5 and w[0] == 'bond':
            res.append((int(w[1]), int(w[2]), float(w[3]), float(w[4])))
    return res
_bond_param = get_bond_param(_bond_param_str.split('\n'))
# connections within residue
CG_CONN = {}
for _bond in _bond_param:
    i1, i2, kb, b0 = _bond
    n1, n2 = (_type_to_name[i1], _type_to_name[i2])
    CG_CONN[(n1, n2)] = b0
    CG_CONN[(n2, n1)] = b0
# connections in backbone
CG_CONN_B = {}
for n1, n2 in CG_CONN:
    if n1 in ('P1', 'P2') and n2 in ('S1', 'S2', 'G1', 'G2'):
        CG_CONN_B[(n1, n2)] = CG_CONN[(n1, n2)]
        CG_CONN_B[(n2, n1)] = CG_CONN[(n1, n2)]

# CG_CONN format
CG_CONN.update({
('P1', 'S1'): 3.9,
})

def get_atypes(line):
    ''' Look up atom type in the definition
    '''
    name = line[12:16].strip()
    resname = line[17:20].strip()
    atype = []
    for idx in name, (resname, name):
        if idx in CG_ATOMS:
            atype = CG_ATOMS[idx]
    return atype

def pair_to_dict(conns):
    conn_dict = defaultdict(list)
    for pair in conns:
        conn_dict[pair[0]].append(pair[1])
        conn_dict[pair[1]].append(pair[0])
    for iatom1 in conn_dict:
        conn_dict[iatom1] = sorted(conn_dict[iatom1])
    return conn_dict

def convert_line(lines, chain: List[List[int]], istart=0, ftype='pdb'):
    ''' Given input lines and list of line numbers, return formatted lines and CONECT lines

    chain: a list of sub-lists of line numbers. Each sub-list defines a residue
    
    connections are searched within residues and between consecutive residues
    based on distance cutoff
    '''
    lines_pdb = []
    lines_pdbconn = []
    iatom = istart
    nres = len(chain)
    residue_last = []
    coord_last = []
    conns = []
    coords_tinker = []
    for ires, residue in enumerate(chain):
        coord = []
        if len(residue) == 0:
            continue
        is_term = ires in (0, nres-1)
        for iline in residue:
            line = lines[iline]
            atype = get_atypes(line)
            d = ires % 2
            if is_term:
                d += 2
            itype = d%len(atype[0])
            iatom += 1
            line1 = "HETATM" + "%5d"%iatom + line[11] + '%4s'%atype[1][itype] + '%4d'%atype[0][itype] + line[20:]
            lines_pdb.append(line1)
            R = np.array(list(map(float, [line[30:38], line[38:46], line[46:54]])))
            coord.append((atype[1][itype], iatom, R))
            coords_tinker.append((iatom, atype[1][itype], atype[0][itype], R))
        for i1 in range(len(coord)):
            atom1 = coord[i1]
            for conn_table, coord2 in ((CG_CONN, coord[i1+1:]), (CG_CONN_B, coord_last)):
                for atom2 in coord2:
                    if not (atom1[0], atom2[0]) in conn_table:
                        continue
                    r = np.linalg.norm(atom1[2] - atom2[2])
                    if r > 1.2 * conn_table[(atom1[0], atom2[0])]:
                        continue
                    pair = atom1[1], atom2[1]
                    if pair[0] > pair[1]:
                        pair = pair[1], pair[0]
                    conns.append(pair)
        residue_last = residue
        coord_last = coord
    conn_dict = pair_to_dict(conns)
    lines_tinker = []
    for (iatom1, name, atype, R) in coords_tinker:
        str_conn = ''
        if iatom1 in conn_dict:
            str_conn = (''.join('%5d'%_ for _ in conn_dict[iatom1]))
        lines_tinker.append('%5d %5s %11.6f %11.6f %11.6f %5d %s\n'%(iatom1, name, R[0], R[1], R[2], atype, str_conn))

    for iatom1 in sorted(conn_dict):
        iatoms = list(_ for _ in conn_dict[iatom1] if _ > iatom1)
        if len(iatoms) > 0:
            lines_pdbconn.append('CONECT%s\n'%(''.join('%5d'%_ for _ in [iatom1]+iatoms)))
    if ftype == 'tinker':
        return lines_tinker, []
    else:
        return lines_pdb, lines_pdbconn

def read_cg_atoms(finp, ftype='pdb'):
    '''READ CG atoms from PDB file
    '''
    lines = []
    if ftype == 'tinker':
        lines = ['0\n']
    else:
        lines = ['HEADER\n', 'COMPND\n', 'SOURCE\n']
    lines_end = []
    iatom = 0
    chain_res_last = ("",-1)
    with open(finp, 'r') as fh:
        lines0 = fh.readlines()
        curr_chain = [[]]
        for iline, line in enumerate(lines0):
            if not (line[:6] in ("ATOM  ", "HETATM") and line[16] in (" ", "A")):
                continue
            if line.startswith('ENDMDL'):
                break
            #name = line[12:16].strip()
            #resname = line[17:20].strip()
            resid = int(line[22:26])
            chainid = line[21]
            chain_res = (chainid, resid)
            if chain_res[0] != chain_res_last[0]:
                newlines, conn_lines = convert_line(lines0, curr_chain, istart=iatom, ftype=ftype)
                iatom += len(newlines)
                lines += newlines
                lines_end += conn_lines
                curr_chain = [[]]
            elif chain_res != chain_res_last:
                curr_chain.append([])
            chain_res_last = chain_res
            atype = get_atypes(line)
            if len(atype) > 0:
                curr_chain[-1].append(iline)
        newlines, conn_lines = convert_line(lines0, curr_chain, istart=iatom, ftype=ftype)
        iatom += len(newlines)
        lines += newlines
        lines_end += conn_lines
    if ftype == 'tinker':
        lines[0] = '%d\n'%(len(lines)-1)
    return lines + lines_end

def main():
    parser = argparse.ArgumentParser(prog='racerize.py', usage='%(prog)s -i in.pdb [-opdb out.pdb] [-oxyz out.xyz]', description=__doc__)
    parser.add_argument('-i', nargs='?', required=True, metavar='INP', help='input pdb file')
    parser.add_argument('-o','-opdb', nargs='?', required=False, metavar='PDB', default='', help='output pdb file')
    parser.add_argument('-oxyz', nargs='?', required=False, metavar='XYZ', default='', help='output txyz file')

    #v = (parser.parse_args(sys.argv[1:]))
    v = parser.parse_args()
    if v.o != '':
        ftype = 'pdb'
        fout = v.o
    elif v.oxyz != '':
        ftype = 'tinker'
        fout = v.oxyz
    else:
        ftype = 'pdb'
        fout = v.i+'.pdb'
        if os.path.isfile(fout):
            print(fout, 'exists. Please specify output file name.')
            return

    with open(fout, 'w') as fh_out:
        lines = read_cg_atoms(v.i, ftype=ftype)
        fh_out.write(''.join(lines))

if __name__ == "__main__":
    main()
