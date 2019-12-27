#!/usr/bin/env python
import numpy
import sys

#Script to convert coarse-grained xyz to all atom pdb. After this script, convert to xyz and minimize before covnerting back to pdb.:

if len(sys.argv) != 1:
        interest = str(sys.argv[1])
else:
        interest=str(raw_input("Please enter the cg xyz filename: "))

#sequence=['G','G','G','C','G','C','A','A','G','C','C','U'];

f=open(interest,"r")
everything=f.readlines()
f.close()
p=0
g=open(interest+".aa.pdb","w+")
#g.write("HEADER"+"\n")
#g.write("COMPND"+"\n")
#g.write("SOURCE"+"\n")
i=2
j=0
k=1
phos=[1,10]
sugar=[2,11]
chain=["A","B","C","D","E","F"]
chain_num = -1
pdb_atom_number=1;
res=0;


## writeatoms subroutine ############
def writeatoms(res_atoms,res_coords,pdb_atom_number,res,chn):
	# write atoms
	for aa in range(len(res_atoms)):

		if len(res_atoms[aa])==2:
 		       	g.write('ATOM' + '%7d' % pdb_atom_number + '  '+ res_atoms[aa] +' ' + '%4s' % sequence[res] + ' '+chain[chn] + '%4s' % str(res+1) + '%12.3f' % res_coords[aa][0] + '%8.3f' % res_coords[aa][1] + '%8.3f' % res_coords[aa][2] + '\n') 
			pdb_atom_number=pdb_atom_number + 1;

		elif len(res_atoms[aa])==1:
			g.write('ATOM' + '%7d' % pdb_atom_number + '  '+ res_atoms[aa] +' ' + '%5s' % sequence[res] + ' '+chain[chn] + '%4s' % str(res+1) + '%12.3f' % res_coords[aa][0] + '%8.3f' % res_coords[aa][1] + '%8.3f' % res_coords[aa][2] + '\n') 
                        pdb_atom_number=pdb_atom_number + 1;

		elif len(res_atoms[aa])==3:
			g.write('ATOM' + '%7d' % pdb_atom_number + '  '+ res_atoms[aa] +' ' + '%3s' % sequence[res] + ' '+chain[chn] + '%4s' % str(res+1) + '%12.3f' % res_coords[aa][0] + '%8.3f' % res_coords[aa][1] + '%8.3f' % res_coords[aa][2] + '\n') 
                        pdb_atom_number=pdb_atom_number + 1;

		elif len(res_atoms[aa])==4:
			g.write('ATOM' + '%7d' % pdb_atom_number + ' '+ res_atoms[aa] +' ' + '%3s' % sequence[res] + ' '+chain[chn] + '%4s' % str(res+1) + '%12.3f' % res_coords[aa][0] + '%8.3f' % res_coords[aa][1] + '%8.3f' % res_coords[aa][2] + '\n') 
                        pdb_atom_number=pdb_atom_number + 1;

	return pdb_atom_number

############################################

## orthogonal unit vector subroutine ############
def orthogvector(atom1,atom2,atom3):
	v1 = [ (atom2[0]-atom1[0]), (atom2[1]-atom1[1]), (atom2[2]-atom1[2])]
	v2 = [ (atom3[0]-atom2[0]), (atom3[1]-atom2[1]), (atom3[2]-atom2[2])]
	crossp = numpy.cross(v1,v2)             ## cross product
	mag_crossp = numpy.sqrt(crossp[0]*crossp[0] + crossp[1]*crossp[1] + crossp[2]*crossp[2])
	unit_crossp = [ crossp[0]/mag_crossp, crossp[1]/mag_crossp, crossp[2]/mag_crossp ]      ## unit vector

	return unit_crossp

############################################



#first get sequence #########################
xyz_atoms=[]; sequence = []; chain_breaks = [];  # chain breaks is residue number before break
for entry in everything:
	if str(entry[8:11])=='':
		continue
	else:
		xyz_atoms.append(str(entry[8:11].strip()))

residues = 0
for i in range(len(xyz_atoms)):
	
	if xyz_atoms[i]=="CG" and xyz_atoms[i+1]=="O6" and xyz_atoms[i+2]=="N2":	#G
		sequence.append("G")
		residues = residues + 1
	elif xyz_atoms[i]=="CG" and xyz_atoms[i+1]=="N6" and xyz_atoms[i+2]=="CA":	#A
		sequence.append("A")
		residues = residues + 1
	elif xyz_atoms[i]=="CU" and xyz_atoms[i-2]=="O2" and xyz_atoms[i-1]=="N6": 	#C
		sequence.append("C")
		residues = residues + 1
	elif xyz_atoms[i]=="CU" and xyz_atoms[i-2]=="O2" and xyz_atoms[i-1]=="O6": 	#U
		sequence.append("U")
		residues = residues + 1
	elif xyz_atoms[i]=="G1" or xyz_atoms[i]=="G2":
		for j in range(len(xyz_atoms)):
			if j>i and j<=i+6 and (xyz_atoms[j]=="G1" or xyz_atoms[j]=="G2"):
				chain_breaks.append(residues+1)
				break
			elif j>i and j<=i+6 and j>len(xyz_atoms)-4:
				chain_breaks.append(residues+1)
				break

print "Sequence: " + "".join(sequence)
print chain_breaks
###############################################


## Now go through every entry and write base and backbone atoms #####

for line in everything:
	if str(line[8:11])=='':  ## atom name
                continue

	linesp = line.split()  ## linesp[0] = atom number, linesp[1] = atom_name, linesp[2] = x, linesp[3] = y, linesp[4] = z

	## First Sugar in Chain ##
	if str(linesp[1][0])=='G' and (res==0 or res in chain_breaks):  #res is zero initialized, so if res in chain_breaks
		chain_num = chain_num + 1				# actually residue after chain_break (new chain)
		sugar_atom_num = int(linesp[0])
		b_c4p = [float(linesp[2]),float(linesp[3]),float(linesp[4])]
                for entry in everything:
			entrysp = entry.split()
			if len(entrysp)>1:
				if sugar_atom_num < int(entrysp[0]) < sugar_atom_num+6 and (str(entrysp[1])=='CU' or str(entrysp[1])=='CG'):
					b_sc = [float(entrysp[2]),float(entrysp[3]),float(entrysp[4])]
				elif sugar_atom_num < int(entrysp[0]) < sugar_atom_num+6 and str(entrysp[1][0])=='P':
					b_p = [float(entrysp[2]),float(entrysp[3]),float(entrysp[4])]
                               		break


		b_o5p = [ b_c4p[0] - 2*(b_p[0]-b_c4p[0])/3 , b_c4p[1] - 2*(b_p[1]-b_c4p[1])/3 , b_c4p[2] - 2*(b_p[2]-b_c4p[2])/3]
		b_c5p = [ b_c4p[0] - (b_p[0]-b_c4p[0])/3 , b_c4p[1] - (b_p[1]-b_c4p[1])/3 , b_c4p[2] - (b_p[2]-b_c4p[2])/3]
		b_o4p = [ (b_sc[0]-b_c4p[0])/2 + b_c4p[0] , (b_sc[1]-b_c4p[1])/2 + b_c4p[1] , (b_sc[2]-b_c4p[2])/2 + b_c4p[2]]
		b_c3p = [ (b_p[0]-b_c4p[0])/3 + b_c4p[0] , (b_p[1]-b_c4p[1])/3 + b_c4p[1] , (b_p[2]-b_c4p[2])/3 + b_c4p[2]]
		b_o3p = [ 2*(b_p[0]-b_c4p[0])/3 + b_c4p[0] , 2*(b_p[1]-b_c4p[1])/3 + b_c4p[1] , 2*(b_p[2]-b_c4p[2])/3 + b_c4p[2]]
		b_c2p = [ (b_sc[0]-b_c3p[0])/3 + b_c3p[0] , (b_sc[1]-b_c3p[1])/3 + b_c3p[1] , (b_sc[2]-b_c3p[2])/3 + b_c3p[2]]
		b_o2p = [ (b_c2p[0]-b_c5p[0])/3 + b_c2p[0] , (b_c2p[1]-b_c5p[1])/2 + b_c2p[1] , (b_c2p[2]-b_c5p[2])/2 + b_c2p[2]]
		b_c1p = [ 2*(b_sc[0]-b_c3p[0])/3 + b_c3p[0] , 2*(b_sc[1]-b_c3p[1])/3 + b_c3p[1] , 2*(b_sc[2]-b_c3p[2])/3 + b_c3p[2]] 

                unit_crossp = orthogvector(b_c4p,b_c2p,b_c1p) # hydrogen on c2' atom

                b_h2p = [ 1.1*unit_crossp[0] + b_c2p[0], 1.1*unit_crossp[1] + b_c2p[1], 1.1*unit_crossp[2] + b_c2p[2]]

		unit_crossp = orthogvector(b_c5p,b_c4p,b_o4p) # hydrogen on c4' atom

		b_h4p = [ 1.1*unit_crossp[0] + b_c4p[0], 1.1*unit_crossp[1] + b_c4p[1], 1.1*unit_crossp[2] + b_c4p[2]]

		unit_crossp = orthogvector(b_o5p,b_c5p,b_o4p) # hydrogens on c5' atom

                b_h5p = [ 1.1*unit_crossp[0] + b_c5p[0], 1.1*unit_crossp[1] + b_c5p[1], 1.1*unit_crossp[2] + b_c5p[2]]
                b_h5pp = [ -1.1*unit_crossp[0] + b_c5p[0], -1.1* unit_crossp[1] + b_c5p[1], -1.1*unit_crossp[2] + b_c5p[2]]

		s_atoms =  ["O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'", "H2'", "H4'", "H5'", "H5''"]
		s_coords = [b_o5p, b_c5p, b_c4p, b_o4p, b_c3p, b_o3p, b_c2p, b_o2p, b_c1p, b_h2p, b_h4p, b_h5p, b_h5pp]

		# write all sugar atoms
                pdb_atom_number = writeatoms(s_atoms,s_coords,pdb_atom_number,res,chain_num)
                continue        # move on to next iteration



	## Non-terminal Sugars and Phosphates, Complete Backbone ####
        if str(linesp[1][0])=='S': 
                sugar_atom_num = int(linesp[0])
                b_c4p = [float(linesp[2]),float(linesp[3]),float(linesp[4])]

                for entry in everything:
			entrysp = entry.split()
			if len(entrysp)>1:
	                        if sugar_atom_num < int(entrysp[0]) < sugar_atom_num+6 and (str(entrysp[1])=='CU' or str(entrysp[1])=='CG'):
        	                        b_sc = [float(entrysp[2]),float(entrysp[3]),float(entrysp[4])]	## sugar connect
				elif sugar_atom_num-3 < int(entrysp[0]) < sugar_atom_num and str(entrysp[1][0])=='P':
					b_p_curr = [float(entrysp[2]),float(entrysp[3]),float(entrysp[4])] # current phosphate
				elif sugar_atom_num-6 < int(entrysp[0]) < sugar_atom_num and (str(entrysp[1][0])=='S' or str(entrysp[1][0])=='G'):
					b_c4p_prev = [float(entrysp[2]),float(entrysp[3]),float(entrysp[4])]  # previous sugar
				elif sugar_atom_num < int(entrysp[0]) < sugar_atom_num+6 and str(entrysp[1][0])=='P':
                        	        b_p_next = [float(entrysp[2]),float(entrysp[3]),float(entrysp[4])]        # next phosphate
					break
	

		unit_crossp = orthogvector(b_c4p_prev,b_p_curr,b_c4p) # oxygens of phosphate group

		b_op1 = [ 1.5*unit_crossp[0] + b_p_curr[0], 1.5*unit_crossp[1] + b_p_curr[1], 1.5*unit_crossp[2] + b_p_curr[2]]
		b_op2 = [ -1.5*unit_crossp[0] + b_p_curr[0], -1.5* unit_crossp[1] + b_p_curr[1], -1.5*unit_crossp[2] + b_p_curr[2]]

                b_o5p = [ (b_c4p[0]-b_p_curr[0])/3 + b_p_curr[0], (b_c4p[1]-b_p_curr[1])/3 + b_p_curr[1], (b_c4p[2]-b_p_curr[2])/3 + b_p_curr[2]]
                b_c5p = [ 2*(b_c4p[0]-b_p_curr[0])/3 + b_p_curr[0], 2*(b_c4p[1]-b_p_curr[1])/3 + b_p_curr[1], 2*(b_c4p[2]-b_p_curr[2])/3 + b_p_curr[2]]
                b_o4p = [ (b_sc[0]-b_c4p[0])/2 + b_c4p[0] , (b_sc[1]-b_c4p[1])/2 + b_c4p[1] , (b_sc[2]-b_c4p[2])/2 + b_c4p[2]]
                b_c3p = [ (b_p_next[0]-b_c4p[0])/3 + b_c4p[0] , (b_p_next[1]-b_c4p[1])/3 + b_c4p[1] , (b_p_next[2]-b_c4p[2])/3 + b_c4p[2]]
                b_o3p = [ 2*(b_p_next[0]-b_c4p[0])/3 + b_c4p[0] , 2*(b_p_next[1]-b_c4p[1])/3 + b_c4p[1] , 2*(b_p_next[2]-b_c4p[2])/3 + b_c4p[2]]
                b_c2p = [ (b_sc[0]-b_c3p[0])/3 + b_c3p[0] , (b_sc[1]-b_c3p[1])/3 + b_c3p[1] , (b_sc[2]-b_c3p[2])/3 + b_c3p[2]]
                b_o2p = [ (b_c2p[0]-b_c5p[0])/3 + b_c2p[0] , (b_c2p[1]-b_c5p[1])/2 + b_c2p[1] , (b_c2p[2]-b_c5p[2])/2 + b_c2p[2]]
                b_c1p = [ 2*(b_sc[0]-b_c3p[0])/3 + b_c3p[0] , 2*(b_sc[1]-b_c3p[1])/3 + b_c3p[1] , 2*(b_sc[2]-b_c3p[2])/3 + b_c3p[2]] 

                unit_crossp = orthogvector(b_c4p,b_c2p,b_c1p) # hydrogen on c2' atom

                b_h2p = [ 1.1*unit_crossp[0] + b_c2p[0], 1.1*unit_crossp[1] + b_c2p[1], 1.1*unit_crossp[2] + b_c2p[2]]

		unit_crossp = orthogvector(b_o5p,b_c5p,b_o4p)

                b_h5p = [ 1.1*unit_crossp[0] + b_c5p[0], 1.1*unit_crossp[1] + b_c5p[1], 1.1*unit_crossp[2] + b_c5p[2]]
                b_h5pp = [ -1.1*unit_crossp[0] + b_c5p[0], -1.1* unit_crossp[1] + b_c5p[1], -1.1*unit_crossp[2] + b_c5p[2]]

                b_atoms =  ["P",      "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'", "H2'", "H5'", "H5''"]
                b_coords = [b_p_curr, b_op1, b_op2, b_o5p, b_c5p, b_c4p, b_o4p, b_c3p, b_o3p, b_c2p, b_o2p, b_c1p, b_h2p, b_h5p, b_h5pp]

                # write all sugar atoms
                pdb_atom_number = writeatoms(b_atoms,b_coords,pdb_atom_number,res,chain_num)
                continue        # move on to next iteration


        ## Last Sugar ##
        if str(linesp[1][0])=='G' and (res+1 in chain_breaks):  # and res!=0:
                sugar_atom_num = int(linesp[0])
                b_c4p = [float(linesp[2]),float(linesp[3]),float(linesp[4])]

                for entry in everything:
			entrysp = entry.split()
			if len(entrysp)>1:
	                        if sugar_atom_num-3 < int(entrysp[0]) < sugar_atom_num and str(entrysp[1][0])=='P':
        	                        b_p_curr = [float(entrysp[2]),float(entrysp[3]),float(entrysp[4])] # current phosphate
                	        elif sugar_atom_num-6 < int(entrysp[0]) < sugar_atom_num and (str(entrysp[1][0])=='S' or str(entrysp[1][0])=='G'):
                        	        b_c4p_prev = [float(entrysp[2]),float(entrysp[3]),float(entrysp[4])]  # previous sugar
				elif sugar_atom_num < int(entrysp[0]) < sugar_atom_num+6 and (str(entrysp[1])=='CU' or str(entrysp[1])=='CG'):
        	                        b_sc = [float(entrysp[2]),float(entrysp[3]),float(entrysp[4])]    ## sugar connect
					break

		unit_crossp = orthogvector(b_c4p_prev,b_p_curr,b_c4p)

                b_op1 = [ 1.5*unit_crossp[0] + b_p_curr[0], 1.5*unit_crossp[1] + b_p_curr[1], 1.5*unit_crossp[2] + b_p_curr[2]]
                b_op2 = [ -1.5*unit_crossp[0] + b_p_curr[0], -1.5* unit_crossp[1] + b_p_curr[1], -1.5*unit_crossp[2] + b_p_curr[2]]

                b_o5p = [ (b_c4p[0]-b_p_curr[0])/3 + b_p_curr[0], (b_c4p[1]-b_p_curr[1])/3 + b_p_curr[1], (b_c4p[2]-b_p_curr[2])/3 + b_p_curr[2]]
                b_c5p = [ 2*(b_c4p[0]-b_p_curr[0])/3 + b_p_curr[0], 2*(b_c4p[1]-b_p_curr[1])/3 + b_p_curr[1], 2*(b_c4p[2]-b_p_curr[2])/3 + b_p_curr[2]]
                b_o4p = [ (b_sc[0]-b_c4p[0])/2 + b_c4p[0] , (b_sc[1]-b_c4p[1])/2 + b_c4p[1] , (b_sc[2]-b_c4p[2])/2 + b_c4p[2]]
                b_c3p = [ (b_c4p[0]-b_p_curr[0])/3 + b_c4p[0] , (b_c4p[1]-b_p_curr[1])/3 + b_c4p[1] , (b_c4p[2]-b_p_curr[2])/3 + b_c4p[2]]
                b_o3p = [ 2*(b_c4p[0]-b_p_curr[0])/3 + b_c4p[0] , 2*(b_c4p[1]-b_p_curr[1])/3 + b_c4p[1] , 2*(b_c4p[2]-b_p_curr[2])/3 + b_c4p[2]]
                b_c2p = [ (b_sc[0]-b_c3p[0])/3 + b_c3p[0] , (b_sc[1]-b_c3p[1])/3 + b_c3p[1] , (b_sc[2]-b_c3p[2])/3 + b_c3p[2]]
                b_o2p = [ (b_c2p[0]-b_c5p[0])/3 + b_c2p[0] , (b_c2p[1]-b_c5p[1])/2 + b_c2p[1] , (b_c2p[2]-b_c5p[2])/2 + b_c2p[2]]
                b_c1p = [ 2*(b_sc[0]-b_c3p[0])/3 + b_c3p[0] , 2*(b_sc[1]-b_c3p[1])/3 + b_c3p[1] , 2*(b_sc[2]-b_c3p[2])/3 + b_c3p[2]]

                unit_crossp = orthogvector(b_c4p,b_c2p,b_c1p) # hydrogen on c2' atom

                b_h2p = [ 1.1*unit_crossp[0] + b_c2p[0], 1.1*unit_crossp[1] + b_c2p[1], 1.1*unit_crossp[2] + b_c2p[2]]

		unit_crossp = orthogvector(b_c5p,b_c4p,b_o4p) # hydrogen on c4' atom

                b_h4p = [ 1.1*unit_crossp[0] + b_c4p[0], 1.1*unit_crossp[1] + b_c4p[1], 1.1*unit_crossp[2] + b_c4p[2]]

		unit_crossp = orthogvector(b_o5p,b_c5p,b_o4p) # hydrogens on c5' atom

                b_h5p = [ 1.1*unit_crossp[0] + b_c5p[0], 1.1*unit_crossp[1] + b_c5p[1], 1.1*unit_crossp[2] + b_c5p[2]]
                b_h5pp = [ -1.1*unit_crossp[0] + b_c5p[0], -1.1* unit_crossp[1] + b_c5p[1], -1.1*unit_crossp[2] + b_c5p[2]]


                b_atoms =  ["P",      "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'", "H2'", "H4'", "H5'", "H5''"]
                b_coords = [b_p_curr, b_op1, b_op2, b_o5p, b_c5p, b_c4p, b_o4p, b_c3p, b_o3p, b_c2p, b_o2p, b_c1p, b_h2p, b_h4p, b_h5p, b_h5pp]

                # write all sugar atoms
                pdb_atom_number = writeatoms(b_atoms,b_coords,pdb_atom_number,res,chain_num)
                continue        # move on to next iteration


	## Guanine ##
	if sequence[res]=='G':
		# First get cg coordinates
                if str(linesp[1])=='CG':
                        g_cg = [float(linesp[2]),float(linesp[3]),float(linesp[4])]
                if str(linesp[1])=='O6':
                        g_o6 = [float(linesp[2]),float(linesp[3]),float(linesp[4])]
                if str(linesp[1])=='N2':
                        g_n2 = [float(linesp[2]),float(linesp[3]),float(linesp[4])] # last base atom entry

                        # set atomic positions along triangle of cg base, then minimize separately
                        g_n7 = [ (g_o6[0]-g_cg[0])/3 + g_cg[0] , (g_o6[1]-g_cg[1])/3 + g_cg[1] , (g_o6[2]-g_cg[2])/3 + g_cg[2]]
                        g_c5 = [ 2*(g_o6[0]-g_cg[0])/3 + g_cg[0] , 2*(g_o6[1]-g_cg[1])/3 + g_cg[1] , 2*(g_o6[2]-g_cg[2])/3 + g_cg[2]]
			g_c6 = [ (g_n2[0]-g_o6[0])/4 + g_o6[0] , (g_n2[1]-g_o6[1])/4 + g_o6[1] , (g_n2[2]-g_o6[2])/4 + g_o6[2]]
                        g_n1 = [ 2*(g_n2[0]-g_o6[0])/4 + g_o6[0] , 2*(g_n2[1]-g_o6[1])/4 + g_o6[1] , 2*(g_n2[2]-g_o6[2])/4 + g_o6[2]]
                        g_c2 = [ 3*(g_n2[0]-g_o6[0])/4 + g_o6[0] , 3*(g_n2[1]-g_o6[1])/4 + g_o6[1] , 3*(g_n2[2]-g_o6[2])/4 + g_o6[2]]
                        g_n3 = [ (g_cg[0]-g_n2[0])/4 + g_n2[0] , (g_cg[1]-g_n2[1])/4 + g_n2[1] , (g_cg[2]-g_n2[2])/4 + g_n2[2]]
			g_c4 = [ 2*(g_cg[0]-g_n2[0])/4 + g_n2[0] , 2*(g_cg[1]-g_n2[1])/4 + g_n2[1] , 2*(g_cg[2]-g_n2[2])/4 + g_n2[2]]
			g_n9 = [ 3*(g_cg[0]-g_n2[0])/4 + g_n2[0] , 3*(g_cg[1]-g_n2[1])/4 + g_n2[1] , 3*(g_cg[2]-g_n2[2])/4 + g_n2[2]]

			g_h21 = [ (g_n2[0]-g_o6[0])/4 + g_n2[0] , (g_n2[1]-g_o6[1])/4 + g_n2[1] , (g_n2[2]-g_o6[2])/4 + g_n2[2]]
			g_h22 = [ (g_n2[0]-g_cg[0])/5 + g_n2[0] , (g_n2[1]-g_cg[1])/5 + g_n2[1] , (g_n2[2]-g_cg[2])/5 + g_n2[2]]

                        # combine for easy writing  (N6 becomes N4, CU becomes C6)
                        g_atoms =  ['N9', 'C8', 'N7', 'C5', 'C6', 'O6', 'N1', 'C2', 'N2', 'N3', 'C4', 'H21', 'H22']
                        g_coords = [g_n9, g_cg, g_n7, g_c5, g_c6, g_o6, g_n1, g_c2, g_n2, g_n3, g_c4, g_h21, g_h22]

                        # write all base atoms
			pdb_atom_number = writeatoms(g_atoms,g_coords,pdb_atom_number,res,chain_num)

                        res = res + 1;  # increment residue counter
                        continue        # move on to next iteration

	## Adenine ##
	if sequence[res]=='A':
		# First get CG coordinates
                if str(linesp[1])=='CG':
                        a_cg = [float(linesp[2]),float(linesp[3]),float(linesp[4])]
                if str(linesp[1])=='N6':
                        a_n6 = [float(linesp[2]),float(linesp[3]),float(linesp[4])]
                if str(linesp[1])=='CA':
                        a_ca = [float(linesp[2]),float(linesp[3]),float(linesp[4])] # last base atom entry

                        # set atomic positions along triangle of cg base, then minimize separately
                        a_n7 = [ (a_n6[0]-a_cg[0])/3 + a_cg[0] , (a_n6[1]-a_cg[1])/3 + a_cg[1] , (a_n6[2]-a_cg[2])/3 + a_cg[2]]
                        a_c5 = [ 2*(a_n6[0]-a_cg[0])/3 + a_cg[0] , 2*(a_n6[1]-a_cg[1])/3 + a_cg[1] , 2*(a_n6[2]-a_cg[2])/3 + a_cg[2]]
                        a_c6 = [ (a_ca[0]-a_n6[0])/3 + a_n6[0] , (a_ca[1]-a_n6[1])/3 + a_n6[1] , (a_ca[2]-a_n6[2])/3 + a_n6[2]]
                        a_n1 = [ 2*(a_ca[0]-a_n6[0])/3 + a_n6[0] , 2*(a_ca[1]-a_n6[1])/3 + a_n6[1] , 2*(a_ca[2]-a_n6[2])/3 + a_n6[2]]
                        a_n3 = [ (a_cg[0]-a_ca[0])/4 + a_ca[0] , (a_cg[1]-a_ca[1])/4 + a_ca[1] , (a_cg[2]-a_ca[2])/4 + a_ca[2]]
                        a_c4 = [ 2*(a_cg[0]-a_ca[0])/4 + a_ca[0] , 2*(a_cg[1]-a_ca[1])/4 + a_ca[1] , 2*(a_cg[2]-a_ca[2])/4 + a_ca[2]]
                        a_n9 = [ 3*(a_cg[0]-a_ca[0])/4 + a_ca[0] , 3*(a_cg[1]-a_ca[1])/4 + a_ca[1] , 3*(a_cg[2]-a_ca[2])/4 + a_ca[2]]

                        # combine for easy writing  (N6 becomes N4, CU becomes C6)
                        a_atoms =  ['N9', 'C8', 'N7', 'C5', 'C6', 'N6', 'N1', 'C2', 'N3', 'C4']
                        a_coords = [a_n9, a_cg, a_n7, a_c5, a_c6, a_n6, a_n1, a_ca, a_n3, a_c4]

                        # write all base atoms
                        pdb_atom_number = writeatoms(a_atoms,a_coords,pdb_atom_number,res,chain_num)

                        res = res + 1;  # increment residue counter
                        continue        # move on to next iteration


	## Cytosine ##
        if sequence[res]=='C':
		# First get cg coordinates
                if str(linesp[1])=='O2':
                        c_o2 = [float(linesp[2]),float(linesp[3]),float(linesp[4])]
                if str(linesp[1])=='N6':
                        c_n6 = [float(linesp[2]),float(linesp[3]),float(linesp[4])]
                if str(linesp[1])=='CU':
                        c_cu = [float(linesp[2]),float(linesp[3]),float(linesp[4])] # last base atom entry

                        # set atomic positions along triangle of cg base, then minimize separately
                        c_n1 = [ (c_o2[0]-c_cu[0])/2 + c_cu[0] , (c_o2[1]-c_cu[1])/2 + c_cu[1] , (c_o2[2]-c_cu[2])/2 + c_cu[2]]
                        c_c2 = [ (c_n6[0]-c_o2[0])/3 + c_o2[0] , (c_n6[1]-c_o2[1])/3 + c_o2[1] , (c_n6[2]-c_o2[2])/3 + c_o2[2]]
                        c_n3 = [ 2*(c_n6[0]-c_o2[0])/3 + c_o2[0] , 2*(c_n6[1]-c_o2[1])/3 + c_o2[1] , 2*(c_n6[2]-c_o2[2])/3 + c_o2[2]]
                        c_c4 = [ (c_cu[0]-c_n6[0])/3 + c_n6[0] , (c_cu[1]-c_n6[1])/3 + c_n6[1] , (c_cu[2]-c_n6[2])/3 + c_n6[2]]
                        c_c5 = [ 2*(c_cu[0]-c_n6[0])/3 + c_n6[0] , 2*(c_cu[1]-c_n6[1])/3 + c_n6[1] , 2*(c_cu[2]-c_n6[2])/3 + c_n6[2]]

                        # combine for easy writing  (N6 becomes N4, CU becomes C6)
                        c_atoms =  ['N1', 'C2', 'O2', 'N3', 'C4', 'N4', 'C5', 'C6']
                        c_coords = [c_n1, c_c2, c_o2, c_n3, c_c4, c_n6, c_c5, c_cu]

                        # write all base atoms
			pdb_atom_number = writeatoms(c_atoms,c_coords,pdb_atom_number,res,chain_num)

                        res = res + 1;  # increment residue counter
			continue 	# move on to next iteration


	## Uracil ##
	if sequence[res]=='U':
		# First get cg coordinates
		if str(linesp[1])=='O2':
			u_o2 = [float(linesp[2]),float(linesp[3]),float(linesp[4])]
		if str(linesp[1])=='O6':
			u_o6 = [float(linesp[2]),float(linesp[3]),float(linesp[4])]
		if str(linesp[1])=='CU':
			u_cu = [float(linesp[2]),float(linesp[3]),float(linesp[4])] # last base atom entry

			# set atomic positions along triangle of cg base, then minimize separately
			u_n1 = [ (u_o2[0]-u_cu[0])/2 + u_cu[0] , (u_o2[1]-u_cu[1])/2 + u_cu[1] , (u_o2[2]-u_cu[2])/2 + u_cu[2]]
			u_c2 = [ (u_o6[0]-u_o2[0])/3 + u_o2[0] , (u_o6[1]-u_o2[1])/3 + u_o2[1] , (u_o6[2]-u_o2[2])/3 + u_o2[2]]
			u_n3 = [ 2*(u_o6[0]-u_o2[0])/3 + u_o2[0] , 2*(u_o6[1]-u_o2[1])/3 + u_o2[1] , 2*(u_o6[2]-u_o2[2])/3 + u_o2[2]]
			u_c4 = [ (u_cu[0]-u_o6[0])/3 + u_o6[0] , (u_cu[1]-u_o6[1])/3 + u_o6[1] , (u_cu[2]-u_o6[2])/3 + u_o6[2]]
			u_c5 = [ 2*(u_cu[0]-u_o6[0])/3 + u_o6[0] , 2*(u_cu[1]-u_o6[1])/3 + u_o6[1] , 2*(u_cu[2]-u_o6[2])/3 + u_o6[2]]

			# combine for easy writing  (O6 becomes O4, CU becomes C6)
			u_atoms =  ['N1', 'C2', 'O2', 'N3', 'C4', 'O4', 'C5', 'C6']
			u_coords = [u_n1, u_c2, u_o2, u_n3, u_c4, u_o6, u_c5, u_cu]

			# write all base atoms
			pdb_atom_number = writeatoms(u_atoms,u_coords,pdb_atom_number,res,chain_num)
			
			res = res + 1;  # increment residue counter
			continue	# move on to next iteration
