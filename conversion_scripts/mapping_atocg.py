#!/usr/bin/env python
import sys

# Red warning = definite problem
# Blue warning = probable problem
# Green warning = possible problem

#Script to convert normal pdb to pdb for coarse-grained

if len(sys.argv) != 1:
	interest = str(sys.argv[1])
else:
	interest=str(raw_input("Please enter the pdb filename: "))

f=open(interest,"r")
everything=f.readlines()
f.close()
p=0
g=open(interest+".pdb","w+")
g.write("HEADER"+"\n")
g.write("COMPND"+"\n")
g.write("SOURCE"+"\n")
h=0
i=0
j=1
k=1
phos=[1,10]
sugar=[2,11]
term_sugar=[12,13]
chain=["A","B","C","D","E","F","G","H"]
term_res = []; term_chain = [];
previous_one = -5;
biomt = 0
#first get the # of chains and the residue numbers of chain ends after the starting point
ok_residues = ['A','C','G','U']
others = ['MG','CA','NA','CL','K','TL','BA','NH4','HOH'] #TL=Thallium; Ions and other heteroatoms in pdb
for line in everything:
	if line[0:4] == "ATOM" or line[0:6] == "HETATM":
		chain_entry = line[20:22].strip()
		resname = line[17:20].strip()
		resid = int(line[22:26])
		if len(chain_entry)!=0 and chain_entry!=chain[i] and resname not in others:
			term_res.append(previous_one)
			term_res.append(resid)
			term_chain.append(chain[i])
			term_chain.append(chain_entry)

			i=chain.index(chain_entry)


		elif resname not in others:
			ending = resid 
			end_chain = chain[i]

		previous_one = resid;
		
		if resname not in ok_residues and resname not in others:
			print "Red Warning: unrecognized nucleotide '"+resname+"' present on chain "+chain_entry+" at residue %3d" % resid + "; will attempt to ignore"
		
	elif line[0:6] == "ENDMDL":
		break
	elif line[13:18] == "BIOMT":
		biomt = 1		

term_res.append(ending)
term_chain.append(end_chain)
#print term_res, term_chain

if i > 1:
	print "Blue Warning: more than two chains, possibly more than one structure; will attempt to ignore"
if biomt == 1:
        print "Green Warning: possible multimeric structure; will attempt to ignore"

## Now, collect and write coarse-grained atom entries ##
for line in everything:
	if (line[0:4] == "ATOM" or line[0:6] == "HETATM") and line[16]==" " or (line[0:4] == "ATOM" or line[0:6] == "HETATM") and line[16]=="A": #if have multiple coordinates for one atom, choose the "A" coordinate
       		if int(line[7:11])>p:
			p = p + 1
			if line[13:16].strip()=="P":
				y=abs((j%2)-2)
				g.write("HETATM"+'%5d' % (k)+" "+'%3s' % "P"+str(y)+" "+'%3d' % (phos[y-1])+line[20:22]+line[22:26]+line[26:54]+"\n")
				k=k+1

			if line[13:16]=="C4'":
				if j==1:  			## beginning sugar
					g.write("HETATM"+'%5d' % (k)+" "+'%4s' % "G1"+" "+'%3d' % 12 +line[20:22]+line[22:26]+line[26:54]+"\n")
				elif int(line[22:26]) not in term_res:    ## line[22:26] is res_id
					## sugar in middle of chain not overlapping terminal res_ids
					d=abs((j%2)-2)
                                	g.write("HETATM"+'%5d' % (k)+" "+'%3s' % "S"+str(d)+" "+'%3d' % (sugar[d-1])+line[20:22]+line[22:26]+line[26:54]+"\n")
				elif int(line[22:26]) in term_res:  ## sugar resid included in terminal res_id
					term_true = 0
					indices = []
					indices = [i for i , x in enumerate(term_res) if x==int(line[22:26])]
					for i in range(len(indices)):  
						if int(line[22:26])==term_res[indices[i]] and line[20:22].strip()==term_chain[indices[i]]:
							term_true = 1 

					if term_true==0:  ## sugar in middle of chain, overlaps terminal res_id
						d=abs((j%2)-2)
                                        	g.write("HETATM"+'%5d' % (k)+" "+'%3s' % "S"+str(d)+" "+'%3d' % (sugar[d-1])+line[20:22]+line[22:26]+line[26:54]+"\n")
						
					elif term_true==1:  ## sugar is terminal residue
						#print int(line[22:26])
                                		d=abs((j%2)-2)
                                		g.write("HETATM"+'%5d' % (k)+" "+'%3s' % "G"+str(d)+" "+'%3d' % (term_sugar[d-1])+line[20:22]+line[22:26]+line[26:54]+"\n")

				
				j=j+1
				k=k+1

			elif line[19]=="A":
				if line[12:15]==" N6":
					g.write("HETATM"+'%5d' % (k)+" "+'%4s' % "N6"+" "+'%3d' % 4 +line[20:22]+line[22:26]+line[26:54]+"\n")
					k=k+1

				if line[12:15]==" C8":
					g.write("HETATM"+'%5d' % (k)+" "+'%4s' % "CG"+" "+'%3d' % 3 +line[20:22]+line[22:26]+line[26:54]+"\n")
					k=k+1
			
				if line[12:16]==" C2 ":
					g.write("HETATM"+'%5d' % (k)+" "+'%4s' % "CA"+" "+'%3d' % 9 +line[20:22]+line[22:26]+line[26:54]+"\n")
					k=k+1
			
			elif line[19]=="G":
                                if line[12:15]==" O6":
					g.write("HETATM"+'%5d' % (k)+" "+'%4s' % "O6"+" "+'%3d' % 6 +line[20:22]+line[22:26]+line[26:54]+"\n")
					k=k+1

				if line[12:15]==" N2":
					g.write("HETATM"+'%5d' % (k)+" "+'%4s' % "N2"+" "+'%3d' % 5 +line[20:22]+line[22:26]+line[26:54]+"\n")
					k=k+1

				if line[12:15]==" C8":
					g.write("HETATM"+'%5d' % (k)+" "+'%4s' % "CG"+" "+'%3d' % 3 +line[20:22]+line[22:26]+line[26:54]+"\n")
					k=k+1

			elif line[19]=="C":
                                if line[12:16]==" O2 ":
					g.write("HETATM"+'%5d' % (k)+" "+'%4s' % "O2"+" "+'%3d' % 7 +line[20:22]+line[22:26]+line[26:54]+"\n")
					k=k+1
		
				if line[12:15]==" C6":
					g.write("HETATM"+'%5d' % (k)+" "+'%4s' % "CU"+" "+'%3d' % 8 +line[20:22]+line[22:26]+line[26:54]+"\n")
					k=k+1
			
				if line[12:15]==" N4":
					g.write("HETATM"+'%5d' % (k)+" "+'%4s' % "N6"+" "+'%3d' % 4 +line[20:22]+line[22:26]+line[26:54]+"\n")
					k=k+1

			elif line[19]=="U":
                               	if line[12:16]==" O2 ":
					g.write("HETATM"+'%5d' % (k)+" "+'%4s' % "O2"+" "+'%3d' % 7 +line[20:22]+line[22:26]+line[26:54]+"\n")
					k=k+1

				if line[12:15]==" C6":
					g.write("HETATM"+'%5d' % (k)+" "+'%4s' % "CU"+" "+'%3d' % 8 +line[20:22]+line[22:26]+line[26:54]+"\n")
					k=k+1

				if line[12:16]==" O4 ":
					g.write("HETATM"+'%5d' % (k)+" "+'%4s' % "O6"+" "+'%3d' % 6 +line[20:22]+line[22:26]+line[26:54]+"\n")	
					k=k+1
	
#	if line[0:6] == "CONECT":
#		g.write(line)

	elif line[0:6]== "ENDMDL":
		break

g.close()

#g=open(interest+".pdb","r+")
#goods=g.readlines()
#goods.reverse()
#g.close()
#g=open(interest+".pdb","w+")
#dd=abs((j-1)%2-2)
#ee=["12","13"]
#q=0
#for line in goods:
#	if line[12:15]==" S"+str(dd) and q==0:
#		g.write("HETATM"+line[6:13]+"G"+str(dd)+" "+line[28:30]+ee[dd-1]+line[26:31]+str(ch)+line[26:53]+"\n")
#		q=1
#		print "Hello there line 127"
#
#	else:
#		g.write(line)
#
#g.close()
#h=open(interest+".pdb","r+")
#newgoods=h.readlines()
#newgoods.reverse()
#h.close()
#h=open(interest+".pdb","w+")
#for line in newgoods:
#	h.write(line)
#
#h.close()
#
#d=open(interest+".pdb","r+")
#books=d.readlines()
#d.close()
#
#k=k-1;
#g=open(interest+".pdb","a+")
#for line in books:
#	if line[0:6] == "HETATM":
#		if line[13:15]=="CG":
#			ua=int(line[7:12])
#			if ua-1==1:
#				g.write("CONECT"+'%5d' % (ua-1)+'%5d' % (ua)+'%5d' % (ua+3)+"\n")
#				g.write("CONECT"+'%5d' % (ua)+'%5d' % (ua-1)+'%5d' % (ua+1)+'%5d' % (ua+2)+"\n")
#                                g.write("CONECT"+'%5d' % (ua+1)+'%5d' % (ua)+'%5d' % (ua+2)+"\n")
#                                g.write("CONECT"+'%5d' % (ua+2)+'%5d' % (ua)+'%5d' % (ua+1)+"\n")
#			
#			elif ua+3>k:
#				g.write("CONECT"+'%5d' % (ua-1)+'%5d' % (ua-2)+'%5d' % (ua)+"\n")
#				g.write("CONECT"+'%5d' % (ua)+'%5d' % (ua-1)+'%5d' % (ua+1)+'%5d' % (ua+2)+"\n")
#				g.write("CONECT"+'%5d' % (ua+1)+'%5d' % (ua)+'%5d' % (ua+2)+"\n")
#				g.write("CONECT"+'%5d' % (ua+2)+'%5d' % (ua)+'%5d' % (ua+1)+"\n")
#
#			else:
#				g.write("CONECT"+'%5d' % (ua-1)+'%5d' % (ua-2)+'%5d' % (ua)+'%5d' % (ua+3)+"\n")
#                                g.write("CONECT"+'%5d' % (ua)+'%5d' % (ua-1)+'%5d' % (ua+1)+'%5d' % (ua+2)+"\n")
#                                g.write("CONECT"+'%5d' % (ua+1)+'%5d' % (ua)+'%5d' % (ua+2)+"\n")
#				g.write("CONECT"+'%5d' % (ua+2)+'%5d' % (ua)+'%5d' % (ua+1)+"\n")
#
#		if line[13:15]=="CU":
#                        ub=int(line[7:11])
#			if ub-3==1:
#				g.write("CONECT"+'%5d' % (ub-3)+'%5d' % (ub)+'%5d' % (ub+1)+"\n")
#                        	g.write("CONECT"+'%5d' % (ub-2)+'%5d' % (ub-1)+'%5d' % (ub)+"\n")
#	                       	g.write("CONECT"+'%5d' % (ub-1)+'%5d' % (ub-2)+'%5d' % (ub)+"\n")
#                        	g.write("CONECT"+'%5d' % (ub)+'%5d' % (ub-3)+'%5d' % (ub-2)+'%5d' % (ub-1)+"\n")
#			elif ub+1>k:
#                        	g.write("CONECT"+'%5d' % (ub-3)+'%5d' % (ub-4)+'%5d' % (ub)+"\n")
#                        	g.write("CONECT"+'%5d' % (ub-2)+'%5d' % (ub-1)+'%5d' % (ub)+"\n")
#                        	g.write("CONECT"+'%5d' % (ub-1)+'%5d' % (ub-2)+'%5d' % (ub)+"\n")
#                        	g.write("CONECT"+'%5d' % (ub)+'%5d' % (ub-3)+'%5d' % (ub-2)+'%5d' % (ub-1)+"\n")
#			else:
#				g.write("CONECT"+'%5d' % (ub-3)+'%5d' % (ub-4)+'%5d' % (ub)+'%5d' % (ub+1)+"\n")
#                        	g.write("CONECT"+'%5d' % (ub-2)+'%5d' % (ub-1)+'%5d' % (ub)+"\n")
#                        	g.write("CONECT"+'%5d' % (ub-1)+'%5d' % (ub-2)+'%5d' % (ub)+"\n")
#                        	g.write("CONECT"+'%5d' % (ub)+'%5d' % (ub-3)+'%5d' % (ub-2)+'%5d' % (ub-1)+"\n")
#		
#		if line[13]=="P":
#			uk=int(line[7:11])
#			if uk-4<=0:
#				g.write("CONECT"+line[6:11]+ '%5d' % (uk+1)+"\n")
#			elif uk+1>k:
#				g.write("CONECT"+line[6:11]+ '%5d' % (uk-4)+"\n")
#			else:
#				g.write("CONECT"+line[6:11]+ '%5d' % (uk-4)+'%5d' % (uk+1)+"\n")
#
#g.close()
#
#g=open(interest+".pdb","a+")
#g.write("END"+"\n")
#g.close()

