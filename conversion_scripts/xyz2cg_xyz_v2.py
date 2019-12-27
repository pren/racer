#!/usr/bin/env python
import sys

if len(sys.argv) != 1:
        interest = str(sys.argv[1])
else:
        interest=str(raw_input("Please enter the pdb filename: "))

f=open(interest,"r")
everything=f.readlines()
f.close()

num_lines = sum(1 for line in open(interest))

k=num_lines-1

###### Terminal sugars/phosphates ###################

beg = []; end = []; dummy = 0; phos_beg = [-5];
for line in everything:
	jj = line.split()
	if jj!=[] and len(jj)==6 and ( jj[1]=="G1" or jj[1]=="G2" ):
		if dummy % 2 == 0:
			beg.append(int(jj[0]))
			dummy = dummy + 1
			if int(jj[0])-1!=0:
                        	for inner in everything:
					aa = inner.split()
                                	if aa!=[] and len(aa)==6 and (int(jj[0])-1)==int(aa[0]) and (aa[1].strip()=="P1" or aa[1].strip()=="P2"):
						phos_beg.append(int(aa[0]))
		elif dummy % 2 == 1:
			end.append(int(jj[0]))
			dummy = dummy + 1

######################################################



g=open(interest+".xyz","w+")
g.write('%6d' % (k)+"\n")

###### Write_atom subroutine #######################

def write_atom(entry,connects):

	if len(connects)==1:
		g.write('%6s' % (entry[0])+'%5s' % (entry[1])+'%12s' % entry[2]+'%12s' % entry[3]+'%12s' % entry[4]+'%6s' % entry[5]+'%6d' % (connects[0])+"\n")

	elif len(connects)==2:
		g.write('%6s' % (entry[0])+'%5s' % (entry[1])+'%12s' % entry[2]+'%12s' % entry[3]+'%12s' % entry[4]+'%6s' % entry[5]+'%6d' % (connects[0])+'%6d' % (connects[1])+"\n")

	elif len(connects)==3:
		g.write('%6s' % (entry[0])+'%5s' % (entry[1])+'%12s' % entry[2]+'%12s' % entry[3]+'%12s' % entry[4]+'%6s' % entry[5]+'%6d' % (connects[0])+'%6d' % (connects[1])+ '%6d' % (connects[2])+ "\n")

	elif len(connects)==4:
		g.write('%6s' % (entry[0])+'%5s' % (entry[1])+'%12s' % entry[2]+'%12s' % entry[3]+'%12s' % entry[4]+'%6s' % entry[5]+'%6d' % (connects[0])+'%6d' % (connects[1])+ '%6d' % (connects[2])+ '%6d' % (connects[3])+ "\n")

####################################################


sugars = [];

for line in everything:
	mm = line.split()

	if mm!=[] and len(mm)==6 and (mm[1][0]=="S" or mm[1][0]=="G"):
		sugars.append(int(mm[0]))
		continue

	## Normal purine order ##
	if mm!=[] and len(mm)==6 and mm[1]=="CG" and int(mm[0])-1 in sugars:
		ua=int(mm[0])
#		g.write('%5d' % (ua)+ '%5d' % (num_lines)+"\n")
		
		## beginning of chain ##
		if ua-1==1 or (ua-1) in beg:
			for lin in everything:
				nn = lin.split()
				if nn!=[] and int(nn[0])==ua-1 and (ua-2) not in phos_beg:
					write_atom(nn[:],[ua,ua+3])

				elif nn!=[] and int(nn[0])==ua-1 and (ua-2) in phos_beg:
					write_atom(nn[:],[ua-2,ua,ua+3])

			write_atom(mm[:],[ua-1,ua+1,ua+2])

			for lin in everything:
				nn = lin.split()
				if nn!=[] and int(nn[0])==ua+1:
					write_atom(nn[:],[ua,ua+2])

                                elif nn!=[] and int(nn[0])==ua+2:
					write_atom(nn[:],[ua,ua+1])

		## end of chain ##
              	elif (ua+3)>k or (ua-1) in end:
			for lin in everything:
				nn = lin.split()
                                if nn!=[] and int(nn[0])==ua-1:
					write_atom(nn[:],[ua-2,ua])

			write_atom(mm[:],[ua-1,ua+1,ua+2])

			for lin in everything:
				nn = lin.split()
                                if nn!=[] and int(nn[0])==ua+1:
					write_atom(nn[:],[ua,ua+2])
                
		                elif nn!=[] and int(nn[0])==ua+2 and len(nn)==6: #lin[8:11].strip()!="": # the empty string negates the first line
					write_atom(nn[:],[ua,ua+1])
	
		## all others ##
		else: 
			for lin in everything:
				nn = lin.split()
                                if nn!=[] and int(nn[0])==ua-1:
					write_atom(nn[:],[ua-2,ua,ua+3])

			write_atom(mm[:],[ua-1,ua+1,ua+2])

			for lin in everything:
				nn = lin.split()
                                if nn!=[] and int(nn[0])==ua+1:
					write_atom(nn[:],[ua,ua+2])

                                elif nn!=[] and int(nn[0])==ua+2:
					write_atom(nn[:],[ua,ua+1])


	## Pyrimidines ##
	if mm!=[] and len(mm)==6 and mm[1]=="CU":
        	ub=int(mm[0])

		## beginning of chain ##
             	if ub-3==1 or (ub-3) in beg:
			for lin in everything:
				nn = lin.split()
                                if nn!=[] and int(nn[0])==ub-3 and (ub-4) not in phos_beg:
					write_atom(nn[:],[ub,ub+1])

				elif nn!=[] and int(nn[0])==ub-3 and (ub-4) in phos_beg:
					write_atom(nn[:],[ub-4,ub,ub+1])
				
				elif nn!=[] and int(nn[0])==ub-2:
					write_atom(nn[:],[ub-1,ub])

                                elif nn!=[] and int(nn[0])==ub-1:
					write_atom(nn[:],[ub-2,ub])

			write_atom(mm[:],[ub-3,ub-2,ub-1])
         	
		## end of chain ##
		elif ub+1>k or (ub-3) in end:
			for lin in everything:
				nn = lin.split()
                                if nn!=[] and int(nn[0])==ub-3:
					write_atom(nn[:],[ub-4,ub])

  			 	elif nn!=[] and int(nn[0])==ub-2:
					write_atom(nn[:],[ub-1,ub])

                                elif nn!=[] and int(nn[0])==ub-1:
					write_atom(nn[:],[ub-2,ub])

			write_atom(mm[:],[ub-3,ub-2,ub-1])

		## all others ##
		else:
			for lin in everything:
				nn = lin.split()
                                if nn!=[] and int(nn[0])==ub-3:
					write_atom(nn[:],[ub-4,ub,ub+1])

				elif nn!=[] and int(nn[0])==ub-2:
					write_atom(nn[:],[ub-1,ub])

                                elif nn!=[] and int(nn[0])==ub-1:
					write_atom(nn[:],[ub-2,ub])

			write_atom(mm[:],[ub-3,ub-2,ub-1])


	## Abnormal purine order ##
	if mm!=[] and len(mm)==6 and mm[1]=="CG" and int(mm[0])-3 in sugars:
                ua=int(mm[0])
#               g.write('%5d' % (ua)+ '%5d' % (num_lines)+"\n")

		guanine = {};
		adenine = {};
		for record in everything:
			qq = record.split()
			if qq!=[] and (int(qq[0]) == ua-2 or int(qq[0]) == ua-1):
				if qq[1]=="O6":
					guanine[qq[1]]=qq[:]

				elif qq[1]=="N2":
					guanine[qq[1]]=qq[:]

				elif qq[1]=="N6":
					adenine[qq[1]]=qq[:]

				elif qq[1]=="CA":
					adenine[qq[1]]=qq[:]

		## beginning of chain ##
                if ua-3==1 or (ua-3) in beg:
                        for lin in everything:
                                nn = lin.split()
                                if nn!=[] and int(nn[0])==ua-3 and (ua-4) not in phos_beg:
					write_atom(nn[:],[ua-2,ua+1])

                                elif nn!=[] and int(nn[0])==ua-3 and (ua-4) in phos_beg:
					write_atom(nn[:],[ua-4,ua-2,ua+1])

			write_atom([str(ua-2),mm[1],mm[2],mm[3],mm[4],mm[5]],[ua-3,ua-1,ua])

                        if adenine=={}:
				write_atom([str(ua-1),guanine["O6"][1],guanine["O6"][2],guanine["O6"][3],guanine["O6"][4],guanine["O6"][5]],[ua-2,ua])
				write_atom([str(ua),guanine["N2"][1],guanine["N2"][2],guanine["N2"][3],guanine["N2"][4],guanine["N2"][5]],[ua-2,ua-1])

			elif guanine=={}:
				write_atom([str(ua-1),adenine["N6"][1],adenine["N6"][2],adenine["N6"][3],adenine["N6"][4],adenine["N6"][5]],[ua-2,ua])
                                write_atom([str(ua),adenine["CA"][1],adenine["CA"][2],adenine["CA"][3],adenine["CA"][4],adenine["CA"][5]],[ua-2,ua-1])

		## end of chain ##
		elif (ua+1)>k or (ua-3) in end:
                        for lin in everything:
                                nn = lin.split()
                                if nn!=[] and int(nn[0])==ua-3:
					write_atom(nn[:],[ua-4,ua-2])

			write_atom([str(ua-2),mm[1],mm[2],mm[3],mm[4],mm[5]],[ua-3,ua-1,ua])

			if adenine=={}:
                                write_atom([str(ua-1),guanine["O6"][1],guanine["O6"][2],guanine["O6"][3],guanine["O6"][4],guanine["O6"][5]],[ua-2,ua])
                                write_atom([str(ua),guanine["N2"][1],guanine["N2"][2],guanine["N2"][3],guanine["N2"][4],guanine["N2"][5]],[ua-2,ua-1])

                        elif guanine=={}:
                                write_atom([str(ua-1),adenine["N6"][1],adenine["N6"][2],adenine["N6"][3],adenine["N6"][4],adenine["N6"][5]],[ua-2,ua])
                                write_atom([str(ua),adenine["CA"][1],adenine["CA"][2],adenine["CA"][3],adenine["CA"][4],adenine["CA"][5]],[ua-2,ua-1])
 
		## all others ##
                else:
                        for lin in everything:
                                nn = lin.split()
                                if nn!=[] and int(nn[0])==ua-3:
					write_atom(nn[:],[ua-4,ua-2,ua+1])
    
			write_atom([str(ua-2),mm[1],mm[2],mm[3],mm[4],mm[5]],[ua-3,ua-1,ua])

			if adenine=={}:
                                write_atom([str(ua-1),guanine["O6"][1],guanine["O6"][2],guanine["O6"][3],guanine["O6"][4],guanine["O6"][5]],[ua-2,ua])
                                write_atom([str(ua),guanine["N2"][1],guanine["N2"][2],guanine["N2"][3],guanine["N2"][4],guanine["N2"][5]],[ua-2,ua-1])

                        elif guanine=={}:
                                write_atom([str(ua-1),adenine["N6"][1],adenine["N6"][2],adenine["N6"][3],adenine["N6"][4],adenine["N6"][5]],[ua-2,ua])
                                write_atom([str(ua),adenine["CA"][1],adenine["CA"][2],adenine["CA"][3],adenine["CA"][4],adenine["CA"][5]],[ua-2,ua-1])



	## Phosphates ##
    	if mm!=[] and len(mm)==6 and (mm[1]=="P1" or mm[1]=="P2"):
         	uk=int(mm[0])
               	if uk not in phos_beg:
			write_atom(mm[:],[uk-4,uk+1])

		elif uk in phos_beg:
			write_atom(mm[:],[uk+1])

g.close()
