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

beg = []; end = []; dummy = 0; phos_beg = [-5];
for liner in everything:
	jj = liner.split()
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

g=open(interest+".xyz","w+")
g.write('%6d' % (k)+"\n")

for line in everything:
	mm = line.split()
	if mm!=[] and len(mm)==6 and mm[1]=="CG":
		ua=int(mm[0])
#		g.write('%5d' % (ua)+ '%5d' % (num_lines)+"\n")
		
		
		if ua-1==1 or (ua-1) in beg:
			for lin in everything:
				nn = lin.split()
				if nn!=[] and int(nn[0])==ua-1 and (ua-2) not in phos_beg:
					g.write('%6s' % (nn[0])+'%5s' % (nn[1])+'%12s' % nn[2]+'%12s' % nn[3]+'%12s' % nn[4]+'%6s' % nn[5]+'%6d' % (ua)+'%6d' % (ua+3)+"\n")
				elif nn!=[] and int(nn[0])==ua-1 and (ua-2) in phos_beg:
					g.write('%6s' % (nn[0])+'%5s' % (nn[1])+'%12s' % nn[2]+'%12s' % nn[3]+'%12s' % nn[4]+'%6s' % nn[5]+'%6d' % (ua-2)+'%6d' % (ua)+'%6d' % (ua+3)+"\n")

	                g.write('%6s' % (mm[0])+'%5s' % (mm[1])+'%12s' % mm[2]+'%12s' % mm[3]+'%12s' % mm[4]+'%6s' % mm[5]+'%6d' % (ua-1)+'%6d' % (ua+1)+'%6d' % (ua+2)+"\n")
			for lin in everything:
				nn = lin.split()
				if nn!=[] and int(nn[0])==ua+1:
					g.write('%6s' % (nn[0])+'%5s' % (nn[1])+'%12s' % nn[2]+'%12s' % nn[3]+'%12s' % nn[4]+'%6s' % nn[5]+'%6d' % (ua)+'%6d' % (ua+2)+"\n")

                                elif nn!=[] and int(nn[0])==ua+2:
                                        g.write('%6s' % (nn[0])+'%5s' % (nn[1])+'%12s' % nn[2]+'%12s' % nn[3]+'%12s' % nn[4]+'%6s' % nn[5]+'%6d' % (ua)+'%6d' % (ua+1)+"\n")

              	elif (ua+3)>k or (ua-1) in end:
			for lin in everything:
				nn = lin.split()
                                if nn!=[] and int(nn[0])==ua-1:
                                        g.write('%6s' % (nn[0])+'%5s' % (nn[1])+'%12s' % nn[2]+'%12s' % nn[3]+'%12s' % nn[4]+'%6s' % nn[5]+'%6d' % (ua-2)+'%6d' % (ua)+"\n")

	                g.write('%6s' % (mm[0])+'%5s' % (mm[1])+'%12s' % mm[2]+'%12s' % mm[3]+'%12s' % mm[4]+'%6s' % mm[5]+'%6d' % (ua-1)+'%6d' % (ua+1)+'%6d' % (ua+2)+"\n")
			for lin in everything:
				nn = lin.split()
                                if nn!=[] and int(nn[0])==ua+1:
                                        g.write('%6s' % (nn[0])+'%5s' % (nn[1])+'%12s' % nn[2]+'%12s' % nn[3]+'%12s' % nn[4]+'%6s' % nn[5]+'%6d' % (ua)+'%6d' % (ua+2)+"\n")
                
		                elif nn!=[] and int(nn[0])==ua+2 and len(nn)==6: #lin[8:11].strip()!="": # the empty string negates the first line
                                        g.write('%6s' % (nn[0])+'%5s' % (nn[1])+'%12s' % nn[2]+'%12s' % nn[3]+'%12s' % nn[4]+'%6s' % nn[5]+'%6d' % (ua)+'%6d' % (ua+1)+"\n")
              	
		else:
			for lin in everything:
				nn = lin.split()
                                if nn!=[] and int(nn[0])==ua-1:
                                        g.write('%6s' % (nn[0])+'%5s' % (nn[1])+'%12s' % nn[2]+'%12s' % nn[3]+'%12s' % nn[4]+'%6s' % nn[5]+'%6d' % (ua-2)+'%6d' % (ua)+'%6d' % (ua+3)+"\n")
        
			g.write('%6s' % (mm[0])+'%5s' % (mm[1])+'%12s' % mm[2]+'%12s' % mm[3]+'%12s' % mm[4]+'%6s' % mm[5]+'%6d' % (ua-1)+'%6d' % (ua+1)+'%6d' % (ua+2)+"\n") 
			for lin in everything:
				nn = lin.split()
                                if nn!=[] and int(nn[0])==ua+1:
                                        g.write('%6s' % (nn[0])+'%5s' % (nn[1])+'%12s' % nn[2]+'%12s' % nn[3]+'%12s' % nn[4]+'%6s' % nn[5]+'%6d' % (ua)+'%6d' % (ua+2)+"\n")

                                elif nn!=[] and int(nn[0])==ua+2:
                                        g.write('%6s' % (nn[0])+'%5s' % (nn[1])+'%12s' % nn[2]+'%12s' % nn[3]+'%12s' % nn[4]+'%6s' % nn[5]+'%6d' % (ua)+'%6d' % (ua+1)+"\n")

	if mm!=[] and len(mm)==6 and mm[1]=="CU":
        	ub=int(mm[0])
             	if ub-3==1 or (ub-3) in beg:
			for lin in everything:
				nn = lin.split()
                                if nn!=[] and int(nn[0])==ub-3 and (ub-4) not in phos_beg:
                                        g.write('%6s' % (nn[0])+'%5s' % (nn[1])+'%12s' % nn[2]+'%12s' % nn[3]+'%12s' % nn[4]+'%6s' % nn[5]+'%6d' % (ub)+'%6d' % (ub+1)+"\n")
				elif nn!=[] and int(nn[0])==ub-3 and (ub-4) in phos_beg:
					g.write('%6s' % (nn[0])+'%5s' % (nn[1])+'%12s' % nn[2]+'%12s' % nn[3]+'%12s' % nn[4]+'%6s' % nn[5]+'%6d' % (ub-4)+'%6d' % (ub)+'%6d' % (ub+1)+"\n")
				
				elif nn!=[] and int(nn[0])==ub-2:
                                        g.write('%6s' % (nn[0])+'%5s' % (nn[1])+'%12s' % nn[2]+'%12s' % nn[3]+'%12s' % nn[4]+'%6s' % nn[5]+'%6d' % (ub-1)+'%6d' % (ub)+"\n")

                                elif nn!=[] and int(nn[0])==ub-1:
                                        g.write('%6s' % (nn[0])+'%5s' % (nn[1])+'%12s' % nn[2]+'%12s' % nn[3]+'%12s' % nn[4]+'%6s' % nn[5]+'%6d' % (ub-2)+'%6d' % (ub)+"\n")

         		g.write('%6s' % (mm[0])+'%5s' % (mm[1])+'%12s' % mm[2]+'%12s' % mm[3]+'%12s' % mm[4]+'%6s' % mm[5]+'%6d' % (ub-3)+'%6d' % (ub-2)+'%6d' % (ub-1)+"\n")
         	
		elif ub+1>k or (ub-3) in end:
			for lin in everything:
				nn = lin.split()
                                if nn!=[] and int(nn[0])==ub-3:
                                        g.write('%6s' % (nn[0])+'%5s' % (nn[1])+'%12s' % nn[2]+'%12s' % nn[3]+'%12s' % nn[4]+'%6s' % nn[5]+'%6d' % (ub-4)+'%6d' % (ub)+"\n")

  			 	elif nn!=[] and int(nn[0])==ub-2:
                                        g.write('%6s' % (nn[0])+'%5s' % (nn[1])+'%12s' % nn[2]+'%12s' % nn[3]+'%12s' % nn[4]+'%6s' % nn[5]+'%6d' % (ub-1)+'%6d' % (ub)+"\n")

                                elif nn!=[] and int(nn[0])==ub-1:
                                        g.write('%6s' % (nn[0])+'%5s' % (nn[1])+'%12s' % nn[2]+'%12s' % nn[3]+'%12s' % nn[4]+'%6s' % nn[5]+'%6d' % (ub-2)+'%6d' % (ub)+"\n")

			g.write('%6s' % (mm[0])+'%5s' % (mm[1])+'%12s' % mm[2]+'%12s' % mm[3]+'%12s' % mm[4]+'%6s' % mm[5]+'%6d' % (ub-3)+'%6d' % (ub-2)+'%6d' % (ub-1)+"\n")
             	
		else:
			for lin in everything:
				nn = lin.split()
                                if nn!=[] and int(nn[0])==ub-3:
                                        g.write('%6s' % (nn[0])+'%5s' % (nn[1])+'%12s' % nn[2]+'%12s' % nn[3]+'%12s' % nn[4]+'%6s' % nn[5]+'%6d' % (ub-4)+'%6d' % (ub)+'%6d' % (ub+1)+"\n")

				elif nn!=[] and int(nn[0])==ub-2:
                                        g.write('%6s' % (nn[0])+'%5s' % (nn[1])+'%12s' % nn[2]+'%12s' % nn[3]+'%12s' % nn[4]+'%6s' % nn[5]+'%6d' % (ub-1)+'%6d' % (ub)+"\n")

                                elif nn!=[] and int(nn[0])==ub-1:
                                        g.write('%6s' % (nn[0])+'%5s' % (nn[1])+'%12s' % nn[2]+'%12s' % nn[3]+'%12s' % nn[4]+'%6s' % nn[5]+'%6d' % (ub-2)+'%6d' % (ub)+"\n")

  	                g.write('%6s' % (mm[0])+'%5s' % (mm[1])+'%12s' % mm[2]+'%12s' % mm[3]+'%12s' % mm[4]+'%6s' % mm[5]+'%6d' % (ub-3)+'%6d' % (ub-2)+'%6d' % (ub-1)+"\n")

    	if mm!=[] and len(mm)==6 and (mm[1]=="P1" or mm[1]=="P2"):
         	uk=int(mm[0])
               	if uk not in phos_beg:
			g.write('%6s' % (mm[0])+'%5s' % (mm[1])+'%12s' % mm[2]+'%12s' % mm[3]+'%12s' % mm[4]+'%6s' % mm[5]+'%6d' % (uk-4)+'%6d' % (uk+1)+"\n")
		elif uk in phos_beg:
			g.write('%6s' % (mm[0])+'%5s' % (mm[1])+'%12s' % mm[2]+'%12s' % mm[3]+'%12s' % mm[4]+'%6s' % mm[5]+'%6d' % (uk+1)+"\n")

g.close()
