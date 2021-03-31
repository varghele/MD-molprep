###-----------------------------------------------------------------------
#Program for MD Prep, used before placer.py
#
#Steps:
# 1 - read in molweights dictionary
# 2 - read .pdb in inp_rotater with name N_LIG.pdb  #--N is the number of copies you want, LIG is the Ligand name you want the molecule to have
# 3 - calculate center of mass
# - translate molecule to COM
# 4 - For Loop:
# 		for K in N:
#			-5- get random rotation quaternion
#			
#			-6- rotate molecule
#			-7- insert LIG and K
#			-8- save to inp_placer
###-----------------------------------------------------------------------

import os
from molweights import elements_dict
import random
from math import cos, sin, pi

# Function to get atom type from the .pdb, which is mostly in the way of X00 (string and number)
def get_atom_type_and_number(muddled_string):
	string_out=""
	number_out=""
	for i in muddled_string:
		try:
			number_out+=str(int(i))
		except ValueError:
			string_out+=i
	return string_out, number_out

#Quaternion functions
##----
def q_mult(q1, q2):
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    w = w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2
    x = w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2
    y = w1 * y2 + y1 * w2 + z1 * x2 - x1 * z2
    z = w1 * z2 + z1 * w2 + x1 * y2 - y1 * x2
    return w, x, y, z
def q_conjugate(q):
    w, x, y, z = q
    return (w, -x, -y, -z)
def qv_mult(q1, v1):
    q2 = (0.0,) + v1
    return q_mult(q_mult(q1, q2), q_conjugate(q1))[1:]

##----
pth="inp_rotater"
mols=os.listdir(pth)

for mol in mols:
	##----STEP 2: Read .pdbs from inp_rotater
	N=int(mol.split("_")[0])
	LIG=mol.split("_")[1][:-5]
	CHAIN=mol.split("_")[1][-5]
	
	pdb=[]
	with open(pth+"\\"+mol) as fh:
		for line in fh:
			pdb.append(line.split())
	
	##----STEP 3: Calculate center of mass
	m,x,y,z=0,0,0,0
	for line in pdb:
		if line[-1] in elements_dict:
			mi=elements_dict[line[-1]]
			x+=float(line[6])*mi
			y+=float(line[7])*mi
			z+=float(line[8])*mi
			m+=mi
	x/=m
	y/=m
	z/=m
	##---- TRANSLATE
	for line in range(len(pdb)):
		if pdb[line][-1] in elements_dict:
			xt=float(pdb[line][6])-x
			yt=float(pdb[line][7])-y
			zt=float(pdb[line][8])-z
			pdb[line][6]=xt
			pdb[line][7]=yt
			pdb[line][8]=zt
	
	##---STEP 4: FOR LOOP
	for k in range(1,N+1):
		##----STEP 5: Get random rotation unit quaternion
		angl=random.uniform(0,2*pi)
		vec=[random.uniform(-1,1) for i in range(3)]
		quat=[cos(angl/2)]+[sin(angl/2)*i for i in vec]
		norm=sum([i**2 for i in quat])**0.5
		quat=tuple([i/norm for i in quat])
		
		for line in range(len(pdb)):
			if pdb[line][-1] in elements_dict:
				xt=pdb[line][6]
				yt=pdb[line][7]
				zt=pdb[line][8]
				##----STEP 7: Rotate molecule
				xr,yr,zr=qv_mult(quat,(xt,yt,zt))
				pdb[line][6]=xr
				pdb[line][7]=yr
				pdb[line][8]=zr
				
				##----STEP 8: INSERT LIG and N
				pdb[line][3]=LIG
				pdb[line][4]=CHAIN
				pdb[line][5]=str(k)
				
		##----STEP 9: Save .pdb to inp_placer
		dat=open("inp_placer/"+str(k)+"_"+LIG+CHAIN+".pdb","w")
		dat.write("CRYST"+str(k)+"\n")
		for line in pdb:
			if line[-1] in elements_dict:
				dat.write("{:6s}".format(line[0])) 		#1 ATOM
				dat.write("{:>5s}".format(line[1])) 	#2 atom serial number
				dat.write("{:1s}".format(" ")) 			#space?
				
				dat.write("{:^4s}".format(line[2]))		#3 atom name
				dat.write("{:1s}".format(" ")) 			#4 alternate location indicator
				dat.write("{:3s}".format(LIG)) 		#5 residue name
				dat.write("{:>2s}".format(CHAIN))	#6 chain identifier
				dat.write("{:4d}".format(int(line[5]))) 		#7 residue sequence number
				
				dat.write("{:3s}".format(" ")) 			#space?
				
				dat.write("{:8.3f}".format(float(line[6])))	#9 orthogonal coordinates for X (in Angstroms)
				dat.write("{:8.3f}".format(float(line[7])))	#10 orthogonal coordinates for Y (in Angstroms) 	
				dat.write("{:8.3f}".format(float(line[8])))	#11 orthogonal coordinates for Z (in Angstroms)
				dat.write("{:6.2f}".format(float(line[9])))		#12	occupancy
				dat.write("{:6.2f}".format(float(line[10])))		#13 temperature factor
				
				dat.write("{:10s}".format(" ")) 			#space?
				
				dat.write("{:>2s}".format(get_atom_type_and_number(line[2])[0]))	#14 element symbol
				
				dat.write("{:2s}".format(" ")) 			#Charge on atom
				
				dat.write("\n")

		dat.write("END")
		dat.close()
		
	
	
	
	

