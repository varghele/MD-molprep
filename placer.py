import os
import random

pth_pdb="inp_placer"
pth_var="var_placer"

upper=pth_pdb+"\\upper\\"
lower=pth_pdb+"\\lower\\"

var_upper=pth_var+"\\pos_upper.txt"
var_lower=pth_var+"\\pos_lower.txt"

Notorious_PDB=[]

##----UPPER LEAFLET
#--Read in positions
pos_uppers=[]
with open(var_upper) as fh:
	for line in fh:
		pos_uppers.append(line.split(","))

#--Shuffling the positions (not for exact placement)
random.shuffle(pos_uppers)

#--Read in .pdbs
i=0
for mol in os.listdir(upper):
	pdb=[]
	with open(upper+"\\"+mol) as fh:
		for line in fh:
			pdb.append(line.split())
	
	pdb=pdb[1:-1]	#Get rid of header and end token
	
	#--Transit
	for line in range(len(pdb)):
		xt=float(pdb[line][5])+float(pos_uppers[i][0])
		yt=float(pdb[line][6])+float(pos_uppers[i][1])
		zt=float(pdb[line][7])+float(pos_uppers[i][2])
		pdb[line][5]=xt
		pdb[line][6]=yt
		pdb[line][7]=zt
		
		Notorious_PDB.append(pdb[line])	#Append all to big PDB
	
	i+=1
	


##----LOWER LEAFLET
#--Read in positions
pos_lowers=[]
with open(var_lower) as fh:
	for line in fh:
		pos_lowers.append(line.split(","))

#--Shuffling the positions (not for exact placement)
random.shuffle(pos_lowers)

#--Read in .pdbs
i=0
for mol in os.listdir(lower):
	pdb=[]
	with open(lower+"\\"+mol) as fh:
		for line in fh:
			pdb.append(line.split())
	
	pdb=pdb[1:-1]	#Get rid of header and end token
	
	#--Transit
	for line in range(len(pdb)):
		xt=float(pdb[line][5])+float(pos_lowers[i][0])
		yt=float(pdb[line][6])+float(pos_lowers[i][1])
		zt=float(pdb[line][7])+float(pos_lowers[i][2])
		pdb[line][5]=xt
		pdb[line][6]=yt
		pdb[line][7]=zt
		
		Notorious_PDB.append(pdb[line])	#Append all to big PDB
	i+=1

for line in Notorious_PDB:
	print(line)

##----WRITEOUT
dat=open(pth_pdb+"\\BigPDB.pdb","w")
dat.write("FakePeptide\n")
for line in Notorious_PDB:
	dat.write("{:6s}".format(line[0])) 		#1 ATOM
	dat.write("{:>5s}".format(line[1])) 	#2 atom serial number
	dat.write("{:2s}".format(" "))
	dat.write("{:4s}".format(line[2]))		#3 atom name
	dat.write("{:1s}".format(" ")) 			#4 alternate location indicator
	dat.write("{:5s}".format(line[3])) 		#5 residue name
	dat.write("{:>2s}".format(line[4])) 		#7 residue sequence number
	dat.write("{:4s}".format(" ")) 			#8 code for insertion of residues
	dat.write("{:8.3f}".format(line[5]))	#9 orthogonal coordinates for X (in Angstroms)
	dat.write("{:8.3f}".format(line[6]))	#10 orthogonal coordinates for Y (in Angstroms) 	
	dat.write("{:8.3f}".format(line[7]))	#11 orthogonal coordinates for Z (in Angstroms)
	dat.write("{:>6s}".format(line[8]))		#12	occupancy
	dat.write("{:>6s}".format(line[9]))		#13 temperature factor
	dat.write("{:>7s}".format("L"))			#7 block?
	dat.write("{:3s}".format(" "))
	dat.write("{:>2s}".format(line[11]))	#14 element symbol
	dat.write("{:2s}".format(" "))			#15 charge on the atom
	dat.write("\n")

dat.write("END")
dat.close()



