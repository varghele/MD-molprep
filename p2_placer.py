import os
import random
import numpy as np
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

#test

pth_pdb="inp_placer"
pth_var="var_placer"

upper=pth_pdb+"\\upper\\"
lower=pth_pdb+"\\lower\\"

pars=pth_var+"\\params.txt"
var_upper=pth_var+"\\pos_upper.txt"
var_lower=pth_var+"\\pos_lower.txt"

##---Read in params
N_u=len(os.listdir(upper))
N_l=len(os.listdir(lower))
crit=0.00001
delta=3 #spatial variance in t-direction of COMs
mindist=1.2  #minimum distance in angstrom between two atoms of different residues

params=[]
with open(pars) as fh:
	for line in fh:
		params.append(line.split()[1])

##---Check if you have to run thomson, if yes, do so
#--Upper leaflet
if not os.path.exists(pth_var+"\\"+"coords"+"\\N"+str(N_u)+"_"+params[0]+"_"+params[1]):
	osstr = "thomson.py" + " " + str(N_u) + " " + params[0] + " " + params[1] + " " + str(crit)
	cmdstr = "cmd /c " + '"' + osstr + '"'
	os.chdir(os.path.dirname(os.path.realpath(__file__)))
	os.system(cmdstr)
#--Lower Leaflet
if not os.path.exists(pth_var+"\\"+"coords"+"\\N"+str(N_l)+"_"+params[0]+"_"+params[1]):
	osstr = "thomson.py" + " " + str(N_l) + " " + params[0] + " " + params[1] + " " + str(crit)
	cmdstr = "cmd /c " + '"' + osstr + '"'
	os.chdir(os.path.dirname(os.path.realpath(__file__)))
	os.system(cmdstr)

##----MAKE COORDINATES
angstrom_to_au=1.8897259886
angstrom_to_au=1 #MD works in angstrom
dims_z=[float(i)*angstrom_to_au for i in params[-2:]]

pos_upper=[]
with open(pth_var+"\\"+"coords"+"\\N"+str(N_u)+"_"+params[0]+"_"+params[1]) as fh:
	for line in fh:
		pos_upper.append([float(i)*angstrom_to_au for i in line.split("\t")])
pos_lower=[]
with open(pth_var+"\\"+"coords"+"\\N"+str(N_l)+"_"+params[0]+"_"+params[1]) as fh:
	for line in fh:
		pos_lower.append([float(i)*angstrom_to_au for i in line.split("\t")])

dat=open(pth_var+"\\pos_upper.txt","w")
for c in pos_upper:
	dat.write(str(c[0]))
	dat.write(",")
	dat.write(str(c[1]))
	dat.write(",")
	dat.write(str(random.uniform(dims_z[0]-delta,dims_z[0]+delta)))
	dat.write("\n")
dat.close()

dat=open(pth_var+"\\pos_lower.txt","w")
for c in pos_lower:
	dat.write(str(c[0]))
	dat.write(",")
	dat.write(str(c[1]))
	dat.write(",")
	dat.write(str(random.uniform(dims_z[1]-delta,dims_z[1]+delta)))
	dat.write("\n")
dat.close()

Notorious_PDB=[]
all_mol_xyz=[]

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
for mol in sorted(os.listdir(upper),key=lambda x: (x[-8:-4],len(x))):
	print("U"+str(i))
	pdb=[]
	with open(upper+"\\"+mol) as fh:
		for line in fh:
			pdb.append(line.split())
	
	pdb=pdb[1:-1]	#Get rid of header and end token
	LEN=len(pdb)

	temp_mol=[[float(pdb[line][6]),float(pdb[line][7]),float(pdb[line][8])] for line in range(len(pdb))]

	rep=False
	while not rep:
		#--Rotate
		angl = random.uniform(0, 2 * pi)
		vec = [random.uniform(-1, 1) for i in range(3)]
		quat = [cos(angl / 2)] + [sin(angl / 2) * i for i in vec]
		norm = sum([i ** 2 for i in quat]) ** 0.5
		quat = tuple([i / norm for i in quat])

		for line in range(len(pdb)):
			xi,yi,zi=float(temp_mol[line][0]),float(temp_mol[line][1]),float(temp_mol[line][2])
			xr,yr,zr=qv_mult(quat,(xi,yi,zi))

			pdb[line][6] = xr
			pdb[line][7] = yr
			pdb[line][8] = zr


		#--Transit and number residues
		for line in range(len(pdb)):
			#pdb[line][5]=str(i+1)

			xt=float(pdb[line][6])+float(pos_uppers[i][0])
			yt=float(pdb[line][7])+float(pos_uppers[i][1])
			zt=float(pdb[line][8])+float(pos_uppers[i][2])
			pdb[line][6]=xt
			pdb[line][7]=yt
			pdb[line][8]=zt
		
		temp_mol_rot=[[float(pdb[line][6]),float(pdb[line][7]),float(pdb[line][8])] for line in range(len(pdb))]

		#Check rep
		rep=True
		for elm0 in all_mol_xyz:
			for elm1 in temp_mol_rot:
				if np.linalg.norm(np.array(elm0)-np.array(elm1))<mindist:
					rep=False

	for line in range(len(pdb)):
		xi,yi,zi=float(pdb[line][6]),float(pdb[line][7]),float(pdb[line][8])
		all_mol_xyz.append([xi,yi,zi])

	for line in range(len(pdb)):
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
j=0
for mol in sorted(os.listdir(lower),key=lambda x: (x[-8:-4],len(x))):
	print("L"+str(j))
	pdb=[]
	with open(lower+"\\"+mol) as fh:
		for line in fh:
			pdb.append(line.split())
	
	pdb=pdb[1:-1]	#Get rid of header and end token

	temp_mol = [[float(pdb[line][6]), float(pdb[line][7]), float(pdb[line][8])] for line in range(len(pdb))]

	rep = False
	while not rep:
		# --Rotate
		angl = random.uniform(0, 2 * pi)
		vec = [random.uniform(-1, 1) for i in range(3)]
		quat = [cos(angl / 2)] + [sin(angl / 2) * i for i in vec]
		norm = sum([i ** 2 for i in quat]) ** 0.5
		quat = tuple([i / norm for i in quat])

		for line in range(len(pdb)):
			xi, yi, zi = float(temp_mol[line][0]), float(temp_mol[line][1]), float(temp_mol[line][2])
			xr, yr, zr = qv_mult(quat, (xi, yi, zi))

			pdb[line][6] = xr
			pdb[line][7] = yr
			pdb[line][8] = zr

		# --Transit and number residues
		for line in range(len(pdb)):
			#pdb[line][5] = str(i+j + 1)

			xt = float(pdb[line][6]) + float(pos_lowers[j][0])
			yt = float(pdb[line][7]) + float(pos_lowers[j][1])
			zt = float(pdb[line][8]) + float(pos_lowers[j][2])
			pdb[line][6] = xt
			pdb[line][7] = yt
			pdb[line][8] = zt
		
		temp_mol_rot=[[float(pdb[line][6]),float(pdb[line][7]),float(pdb[line][8])] for line in range(len(pdb))]

		# Check rep
		rep = True
		for elm0 in all_mol_xyz:
			for elm1 in temp_mol_rot:
				if np.linalg.norm(np.array(elm0) - np.array(elm1)) < mindist:
					rep = False

	for line in range(len(pdb)):
		xi, yi, zi = float(pdb[line][6]), float(pdb[line][7]), float(pdb[line][8])
		all_mol_xyz.append([xi, yi, zi])

	for line in range(len(pdb)):
		Notorious_PDB.append(pdb[line])  # Append all to big PDB
	j+=1

for line in Notorious_PDB:
	print(line)

##----WRITEOUT
dat=open(pth_pdb+"\\BigPDB.pdb","w")
dat.write("FakePeptide\n")
ter=1
sn=1
for line in Notorious_PDB:
	
	#print(line[5])
	#break
	
	if ter!=int(line[5]):
		dat.write("TER\n")
		ter=int(line[5])
	
	dat.write("{:6s}".format(line[0])) 		#1 ATOM
	dat.write("{:>5s}".format(str(sn))) 	#2 atom serial number
	dat.write("{:1s}".format(" ")) 			#space?
	
	dat.write("{:^4s}".format(line[2]))		#3 atom name
	dat.write("{:1s}".format(" ")) 			#4 alternate location indicator
	dat.write("{:3s}".format(line[3])) 		#5 residue name
	dat.write("{:>2s}".format(line[4]))	#6 chain identifier
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
	
	
	sn+=1

dat.write("END")
dat.close()



