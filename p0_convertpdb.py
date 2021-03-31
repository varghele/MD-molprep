#----------------------------------
# Program to take Charrm-GUIs shitty excuse of a .pdb file and fix it
# .pdb format taken from https://cupnet.net/pdb-format/
#----------------------------------

#####------------------------------IMPORTS
import os
import sys

# Set the residue name and chain identifier
residue_name=sys.argv[2]
chain_identifier=sys.argv[3]

# Import the false .pdb
pdb_input=sys.argv[1]

# Read in the false pdb:
pdb_false=[]
with open(pdb_input) as fh:
	for line in fh:
		pdb_false.append(line.split())

# Get rid of HEADER, END, and TER
pdb_false=pdb_false[2:-2]

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

# Write correct .pdb
dat=open(residue_name+chain_identifier+".pdb","w")
dat.write(residue_name+chain_identifier+"\n")
for line in pdb_false:
	dat.write("{:6s}".format(line[0])) 		#1 ATOM
	dat.write("{:>5s}".format(line[1])) 	#2 atom serial number
	dat.write("{:1s}".format(" ")) 			#space?
	
	dat.write("{:^4s}".format(get_atom_type_and_number(line[2])[0]+get_atom_type_and_number(line[2])[1]))		#3 atom name
	dat.write("{:1s}".format(" ")) 			#4 alternate location indicator
	dat.write("{:3s}".format(residue_name)) 		#5 residue name
	dat.write("{:>2s}".format(chain_identifier))	#6 chain identifier
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

for line in pdb_false:
	print(get_atom_type_and_number(line[2]))







