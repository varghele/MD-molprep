import itertools
import matplotlib.pyplot as plt

angstrom_to_au=1.8897259886

dim_xy=74*angstrom_to_au
dims_z=[10*angstrom_to_au,65*angstrom_to_au]
n1=5
n2=4

x=[(dim_xy/n1)*i-(dim_xy/n1)/2 for i in range(1,n1+1)]
print(x)
y=[(dim_xy/n2)*i-(dim_xy/n2)/2 for i in range(1,n2+1)]
print(y)
combs=list(itertools.product(x, y))
#print(combs)

dat=open("pos_lower.txt","w")
for c in combs:
	dat.write(str(c[0]))
	dat.write(",")
	dat.write(str(c[1]))
	dat.write(",")
	dat.write(str(dims_z[0]))
	dat.write("\n")
dat.close()

dat=open("pos_upper.txt","w")
for c in combs:
	dat.write(str(c[0]))
	dat.write(",")
	dat.write(str(c[1]))
	dat.write(",")
	dat.write(str(dims_z[1]))
	dat.write("\n")
dat.close()

for c in combs:
	plt.scatter(c[0],c[1])
plt.savefig("test.png")