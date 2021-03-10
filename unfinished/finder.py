from random import uniform
import matplotlib.pyplot as plt
import numpy as np

class Particle(object):
	box_size=None
	cutoff=5
	dt=0.01
	lamda=0.3 #damping term
	
	def __init__(self,x,y,vx,vy):	#initialize paricle position and velocity, don't need acceleration
		self.x=x
		self.y=y
		self.vx=vx
		self.vy=vy
		self.ax=0
		self.ay=0
	
	def get_force_parts(self,other):
		dx=other.x-self.x	#x-distance between particles
		dy=other.y-self.y	#y-distance between particles
		r2=(dx**2+dy**2)	#euclidean r**2 between particles
		r=r2**0.5			#distance r between particles
		
		if dx==dy==0:
			return	#same particle(or particles at same position)
		
		#self.ax=dx/r2+dy/r2+bxl/bxyd2-bxr/byxd2
		self.ax+=1/dx
		self.ay+=1/dy
	
	def get_force_walls(self):
		bxl=self.x						#distance to left edge of box
		bxr=self.box_size[0]-self.x		#distance to right edge of box
		byb=self.y						#distance to bottom edge of box
		byt=self.box_size[1]-self.y		#distance to top edge of box
		bxyd2=(bxl**2+bxr**2+byb**2+byt**2)
		
		self.ax+=0.5/(bxl)-0.5/(bxr)
		self.ay+=0.5/(byb)-0.5/(byt)
	
	def update(self):
		self.vx+=self.ax*self.dt
		self.vy+=self.ay*self.dt
		self.vx*=self.lamda
		self.vy*=self.lamda
		dv=0.0001
		if self.vx>dv:
			self.vx=dv
		if self.vy>dv:
			self.vy=dv
		if self.vx<-dv:
			self.vx=-dv
		if self.vy<-dv:
			self.vx=-dv
		self.x+=self.vx*self.dt
		self.y+=self.vy*self.dt



#Set box size
Particle.box_size=[70,70]

#Initialize particles
parts=[]
for i in range(4):
	ux=uniform(0,70)
	uy=uniform(0,70)
	p=Particle(ux,uy,0,0)
	parts.append(p)

#Run simulation
T=100000
delta=0.01
t=0

xes=[]
yes=[]

while t<T:
	for p in parts:
		#Set acceleration to 0
		p.ax,p.ay=0,0
		#get particle forces
		for i in parts:
			p.get_force_parts(i)
		#get wall forces
		p.get_force_walls()
	sumv=[]
	xe=[]
	ye=[]
	for p in parts:
		#update
		p.update()
		xe.append(p.x)
		ye.append(p.y)
		#
		sumv.append((p.vx**2+p.vy**2)**0.5)
	if t%10<0.01:
		print(sum(sumv))
		xes.append(xe)
		yes.append(ye)
	t+=delta

plt.scatter(xes[-1],yes[-1])
plt.savefig("end.png")
plt.clf()

xes=np.transpose(xes)
yes=np.transpose(yes)

for p in range(len(parts)):
	plt.scatter(xes[p],yes[p])

plt.savefig("trajectories.png")
