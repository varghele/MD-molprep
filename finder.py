class Particle(object):
	box_size=None
	cutoff=5
	dt=0.01
	
	def __init__(self,x,y,vx,vy):	#initialize paricle position and velocity, don't need acceleration
		self.x=x
		self.y=y
		self.vx=vx
		self.vy=vy
		self.ax=0
		self.ay=0
	
	def get_force(self,other):
		dx=other.x-self.x	#x-distance between particles
		dy=other.y-self.y	#y-distance between particles
		r2=(dx**2+dy**2)	#euclidean r**2 between particles
		r=r2**0.5			#distance r between particles
		
		bxl=self.x					#distance to left edge of box
		bxr=box_size[0]-self.x		#distance to right edge of box
		byb=self.y					#distance to bottom edge of box
		byt=box_size[1]-self.y		#distance to top edge of box
		bxyd2=(bxl**2+bxr**2+byb**2+byt**2)
		
		#self.ax=dx/r2+dy/r2+bxl/bxyd2-bxr/byxd2
		self.ax=dx+bxl-bxr
		self.ay=dy+db-dyt
		
	def update(self):
		self.vx+=self.ax*dt
		self.vy+=self.ay*dt
		self.x+=self.vx+dt+(self.ax/2)*dt**2
		self.y+=