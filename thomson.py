##----Takes input dimensions and computes stable equilibrium positions of N points within a box that repels particles at half the distance they repel each other

import sys
import os
import numpy as np
import tensorflow as tf
import matplotlib.pyplot as plt

pth="var_placer\\coords\\"
if not os.path.exists(pth):
    os.mkdir(pth)

#Get rid of specific TF error, that sometimes doesn't allow matmul
physical_devices = tf.config.list_physical_devices('GPU') 
tf.config.experimental.set_memory_growth(physical_devices[0], True)

N=int(sys.argv[1])
Bx=float(sys.argv[2])
By=float(sys.argv[3])
crit=float(sys.argv[4])

box=np.array([Bx,By],dtype="float64")

r0=np.random.uniform(box,size=[N,2])            #Setting initial particle positions
c=tf.Variable(r0,name="coord",trainable=True)   #trainable coordinates

def squared_diff(A):            #Taken from https://github.com/thoppe/tf_thomson_charges
    r = tf.reduce_sum(A*A, 1)
    r = tf.reshape(r, [-1, 1])
    return r - 2*tf.matmul(A, tf.transpose(A)) + tf.transpose(r)

def modloss():      #Loss Function for the TF optimizer
    ##----Get pairwise distances of all particles
    RR = squared_diff(c)
    mask = np.triu(np.ones((N, N), dtype=np.bool_), 1)
    R = tf.sqrt(tf.boolean_mask(RR, mask))
    ##----Get distance of left and right boundaries
    p = tf.math.square(c)
    m = tf.math.square((tf.repeat([box], N, axis=0) - c))
    ##----Enforce boundaries with OOB Loss
    mi, ni = 2, 100  # y=mx+n
    poob = tf.math.abs(tf.where(c > 0, tf.zeros_like(c), c)) * mi
    moob = tf.math.abs(tf.where(c - box < 0, tf.zeros_like(c), c - box)) * mi
    # -
    poob += tf.where(poob != 0, tf.ones_like(poob) * ni, poob)
    moob += tf.where(moob != 0, tf.ones_like(moob) * ni, moob)

    U = tf.reduce_sum(1 / R) + tf.reduce_sum(tf.sqrt(tf.reduce_sum(1 / p, axis=1, keepdims=True))) + tf.reduce_sum(
        tf.sqrt(tf.reduce_sum(1 / m, axis=1, keepdims=True))) + tf.reduce_sum(poob) + tf.reduce_sum(moob)
    return U

def calcloss(c):        #Exactly the same loss function, only for plotting...and because I can't figure out how to extract the loss from the optimizer directly
    ##----Get pairwise distances of all particles
    RR = squared_diff(c)
    mask = np.triu(np.ones((N, N), dtype=np.bool_), 1)
    R = tf.sqrt(tf.boolean_mask(RR, mask))
    ##----Get distance of left and right boundaries
    p = tf.math.square(c)
    m = tf.math.square((tf.repeat([box], N, axis=0) - c))
    ##----Enforce boundaries with OOB Loss
    mi, ni = 2, 100  # y=mx+n
    poob = tf.math.abs(tf.where(c > 0, tf.zeros_like(c), c)) * mi
    moob = tf.math.abs(tf.where(c - box < 0, tf.zeros_like(c), c - box)) * mi
    # -
    poob += tf.where(poob != 0, tf.ones_like(poob) * ni, poob)
    moob += tf.where(moob != 0, tf.ones_like(moob) * ni, moob)

    U = tf.reduce_sum(1 / R) + tf.reduce_sum(tf.sqrt(tf.reduce_sum(1 / p, axis=1, keepdims=True))) + tf.reduce_sum(
        tf.sqrt(tf.reduce_sum(1 / m, axis=1, keepdims=True))) + tf.reduce_sum(poob) + tf.reduce_sum(moob)
    return U

##----START
opt=tf.keras.optimizers.Nadam(learning_rate=0.1)    #Yes, fancy Nadam instead of Adam. Seems to perform slightly better
loss=modloss

prev_loss=N**3
cur_loss=[]
cur_loss.append(calcloss(c).numpy())

while abs(prev_loss-cur_loss[-1])>crit:
    s=opt.minimize(loss,var_list=[c]).numpy()
    l=calcloss(c).numpy()

    prev_loss=cur_loss[-1]
    cur_loss.append(l)

    if s%100==0:
        print("Cur.It.: "+str(s)+"\t \t"+"Cur.Loss: "+str(l))


fig,axs=plt.subplots(1,2)
axs[0].plot(range(len(cur_loss)),cur_loss)
axs[1].scatter(c.numpy()[:,0],c.numpy()[:,1])

axs[0].set_xlabel("Epochs")
axs[0].set_ylabel("Loss")
axs[1].set_xlabel("X")
axs[1].set_ylabel("Y")

axs[1].set_xlim([0,box[0]])
axs[1].set_xlim([0,box[1]])

fig.set_size_inches(16,8)
plt.savefig(pth+"N"+str(N)+"_"+str(box[0])+"_"+str(box[1])+".png")

dat=open(pth+"N"+str(N)+"_"+str(box[0])+"_"+str(box[1]),"w")
for co in c.numpy():
    dat.write(str(co[0])+"\t")
    dat.write(str(co[1])+"\n")
dat.close()


