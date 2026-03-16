import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

acolors = plt.get_cmap('tab20',24)
ccolors = plt.get_cmap('tab20',24)

k=-2*np.pi/24.
offset = 6*k + 3*k #-pi/2
xarra_1 = np.array([37*np.cos(k*i+offset) for i in np.arange(0,24)])
yarra_1 = np.array([37*np.sin(k*i+offset) for i in np.arange(0,24)])
labelsa_1 = np.array([i for i in np.arange(0,24)])

fig,ax = plt.subplots(figsize=(10,10))
ax.invert_yaxis()
ax.plot(xarra_1,yarra_1,"x",label="anode, z=-L/2")

for x,y,label in zip(xarra_1,yarra_1,labelsa_1):
    ax.text(x,y,label)

kc=2*np.pi/24.
offsetc = -4*kc + 2*kc -  np.pi/24 #-pi/4
xarrc_1 = np.array([42*np.cos(kc*i+offsetc) for i in np.arange(0,24)])
yarrc_1 = np.array([42*np.sin(kc*i+offsetc) for i in np.arange(0,24)])
labelsc_1 = np.array([i for i in np.arange(0,24)])

ax.plot(xarrc_1,yarrc_1,"o",label="cathode, z=-L/2, where they are picked up")

for x,y,label in zip(xarrc_1,yarrc_1,labelsc_1):
    ax.text(x,y,label)
plt.title("z=-L/2 plane, beam going into the plane along +z, (+x right, +y down)")
plt.grid()
plt.legend()
plt.savefig("plane1.png")
plt.show()

fig,ax = plt.subplots(figsize=(10,10))
ax.invert_yaxis()

offset = offset-3*k
xarra_2 = np.array([37*np.cos(k*i+offset) for i in np.arange(0,24)])
yarra_2 = np.array([37*np.sin(k*i+offset) for i in np.arange(0,24)])
labelsa_2 = np.array([i for i in np.arange(0,24)])

ax.plot(xarra_2,yarra_2,"x",label="anode, z=+L/2, where they are picked up")

for x,y,label in zip(xarra_2,yarra_2,labelsa_2):
    ax.text(x,y,label)

offsetc = offsetc-3*kc
xarrc_2 = np.array([42*np.cos(kc*i+offsetc) for i in np.arange(0,24)])
yarrc_2 = np.array([42*np.sin(kc*i+offsetc) for i in np.arange(0,24)])
labelsc_2 = np.array([i for i in np.arange(0,24)])

ax.plot(xarrc_2,yarrc_2,"o",label="cathode, z=+L/2")

for x,y,label in zip(xarrc_2,yarrc_2,labelsc_2):
    ax.text(x,y,label)
plt.title("z=+L/2 plane, beam going into the plane along +z, (+x right, +y down)")
plt.grid()
plt.legend()
plt.savefig("plane2.png")
plt.show()

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111,projection='3d')
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")

xx_a = np.array([[x1,x2] for x1,x2 in zip(xarra_1,xarra_2)])
yy_a = np.array([[y1,y2] for y1,y2 in zip(yarra_1,yarra_2)])
zz_a = np.array([[-173.5,173.5] for x1,x2 in zip(xarra_1,xarra_2)])
for i,[xx,yy,zz] in enumerate(zip(xx_a,yy_a,zz_a)):
    ax.plot(xx,yy,zz,'-',color=acolors(i/24))

for i,[x,y,label] in enumerate(zip(xarra_1,yarra_1,labelsa_1)):
    ax.text(x,y,-173,"a"+str(label),color=acolors(i/24))
for i,[x,y,label] in enumerate(zip(xarra_2,yarra_2,labelsa_2)):
    ax.text(x,y,+173,"a"+str(label),color=acolors(i/24))

xx_c = np.array([[x1,x2] for x1,x2 in zip(xarrc_1,xarrc_2)])
yy_c = np.array([[y1,y2] for y1,y2 in zip(yarrc_1,yarrc_2)])
zz_c = np.array([[-173.5,173.5] for x1,x2 in zip(xarrc_1,xarrc_2)])
for i,[xx,yy,zz] in enumerate(zip(xx_c,yy_c,zz_c)):
    ax.plot(xx,yy,zz,'--',color=ccolors(((47-i)%24)/24))

for i,[x,y,label] in enumerate(zip(xarrc_1,yarrc_1,labelsc_1)):
    ax.text(x,y,-173,"c"+str(label),color=ccolors(((25-i)%24)/24))
for i,[x,y,label] in enumerate(zip(xarrc_2,yarrc_2,labelsc_2)):
    ax.text(x,y,+173,"c"+str(label),color=ccolors(((47-i)%24)/24))
ax.view_init(elev=-53, azim=-106, roll=18)
plt.tight_layout()
plt.show()
phi_qqq = np.array([[2*np.pi*(-i*16+j+0.5)/(16*4) for i in range(4)] for j in range(16)])
print(phi_qqq)
#'''
for i,[phi1,phi2,phi3,phi4] in enumerate(phi_qqq):
    ax.plot([50*np.cos(phi1),100*np.cos(phi1)],[50*np.sin(phi1),100*np.sin(phi1)],[100,100],'-',color='red')
    ax.text(104*np.cos(phi1),104*np.sin(phi1),100,"0_%d"%(i),color="red")

    ax.plot([50*np.cos(phi2),100*np.cos(phi2)],[50*np.sin(phi2),100*np.sin(phi2)],[100,100],'-',color='green')
    ax.text(104*np.cos(phi2),104*np.sin(phi2),100,"1_%d"%(i),color="green")

    ax.plot([50*np.cos(phi3),100*np.cos(phi3)],[50*np.sin(phi3),100*np.sin(phi3)],[100,100],'-',color='blue')
    ax.text(104*np.cos(phi3),104*np.sin(phi3),100,"2_%d"%(i),color="blue")

    ax.plot([50*np.cos(phi4),100*np.cos(phi4)],[50*np.sin(phi4),100*np.sin(phi4)],[100,100],'-',color='brown')
    ax.text(104*np.cos(phi4),104*np.sin(phi4),100,"3_%d"%(i),color="brown")
#'''
#coords_qqq = np.array([[50*np.cos(phi),100*np.sin(phi),100] for phi in phi_qqq]).T
plt.tight_layout()
plt.show()

