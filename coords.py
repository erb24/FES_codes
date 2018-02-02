# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import Markov_Models as mm
import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
path="/home/ebeyerle/Desktop/PCA/"
mode=10
fh=open("mode","w")
fh.write(str(mode)+"\n")
fh.close()

#files = [[os.path.join(path,str(i+1),'Distance'), os.path.join(path,str(i+1),'Twist') ] for i in range(10)]
#path="/home/ebeyerle/Desktop/ApA_test/Data/"
counter=1
with open(path+'anly_'+str(mode)+'.dat','r') as data:
    x=[]
    #x2=[]
    y=[]
    #y2=[]
    for line in data:
        p=line.split()
        #if (counter%2) == 0:
        #    x2.append(float(p[0]))
        #    y2.append(float(p[1]))
        #else:
        x.append(float(p[0]))
        y.append(float(p[1]))
        counter=counter+1
data=[np.column_stack([x,y])]
#data2=[np.column_stack([x2,y2])]
model = mm.MSM(data)
his0 = model.histogram(0, bins=50)
his1 = model.histogram(1, bins=50)
his = model.histogram(bins=50)
ext=model.extent
bins=50
x_axis=np.linspace(ext[0],ext[1],bins)
y_axis=np.linspace(ext[2],ext[3],bins)
#plt.figure(figsize=(12,6))
#plt.subplot(1,2,1)
#plt.plot(x_axis,his0/his0.sum(), linewidth=3, alpha=0.5)
#plt.xlabel(r'$R (nm)$', fontsize=25)
#plt.ylabel(r'$P$', fontsize=25)

#plt.subplot(1,2,2)
#plt.plot(y_axis,his1/his1.sum(), linewidth=3, alpha=0.5)
#plt.xlabel(r'$\varphi \ (Degrees)$', fontsize=25)
fig=plt.figure(figsize=(10,8))
plt.contourf(-np.ma.log(his.T/his.sum()), 25, cmap='gnuplot', extent=model.extent)
plt.xlabel(r'$X$')
plt.ylabel(r'$Y$')
plt.colorbar()
fig.show()
coords = []

#def onclick(event):
#    global ix, iy
#    ix, iy = event.xdata, event.ydata
#    print('x = %d, y = %d'%(
#        ix, iy))
#
#    global coords
#    coords.append((ix, iy))

#    return coords
#cid = fig.canvas.mpl_connect('button_press_event', onclick)
coords=plt.ginput(-1,show_clicks=True)
#plt.savefig("fes_"+str(mode)+".pdf")
print(coords)
np.savetxt(path+'coords_'+str(mode)+'.dat',coords)
