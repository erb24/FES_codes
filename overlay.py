#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 19:24:10 2018

@author: ebeyerle
"""

import Markov_Models as mm
import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import os
path="/home/ebeyerle/Desktop/PCA/"
os.chdir(path)
mode=1
for mode in range(1,11):
    fh=open("mode","w")
    fh.write(str(mode)+"\n")
    fh.close()
    
    with open(path+'anly_'+str(mode)+'.dat','r') as data:
        x=[]
        y=[]
        for line in data:
            p=line.split()
            x.append(float(p[0]))
            y.append(float(p[1]))
    data=[np.column_stack([x,y])]
    model = mm.MSM(data)
    his0 = model.histogram(0, bins=50)
    his1 = model.histogram(1, bins=50)
    his = model.histogram(bins=50)
    ext=model.extent
    #Find the global minimum
    max0=max(his0)
    max1=max(his1)
    min0=[i for i, j in enumerate(his0) if j ==max0][0]
    min1=[i for i, j in enumerate(his1) if j ==max1][0]
    bins=50
    minx=(min0/bins)*180 #Theta min in degrees
    miny=(min1/bins)*360 #Phi min in degrees
    x_axis=np.linspace(ext[0],ext[1],bins)
    y_axis=np.linspace(ext[2],ext[3],bins)
    fig=plt.figure(figsize=(10,8))
    plt.contourf(-np.ma.log(his.T/his.sum()), 25, cmap='gnuplot', extent=model.extent)
    plt.xlabel(r'$\theta$')
    plt.ylabel(r'$\phi$')
    plt.colorbar()
    #fig.show()
    with open(path+"interpol"+str(mode)+"_traj.dat") as points:
        x=[]
        y=[]
        for line in points:
            p=line.split()
            x.append(float(p[0]))
            y.append(float(p[1]))
    #coords=plt.ginput(-1,show_clicks=True)
    for i in range(len(x)):
    	plt.scatter([x[i]],[y[i]],c='w',marker='*')
    	#plt.scatter([minx],[miny],c='w',marker='*')
    plt.scatter([x[0]],[y[0]],c='b',marker='x') #Min proj
    plt.scatter([x[1]],[y[1]],c='r',marker='x') #Max proj
    plt.savefig("fes_points_"+str(mode)+".pdf")
