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
path="/home/ebeyerle/Desktop/PCA/UBQN/400ns/"
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
    cbar=plt.colorbar()
    cbar.set_label("Free Energy (kcal/mol)")
    #fig.show()
    plt.savefig("fes_"+str(mode)+".pdf")
