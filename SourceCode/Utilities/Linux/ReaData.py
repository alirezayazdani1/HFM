# -*- coding: utf-8 -*-
import glob , os , shutil
from math import cos,pi,exp
from numpy import *
import scitools.filetable as ft
#from scitools.easyviz import *
import sys
file = open('result.out' , 'r')
file.readline()
X,Y,Z,U,V,W,P,T =ft.read_columns(file); 
for i in xrange(X.shape[0]):
    print X[i],T[i]