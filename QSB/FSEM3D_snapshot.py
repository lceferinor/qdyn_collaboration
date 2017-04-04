import sys
#P. Galvez, AECOM, Sep. 2016
#Based on FSEM3D_snapshot.m
#sys.path.append('/Users/pgalvez/python_lib')
import struct
import array  
import numpy as N
import os
from numpy import *

class dFault : 
      def _init_(self,frame,xcoord,ycoord,zcoord,d1,d2,v1,v2,t1,t2,t3):
          self.frame=frame
          self.xcoord = xcoord 
          self.ycoord = ycoord
          self.zcoord = zcoord
          self.d1 = d1
          self.d2 = d2
          self.v1 = v1
          self.v2 = v2
          self.t1 = t1
          self.t2 = t2
          self.t3 = t3


def SNAPSHOT(dA):
# Read a binary file Snapshotxxxx_Fx.bin
# Number of variables. 
   NVAR = 14
# Open the file first time to get number of nodes (nglob)
   filename='Snapshot'+str(dA.frame)+'_F1.bin'
   fileobj_f=open(filename,"rb")
   datall = N.fromfile(fileobj_f,dtype=N.int32,count=-1)
   nglob = (len(datall)/NVAR)-2
   fileobj_f.close()
# Open the file again to read the array.
   filename='Snapshot'+str(dA.frame)+'_F1.bin'
   fileobj=open(filename,"rb")
# Format:
# write(IOUT,*) dataXZ%xcoord           !4*dummy,4*nglob,4*dummy                         xcoord(nglob)
# write(IOUT,*) dataXZ%ycoord           !4*dummy,4*nglob,4*dummy 
# write(IOUT,*) dataXZ%zcoord           !4*dummy,4*nglob,4*dummy                 
# write(IOUT,*) dataXZ%d1               !4*dummy,4*nglob,4*dummy        
# write(IOUT,*) dataXZ%d2               !4*dummy,4*nglob,4*dummy        
# write(IOUT,*) dataXZ%v1               !4*dummy,4*nglob,4*dummy                         
# write(IOUT,*) dataXZ%v2               !4*dummy,4*nglob,4*dummy                         
# write(IOUT,*) dataXZ%t1               !4*dummy,4*nglob,4*dummy                         
# write(IOUT,*) dataXZ%t2
# write(IOUT,*) dataXZ%t3
# write(IOUT,*) dataXZ%sta
# write(IOUT,*) dataXZ%stg
# write(IOUT,*) dataXZ%tRUP
# write(IOUT,*) dataXZ%tPZ
   xcoord  = N.fromfile(fileobj,dtype=N.float32,count=1+nglob+1)
   ycoord  = N.fromfile(fileobj,dtype=N.float32,count=1+nglob+1)
   zcoord  = N.fromfile(fileobj,dtype=N.float32,count=1+nglob+1)
   d1      = N.fromfile(fileobj,dtype=N.float32,count=1+nglob+1)
   d2      = N.fromfile(fileobj,dtype=N.float32,count=1+nglob+1)
   v1      = N.fromfile(fileobj,dtype=N.float32,count=1+nglob+1)
   v2      = N.fromfile(fileobj,dtype=N.float32,count=1+nglob+1)
   t1      = N.fromfile(fileobj,dtype=N.float32,count=1+nglob+1)
   t2      = N.fromfile(fileobj,dtype=N.float32,count=1+nglob+1)
   t3      = N.fromfile(fileobj,dtype=N.float32,count=1+nglob+1)
   sta     = N.fromfile(fileobj,dtype=N.float32,count=1+nglob+1)
   stg     = N.fromfile(fileobj,dtype=N.float32,count=1+nglob+1)
   tRUP    = N.fromfile(fileobj,dtype=N.float32,count=1+nglob+1)
   tPZ     = N.fromfile(fileobj,dtype=N.float32,count=1+nglob+1)
   dA.xcoord = xcoord[1:nglob+1]
   dA.ycoord = ycoord[1:nglob+1]
   dA.zcoord = zcoord[1:nglob+1]
   dA.d1     = d1[1:nglob+1]
   dA.d2     = d2[1:nglob+1]
   dA.v1     = v1[1:nglob+1]
   dA.v2     = v2[1:nglob+1]
   dA.t1     = t1[1:nglob+1]
   dA.t2     = t2[1:nglob+1]
   dA.t3     = t3[1:nglob+1]
   fileobj.close()
   return dA

def CSV(dA,xsh,ysh):
# Converting Snapshots to .CSV format to make a table in paraview.
    filename='SnapshotF1'+str(dA.frame)+'.csv'
    fileobj = open(filename,'w')
    txt = ''                   # Saving slip and displacement.
    for i in range(len(dA.xcoord)):
        txt = '%f,%f,%f,%f,%f\n' % (dA.xcoord[i]+xsh,dA.ycoord[i]+ysh,dA.zcoord[i],-dA.v2[i],-dA.d2[i])      
        fileobj.write(txt)
    fileobj.close() 

def CSV_init(iframe):
    dF = dFault()
    dF.frame = iframe
    dF = SNAPSHOT(dF)
    xshi = 43832.3207+55000
    yshi = -23381.37+400000 
    CSV(dF,xshi,yshi)


def CSV_all(val):
    for i in range(0,val):
        CSV_init((i+1)*1000)
        print 'Snap:',(i+1)*1000
