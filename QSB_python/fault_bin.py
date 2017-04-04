#P. Galvez, AECOM, Sep. 2016

import struct
import array  
import numpy as N
import os
from numpy import *

def fault_bin(dA):
# Read a binary file Snapshotxxxx_Fx.bin
# Number of variables. 
   NVAR = 8
# Open the file first time to get number of nodes (nglob)
   filename='proc'+str(int(dA.frame)).zfill(6)+'_fault_db.bin'
   print filename
   fileobj=open(filename,"rb")
   nspec_nglob = N.fromfile(fileobj,dtype=N.int32,count=1+2+1)
   print nspec_nglob
   nspec = nspec_nglob[0]
   nglob = nspec_nglob[1]
   nspec_nglob = N.fromfile(fileobj,dtype=N.int32,count=1)
   print nspec_nglob
   NGLLX = 5
   NGLLY = 5
   NDIM  = 3
   NGLLSQUARE   = NGLLX*NGLLY
   nbool        = nspec*NGLLSQUARE
   njacobian2Dw = nbool
   nnormal      = nbool*NDIM
   print nspec,nglob
# Format:
# write(IOUT)   f%nspec,f%nglob    !4*dummy,4*1,4*1,4*dummy
# write(IOUT,*) f%ibool1           !4*dummy,4*NGLLSQUARE*nspec,4*dummy 
# write(IOUT,*) f%jacobian2Dw      !4*dummy,4*NGLLSQUARE*nspec,4*dummy 
# write(IOUT,*) f%normal           !4*dummy,4*nglob,4*dummy                 
# write(IOUT,*) f%ibulk1           !4*dummy,4*nglob,4*dummy        
# write(IOUT,*) f%ibulk2           !4*dummy,4*nglob,4*dummy        
# write(IOUT,*) f%xcoordbulk1      !4*dummy,4*nglob,4*dummy                         
# write(IOUT,*) f%ycoordbulk1      !4*dummy,4*nglob,4*dummy                         
# write(IOUT,*) f%zcoordbulk1      !4*dummy,4*nglob,4*dummy

   if (nspec ==0): return dA                         
   ibool            = N.fromfile(fileobj,dtype=N.int32,count=1+nbool+1)
   jacobian2Dw      = N.fromfile(fileobj,dtype=N.float32,count=1+njacobian2Dw+1)
   normal           = N.fromfile(fileobj,dtype=N.float32,count=1+nnormal+1)
   ibulk1           = N.fromfile(fileobj,dtype=N.int32,count=1+nglob+1)
   ibulk2           = N.fromfile(fileobj,dtype=N.int32,count=1+nglob+1)
   xcoordbulk1      = N.fromfile(fileobj,dtype=N.float32,count=1+nglob+1)
   ycoordbulk1      = N.fromfile(fileobj,dtype=N.float32,count=1+nglob+1)
   zcoordbulk1      = N.fromfile(fileobj,dtype=N.float32,count=1+nglob+1)
   dA.ibool         = ibool[1:nbool+1]
   dA.jacobian2Dw   = jacobian2Dw[1:njacobian2Dw+1]
   dA.normal        = normal[1:nnormal+1]
   dA.ibulk1        = ibulk1[1:nglob+1]
   dA.ibulk2        = ibulk2[1:nglob+1]
   dA.xcoordbulk1   = xcoordbulk1[1:nglob+1]
   dA.ycoordbulk1   = ycoordbulk1[1:nglob+1]
   dA.zcoordbulk1   = zcoordbulk1[1:nglob+1]
   dA.nspec = nspec
   dA.nglob = nglob
   fileobj.close()
   return dA

def fault_txt(dA):
# Read a binary file Snapshotxxxx_Fx.bin
# Number of variables. 
# Open the file first time to get number of nodes (nglob)
   filename='proc'+str(int(dA.frame)).zfill(6)+'_fault_xyz.dat'
   print filename
   fileobj=open(filename,"r")
   rline=fileobj.readline()
   nglob = int(rline)
   dA.nglob = nglob
   if (nglob==0) : 
      fileobj.close()
      return dA
   xcoordbulk1=N.zeros(nglob)
   ycoordbulk1=N.zeros(nglob)
   zcoordbulk1=N.zeros(nglob)
   for i in (range(nglob)):
       rline=fileobj.readline().split()
       xcoordbulk1[i]=float(rline[0])
       ycoordbulk1[i]=float(rline[1])
       zcoordbulk1[i]=float(rline[2])
# Format :
# write(IOUT,*) f%xcoordbulk1, f%ycoordbulk1, f%zcoordbulk1
   dA.xcoordbulk1   = xcoordbulk1
   dA.ycoordbulk1   = ycoordbulk1
   dA.zcoordbulk1   = zcoordbulk1
   fileobj.close()
   return dA



class dFault_bin : 
      def _init_(self,frame,ibool,jacobian2Dw,normal,ibulk1,ibulk2,xcoordbulk1,ycoordbulk1,zcoordbulk1,nspec,nglob):
          self.frame=frame
          self.nspec = nspec
          self.nglob = nglob
          self.ibool = ibool 
          self.jacobian2Dw = jacobian2Dw
          self.normal = normal
          self.ibulk1 = ibulk1
          self.ibulk2 = ibulk2
          self.xcoordbulk1 = xcoordbulk1
          self.ycoordbulk1 = ycoordbulk1
          self.zcoordbulk1 = zcoordbulk1

class dFault_txt : 
      def _init_(self,nglob,xcoordbulk1,ycoordbulk1,zcoordbulk1):
          self.frame=frame
          self.nglob = nglob
          self.xcoordbulk1 = xcoordbulk1
          self.ycoordbulk1 = ycoordbulk1
          self.zcoordbulk1 = zcoordbulk1


def CSV(dA):
# Converting Snapshots to .CSV format to make a table in paraview.
    filename='xyz_fault_'+str(int(dA.frame)).zfill(6)+'.csv'
    print filename
    fileobj = open(filename,'w')
    txt = ''                   # Saving slip and displacement.
    if (dA.nglob)>0:
       for i in range(dA.nglob):
           txt = '%f,%f,%f\n' % (dA.xcoordbulk1[i],dA.ycoordbulk1[i],dA.zcoordbulk1[i])      
           fileobj.write(txt)
    fileobj.close() 

def CSV_frame(iframe):
#   dF = dFault_bin()
    dF = dFault_txt()
    dF.frame = iframe
    dF = fault_txt(dF)
    CSV(dF)
    return dF

nproc=6
dFtot=dFault_txt()
fid=open('nodesonfault.csv','w')
ikglob=0
for i in range(nproc):
    dFtot=CSV_frame(i)
    print 'iframe:',str(i).zfill(6)
    print 'nodes:',dFtot.nglob
    if dFtot.nglob>0:
       for k in range(dFtot.nglob):
           ikglob=ikglob+1
           txt='%6u %6u %15g %15g %15g\n'%(i+1,ikglob,dFtot.xcoordbulk1[k],dFtot.ycoordbulk1[k],dFtot.zcoordbulk1[k])
           fid.write(txt)
fid.close()
        
