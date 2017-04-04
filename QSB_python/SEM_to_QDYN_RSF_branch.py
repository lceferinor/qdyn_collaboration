#P. Galvez, AECOM, Sep. 2016
#Input (SEM to QDYN) .
import numpy as np
import math
import sys

sys.path.append('/Users/percy.galvez/Dropbox/qdyn_developer/QSB_python')

import FSEM3D_snapshot as fsnap
import Qdyn_read_branch_in as qdyn_branch
import Qdyn_read_ox_seq as qdyn_ox

x_off=-64000.0
y_off=0
z_off=0

isnap=60000
dt=0.001
ft_id=1
i_v_corr = 0 # = 1 to set Vmin = v_corr
v_corr = 1e-20 


print '[Processing Snapshot #',isnap,'...]'
dF=fsnap.dFault()
dF.frame=isnap
ifault=2 # for branch
d=fsnap.SNAPSHOT(dF,ifault)

#if i_v_corr==1:
#   d.V1 = max(d.V1,v_corr)

iiv = np.where(d.v1>=0)
d.Vx = d.v1[iiv]
d.Vz = d.v2[iiv]
d.Dx = d.d1[iiv]
d.Dz = d.d2[iiv]
d.Tx = d.t1[iiv]
d.Ty = d.t2[iiv]
d.Tz = d.t3[iiv]
d.X = d.xcoord[iiv]
d.Y = d.ycoord[iiv]
d.Z = d.zcoord[iiv]
d.S = d.sta[iiv]
d.Sg = d.stg[iiv]
d.Trup = d.tRUP[iiv]
d.Tpz = d.tPZ[iiv]

print 'Done'

print 'Loading qdyn_branch.in ...'
N_branch,DIP0,lam,mu,X_b,Y_b,Z_b,SIGMA_b,V0_b,THSS_b,A_b,B_b,DC_b,V1_b,V2_b,MUSS_b,VSS_b,CO_b,TAU_b=qdyn_branch.read_branch_in()

print 'Done'

xoff = -64000.00
yoff = 0.0
zoff = 0.0
X_b = X_b+x_off
Y_b = Y_b+y_off
Z_b = Z_b+z_off

print 'Matching Values ...'

for i in range(N_branch):
    idex=np.argmin(np.square(d.X-X_b[i]) + np.square(d.Z-Z_b[i]))
    THSS_b[i]  = d.S[idex]
    V0_b[i]    = d.Vx[idex]*0.01
    SIGMA_b[i] = -d.Tz[idex]
    if (i%math.ceil(N_branch/100)==0):
        print(i/math.ceil(N_branch/100), '% Complete')

print 'Done'

X_b = X_b - x_off
Y_b = Y_b - y_off
Z_b = Z_b - z_off 

print 'Writting qdyn_branch.in ...'
fid = open('qdyn_branch.in','w')
txt='%u\n'%(N_branch)
fid.write(txt)
txt='%15.6f\n'%(DIP0)
fid.write(txt)
txt='%20.6f %20.6f\n'%(lam,mu)
fid.write(txt)

for l in range(N_branch):
      txt='%.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g\n'%(X_b[l],Y_b[l],Z_b[l],SIGMA_b[l],V0_b[l],\
            THSS_b[l],A_b[l],B_b[l],DC_b[l],V1_b[l],V2_b[l],MUSS_b[l],VSS_b[l],CO_b[l])
      fid.write(txt)
fid.close()
