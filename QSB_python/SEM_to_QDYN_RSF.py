#P. Galvez, AECOM, Sep. 2016
#Input (SEM to QDYN) .
import numpy as np
import math
import sys

sys.path.append('/Users/percy.galvez/Dropbox/qdyn_developer/QSB_python')

import FSEM3D_snapshot as fsnap
import Qdyn_read_in as qdyn
import Qdyn_read_ox_seq as qdyn_ox

nprocs = 6
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
d=fsnap.SNAPSHOT(dF)

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

## Reading mesh chunk
print 'Loading qdynXXXXXX.in ...'

for iproc in range(nprocs):
  qdynin='qdyn%06d.in'%(iproc) 
  model='./'
  pdyn=qdyn.dqyn()
  p=qdyn.Qdyn_read_in(pdyn,model,qdynin)
  print p.X

  print 'Done'

#  rname = 'fort.'+str(20001)
#  print 'Loading QDYN output : ',rname,'...'
#    
#  osnap = qdyn_ox.dsnap()
#  o = qdyn_ox.Qdyn_read_ox_seq(osnap,rname)

  p.X = p.X+x_off
  p.Y = p.Y+y_off
  p.Z = p.Z+z_off

  print 'Matching Values ...'

  p.N = p.NX*p.NW

  for i in range(p.N):
      idex=np.argmin(np.square(d.X-p.X[i]) + np.square(d.Z-p.Z[i]))
      p.TH_0[i] = d.S[idex]
      p.V_0[i]  = d.Vx[idex]
      p.SIGMA[i] = -d.Tz[idex]
      if (i%math.ceil(p.N/100)==0):
          print(i/math.ceil(p.N/100), '% Complete')

  p.X = p.X - x_off
  p.Y = p.Y - y_off
  p.Z = p.Z - z_off

  print 'Done' 

  print 'Writting qdyn.in ...'
  print p.MESHDIM

  qdynout='qdyn%06d.in'%(iproc)
  print qdynout 
  fid = open(qdynout,'w')
  txt='%u		meshdim\n'%(float(p.MESHDIM))
  fid.write(txt)

  if p.SIGMA_CPL==1:
     p.NEQS=3

  if ((float(p.MESHDIM)==2) or (float(p.MESHDIM)==4)):
     txt='%u %u	NX, NW\n'%(p.NX,p.NW)
     fid.write(txt)
     txt='%.15g %.15g %.15g	L,W,Z_CORNER\n'%(p.L,p.W,p.Z_CORNER)
     fid.write(txt)
     for i in range(len(p.DW)):
         txt='%.15g %.15g\n'%(p.DW[i],p.DIP_W[i])
         fid.write(txt)
  else:
     txt='%u	NN\n'%(p.N)
     fid.write(txt)
     txt='%.15g %.15g	L,W\n'%(p.L,p.W)
     fid.write(txt)

  if (p.MESHDIM==1):
     txt='%u	finite\n'%(p.FINITE)
     fid.write(txt)

  txt='%u   itheta_law\n'%(int(p.THETA_LAW))
  fid.write(txt)
  txt='%u   i_rns_law\n'%(int(p.RNS_LAW))
  fid.write(txt)
  txt='%u   i_sigma_cpl\n'%(int(p.SIGMA_CPL))
  fid.write(txt)
  txt='%u   n_equations\n'%(int(p.NEQS))
  fid.write(txt)
  txt='%u %u %u %u %u %u  ntout, nt_coord, nxout, nxout_DYN, ox_SEQ, ox_DYN\n'%(int(p.NTOUT),int(p.IC),int(p.NXOUT),int(p.NXOUT_DYN),int(p.OX_SEQ),int(p.OX_DYN))
  fid.write(txt)    
  txt='%.15g %.15g %.15g %.15g    beta, smu, lambda, v_th\n'%(float(p.VS),float(p.MU),float(p.LAM),float(p.V_TH))
  fid.write(txt)
  txt='%.15g %.15g    Tper, Aper\n'%(float(p.TPER),float(p.APER))
  fid.write(txt)
  txt='%.15g %.15g %.15g %.15g    dt_try, dtmax, tmax, accuracy\n'%(float(p.DTTRY),float(p.DTMAX),float(p.TMAX),float(p.ACC))
  fid.write(txt)
  txt='%u   nstop\n'%(int(p.NSTOP))
  fid.write(txt)
  txt='%u %u  DYN_FLAG, DYN_SKIP\n'%(int(p.DYN_FLAG),int(p.DYN_SKIP))
  fid.write(txt)
  txt='%.15g %.15g %.15g    M0, DYN_th_on, DYN_th_off\n'%(float(p.DYN_M),float(p.DYN_TH_ON),float(p.DYN_TH_OFF))
  fid.write(txt)
  for l in range(p.N):
      txt='%.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g %.15g\n'%(p.SIGMA[l],p.V_0[l],\
            p.TH_0[l],p.A[l],p.B[l],p.DC[l],p.V1[l],p.V2[l],p.MU_SS[l],p.V_SS[l],p.IOT[l],p.IASP[l],p.CO[l])
      fid.write(txt)
  for k in range(p.N):
      txt='%.15g %.15g %.15g %15.g\n'%(p.X[k],p.Y[k],p.Z[k],p.DIP[k])
      fid.write(txt)
  fid.close()
