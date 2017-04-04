#P. Galvez, AECOM, Sep. 2016
#Reading qdyn.in format. 
import numpy as np

class dqyn :
   def _init_(self,NX,NW,L,W,Z_CORNER,DW,DIP_W,THETA_LAW,RNS_LAW,SIGMA_CPL,NEQS,NTOUT,NXOUT_DYN,OX_SEQ,OX_DYN,VS,MU,LAM,V_TH,TPER,APER,DTTRY,DTMAX,TMAX,ACC,NSTOP,DYN_FLAG,DYN_SKIP,DYN_M,DYN_TH_ON,DYN_TH_OFF,SIGMA,V_0,TH_0,A,B,DC,V1,V2,MU_SS,V_SS,IOT,IASP,CO,X,Y,Z,DIP):
       self.NX = NX
       self.NW = NW
       self.L  = L
       self.W  = W
       self.Z_CORNER = Z_CORNER
       self.DW = DW
       self.DIP_W = DIP_W
       self.THETA_LAW = THETA_LAW
       self.RNS_LAW = RNS_LAW
       self.SIGMA_CPL = SIGMA_CPL
       self.NEQS = NEQS
       self.NTOUT = NTOUT
       self.IC = IC
       self.NXOUT = NXOUT
       self.NXOUT_DYN = NXOUT_DYN
       self.OX_SEQ = OX_SEQ
       self.OX_DYN = OX_DYN
       self.VS = VS
       self.MU = MU
       self.LAM = LAM
       self.V_TH = V_TH
       self.TPER = TPER 
       self.APER = APER
       self.DTTRY = DTTRY
       self.DTMAX = DTMAX
       self.TMAX = TMAX
       self.ACC = ACC
       self.NSTOP = NSTOP
       self.DYN_FLAG = DYN_FLAG
       self.DYN_SKIP = DYN_SKIP 
       self.DYN_M = DYN_M
       self.DYN_TH_ON = DYN_TH_ON
       self.DYN_TH_OFF = DYN_TH_OFF
       self.SIGMA = SIGMA
       self.V_0 = V_0
       self.TH_0 = TH_0
       self.A = A
       self.B = B
       self.DC = DC
       self.V1 = V1
       self.V2 = V2
       self.MU_SS = MU_SS
       self.V_SS = V_SS
       self.IOT = IOT
       self.IASP = IASP
       self.CO = CO
### Adding for parallel applications
       self.X = X
       self.Y = Y
       self.Z = Z
       self.DIP = DIP
      
   
def Qdyn_read_in(d,model,iname=None):
  if iname is None:
    filename=model+'qdyn.in'
  else:
    filename=model+iname
  print filename
  fid=open(filename)
  rline=fid.readline().split()
  d.MESHDIM=rline[0]
  print d.MESHDIM
  rline=fid.readline().split()
  d.NX=int(rline[0])
  d.NW=int(rline[1])
  print d.NX, d.NW
  rline=fid.readline().split()
  d.L = float(rline[0])
  d.W = float(rline[1])
  d.Z_CORNER = float(rline[2])
  d.DW = np.zeros(d.NW) 
  d.DIP_W = np.zeros(d.NW) 
  for i in range(d.NW):
    rline=fid.readline().split()
    d.DW[i]=rline[0] 
    d.DIP_W[i]=rline[1] 
  rline=fid.readline().split()
  d.THETA_LAW = rline[0] 
  rline=fid.readline().split()
  d.RNS_LAW = rline[0]
  rline=fid.readline().split()
  d.SIGMA_CPL = rline[0]
  rline=fid.readline().split()
  d.NEQS = rline[0]
  rline=fid.readline().split()
  d.NTOUT = rline[0] 
  d.IC = rline[1] 
  d.NXOUT = rline[2] 
  d.NXOUT_DYN = rline[3]
  d.OX_SEQ = rline[4]
  d.OX_DYN = rline[5]
  rline=fid.readline().split()
  d.VS = rline[0]
  d.MU = rline[1]
  d.LAM = rline[2]
  d.V_TH = rline[3]
  rline=fid.readline().split()
  d.TPER = rline[0] 
  d.APER = rline[1]
  rline=fid.readline().split()
  d.DTTRY = rline[0]
  d.DTMAX = rline[1]
  d.TMAX  = rline[2]
  d.ACC = rline[3]
  rline=fid.readline().split()
  d.NSTOP = rline[0] 
  rline=fid.readline().split()
  d.DYN_FLAG = rline[0]
  d.DYN_SKIP = rline[1]
  rline=fid.readline().split()
  d.DYN_M = rline[0]
  d.DYN_TH_ON = rline[1]
  d.DYN_TH_OFF = rline[2]
  Nxw = int(d.NX*d.NW)
  d.SIGMA=np.zeros(Nxw)
  d.V_0=np.zeros(Nxw)
  d.TH_0=np.zeros(Nxw)
  d.A=np.zeros(Nxw)
  d.B=np.zeros(Nxw)
  d.DC=np.zeros(Nxw)
  d.V1=np.zeros(Nxw)
  d.V2=np.zeros(Nxw)
  d.MU_SS=np.zeros(Nxw)
  d.V_SS=np.zeros(Nxw)
  d.IOT=np.zeros(Nxw)
  d.IASP=np.zeros(Nxw)
  d.CO=np.zeros(Nxw)
  d.X = np.zeros(Nxw)
  d.Y = np.zeros(Nxw)
  d.Z = np.zeros(Nxw)
  d.DIP = np.zeros(Nxw)
  for i in range(Nxw):
     rline = fid.readline().split()
     d.SIGMA[i] = float(rline[0])
     d.V_0[i]   = float(rline[1])
     d.TH_0[i]  = float(rline[2])
     d.A[i]     = float(rline[3])
     d.B[i]     = float(rline[4])
     d.DC[i]    = float(rline[5])
     d.V1[i]    = float(rline[6])
     d.V2[i]    = float(rline[7])
     d.MU_SS[i] = float(rline[8])
     d.V_SS[i]  = float(rline[9])
     d.IOT[i]   = float(rline[10])
     d.IASP[i]  = float(rline[11])
     d.CO[i]    = float(rline[12])
  if iname is not None:
     for ik in range(Nxw):
         rline = fid.readline().split()
         d.X[ik]   = float(rline[0])
         d.Y[ik]   = float(rline[1])
         d.Z[ik]   = float(rline[2])
         d.DIP[ik] = float(rline[3])
  fid.close()
#  print d.NX,d.NW,d.L,d.W,d.Z_CORNER,d.THETA_LAW,d.RNS_LAW,d.SIGMA_CPL,d.NEQS,d.NTOUT,d.IC,d.NXOUT,d.NXOUT_DYN,d.OX_SEQ,d.OX_DYN,d.VS,d.MU
#  print d.LAM,d.V_TH,d.TPER,d.APER,d.DTTRY,d.DTMAX,d.TMAX,d.ACC,d.NSTOP,d.DYN_FLAG,d.DYN_SKIP,d.DYN_M,d.DYN_TH_ON,d.DYN_TH_OFF
#  print d.SIGMA
  return d

#d=dqyn()
#dq=Qdyn_read_in(d)
#print dq.SIGMA
