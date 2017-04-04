#P. Galvez, AECOM, Sep. 2016
#Input (QDYN to SEM) .
import numpy as np

class dsnap:
    def _init_(self,X,Y,Z,time,v,th,vd,dtau,dtaud,d,sigma,t):
        self.X = X 
        self.Y = Y 
        self.Z = Z 
        self.time = time 
        self.v = v
        self.t = t
        self.th = th 
        self.vd = vd
        self.dtau = dtau 
        self.dtaud = dtaud
        self.d = d
        self.sigma = sigma


def Qdyn_read_ox_seq(d,filename):
    fid = open(filename)
    rline=fid.readline().split()
    d.it = int(rline[0])
    d.ivmax = int(rline[1])
    d.n = int(rline[2])
    d.t = float(rline[3])
    rline=fid.readline().split()
    rlines=fid.readlines()
    xl=[]
    yl=[]
    zl=[]
    timel=[]
    vl=[]
    thl=[]
    vdl=[]
    dtaul=[]
    dtaudl=[]
    ddl=[]
    sigmal=[]
    for rline in rlines:
        line=rline.split()
        xl.append(float(line[0])) 
        yl.append(float(line[1])) 
        zl.append(float(line[2])) 
        timel.append(float(line[3])) 
        vl.append(float(line[4])) 
        thl.append(float(line[5])) 
        vdl.append(float(line[6])) 
        dtaul.append(float(line[7])) 
        dtaudl.append(float(line[8])) 
        ddl.append(float(line[9])) 
        sigmal.append(float(line[10]))
    d.X = np.array(xl) 
    d.Y = np.array(yl) 
    d.Z = np.array(zl) 
    d.time = np.array(timel) 
    d.v = np.array(vl) 
    d.th = np.array(thl) 
    d.vd = np.array(vdl) 
    d.dtau = np.array(dtaul) 
    d.dtaud = np.array(dtaudl) 
    d.d = np.array(ddl) 
    d.sigma = np.array(sigmal)
    return d

#ds=dsnap()
#d=Qdyn_read_ox_seq(ds,'fort.20001')
