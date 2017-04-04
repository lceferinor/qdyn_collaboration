#P. Galvez, AECOM, Sep. 2016
#Input (QDYN to SEM) . 
import numpy as np
import math 
import Qdyn_read_in as qdyn
import Qdyn_read_ox_seq as qdyn_ox
import os.path

t_dyn = 120;

x_off=0
y_off=0
z_off=0
istart=1

mname='nodesonfault_branch.csv'
NGLL = 5
itp =1

idummy =1
fdummy =1


# Loading mesh
print 'Loading SEM mesh:',mname 

fid = open(mname,'r')
rlines=fid.readlines()
qXl =[]
qYl =[]
qZl =[] 
for rline in rlines:
    line=rline.split()
    qXl.append(float(line[2]))
    qYl.append(float(line[3]))
    qZl.append(float(line[4]))
fid.close()
qX=np.array(qXl)
qY=np.array(qYl)
qZ=np.array(qZl)
qN = len(qXl)
#qNX = 0
#for i in range(len(qZl)):
#    if abs(qZ[i]-qZ[0])<1:
#       qNX = qNX+1
#qNZ=qN/qNX

#print 'qNX:',qNX,'qNZ:',qNZ
print 'qN:',qN

print 'Loading qdyn.in ...'

model = './'
pdyn=qdyn.dqyn()
p=qdyn.Qdyn_read_in(pdyn,model)
print p.SIGMA

print 'Done'

#Loading QDYN output
rname = model+'fort.'+str(20000)
print 'Loading QDYN output : ',rname,'...'
    
osnap = qdyn_ox.dsnap()
o = qdyn_ox.Qdyn_read_ox_seq(osnap,rname)
print o.d

p.X = o.X+x_off
p.Y = o.Y+y_off
p.Z = o.Z+z_off

# Hypocenter
loc_i=np.argmax(o.v)
v_max_i=max(o.v)
hx=p.X[loc_i]
hy=p.Y[loc_i]
hz=p.Z[loc_i]

co = np.zeros(qN)
tau = np.zeros(qN)

qSIGMA=np.zeros(qN)
qA=np.zeros(qN)
qB=np.zeros(qN)
qDC=np.zeros(qN)
qV0=np.zeros(qN)
qf0=np.zeros(qN)
qV_ini=np.zeros(qN)
qTH_ini=np.zeros(qN)

print 'Matching Values ...'

for i in range(qN):
    idexmin=np.argmin(np.square(p.X-qX[i]) + np.square(p.Z-qZ[i])) 
    idex=idexmin
    qSIGMA[i]=p.SIGMA[idex]
    qA[i]=p.A[idex]
    qB[i]=p.B[idex]
    qDC[i]=p.DC[idex]
    qV0[i]=p.V_SS[idex]
    qf0[i]=p.MU_SS[idex]
    qV_ini[i]=o.v[idex]
    qTH_ini[i]=o.th[idex]
    co[i]=p.CO[idex]
    tau[i]=qSIGMA[i]*(qf0[i]+qA[i]*np.log(qV_ini[i]/qV0[i])+qB[i]*np.log(qTH_ini[i]*qV0[i]/qDC[i]))
    if (i%math.ceil(qN/100)==0):
        print(i/math.ceil(qN/100), '% Complete')

print 'Done'

fid=open('rsf_hete_input_file.txt','w')
#txt='%u %u %E %E\n' % (qNX,qNZ,max(qX)-min(qX),max(abs(qZ))-min(abs(qZ)))
#fid.write(txt)
txt='%u\n' % (qN)
fid.write(txt)

print 'Start write SEM input:...'
for i in range(qN):
    txt='%E %E %E %E %E %E %E %E %E %E %E %E %E\n'%(qX[i],qZ[i],qSIGMA[i],tau[i],0.0,qV0[i],qf0[i],qA[i],qB[i],qDC[i],qV_ini[i],qTH_ini[i],co[i])
    fid.write(txt)
    if (i%math.ceil(qN/100)==0):
        print(i/math.ceil(qN/100), '% Complete')

print 'Generated SEM input: rsf_hete_input_file.txt'
fid.close()

if (os.path.isfile('timestamp.txt')):
    fid=open('timestamp.txt','r')
    rline=fid.readline().split()
    t0=float(rline[0])
    t1=float(rline[1])
    print '[t0 =',t0,'	t1=',t1,']'
    fid.close()
    fid=open('timestamp.txt','w')
    txt='%E %E'%(t1+t_dyn,t1+t_dyn+o.t)
    print '[t0 = ',t1+t_dyn,'	t1 = ',t1+t_dyn+o.t
    fid.close()
else:
    fid=open('timestamp.txt','w')
    txt='%E %E'%(0,o.t)
    fid.write(txt)
    fid.close()
    print '[t0 =',0,'	t1=',o.t,']'
