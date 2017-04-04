#Importing dynamic parameters from the main fault loading to the branch fault
# Loading branch fault parameters qdyn_branch

#P. Galvez, AECOM, Oct. 2016
#Input (QDYN to SEM). 
import numpy as np
import math
import sys
import os 

sys.path.append('/Users/percy.galvez/Dropbox/qdyn_developer/utils/kernel_calc')
sys.path.append('/Users/percy.galvez/Dropbox/qdyn_developer/QSB_python')
 
import Qdyn_read_ox_seq as qdyn_ox
import Qdyn_read_in as qdyn
import Qdyn_read_branch_in as qdyn_branch
import okada_gen as ok

## Calling okada subroutine in fortran.
## and generte okada_gen.so python module 
## f2py okada_gen.f90 -m okada_gen -h okada_gen.pyf
## f2py -c -m okada_gen okada_gen.f90

N_branch,DIP0,lam,mu,X_b,Y_b,Z_b,SIGMA_b,V0_b,THSS_b,A_b,B_b,DC_b,V1_b,V2_b,MUSS_b,VSS_b,CO_b,TAU_b=qdyn_branch.read_branch_in() 
xoff = -64000.00
yoff = 0.0
zoff = 0.0
print SIGMA_b
#swifting 64km
#Loading QDYN output
X_b = X_b + xoff
Y_b = Y_b + yoff
Z_b = Z_b + zoff

model='./'
rname = model+'fort.'+str(20001) 
print 'Loading QDYN output : ',rname,'...'
osnap = qdyn_ox.dsnap()
o = qdyn_ox.Qdyn_read_ox_seq(osnap,rname)

print 'Loading qdyn.in ...' 
model = './'
pdyn=qdyn.dqyn()
p=qdyn.Qdyn_read_in(pdyn,model)
print p.SIGMA

rname0 = model+'fort.'+str(1000) 
print 'Loading QDYN output : ',rname0,'...'
osnap0 = qdyn_ox.dsnap()
o0 = qdyn_ox.Qdyn_read_ox_seq(osnap0,rname)
Ustrike = o.d-o0.d
Udip = 0
mode = 1
print 'Ustrike:',Ustrike

DIPM = p.DIP_W[0] 
DWM  = p.DW[0]
DLM  = DWM
print DIPM, DWM,DLM
N_main = o.n
print 'number of nodes main: ',N_main
X_m=np.zeros([N_main])    
Y_m=np.zeros([N_main])    
Z_m=np.zeros([N_main])    
X_m=o.X+xoff
Y_m=o.Y+yoff
Z_m=o.Z+zoff
print Z_m
print 'number of nodes branch: ',N_branch
V_new=np.zeros([N_branch])
th_new=np.zeros([N_branch])
sigma_new=np.zeros([N_branch])
tau_new=np.zeros([N_branch])


#   sigma=np.zeros([N_main,N_branch])
#   Calling okada to compute initial stresses on the branch 
#   as a result of the slip on the main fault
for j in range(N_branch):
      itau_new   = 0.0
      isigma_new = 0.0
      if (j%math.ceil(N_branch/100)==0):
         print(j/math.ceil(N_branch/100), '% Complete')
      for i in range(N_main):
         [itau,isigma,IRET]=ok.okada.compute_kernel_gen(lam,mu,X_m[i],Y_m[i],Z_m[i],DIPM,DLM,DWM,\
                            X_b[j],Y_b[j],Z_b[j],DIP0,Ustrike[j],Udip,mode)
         itau_new      = itau   + itau_new
         isigma_new    = isigma + isigma_new
      sigma_new[j] = SIGMA_b[j] + isigma_new
      tau_new[j] = TAU_b[j]+ itau_new
      if (abs(A_b[j]-B_b[j])<0.25*B_b[j]):
         V_new[j] = V0_b[j]
      else : 
         V_new[j] = V0_b[j]*np.exp((tau_new[j]/sigma_new[j] - MUSS_b[j])/(A_b[j]-B_b[j]))
      th_new[j] = DC_b[j]/V_new[j]
#      tau_new[j]=sigma_new[j]*(MUSS_b[j]+A_b[j]*np.log(V_new[j]/V0_b[j])+B_b[j]*np.log(th_new[j]*V0_b[j]/DC_b[j]))
print 'max_tau:',max(tau_new),'min_tau:',min(tau_new)
print 'max_sigma:',max(sigma_new),'min_tau:',min(sigma_new)
    
print 'Start write SEM branch input: rsf_hete_input_file2.txt'
fid_branch=open('rsf_hete_input_file2.txt','w')
txt='%u\n'%(N_branch)
fid_branch.write(txt)
for i in range(N_branch):
      txt='%E %E %E %E %E %E %E %E %E %E %E %E %E\n'%(X_b[i],Z_b[i],sigma_new[i],tau_new[i],0.0,V0_b[i],\
                                                    MUSS_b[i],A_b[i],B_b[i],DC_b[i],V_new[i],th_new[i],CO_b[i])
      fid_branch.write(txt)
      if (i%math.ceil(N_branch/100)==0):
         print(i/math.ceil(N_branch/100), '% Complete')
fid_branch.close()

