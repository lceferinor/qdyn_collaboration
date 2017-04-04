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
import okada_gen as ok

## Calling okada subroutine in fortran.
## and generte okada_gen.so python module 
## f2py okada_gen.f90 -m okada_gen -h okada_gen.pyf
## f2py -c -m okada_gen okada_gen.f90

def read_branch_in():
    xoff = -64000.00
    yoff = 0.0
    zoff = 0.0
    phi = 90.0
    model = './'
    fileb = 'qdyn_branch.in'
    nprocs=6
    fid = open(fileb)
    rline = fid.readline()
    line=rline.split()
    N_branch=int(line[0])
    rline = fid.readline()
    line=rline.split()
    DIP0=float(line[0])
    rline = fid.readline()
    line=rline.split()
    lam=float(line[0])
    mu=float(line[1])
    rlines = fid.readlines()
    print N_branch,DIP0,lam,mu
    X_b = np.zeros([N_branch])
    Y_b = np.zeros([N_branch])
    Z_b = np.zeros([N_branch])
    SIGMA_b = np.zeros([N_branch])
    V0_b  = np.zeros([N_branch])
    THSS_b = np.zeros([N_branch])
    A_b   = np.zeros([N_branch])
    B_b   = np.zeros([N_branch])
    DC_b  = np.zeros([N_branch])
    V1_b  = np.zeros([N_branch])
    V2_b  = np.zeros([N_branch])
    MUSS_b= np.zeros([N_branch])
    VSS_b = np.zeros([N_branch])
    C0_b  = np.zeros([N_branch])
    TAU_b  = np.zeros([N_branch])
    j=0
    for rline in rlines:
        line=rline.split()
        X_b[j]=float(line[0]) 
        Y_b[j]=float(line[1]) 
        Z_b[j]=float(line[2]) 
        SIGMA_b[j]=float(line[3]) 
        V0_b[j]=float(line[4]) 
        THSS_b[j]=float(line[5]) 
        A_b[j]=float(line[6]) 
        B_b[j]=float(line[7]) 
        DC_b[j]=float(line[8]) 
        V1_b[j]=float(line[9])
        V2_b[j]=float(line[10])
        MUSS_b[j]=float(line[11])
        VSS_b[j]=float(line[12])
        C0_b[j]=float(line[13])
        TAU_b[j]=SIGMA_b[j]*(MUSS_b[j]+A_b[j]*np.log(VSS_b[j]/V0_b[j])+B_b[j]*np.log(THSS_b[j]*V0_b[j]/DC_b[j]))
        j = j+1
#    SPECFEM3D already seek for the neareast node on the mesh. 
#    nameb = 'nodesonfault_branch.csv'
    # Loading mesh
#    print 'Loading SEM mesh:',nameb
#    fid = open(nameb,'r')
#    rlines=fid.readlines()
#    qXl =np.zeros([N_branch])
#    qYl =np.zeros([N_branch])
#    qZl =[] 
#    for rline in rlines:
#      line=rline.split()
#      qXl.append(float(line[2]))
#      qYl.append(float(line[3]))
#      qZl.append(float(line[4]))
#    fid.close()
#    qX=np.array(qXl)
#    qY=np.array(qYl)
#    qZ=np.array(qZl)
#    qN = len(qXl)
#    fid.close()
    print 'Branching fault'
    print 'mintaub',min(TAU_b) 
    print 'maxtaub',max(TAU_b) 
    #Translation swifting 64km
#    qXl = X_b + xoff
#    qYl = Y_b + yoff
    X_b = X_b + xoff
    Y_b = Y_b + yoff
    Z_b = Z_b + zoff
#   Rotation of coordiantes the rotation after the computation of Okada
#    for i in range(N_branch):
#       X_b[i] =  qXl[i]*np.cos(phi*np.pi/180)+qYl[i]*np.sin(phi*np.pi/180)
#       Y_b[i] = -qXl[i]*np.sin(phi*np.pi/180)+qYl[i]*np.cos(phi*np.pi/180)
    
#    rname = model+'fort.'+str(20001) 
    rname = model+'fort.'+str(20000) 
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
#   Translation    
#    o.X=o.X+xoff
#    o.Y=o.Y+yoff
    X_m=o.X+xoff
    Y_m=o.Y+yoff
    Z_m=o.Z+zoff
#   Rotation
#    for i in range(N_branch):
#       X_m[i] =  o.X[i]*np.cos(phi*np.pi/180)+o.Y[i]*np.sin(phi*np.pi/180)
#       Y_m[i] = -o.X[i]*np.sin(phi*np.pi/180)+o.Y[i]*np.cos(phi*np.pi/180)

    print Y_m
    print 'number of nodes branch: ',N_branch
    V_new=np.zeros([N_branch])
    th_new=np.zeros([N_branch])
    sigma_new=np.zeros([N_branch])
    tau_new=np.zeros([N_branch])
    N_branch_loc=np.zeros([nprocs])
    n0=np.zeros([nprocs+1])
    j0=0
    if (N_branch%(nprocs)==0):
      for iproc in (range(nprocs)):
        N_branch_loc[iproc] = int(math.floor(N_branch/(nprocs)))
        j0=N_branch_loc[iproc]+j0
        n0[iproc+1]=j0
    else:
      for iproc in (range(nprocs)):
        if (iproc==0):
          N_branch_loc[iproc] = int(N_branch - (nprocs-1)*math.floor(N_branch/(nprocs-1)))
          j0=N_branch_loc[iproc]+j0         
        else:
          N_branch_loc[iproc] = int(math.floor(N_branch/(nprocs-1)))
          j0=N_branch_loc[iproc]+j0
        n0[iproc+1]=j0
    print n0,N_branch_loc,sum(N_branch_loc)
    for iproc in (range(nprocs)):
      filein='dynamic_branch'+str(iproc)+'.par'
      fid_out=open(filein,'w')
      txt = '%d\n'%(N_main)
      fid_out.write(txt)
      txt = '%20.6f %20.6f %20.6f %20.6f %20.6f\n'%(lam,mu,DIPM,DLM,DWM)
      fid_out.write(txt)
      for i in (range(N_main)):
        txt = '%20.6f %20.6f %20.6f\n'%(X_m[i],Y_m[i],Z_m[i])
        fid_out.write(txt)
      txt = '%d\n'%(N_branch_loc[iproc])
      fid_out.write(txt)
      txt = '%d\n'%(mode)
      fid_out.write(txt)
      txt = '%20.6f\n'%(DIP0)
      fid_out.write(txt)
      print N_branch_loc[iproc]
      for jloc in (range(int(N_branch_loc[iproc]))):
        j = jloc+n0[iproc]
        txt = '%20.6E %20.6E %20.6E %20.6E %20.6E %20.6E %20.6E %20.6E %20.6E %20.6E %20.6E %20.6E %20.6E\n'%\
               (X_b[j],Y_b[j],Z_b[j],SIGMA_b[j],TAU_b[j],V0_b[j],\
               A_b[j],B_b[j],DC_b[j],MUSS_b[j],Ustrike[j],Udip,C0_b[j]) 
        fid_out.write(txt)
     #   print j,n0[iproc],jloc,iproc
      fid_out.close()
         
read_branch_in()
#   sigma=np.zeros([N_main,N_branch])
#   Calling okada to compute initial stresses on the branch 
#   as a result of the slip on the main fault
#    for j in range(N_branch):
#      itau_new   = 0.0
#      isigma_new = 0.0
#      if (j%math.ceil(N_branch/100)==0):
#         print(j/math.ceil(N_branch/100), '% Complete')
#      for i in range(N_main):
#         [itau,isigma,IRET]=ok.okada.compute_kernel_gen(lam,mu,X_m[i],Y_m[i],Z_m[i],DIPM,DLM,DWM,\
#                            X_b[j],Y_b[j],Z_b[j],DIP0,Ustrike[j],Udip,mode)
#         itau_new      = itau   + itau_new
#         isigma_new    = isigma + isigma_new
#      sigma_new[j] = SIGMA_b[j] + isigma_new
#      tau_new[j] = TAU_b[j]+ itau_new
#      V_new[j]  = V0_b[j]*np.exp((tau_new[j]/sigma_new[j] - MUSS_b[j])/(A_b[j]-B_b[j]))
#      th_new[j] = DC_b[j]/V_new[j]
#      tau_new[j]=sigma_new[j]*(MUSS_b[j]+A_b[j]*np.log(V_new[j]/V0_b[j])+B_b[j]*np.log(th_new[j]*V0_b[j]/DC_b[j]))
#    print 'max_tau:',max(tau_new),'min_tau:',min(tau_new)
#    print 'max_sigma:',max(sigma_new),'min_tau:',min(sigma_new)
#    
#    print 'Start write SEM branch input: rsf_hete_input_file2.txt'
#    fid_branch=open('rsf_hete_input_file2.txt','w')
#    txt=''
#    for i in range(N_branch):
#      txt='%E %E %E %E %E %E %E %E %E %E %E %E %E\n'%(X_b[i],Z_b[i],sigma_new[i],tau_new[i],0.0,V0_b[i],\
#                                                    MUSS_b[i],A_b[i],B_b[i],DC_b[i],V_new[i],th_new[i],CO_b[i])
#      fid_branch.write(txt)
#      if (i%math.ceil(N_branch/100)==0):
#         print(i/math.ceil(N_branch/100), '% Complete')
#
#read_branch_in()
