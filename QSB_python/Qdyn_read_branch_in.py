#Importing dynamic parameters from the main fault loading to the branch fault
# Loading branch fault parameters qdyn_branch

#P. Galvez, AECOM, Oct. 2016
#Input (QDYN to SEM). 
import numpy as np
import os 

## Calling okada subroutine in fortran.
## and generte okada_gen.so python module 
## f2py okada_gen.f90 -m okada_gen -h okada_gen.pyf
## f2py -c -m okada_gen okada_gen.f90

def read_branch_in():
    model = './'
    fileb = 'qdyn_branch.in'
    fid = open(fileb,'r')
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
    CO_b  = np.zeros([N_branch])
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
        CO_b[j]=float(line[13])
        TAU_b[j]=SIGMA_b[j]*(MUSS_b[j]+A_b[j]*np.log(VSS_b[j]/V0_b[j])+B_b[j]*np.log(THSS_b[j]*V0_b[j]/DC_b[j]))
        j = j+1
    fid.close()
    return N_branch,DIP0,lam,mu,X_b,Y_b,Z_b,SIGMA_b,V0_b,THSS_b,A_b,B_b,DC_b,V1_b,V2_b,MUSS_b,VSS_b,CO_b,TAU_b
