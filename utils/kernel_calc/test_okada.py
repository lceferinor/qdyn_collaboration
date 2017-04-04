#Test 
import okada_gen as ok

fid=open('kernel.in','r')
fidout=open('kernel_py.in','w')
rline = fid.readline()
mode=int(rline)
print mode
rline = fid.readline()
N = int(rline)
print N
rline = fid.readline()
line=rline.split()
lam = float(line[0])
mu  = float(line[1])
print lam, mu
rlines=fid.readlines()
X=[]
Y=[]
Z=[]
XX=[]
WW=[]
DIP=[]
tau=0.0
sigma=0.0
IRET=0
for rline in rlines:
    line=rline.split()
    X.append(float(line[0]))
    Y.append(float(line[1]))
    Z.append(float(line[2]))
    DIP.append(float(line[3]))
    XX.append(float(line[4]))
    WW.append(float(line[5]))

Ustrike=1.0
Udip = 0.0
for i in range(N):
  for j in range(N):
    [tau,sigma,IRET]=ok.okada.compute_kernel_gen(lam,mu,X[i],Y[i],Z[i],DIP[i],XX[i],WW[i],\
                         X[j],Y[j],Z[j],DIP[j],Ustrike,Udip,mode)
    if (IRET!=0) :
       print 'WARNING: Kernel singular, set value to 0, (i,j)=',i,j
       tau=0.0
    txt='%20.6f\n'%(float(tau))
    fidout.write(txt)
fidout.close()
fid.close()
