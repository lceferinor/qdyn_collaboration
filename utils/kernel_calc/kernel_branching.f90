! Computes Okada shear stress kernel 
! To compile: ifort -o kernel.exe  -I../../src/ kernel.f90 ../../src/okada.o

program kernel

  use okada, only : compute_kernel

  integer, parameter :: iin =15, iout=16
  integer, parameter :: FAULT_TYPE = 1 ! 1= strike-slip, 2=thrust

  double precision, dimension(:), allocatable :: x,y,z,dip,xx,ww 
  double precision, dimension(:), allocatable :: xb,yb,zb,dipb 
  double precision :: mu, lam, tau, sigma
  integer :: iret, i,j,ii, mode, nn, nw, nx

  open(unit=iin,FILE= 'kernel_main.in') 
  open(unit=iin_branch,FILE= 'kernel_branch.in') 
  open(unit=iout,FILE= 'kernel.out') 

  write(6,*) 'Calculate Okada kernel for general mesh'
  read(iin,*) nn
  read(iin,*) lam, mu
  allocate(x(nn),y(nn),z(nn),dip(nn),xx(nn),ww(nn))
  do i =1,nn
    read(iin,*) x(i),y(i),z(i),dip(i),xx(i),ww(i)
  end do
  !Reading branching mesh 
  write(6,*) 'Reading branch mesh'
  read(iin,*) nb
  read(iin,*) lam, mu
  allocate(xb(nb),yb(nb),zb(nb),dipb(nb))
  do i =1,nn
    read(iin,*) xb(i),yb(i),zb(i),dipb(i)
  end do
  do i=1,nn
    do j=1,nb
      call compute_kernel(lam,mu,x(i),y(i),z(i),dip(i),xx(i),ww(i), &
             xb(j),yb(j),zb(j),dipb(j),iret,tau,sigma,FAULT_TYPE)
      if (iret /= 0) then
        write(6,*) 'WARNING : Kernel singular, set value to 0, (i,j)=',i,j
        tau = 0d0
      endif
      write(iout,*) xb(j),yb(j),zb(j),tau,sigma
    end do
  end do
        
  write(6,*) 'Kernel calculation completed and stored in kernel.out'    
  close(iin)
  close(iout)
  
end program kernel
