!----------------------------------------------------------------------------
! MPI parallel routines for QDYN. 
! This file containes all the MPI routines used in QDYN to run in parallel.
! Modified from Parallel.f90 SPECFEM3D.
module my_mpi

! main parameter module for specfem simulations

  use mpi

  implicit none

  integer :: my_local_mpi_comm_world, my_local_mpi_comm_for_bcast

end module my_mpi

!----------------------------------------------------------------------------
! Subroutine to initialize MPI 
subroutine init_mpi()

  use my_mpi
  
  implicit none
  
  integer :: ier 

  call MPI_INIT(ier)
  if (ier /= 0 ) stop 'Error initializing MPI'

end subroutine init_mpi 

!-------------------------------------------------------------------------------------------------
subroutine world_size(sizeval)

  use my_mpi

  implicit none

  integer,intent(out) :: sizeval

  ! local parameters
  integer :: ier

  call MPI_COMM_SIZE(my_local_mpi_comm_world,sizeval,ier)
  if (ier /= 0 ) stop 'Error getting MPI world size'

end subroutine world_size
!----------------------------------------------------------------------------
! Retrieve processor number
subroutine world_rank(rank)

  use my_mpi

  implicit none

  integer,intent(out) :: rank

  integer :: ier

  call MPI_COMM_RANK(my_local_mpi_comm_world,rank,ier)
  if (ier /= 0 ) stop 'Error getting MPI rank'

end subroutine world_rank

!-------------------------------------------------------------------------------------------------
!Gather all MPI.
subroutine gather_allv(sendbuf, scounts, recvbufall, recvcountsall, recvoffsetall,recvcountstotal, NPROC)

  use my_mpi
  use constants, only: CUSTOM_REAL,CUSTOM_MPI_TYPE

  implicit none
  
  integer :: NPROC,myrank
  integer :: scounts,recvcountstotal
  real(kind=CUSTOM_REAL), dimension(scounts) :: sendbuf
  real(kind=CUSTOM_REAL), dimension(recvcountstotal) :: recvbufall
  integer, dimension(0:NPROC-1) :: recvcountsall,recvoffsetall

  integer ier,iproc

  !PG: sending the each processor data to the corresponding index in the recvbufall array.
  
  call MPI_ALLGATHERV(sendbuf,scounts,CUSTOM_MPI_TYPE,recvbufall,recvcountsall,&
                      recvoffsetall,CUSTOM_MPI_TYPE,my_local_mpi_comm_world,ier)

end subroutine gather_allv

!-------------------------------------------------------------------------------------------------
!Gather all MPI.
subroutine gather_allvdouble(sendbuf, scounts, recvbufall, recvcountsall, recvoffsetall,recvcountstotal, NPROC)

  use my_mpi
  use constants

  implicit none
  
  integer :: NPROC,myrank
  integer :: scounts,recvcountstotal
  double precision, dimension(scounts) :: sendbuf
  double precision, dimension(recvcountstotal) :: recvbufall
  integer, dimension(0:NPROC-1) :: recvcountsall,recvoffsetall

  integer ier,iproc

  !PG: sending the each processor data to the corresponding index in the recvbufall array.
  
  call MPI_ALLGATHERV(sendbuf,scounts,MPI_DOUBLE_PRECISION,recvbufall,recvcountsall,&
                      recvoffsetall,MPI_DOUBLE_PRECISION,my_local_mpi_comm_world,ier)

end subroutine gather_allvdouble

!-------------------------------------------------------------------------------------------------
subroutine gather_alli(sendbuf, recvbuf, myrank, NPROC)

  use my_mpi

  implicit none

  integer :: NPROC,myrank
  integer :: sendbuf
  integer, dimension(0:NPROC-1) :: recvbuf

  integer :: ier

  call MPI_ALLGATHER(sendbuf,1,MPI_INTEGER, &
                  recvbuf,1,MPI_INTEGER, &
                  my_local_mpi_comm_world,ier)

end subroutine gather_alli

!-------------------------------------------------------------------------------------------------
subroutine synchronize_all()

  use my_mpi

  implicit none

  integer :: ier

  ! synchronizes MPI processes
  call MPI_BARRIER(my_local_mpi_comm_world,ier)
  if (ier /= 0 ) stop 'Error synchronize MPI processes'

end subroutine synchronize_all

!-------------------------------------------------------------------------------------------------

