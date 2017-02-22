!**************************************************************************************************
!* Parallel Code
!* For Grid Spacing, dx =0.001
!*
!* Created on : 11-02-2016
!*
!* by : Akshat Srivastava
!* @ Cranfield University (2015-2016)
!***************************************************************************************************
!* 
!*****************************************Main Program Begin******************************************
!*
program OneDHeatTransfer_MSc_HPC_Assignment
implicit none
Include 'mpif.h'
  
integer 					  :: i, iNx, iNxg, myid, ierr, nprocs, ChunkSize, request, status(MPI_STATUS_SIZE), tag1, tag2
double precision 			  :: dx, dt, t, Time, pi
double precision, allocatable :: rF_T_new(:), rF_T(:),Xb(:),Xc(:)
character*50				  :: file_num
character*50				  :: filename
double precision 			  :: t1,t2

call cpu_time(t1) ! Starts Computation time
call MPI_INIT(ierr) ! Initialize MPI
call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)

!************************************Initialization of Parameters**********************************************
pi = 3.146
dx = 0.001 				! Grid Spacing
iNx = 1/dx 				! Number of Cell centres
iNxg = iNx+1 			! Number of Cell Boundaries
Time = 0.002 			! Total time for which program run
t = 0.0 				! Initialization of time
dt = (0.5*dx**2)/1.001  ! Time step
ChunkSize = 333 		! Domain discretised for three processors
allocate (rF_T_new(iNxg),rF_T(iNxg),Xb(iNxg),Xc(iNx))
!*
!**************************************************************************************************************
!*
!****************************** Initial condition *************************************************************
Xb(1)=0.0 
Xb(iNxg)=1.0
Xc(1)=(Xb(1))+(0.5*dx)
!*
!***********************************************************************************************************
!*
do i=2, iNxg-1
Xb(i)=(Xb(i-1))+dx
Xc(i)=(Xc(i-1))+dx
end do
!****************************************************************************************************************
!*********************************** Time Loop Begin ************************************************************
!****************************************************************************************************************
do while(t.le.Time)
tag1=1
tag2=2
!********************************************************************************************************************
!****                                                Processor_0                                                *****
!********************************************************************************************************************
if (myid .eq. 0) then 

rF_T(1)= 1.0
rF_T_new(1)=1.0

if (t==0.0) then
F_T(i)+(dt*(1.0 + 0.001 * rF_T(i)**pi)/(dx**2))*((rF_T(i+1))-(2*rF_T(i))+(rF_T(i-1)))
end do

!******************************Time Filter for Specific files at Specific Time Only**************************************
if((0.0<=t - 2e-04 .AND. t - 2e-04 <= 1e-06).OR.(0.0 <= t - 1.9995e-03)) then
write (file_num,*) t
filename= 'solution_P0'//trim(file_num)//'.dat'
open (unit = 12, file = filename ,status='replace')
write(12,*) 'TITLE="SOULTION FOR THE 1D HEAT EQUATION"'
write(12,*) 'VARIABLES= "Cell Centers" "Temperature Flux" '  
write(12,*) "ZONE I=",ChunkSize,"F=Point"

do i = 1, chunksize  
write(12,"(7E15.7)") Xc(i), rF_T_new(i)	
end do
close (12)								
end if
!*********************************************Time Filter Ends*************************************************************
rF_T= rF_T_new
call MPI_issend (rF_T(ChunkSize), 1, MPI_DOUBLE_PRECISION, myid+1, tag1, MPI_COMM_WORLD, request, ierr)
Call MPI_request_free (request, ierr)
 
end if
!********************************************************************************************************************
!*****                                                   Processor_1                                            *****
!********************************************************************************************************************
if (myid .eq. 1) then
  
if (t==0.0) then
do i=1,ChunkSize
rF_T(i)=0.0
end do
end if
  
call MPI_irecv(rF_T(1), 1, MPI_DOUBLE_PRECISION, myid-1, tag1, MPI_COMM_WORLD, request, ierr)
Call MPI_wait (request, status, ierr)  

do i=1,ChunkSize
rF_T_new(i)=rF_T(i)+(dt*(1.0 + 0.001 * rF_T(i)**pi)/(dx**2))*((rF_T(i+1))-(2*rF_T(i))+(rF_T(i-1)))
rF_T= rF_T_new
end do
  
!******************************Time Filter for Specific files at Specific Time Only**************************************
if((0.0<=t - 2e-04 .AND. t - 2e-04 <= 1e-06).OR.(0.0 <= t - 1.9995e-03)) then
write (file_num,*) t
filename= 'solution_P1'//trim(file_num)//'.dat'
open (unit = 13, file = filename ,status='replace')
write(13,*) 'TITLE="SOULTION FOR THE 1D HEAT EQUATION"'
write(13,*) 'VARIABLES= "Cell Centers" "Temperature Flux" '  
write(13,*) "ZONE I=",ChunkSize,"F=Point"

do i = 1, chunksize  
write(13,"(7E15.7)") Xc(i+333), rF_T_new(i)	
end do
close (13)								
end if
!*********************************************Time Filter Ends*************************************************************

Call MPI_issend (rF_T(ChunkSize), 1, MPI_DOUBLE_PRECISION, myid+1, tag2, MPI_COMM_WORLD, request, ierr)
Call MPI_request_free (request, ierr)

end if
!**************************************************************************************************************************
!*****                                                   Processor_2                                                  *****
!**************************************************************************************************************************
if(myid.eq.2)then !Processor_2
  
if (t==0.0) then
do i=1,ChunkSize
rF_T(i)=0.0
end do
end if
  
rF_T(ChunkSize)= 0.0 !Boundary Conditions

call MPI_irecv(rF_T(1), 1, MPI_DOUBLE_PRECISION, myid-1, tag2, MPI_COMM_WORLD, request, ierr)
Call MPI_wait (request, status, ierr)  

do i=1,ChunkSize
rF_T_new(i)=rF_T(i)+(dt*(1.0 + 0.001 * rF_T(i)**pi)/(dx**2))*((rF_T(i+1))-(2*rF_T(i))+(rF_T(i-1)))
rF_T=rF_T_new
end do

!******************************Time Filter for Specific files at Specific Time Only**************************************
if((0.0<=t - 2e-04 .AND. t - 2e-04 <= 1e-06).OR.(0.0 <= t - 1.9995e-03)) then
write (file_num,*) t
filename= 'solution_P2'//trim(file_num)//'.dat'
open (unit = 14, file = filename ,status='replace')
write(14,*) 'TITLE="SOULTION FOR THE 1D HEAT EQUATION"'
write(14,*) 'VARIABLES= "Cell Centers" "Temperature Flux" '  
write(14,*) "ZONE I=",ChunkSize,"F=Point"

do i = 1, chunksize  
write(14,"(7E15.7)") Xc(i+666), rF_T_new(i)
end do	
close (14)								
end if
!*********************************************Time Filter Ends*************************************************************
end if
t = t + dt
print*, time, nprocs, myid ! Validation for Parallel Coding
end do
!**************************************************************************************************************************
!************************************************Time Loop Ends************************************************************  
!**************************************************************************************************************************
call MPI_FINALIZE(ierr)
call cpu_time(t2) ! End Computation Time
print*, 'Start time', t1, 'End time', t2, 'time Taken', t2-t1, myid ! CPU time of all 3 Processors
end program OneDHeatTransfer_MSc_HPC_Assignment
!*
!***********************************************Program End******************************************************************
!****************************************************************************************************************************
!*****                Time Filters need to be deactivated before running this code for dx = 0.0001                      *****
!****************************************************************************************************************************