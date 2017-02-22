!**************************************************************************************************
!* Sequential Code
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
Program HeatEquation

IMPLICIT NONE

integer							:: i, iNx, iNxg
double precision				:: dx, dt, t , Time, pi
double precision, allocatable	:: rF_T(:), rF_T_new(:) , Xc(:) , Xb(:)						
character*50					:: file_num
character*50					:: filename
double precision 			  	:: t1,t2

call cpu_time(t1) ! Starts Computation time

!************************************Initialization of Parameters**********************************************
pi = 3.1415927
dx = 0.001				! Grid Spacing
dt = (0.5*dx**2)/1.001	! Time Step
iNx = 1/dx 				! Number of Cell centers
iNxg = iNx+1 			! Number of Cell Boundaries
t = 0.0					! Initialization of time
Time = 2e-03			! Total time for which program run
allocate (rF_T(iNxg),rF_T_new(iNxg),Xb(iNxg),Xc(iNx))
!**************************************************************************************************************
!*
!****************************** Initial and Boundary conditions ***********************************************
Xb(1) = 0.0
Xc(1) = Xb(1) + 0.5*dx
rF_T(1) = 1.0
Xb(iNxg) = 1.0
rF_T_new (1) = 1.0
rF_T_new (iNx) = 0.0
!***********************************************************************************************************
Do i=2,iNx 
  Xb(i) = Xb(i-1)+dx
  Xc(i) = Xc(i-1)+dx
    rF_T(i) = 0.0		 				
End do
!****************************************************************************************************************
!*********************************** Time Loop Begin ************************************************************
!****************************************************************************************************************   
DO while (t <= Time)
!******************************Time Filter for Specific files at Specific Time Only (1 of 3)*********************
if((0.0<=t - 2e-04 .AND. t - 2e-04 <= 1e-06).OR.(0.0 <= t - 1.9995e-03)) then  
write (file_num,*) t
filename= 'solution_'//trim(file_num)//'.dat'
open (unit = 12, file = filename ,status='replace')
write(12,*) 'TITLE="NUMERICAL SOULTION FOR THE 1D HEAT ADVECTIVE EQUATION"'
write(12,*) 'VARIABLES= "Xc" "F_Central" '  
write(12,*) "ZONE I=",iNx,"F=Point"
write(12,"(7E15.7)")  0.0,rF_T_new(1)				
write(12,"(7E15.7)")  Xc(1),rF_T_new(1)				
end if
!*********************************************Time Filter Ends (1 End)*********************************************
do i=2,(iNx-1)
rF_T_new(i) = rF_T(i)+ (((1.0+0.001*rF_T(i)**pi)*dt)/(dx**2))*(rF_T(i+1)-(2*rF_T(i))+rF_T(i-1))												
!******************************Time Filter for Specific files at Specific Time Only (2 of 3)************************
if((0.0<=t - 2e-04 .AND. t - 2e-04 <= 1e-06).OR.(0.0 <= t - 1.9995e-03)) then
  write(12,"(7E15.7)") Xc(i), rF_T_new(i)				
end if  
!*********************************************Time Filter Ends (2 End)**********************************************
end do
!******************************Time Filter for Specific files at Specific Time Only (3 of 3)************************
if((0.0<=t - 2e-04 .AND. t - 2e-04 <= 1e-06).OR.(0.0 <= t - 1.9995e-03)) then
write(12,"(7E15.7)") Xc(iNx),rF_T_new(iNx)
write(12,"(7E15.7)") 1.0, rF_T_new(iNx)
close(12)
end if
!*********************************************Time Filter Ends (3 End)**********************************************
rF_T = rF_T_new
t = t + dt

END DO
call cpu_time(t2) ! End Computation Time
print*, 'Start time', t1, 'End time', t2, 'time Taken', t2-t1
end program HeatEquation
!*
!***********************************************Program End******************************************************************
!****************************************************************************************************************************
!*****                Time Filters need to be deactivated before running this code for dx = 0.0001                      *****
!****************************************************************************************************************************