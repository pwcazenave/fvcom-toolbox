! a program to make a binary file to force unstructured swan4082 with idealized sea level
! z is positive up

!INPGRID WLEV UNSTRUCTURED EXC -99. NONSTAT 20070101_000000 1 HR 20070103_000000
!READINP WLEV 1. 'sealevel.bin' 1 0 UNF  
!
program fake_sealevel 
implicit none
real*4  :: amp = 2. !tidal amplitude in meters 
real*4  :: period = 44712.0 !period of M2 tide in seconds 
real*4  :: slev,time
integer :: Nverts  = 845 !number of nodes in the mesh
integer :: Nframes = 241 !number of hourly frames 
integer :: i,ii 



open(unit=33,file='sealevel.bin',form='unformatted')
time = 0.0
do i=1,nframes
  slev = amp*sin(time*2*3.14159/period) 
  write(*,*)i,float(i)/24.,time,slev 
  write(33)(slev,ii=1,Nverts)
  time = time + 3600
end do
close(33)

  
end program fake_sealevel
