! a program to make a binary file to force unstructured swan4082 with idealized wind 
!
!INPGRID WIND UNSTRUCTURED EXC -99. NONSTAT 20070101_000000 24 HR 20070108_000000
!READINP WIND 1. 'wind.bin' 1 0 UNF

program const_wind
implicit none
real*4  :: mag(1) = (/5/) !wind mag in m/s
real*4  :: dang,u,v,unitx,unity,ang,cang
integer :: Nverts,i,ii,j,k,cnt,Switch
Switch = 48   !switch every Switch hours
Nverts = 845  !number of nodes in the mesh
dang   = 45.0 


open(unit=33,file='wind.bin',form='unformatted')
ang = 0
cnt = 0
do i=1,9 
  cang  = 3.14159*(90-ang)/180.
  unitx = cos(cang)
  unity = sin(cang)
  do j=1,1
    u   = -unitx*mag(j)
    v   = -unity*mag(j)
    do k=1,switch
      write(*,*)cnt,u,v,cang*180/3.14159,mag(j)
!      write(33)float(cnt)
!      write(33)(u,v,ii=1,Nverts)
      write(33)(u,ii=1,Nverts)
      write(33)(v,ii=1,Nverts)
      cnt = cnt + 1
    end do
  end do
  ang = ang + dang
end do
write(*,*)'final cnt',cnt
close(33)
end
