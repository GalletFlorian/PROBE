! This is with the new data files from Florian Gallet, so the normalization 
! was done with respect to the Sun.

module Msun_data


use constants

implicit none


contains


subroutine readdata(Lum_tabltmp,Teff_tabltmp,toto_tabltmp,radius_tabltmp,r2conv_tabltmp,r2rad_tabltmp, &
				Mrad_tabltmp,Rrad_tabltmp,Qs_inv_tabltmp,sizetot,timedec)

use constants

implicit none

integer :: nlineheader = 1
character, dimension(:), allocatable :: header
integer :: nline(3),sizetot(3)
character :: filename*100, head*130
integer :: nt,nmax,i,j,k,dumm,l

real(dp),dimension(3,100000) :: toto_tabltmp,radius_tabltmp,Qs_inv_tabltmp,rots_tabl

real(dp) :: timedec

real(dp),dimension(3,100000) :: Lum_tabltmp,Teff_tabltmp, tstar_tabltmp,Rstar_tabltmp,r2conv_tabltmp
real(dp),dimension(3,100000) :: r2rad_tabltmp,Mrad_tabltmp,Rrad_tabltmp

real(dp), dimension(:,:), allocatable :: read_array2 

do i= 1,3
    if (i .eq. 1)  filename = './Model/mod10moa17short2.dat'
    if (i .eq. 2)  filename = './Model/mod10moa17short2.dat'
    if (i .eq. 3)  filename = './Model/mod10moa17short2.dat'
    
            open(unit=1,file=filename,status='old')
            read(1,*) head
            nt = 0
            do
            	read(1,*,end=10)
            	nt = nt + 1
            enddo	
            10 close(1)    
    nline(i) = nt  
enddo

	

nmax = int(maxval(nline))

allocate( header(nlineheader) )


do i = 1,3 

    if (i .eq. 1)  filename = './Model/mod10moa17short2.dat'
    if (i .eq. 2)  filename = './Model/mod10moa17short2.dat'
    if (i .eq. 3)  filename = './Model/mod10moa17short2.dat'
    nt = nline(i)
    
    allocate( read_array2(9, nt))    

    open(unit=1,file=filename,status='old')
    read(1,*) head
    do j =1,nline(i)
    	read(1,*) dumm,(read_array2(k,j),k=1,9)
    
    	
    enddo
            
    close(1)

    do j=1,nt
     	
    	Lum_tabltmp(i,j)	= read_array2(1,j)  
    	Teff_tabltmp(i,j)	= read_array2(2,j) 
    	toto_tabltmp(i,j)	= read_array2(3,j) 
    	radius_tabltmp(i,j)	= read_array2(4,j) 
    	r2conv_tabltmp(i,j)	= read_array2(5,j) 
    	r2rad_tabltmp(i,j)	= read_array2(6,j) 
    	Mrad_tabltmp(i,j)	= read_array2(7,j) 
    	Rrad_tabltmp(i,j)	= read_array2(8,j)
		Qs_inv_tabltmp(i,j) = read_array2(9,j)
    	   
    enddo
    
    
    l=1
    do while (l .le. nt)
    	if (r2rad_tabltmp(1,l) .ne. 0.0) then
    		timedec =  toto_tabltmp(1,l)
    		l = nt + 1
    	else 
    		l = l+1
    	endif
    
    enddo
    
    sizetot(i) = nt
    
	deallocate(read_array2)
enddo




end subroutine readdata

end module Msun_data


