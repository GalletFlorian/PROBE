program condinit

    use fonctions 
    use Msun_data
    use constants
    use tide
    implicit none
    
    
real(dp) :: mass_planet(4),sma_tabl(2),Time_fin(3)

real(dp) :: dttry,Ms,t_init,Kwind,spin_sat
integer :: cst_diss,i

integer :: sizem,sizesmall,numberage

integer :: only_wind,ori,tideop

real(dp) :: t_init_tabl(3),Kwind_tabl(3),spin_sat_tabl(3),findgen(18),findgen2(45),deltat
real(dp) :: mass_star(3),Ms_star(3),ss_eq(3)
integer :: ims,sma,massp,pseudosync,indice0,indicend
real(dp) :: Time ,ap0,ep0,obls0,oblp0,delta,Rs,n_orb
real(dp) :: T_start,rotst,T_star0,wconvi,wradi,Sun_dyn_freq,Q_str_inv,epsilon_squared
real(dp) :: lag_angle,rg2s,delta_t
real(dp) :: Mp,rmf,Rplanr,Rp
real(dp) :: rotp0,k2deltat_plan
real(dp) :: sigmas,sigmas_eq,sigmap,rg2p,a_roche
character :: sonly_wind*10,scst_diss*10,smass_planet*10,ssma_tabl*10,smass_star*10
character :: head *170
real(dp),dimension(:),allocatable :: age_output,age_outputtmp,header
real(dp),dimension(:,:),allocatable :: age_output_final,read_array
real(dp),dimension(:),allocatable :: td,ad,ed,oblpd,oblsd,rotpd,rotsd,Rpd,Rsd


real(dp),dimension(3,100000) :: age_tabltmp,radius_tabltmp,Qs_inv_tabltmp
real(dp),dimension(3,100000) :: Lum_tabltmp,Teff_tabltmp, tstar_tabltmp,Rstar_tabltmp,r2conv_tabltmp
real(dp),dimension(3,100000) :: r2rad_tabltmp,Mrad_tabltmp,Rrad_tabltmp


real(dp), dimension(:,:), allocatable :: age_tabl,radius_tabl,Qs_inv_tabl
real(dp), dimension(:,:), allocatable :: Lum_tabl,Teff_tabl, tstar_tabl,Rstar_tabl,r2conv_tabl
real(dp), dimension(:,:), allocatable :: r2rad_tabl,Mrad_tabl,Rrad_tabl
real(dp), dimension(:), allocatable :: Lum,Teff, tstar,Rstar,r2conv,r2rad,Mrad,Rrad
real(dp), dimension(:), allocatable :: maxv,age,radius,Qs_inv,d2QSdt,d2RSdt
real(dp), dimension(:), allocatable :: d2Ldt,d2Teffdt,d2r2convdt,d2r2raddt,d2Mraddt,d2Rraddt

real(dp), dimension(:), allocatable :: agesmall,radiussmall,Qs_invsmall,Lumsmall,Teffsmall,r2convsmall,r2radsmall
real(dp), dimension(:), allocatable :: Mradsmall,Rradsmall

real(dp), dimension(2) :: Lum2,Teff2,tstar2,Rstar2,r2conv2,r2rad2,Mrad2,Rrad2
real(dp) :: wradf,Wuprad,juprad,deltaJ,Iconvtmp,rotstmp,intmstar
real(dp) :: DAMbrake, DAMtide,Wtide,Wupconv,Wdownconv,Wjuprad,deltaWc

integer :: kk,nlineheader,nt,j,sizetot(3),sizemax,k
character :: filename1*100,filename2*100,filename3*100,arg*30,diff*30

real(dp) :: orbitalper,timedec,Porb, corot, corotm, corotp,  corotm2, corotp2, corotm3, corotp3,  corotm4, corotp4
real(dp) :: corotm5, corotp5,corotm6,corotm7,corotm8,corotm9,corotm10,corotm11,corotm12,corotm13,corotm14,corotm15
real(dp) :: corotm16,corotm17,corotp6,corotp7
    
    !Read data 
    !call readdata(age_tabltmp,radius_tabltmp,Qs_inv_tabltmp,sizetot)
    
    call readdata(Lum_tabltmp,Teff_tabltmp,age_tabltmp,radius_tabltmp,r2conv_tabltmp,r2rad_tabltmp, &
                Mrad_tabltmp,Rrad_tabltmp,Qs_inv_tabltmp,sizetot,timedec)
    
    
    sizemax = max(maxval(sizetot),-1)
        
    allocate(age_tabl(3,sizemax))
    allocate(radius_tabl(3,sizemax))
    allocate(Qs_inv_tabl(3,sizemax))    
    allocate(Lum_tabl(3,sizemax))
    allocate(Teff_tabl(3,sizemax))
    allocate(r2conv_tabl(3,sizemax))
    allocate(r2rad_tabl(3,sizemax))
    allocate(Mrad_tabl(3,sizemax))
    allocate(Rrad_tabl(3,sizemax))
    
    
    
    do ims = 1,3
        do i = 1,sizetot(ims)
            Lum_tabl(ims,i) = Lum_tabltmp(ims,i)
            Teff_tabl(ims,i) = Teff_tabltmp(ims,i)
            age_tabl(ims,i) = age_tabltmp(ims,i)
            radius_tabl(ims,i) = radius_tabltmp(ims,i)
            r2conv_tabl(ims,i) = r2conv_tabltmp(ims,i)
            r2rad_tabl(ims,i) = r2rad_tabltmp(ims,i)
            Mrad_tabl(ims,i) = Mrad_tabltmp(ims,i)
            Rrad_tabl(ims,i) = Rrad_tabltmp(ims,i)
            Qs_inv_tabl(ims,i) = Qs_inv_tabltmp(ims,i)
        enddo
    
    
    enddo
    
    
!***********************************************************************
! Original code? 

    ori = 0 ! 1 = yes;  0 = no
        
    if (ori .eq. 1) then
        write(6,*) "Original code (without decoupling) selected."    
    else 
        ori = 0
        write(6,*) "Code with core/envelope decoupling selected."    
    endif
    
  
! Tides?

	tideop = 1  ! 1 = yes, 0 = No  
	
	if (tideop .eq. 0) then
	
        do i = 1,sizetot(ims)
            Qs_inv_tabl(ims,i) = 0.0
        enddo
	
	endif  
    
!***********************************************************************
    

	!Mass planet
	 mass_planet = (/1.0,1.0,1.0,10./)*318.0	!1Jup

    orbitalper = ((G*(mstar*Msun)*(3.*24.*3600.)**2./(4.*pi**2.))**(1./3.))/ AU

    !To suppress effect of planet
    !orbitalper = 100000
     
    corot = ((G*(mstar*Msun)*(initper*24.*3600.)**2./(4.*pi**2.))**(1./3.))/ AU

	orbitalper = 100000.0
	

    sma_tabl = (/orbitalper, orbitalper/) !choice of initial SMA

	! Integration parameters
	!only_wind = 0: tides and wind
    only_wind = 0


    !MAXVAL
    !Determines the maximum value of the elements in an array value
    !DIM = 2 here : determines the maximum value along y axes. Given an array
    !    3    6    9    12
    !    2    5    8    11
    !    1    4    7    10
    !MAXVAL (array,2) will produce a vector (if dim array =2) with max of each column => 3,6,9,12

    sizem = size(MAXVAL(age_tabl,2))
    
    allocate( maxv(sizem))
    
    ! Maxv (vector) contain the max of each collumn. With above example
    ! maxv(1) = 3, maxv(2) = 6, maxv(3) = 9, maxv(4) = 12
    maxv = MAXVAL(age_tabl,2)

    
    Time_fin =  (/maxv(1),maxv(2), maxv(3)/) * yr     !s
    

	!write(6,*) Time_fin/yr
    dttry    =  1.0d-1 * yr         !first timestep in s (1.0d-3 * yr)
    cst_diss =  0
    
    !t_init_tabl   =  (/max(5d6, age_tabl(1,1)), max(5d6, age_tabl(2, 1)), max(5d6, age_tabl(3, 1))/) * yr        !8.5d6  * yr
    Kwind_tabl    =  (/1.7d40, 1.7d40, 1.7d40/)          !47
    spin_sat_tabl =  (/9.0, 9.0, 9.0/)*spinsun              !s-1 (Bolmont 2012: wsat = 14 spinsun)

    mass_star =  (/mstar*10., mstar*10., mstar*10./)
    Ms_star   =  mass_star/10.d0 * Msun    !kg
    
    
    intmstar = nint(mstar*10.)

! Dissipation equilibrium tide, between (//) is the \overline{sigma_star} of Hansen 2012
    ss_eq = 6.4d-59*(/3.d-7, 3.d-7, 3.d-7/)*1d7
    
    
	ss_eq = (/0.0,0.0,0.0/)

	if (intmstar .eq. 3.0)	 ss_eq = 6.4d-59*(/4.00E-05,4.00E-05,4.00E-05/)*1d7
	if (intmstar .eq. 4.0)   ss_eq = 6.4d-59*(/2.00E-05,2.00E-05,2.00E-05/)*1d7
	if (intmstar .eq. 5.0)	 ss_eq = 6.4d-59*(/5.00E-06,5.00E-06,5.00E-06/)*1d7
	if (intmstar .eq. 6.0)	 ss_eq = 6.4d-59*(/2.50E-06,2.50E-06,2.50E-06/)*1d7
	if (intmstar .eq. 7.0)	 ss_eq = 6.4d-59*(/2.00E-06,2.00E-06,2.00E-06/)*1d7
	if (intmstar .eq. 8.0)	 ss_eq = 6.4d-59*(/1.00E-06,1.00E-06,1.00E-06/)*1d7
	if (intmstar .eq. 9.0)	 ss_eq = 6.4d-59*(/6.00E-07,6.00E-07,6.00E-07/)*1d7
	if (intmstar .eq. 10.)   ss_eq = 6.4d-59*(/3.00E-07,3.00E-07,3.00E-07/)*1d7
	if (intmstar .eq. 11.)	 ss_eq = 6.4d-59*(/2.00E-07,2.00E-07,2.00E-07/)*1d7
	if (intmstar .eq. 12.)	 ss_eq = 6.4d-59*(/7.00E-08,7.00E-08,7.00E-08/)*1d7
	if (intmstar .eq. 13.)	 ss_eq = 6.4d-59*(/2.00E-08,2.00E-08,2.00E-08/)*1d7
	if (intmstar .eq. 14.)	 ss_eq = 6.4d-59*(/3.00E-09,3.00E-09,3.00E-09/)*1d7


	if ( ss_eq(1) .eq. 0.0) then

		write(6,*) "Equilibrium tide not define"
		write(6,*) "It could mean that the required stellar mass is not known" 
		write(6,*) intmstar 
		stop

	endif

    
    !To force the planet to appears at the end of the disk lifetime
    t_init_tabl = (/tdiskparam,tdiskparam,tdiskparam/) * yr
    
    write(6,*) "Tdiskparam = ", tdiskparam
    
    
    !t_init_tabl   =  (/max(5d6, age_tabl(1,1)), max(5d6, age_tabl(2, 1)), max(5d6, age_tabl(3, 1))/) * yr

    
    !do ims = 1, size(Ms_star)
    do ims = 1,1
        Time = Time_fin(ims)
        !do sma = 1,size(sma_tabl)
        do sma = 1,1
            
            !do massp = 1,size(mass_planet)
            do massp = 1,1

            ! Orbital initial conditions:
            ap0   = sma_tabl(sma) * AU  !m
            ep0   = 0.d0
            obls0 = 0.d0
            oblp0 = 0.d0
            
                        
            filename1 = 'infotides.dat'
            filename2 =	'datatides.dat'
            filename3 = 'wconvevol.dat'
                        
            
        
            !***********************************************************************
            ! Two modes : 
            ! - or rotp is taken independently to be equal to
            !   pseudo-synchronisation value : pseudosync = 1.
            ! - or rotp evolves according to its equation, in this
            !   case, it will evolve quickly to the pseudo-synchronisation
            !   value : pseudosync = 2.
            
            

            pseudosync = 1

            !***********************************************************************
            ! Properties of the considered bodies
            !*******************      STAR  PARAMETERS   ********************

            Ms = Ms_star(ims)
            write(6,*)'Ms = ',Ms/Msun,' Msun'

            t_init = t_init_tabl(ims)
                        
            
            write(6,*)'t_init',t_init/yr
            Kwind = Kwind_tabl(ims)
            spin_sat = spin_sat_tabl(ims)
            
            
            allocate (age(sizetot(ims)))
            
            !age = age_tabl(ims, *)
            do i=1,sizetot(ims)
                age(i) = age_tabl(ims,i)
            enddo
            
            
            
            allocate (radius(sizemax))
        
            !radius = radius_tabl(ims, *)
            do i=1,sizetot(ims)
                radius(i) = radius_tabl(ims,i)
            enddo
            
            allocate (Qs_inv(sizemax))
            
            !Qs_inv = Qs_inv_tabl(ims, *)
            do i=1,sizetot(ims)
                 Qs_inv(i) = Qs_inv_tabl(ims,i)
            enddo


            
            allocate (Lum(size(Lum_tabl,2)))

            !Lum = Lum_tabl(ims, *)
            do i=1,size(Lum_tabl,2)
                 Lum(i) = Lum_tabl(ims,i)
            enddo
            
            
            allocate (Teff(size(Teff_tabl,2)))

            !Teff = Teff_tabl(ims, *)
            do i=1,size(Teff_tabl,2)
                 Teff(i) = Teff_tabl(ims,i)
            enddo
            
            
            do i=1,size(Lum_tabl,2)
                !write(6,*) age(i),radius(i),Qs_inv(i)
            enddo
            
            
            allocate (r2conv(size(r2conv_tabl,2)))

            !r2conv = r2conv_tabl(ims, *)
            do i=1,size(r2conv_tabl,2)
                 r2conv(i) = r2conv_tabl(ims,i)
            enddo
            
            
            allocate (r2rad(size(r2rad_tabl,2)))
            
            !r2rad = r2rad_tabl(ims, *)
            do i=1,size(r2rad_tabl,2)
                 r2rad(i) = r2rad_tabl(ims,i)
            enddo
            
            
            allocate (Mrad(size(Mrad_tabl,2)))
            
            !Mrad = Mrad_tabl(ims, *)
            do i=1,size(Mrad_tabl,2)
                 Mrad(i) = Mrad_tabl(ims,i)
            enddo
            
            
            allocate (Rrad(size(Rrad_tabl,2)))
            
            !Mrad = Mrad_tabl(ims, *)
            do i=1,size(Rrad_tabl,2)
                 Rrad(i) = Rrad_tabl(ims,i)
            enddo
            
        
            ! Indice corresponding to t_init
            do i = 1,size(age)
               if (age(i) .ge. (t_init/yr-1d0)) then 
                  indice0 = max(1, i-1)
                  exit
               endif   
            enddo
            ! Indice corresponding to Time (end simu)
            do i = 1,size(age)
               if (age(i) .ge. (Time/yr-1d0)) then
                  indicend = i
                  exit
               endif   
            enddo
            
         

            write(6,*)"Verification time zero", age(indice0), age(1)
            write(6,*)"Verification end time", age(indicend) 
            Rs     =  radius(indice0) * Rsun    !m
            write(6,*)'Rs =',Rs,' m  = ',Rs/Rsun,' Rsun = ',Rs/AU,' AU'
                

            ! Initialization of spin star
            T_start   =  initper* day          !s
            rotst     =  2.*Pi/T_start     !rad.s**-1
            T_star0   =  T_start
            wconvi     =  rotst
            wradi = 0.0 !Initial rotation for the core, if t_init < 1.5 Myr
            
            Porb = (G * (Ms_star(1)+ mass_planet(1)*Mearth))**(1./3.) *  &
            rotst**(-2./3.) / 1.49598d11 
            
            
           !wradi = 51.15 * 2.87e-6
            
            write(6,*) "Prob vs Obitalper= " , Porb,orbitalper
            write(6,*) " Rot init", rotst, T_start, initper    
            write(6,*) "Mass star and planet" , Ms_star(1) ,  mass_planet(1), mass_planet(1)*Mearth
    			
	

            ! Sun dynamical frequency squared
            Sun_dyn_freq = G*Msun/(Rsun*Rsun*Rsun)
            ! Inverse of Q structure
            Q_str_inv = Qs_inv(indice0)
            epsilon_squared = wconvi*wconvi/Sun_dyn_freq


            ! lag angle 
            ! Normal formula = 3.d0*epsilon_squared*Q_str_inv/(4.d0*k2s)
            ! but there is a k2* later, so I just don't put it here
            lag_angle = 3.d0*epsilon_squared*Q_str_inv/4.d0
            ! Dissipation needs Mp to get calculated so it's defined later

            rg2s        = rg2_sun

            !*******************      PLANET  PARAMETERS   ********************

            Mp      =   mass_planet(massp) * Mearth    !kg
            ! Choose what planet composition : (Fortney et al, 2007)

            if (Mp/Mearth.le.20.) then 
                ! Rock mass fraction for rock/iron planet: 
                rmf     =   0.71       ! earth proportion  
                Rplanr   =  Rearth * ((0.0592d0*rmf+0.0975d0)*(log10(mass_planet(massp)))**2 &
                    + (0.2337d0*rmf + 0.4938d0)*log10(mass_planet(massp)) &
                    + (0.3102d0*rmf + 0.7932d0))
                Rp  =  Rplanr
            else 
                Rp = 1.d0 * Rjup
            endif
            
            ! Here value of dissipation star that needs Mp and ap0
            n_orb = sqrt(G*(Ms+Mp)/(ap0*ap0*ap0))
            delta_t = lag_angle/abs(n_orb-wconvi)
            if (cst_diss.eq.0) then 
                sigmas = 2.d0*G* delta_t /(3.d0*Rs*Rs*Rs*Rs*Rs) + ss_eq(ims)
                sigmas_eq = ss_eq(ims)
            endif
            if (cst_diss.eq.1) then
                sigmas_eq = ss_eq(ims)
                sigmas = sigmas_eq
            endif
            write(6,*)'sigmas = ',sigmas
            write(6,*)'wconvi =',wconvi,' s-1'
            write(6,*)'initial rot period = ',2.*Pi/(wconvi*day),'  day'

            ! Planet rotation, here pseudo-synch
            rotp0   = (1+15./2.*ep0**2+45./8.*ep0**4+5./16.*ep0**6)* &
                 1./(1+3*ep0**2+3./8.*ep0**4)*1./(1-ep0**2)**1.5 * n_orb !rad.s**-1
            write(6,*)'rotp0 =',rotp0,' s-1'
            write(6,*)'initial rot period = ',2.*Pi/(rotp0*day),'  day'

            ! Sigma calculation : 
            k2deltat_plan =  213.0                      !s
            !sigmap  =  (2.*G*k2deltat_plan/(3.*Rp**5))   !kg**-1.m**-2.s**-1   
            sigmap  =  0.0d0                            !kg**-1.m**-2.s**-1   

            rg2p    =  rg2_earth

            ! Roche limit
            a_roche = (Rp/0.462)* &
                      (Mp/Ms)**(-1./3.)   !(Faber et al, 2005! Pacynski,1971)
                      
                      
            !***********************************************************************
            ! 		Estimation of wrad during disk lifetime
            !*********************************************************************** 
			open(unit=23,file=filename3,status='new')
			close(23)                
			if (ori .eq. 0) then
				deltat = age(2)*yr
				numberage = 2		
			endif


			do while (deltat .lt. t_init)
		
        		Lum2(1) = Lum(numberage-1)
        		Lum2(2) = Lum(numberage)

        		Teff2(1) = Teff(numberage-1)
        		Teff2(2) = Teff(numberage)

        		tstar2(1) = age(numberage-1)
        		tstar2(2) = age(numberage)

        		Rstar2(1) = radius(numberage-1)
        		Rstar2(2) = radius(numberage)

        		r2conv2(1) = r2conv(numberage-1)
        		r2conv2(2) = r2conv(numberage)

        		r2rad2(1) = r2rad(numberage-1)
        		r2rad2(2) = r2rad(numberage)

        		Mrad2(1) = Mrad(numberage-1)
        		Mrad2(2) = Mrad(numberage)
	
       			Rrad2(1) = Rrad(numberage-1)
        		Rrad2(2) = Rrad(numberage)

        		!call rotevol with initial velocity wconvi and wradi
        		!calculates new rotation rate rotstmp and wradf
    
        		call rotevol(taudecinit,tdiskparam,wconvi,rotstmp,wradi,wradf,Lum2,Teff2,tstar2,Rstar2,r2conv2,r2rad2,Mrad2,Rrad2 &
            	,Wuprad,juprad,deltaJ,timedec,DAMbrake,DAMtide,Wtide,Wupconv,Wdownconv,Wjuprad,deltaWc)

            
        		open(unit=23,file=filename3,access = 'append',status='old')
        		write(23,*) (deltat-t_init)/yr/1.e6, rotstmp/2.87e-6,wradf/2.87e-6,0.00,DAMbrake, DAMtide, &
        		Wtide,Wupconv,Wdownconv,Wjuprad,deltaWc
        		close(23)
        
        		numberage = numberage +1
        		deltat = age(numberage)*yr
        		wconvi   = wconvi
        		wradi    = wradf
			enddo

            !***********************************************************************
            !***********************************************************************
            
            !Produce array of 10 element from 0 to 9  
            do i=1,18
                !findgen(i) = 1.0 * i*0.5
                findgen(i) = 0.5*i+0.5
            enddo
            
            do i=1,45
                !findgen(i) = 1.0 * i*0.5
                findgen2(i) = 0.2*i+0.8
            enddo
        
        	! Allocate array age_output of 10 elements        	
            allocate (age_output(size(findgen2)))
            age_output  =  findgen2*1+t_init/yr
            

            
            do kk = 1,239 
                allocate (age_outputtmp(size(age_output)))
                age_outputtmp = age_output 
                
                if (kk .lt. 6 ) then
                
                    deallocate (age_output)
                    allocate (age_output(size(findgen)+size(age_outputtmp)))
                    age_output = (/age_outputtmp,t_init/yr+findgen*10**(float(kk))/)  
                endif
                 
                if ( (kk .ge. 6) .and. (kk .le. 8)) then
                
                    deallocate (age_output)
                    allocate (age_output(size(findgen2)+size(age_outputtmp)))
                    age_output = (/age_outputtmp,t_init/yr+findgen2*10**(float(kk))/)  
                endif    
                if ( (kk .gt. 8) .and. (kk .le. 189) ) then
                    deallocate (age_output)
                    allocate (age_output(size(age_outputtmp)+1))
                	age_output = (/age_outputtmp,t_init/yr+1d9+(kk-9)*5.d7/) 
                endif 
                if ( kk .gt. 189) then
                    deallocate (age_output)
                    allocate (age_output(size(age_outputtmp)+1))
                	age_output = (/age_outputtmp,t_init/yr+1d10+(kk-189)*1.d8/) 
                endif 
                                
                deallocate (age_outputtmp)
            enddo
            !write(6,*)age_output-t_init/yr
            !write(6,*) "size=", size(age_output)
            !stop


        
            allocate( age_output_final(3,size(age_output)))

            ! Second derivative in Rsun/yr**2
            
            
            allocate(agesmall(indicend-(indice0)+1))
            allocate(radiussmall(indicend-(indice0)+1))
            allocate(Qs_invsmall(indicend-(indice0)+1))
            allocate(Lumsmall(indicend-(indice0)+1))
            allocate(Teffsmall(indicend-(indice0)+1))
            allocate(r2convsmall(indicend-(indice0)+1))
            allocate(r2radsmall(indicend-(indice0)+1))
            allocate(Mradsmall(indicend-(indice0)+1))
            allocate(Rradsmall(indicend-(indice0)+1))
            
            allocate(d2RSdt(indicend-(indice0)+1))
            allocate (d2QSdt(indicend-(indice0)+1) )
            allocate (d2Ldt(indicend-(indice0)+1) )
            allocate (d2Teffdt(indicend-(indice0)+1) )
            allocate (d2r2convdt(indicend-(indice0)+1) )
            allocate (d2r2raddt(indicend-(indice0)+1) )
            allocate (d2Mraddt(indicend-(indice0)+1) )
            allocate (d2Rraddt(indicend-(indice0)+1) )

            
            
            do k=1,indicend-indice0+1
                agesmall(k) = age(k+indice0-1)
                radiussmall(k) = radius(k+indice0-1)
                Qs_invsmall(k) = Qs_inv(k+indice0-1)
                Lumsmall(k) = Lum(k+indice0-1)
                Teffsmall(k) = Teff(k+indice0-1)
                r2convsmall(k) = r2conv(k+indice0-1)
                r2radsmall(k) = r2rad(k+indice0-1)
                Mradsmall(k) = Mrad(k+indice0-1)
                Rradsmall(k) = Rrad(k+indice0-1)                
            enddo
            
                         
            sizesmall = int(indicend-indice0)
            
            
            call splineme(agesmall,radiussmall,indicend-indice0+1,d2RSdt)
            
            
            do i = indice0, indicend-1
                if (age(i+1).eq.age(i)) then 
                    write(6,*) 'two times are identical, spline is not going to work'
                    write(6,*) i, age(i), age(i+1)
                    stop
                endif
            enddo

            ! Second derivative in Rsun/yr**2
            
            call splineme(agesmall,Qs_invsmall,indicend-indice0+1,d2QSdt)    
            call splineme(agesmall,Lumsmall,indicend-indice0+1,d2Ldt)
            call splineme(agesmall,Teffsmall,indicend-indice0+1,d2Teffdt)
            call splineme(agesmall,r2convsmall,indicend-indice0+1,d2r2convdt)
            call splineme(agesmall,r2radsmall,indicend-indice0+1,d2r2raddt)
            call splineme(agesmall,Mradsmall,indicend-indice0+1,d2Mraddt)
            call splineme(agesmall,Rradsmall,indicend-indice0+1,d2Rraddt)
            
                        
            
            
            ! Call function that calculates integration
!            if (only_wind.eq.1) then 
!                !??? Call?
!                call TIDESTEST_NO_EVOLUTION(age_output,age_output_final,t_init,ap0,ep0,obls0,oblp0,wconvi,rotp0 &
!                    ,Mp,Ms,Rp,Rs,d2RSdt,d2QSdt,epsilon_squared,sigmap,sigmas,dttry,rg2s,rg2p &
!                    ,filename1,filename2,Time,indice0,indicend &
!                    ,age,radius,Qs_inv,a_roche,pseudosync,spin_sat,Kwind,cst_diss,sigmas_eq, only_wind)
!                    
!            endif   

     
             !sizearray = size(Lum)
             !sizederiv = size(d2RSdt)
             
             
             write(6,*) "dt try = " ,dttry, yr
             
                                 
     
            if (only_wind.eq.0) then 
                call TIDESTEST(age_output,age_output_final,t_init,ap0,ep0,obls0,oblp0,wconvi,wradi,rotp0 &
                    ,Mp,Ms,Rp,Rs,d2RSdt,d2QSdt,epsilon_squared,sigmap,sigmas,dttry,rg2s,rg2p &
                    ,filename1,filename2,filename3,Time,indice0,indicend &
                    ,agesmall,radiussmall,Qs_invsmall,a_roche,pseudosync,spin_sat,Kwind,cst_diss,sigmas_eq,Sun_dyn_freq, &
                    Lumsmall,d2Ldt,Teffsmall,d2Teffdt,r2convsmall,d2r2convdt,r2radsmall,d2r2raddt,Mradsmall,d2Mraddt, &
                    Rradsmall,d2Rraddt,indicend-(indice0)+1, &
                    indicend-(indice0)+1,ori,timedec,tideop)
                    
            endif
            
            
            write(6,*) "On sort là?"

            nlineheader=1
            
            allocate( header(nlineheader)) 
            
            open(unit=1,file=filename2,status='old')
            nt = 0
            read(1,*)
            do
                read(1,*,end=10)
                nt = nt + 1
            enddo    
            10 close(1)

            allocate( td(nt) )
            allocate( ad(nt) )
            allocate( ed(nt) )
            allocate( oblpd(nt) )
            allocate( oblsd(nt) )
            allocate( rotpd(nt) )
            allocate( rotsd(nt) )
            allocate( Rpd(nt) )
            allocate( Rsd(nt) )
            
            
            allocate( read_array(15,nt) )
                        
            open(unit=1,file=filename2,status='old')
            read(1,*) head
            do j =1,nt
                    read(1,*) (read_array(i,j),i=1,15)
            enddo
            
            close(1)
            
                    
            !do j=1,nt
            !    td(j)     = read_array(1,j)
            !    ad(j)    = read_array(2,j)
            !       ed(j)    = read_array(3,j)
            !    oblpd(j) = read_array(4,j)
            !    oblsd(j) = read_array(5,j)
            !    rotpd(j) = read_array(6,j)
            !    rotsd(j) = read_array(7,j)
            !    Rpd(j)   = read_array(8,j)
            !    Rsd(j)   = read_array(9,j)
            
            !enddo

        enddo
    enddo
enddo

end program condinit
