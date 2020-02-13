program condinit

    use fonctions 
    use Msun_data
    use constants
    use tide
    implicit none
    
    
real(dp) :: mass_planet(4),sma_tabl(2),Time_fin(3)

real(dp) :: dttry,Ms,t_init,Kwind,spin_sat
integer :: cst_diss,i

integer :: sizem,sizesmall,numbertoto

integer :: only_wind,ori,tideop

real(dp) :: t_init_tabl(3),Kwind_tabl(3),spin_sat_tabl(3),findgen(18),findgen2(45),deltat
real(dp) :: intmstar_arr(13),ss_eq_arr(13),d2ss_eq(13)
real(dp) :: ss_eq_val
real(dp) :: mass_star(3),Ms_star(3),ss_eq(3)
integer :: ims,sma,massp,babaorum,indice0,indicend
real(dp) :: Time ,ap0,ep0,obls0,oblp0,delta,Rs,n_orb
real(dp) :: T_start,rotst,T_star0,wconvi,wradi,Sun_dyn_freq,Q_str_inv,epsilon_squared
real(dp) :: lag_angle,rg2s,delta_t
real(dp) :: Mp,rmf,Rplanr,Rp
real(dp) :: grav,Jrad,Jconv,Irad,Iconv,DeltaJ_t,dt_t,dIconv,dMrad,term1,term2,term3,term4
real(dp) :: rotp0,k2deltat_plan
real(dp) :: sigmas,sigmas_eq,sigmap,rg2p,a_roche
character :: sonly_wind*10,scst_diss*10,smass_planet*10,ssma_tabl*10,smass_star*10
character :: head *170
real(dp),dimension(:),allocatable :: hamlet,hamlettmp,header
real(dp),dimension(:,:),allocatable :: horatio,read_array
real(dp),dimension(:),allocatable :: td,ad,ed,oblpd,oblsd,rotpd,rotsd,Rpd,Rsd


real(dp),dimension(3,100000) :: toto_tabltmp,radius_tabltmp,Qs_inv_tabltmp, Mass_startmp
real(dp),dimension(3,100000) :: Lum_tabltmp,Teff_tabltmp, tstar_tabltmp,Rstar_tabltmp,r2conv_tabltmp
real(dp),dimension(3,100000) :: r2rad_tabltmp,Mrad_tabltmp,Rrad_tabltmp


real(dp), dimension(:,:), allocatable :: toto_tabl,radius_tabl,Qs_inv_tabl
real(dp), dimension(:,:), allocatable :: Lum_tabl,Teff_tabl, tstar_tabl,Rstar_tabl,r2conv_tabl
real(dp), dimension(:,:), allocatable :: r2rad_tabl,Mrad_tabl,Rrad_tabl, Mass_star_tabl
real(dp), dimension(:), allocatable :: Lum,Teff, tstar,Rstar,r2conv,r2rad,Mrad,Rrad,Mass_star_evol
real(dp), dimension(:), allocatable :: maxv,toto,radius,Qs_inv,d2QSdt,d2RSdt
real(dp), dimension(:), allocatable :: d2Ldt,d2Teffdt,d2r2convdt,d2r2raddt,d2Mraddt,d2Rraddt,d2Massstardt

real(dp), dimension(:), allocatable :: totosmall,radiussmall,Qs_invsmall,Lumsmall,Teffsmall,r2convsmall,r2radsmall
real(dp), dimension(:), allocatable :: Mradsmall,Rradsmall,Massstarsmall

real(dp), dimension(2) :: Lum2,Teff2,tstar2,Rstar2,r2conv2,r2rad2,Mrad2,Rrad2,Mstarin2
real(dp) :: wradf,Wuprad,juprad,deltaJ,Iconvtmp,rotstmp,intmstar
real(dp) :: DAMbrake, DAMtide,Wtide,Wupconv,Wdownconv,Wjuprad,deltaWc

integer :: kk,nlineheader,nt,j,sizetot(3),sizemax,k,indice_BGM
character :: filename1*100,filename2*100,filename3*100,arg*30,diff*30,filename4*100,filename5*100

real(dp) :: orbitalper,timedec,Porb, corot, corotm, corotp,  corotm2, corotp2, corotm3, corotp3,  corotm4, corotp4
real(dp) :: corotm5, corotp5,corotm6,corotm7,corotm8,corotm9,corotm10,corotm11,corotm12,corotm13,corotm14,corotm15
real(dp) :: corotm16,corotm17,corotp6,corotp7


           
    !call getarg(1,arg)
    !read(arg,*) diff
    
    
    
    !Read data 

    !call readdata(toto_tabltmp,radius_tabltmp,Qs_inv_tabltmp,sizetot)
    
    call readdata(Lum_tabltmp,Teff_tabltmp,toto_tabltmp,radius_tabltmp,r2conv_tabltmp,r2rad_tabltmp, &
                Mrad_tabltmp,Rrad_tabltmp,Qs_inv_tabltmp,sizetot,timedec,Mass_startmp,indice_BGM)
                

    
    
    write(6,*) sizetot
    
    sizemax = max(maxval(sizetot),-1)
    
    write(6,*) sizemax
    
    
    !write(6,*) "Timedec = ", timedec, timedec*yr
    !stop
    
    
    allocate(toto_tabl(3,sizemax))
    allocate(radius_tabl(3,sizemax))
    allocate(Qs_inv_tabl(3,sizemax))    
    allocate(Lum_tabl(3,sizemax))
    allocate(Teff_tabl(3,sizemax))
    allocate(r2conv_tabl(3,sizemax))
    allocate(r2rad_tabl(3,sizemax))
    allocate(Mrad_tabl(3,sizemax))
    allocate(Rrad_tabl(3,sizemax))
    allocate(Mass_star_tabl(3,sizemax))
    
    
    
    
    do ims = 1,3
        do i = 1,sizetot(ims)
            Lum_tabl(ims,i) = Lum_tabltmp(ims,i)
            Teff_tabl(ims,i) = Teff_tabltmp(ims,i)
            toto_tabl(ims,i) = toto_tabltmp(ims,i)
            radius_tabl(ims,i) = radius_tabltmp(ims,i)
            r2conv_tabl(ims,i) = r2conv_tabltmp(ims,i)
            r2rad_tabl(ims,i) = r2rad_tabltmp(ims,i)
            Mrad_tabl(ims,i) = Mrad_tabltmp(ims,i)
            Rrad_tabl(ims,i) = Rrad_tabltmp(ims,i)
            Qs_inv_tabl(ims,i) = Qs_inv_tabltmp(ims,i)
            Mass_star_tabl(ims,i) = Mass_startmp(ims,i)
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
    
    
    


! Mass planet
!	 mass_planet = (/11130.0,1.,2.,10./)!35Jup
!	 mass_planet = (/4134.0,1.,2.,10./)	!13Jup
!	 mass_planet = (/1590.0,1.,2.,10./)	!5Jup
!	 mass_planet = (/1272.0,1.,2.,10./)	!4Jup
!	 mass_planet = (/954.0,1.,2.,10./)	!3Jup
!	 mass_planet = (/652.0,1.,2.,10./)	!WASP-43
!	 mass_planet = (/858.0,1.,2.,10./)	!2.7 Mjup
	 mass_planet = (/636.0,1.,2.,10./)	!2Jup
!	 mass_planet = (/318.0,1.,2.,10./)	!1Jup
!	 mass_planet = (/10.0,1.,2.,10./)	!10Earth
!	 mass_planet = (/1.0,1.,2.,10./)	!1Earth




    
    
    orbitalper = ((G*(mstar*Msun)*(3.*24.*3600.)**2./(4.*pi**2.))**(1./3.))/ AU
    
    
    !orbitalper = 2.e-2
    !To suppress effect of planet
    !orbitalper = 100000
    
    
    
    corot = ((G*(mstar*Msun)*(initper*24.*3600.)**2./(4.*pi**2.))**(1./3.))/ AU
	
	corotm = corot - 2*corot/100 
	corotp = corot + 2*corot/100 


	corotm2 = corot - 4*corot/100 
	corotp2 = corot + 4*corot/100 
	
	corotm3 = corot - 6*corot/100 
	corotp3 = corot + 6*corot/100 


	corotm4 = corot - 8*corot/100 
	corotp4 = corot + 8*corot/100 
	
	corotm5 = corot - 20*corot/100 
	corotp5 = corot + 20*corot/100 
	
	corotp6 = 1.1*corot
	
	corotp7 = 1.05*corot
	
	corotm6 = 0.1*corot
	corotm7 = 0.3*corot 
	corotm8 = 0.5*corot 
    corotm9 = 0.7*corot 
    corotm10 = 0.45*corot 
    corotm11 = 0.55*corot 
    corotm12 = 0.4*corot 
    corotm13 = 0.6*corot
    corotm14 = 0.65*corot
    corotm15 = 0.2*corot
    corotm16 = 0.8*corot
    corotm17 = 0.9*corot




	orbitalper = 100000.0
	
	


    sma_tabl = (/orbitalper, orbitalper/)

! Integration parameters
!only_wind = 0: tides and wind
    only_wind = 0

!Time_fin =  (/2.d10, 1.d10, 6.5d9/) * yr    !s


    !MAXVAL
    !Determines the maximum value of the elements in an array value
    !DIM = 2 here : determines the maximum value along y axes. Given an array
    !    3    6    9    12
    !    2    5    8    11
    !    1    4    7    10
    !MAXVAL (array,2) will produce a vector (if dim array =2) with max of each column => 3,6,9,12

    sizem = size(MAXVAL(toto_tabl,2))
    
    allocate( maxv(sizem))
    
    ! Maxv (vector) contain the max of each collumn. With above example
    ! maxv(1) = 3, maxv(2) = 6, maxv(3) = 9, maxv(4) = 12
    maxv = MAXVAL(toto_tabl,2)

    write(6,*) "Max? ",maxv(1),maxv(2), maxv(3)
    


    Time_fin =  (/maxv(1),maxv(2), maxv(3)/) * yr     !s
    
    
    

    
!write(6,*) Time_fin/yr
    dttry    =  1.0d-1 * yr         !first timestep in s (1.0d-3 * yr)
    cst_diss =  0
    
    !t_init_tabl   =  (/max(5d6, toto_tabl(1,1)), max(5d6, toto_tabl(2, 1)), max(5d6, toto_tabl(3, 1))/) * yr        !8.5d6  * yr
    Kwind_tabl    =  (/1.7d40, 1.7d40, 1.7d40/)          !47
    spin_sat_tabl =  (/9.0, 9.0, 9.0/)*spinsun              !s-1 (Bolmont 2012: wsat = 14 spinsun)

    mass_star =  (/mstar*10., mstar*10., mstar*10./)
    Ms_star   =  mass_star/10.d0 * Msun    !kg
    
    
    intmstar = nint(mstar*10.)

	!Modif pour présentation Grenoble/que pour ori=1
	!t_init_tabl   =  (/max(2.0d6, toto_tabl(1,1)), max(2.0d6, toto_tabl(2, 1)), max(2.0d6, toto_tabl(3, 1))/) * yr        !8.5d6  * yr
	
! Dissipation equilibrium tide, between (//) is the \overline{sigma_star} of Hansen 2012
    ss_eq = 6.4d-59*(/3.d-7, 3.d-7, 3.d-7/)*1d7
    
    
	ss_eq = (/0.0,0.0,0.0/)

	if (intmstar .eq. 3.0)	 ss_eq = 6.4d-59*(/4.00E-05,4.00E-05,4.00E-05/)*1d7
	if (intmstar .eq. 4.0)  ss_eq = 6.4d-59*(/2.00E-05,2.00E-05,2.00E-05/)*1d7
	if (intmstar .eq. 5.0)	 ss_eq = 6.4d-59*(/5.00E-06,5.00E-06,5.00E-06/)*1d7
	if (intmstar .eq. 6.0)	 ss_eq = 6.4d-59*(/2.50E-06,2.50E-06,2.50E-06/)*1d7
	if (intmstar .eq. 7.0)	 ss_eq = 6.4d-59*(/2.00E-06,2.00E-06,2.00E-06/)*1d7
	if (intmstar .eq. 8.0)	 ss_eq = 6.4d-59*(/1.00E-06,1.00E-06,1.00E-06/)*1d7
	if (intmstar .eq. 9.0)	 ss_eq = 6.4d-59*(/6.00E-07,6.00E-07,6.00E-07/)*1d7
	if (intmstar .eq. 10.)  ss_eq = 6.4d-59*(/3.00E-07,3.00E-07,3.00E-07/)*1d7
	if (intmstar .eq. 11.)	 ss_eq = 6.4d-59*(/2.00E-07,2.00E-07,2.00E-07/)*1d7
	if (intmstar .eq. 12.)	 ss_eq = 6.4d-59*(/7.00E-08,7.00E-08,7.00E-08/)*1d7
	if (intmstar .eq. 13.)	 ss_eq = 6.4d-59*(/2.00E-08,2.00E-08,2.00E-08/)*1d7
	if (intmstar .eq. 14.)	 ss_eq = 6.4d-59*(/3.00E-09,3.00E-09,3.00E-09/)*1d7
	if (intmstar .eq. 15.)   ss_eq = 6.4d-59*(/6.00E-13,3.00E-09,3.00E-09/)*1d7
	
	intmstar_arr(1) = 3.0 
	intmstar_arr(2) = 4.0
        intmstar_arr(3) = 5.0
        intmstar_arr(4) = 6.0
        intmstar_arr(5) = 7.0
        intmstar_arr(6) = 8.0
        intmstar_arr(7) = 9.0
        intmstar_arr(8) = 10.0
        intmstar_arr(9) = 11.0
        intmstar_arr(10) = 12.0
        intmstar_arr(11) = 13.0
        intmstar_arr(12) = 14.0
        intmstar_arr(13) = 15.0
 
	ss_eq_arr(1) = 4.00E-05
        ss_eq_arr(2) = 2.00E-05
 	ss_eq_arr(3) = 5.00E-06
 	ss_eq_arr(4) = 2.50E-06
 	ss_eq_arr(5) = 2.00E-06
 	ss_eq_arr(6) = 1.00E-06
 	ss_eq_arr(7) = 6.00E-07
	ss_eq_arr(8) = 3.00E-07
	ss_eq_arr(9) = 2.00E-07
 	ss_eq_arr(10) = 7.00E-08
 	ss_eq_arr(11) = 2.00E-08
 	ss_eq_arr(12) = 3.00E-09
 	ss_eq_arr(13) = 6.00E-13	

	do i = 1,13
		d2ss_eq(i) = 0.0
	enddo
	

	!intmstar_arr(13),ss_eq_arr(13),d2ss_eq(13)
	call splineme(intmstar_arr,ss_eq_arr,13,d2ss_eq)	
	call splint(intmstar_arr,ss_eq_arr,d2ss_eq,13,mstar*10.,ss_eq_val)

	
	ss_eq = 6.4d-59*(/ss_eq_val,ss_eq_val,ss_eq_val/)*1d7

	write(6,*) "Mstar = ", mstar, "SS_eq= ", ss_eq_val

	if ( ss_eq(1) .eq. 0.0) then

		write(6,*) "Equilibrium tide not define"
		write(6,*) "It could mean that the required stellar mass is not known" 
		write(6,*) intmstar 
		stop

	endif

    
    !To force the planet to appears at the end of the disk lifetime
    t_init_tabl = (/tdiskparam,tdiskparam,tdiskparam/) * yr
    
    write(6,*) "Tdiskparam = ", tdiskparam
    
    
    !t_init_tabl   =  (/max(5d6, toto_tabl(1,1)), max(5d6, toto_tabl(2, 1)), max(5d6, toto_tabl(3, 1))/) * yr

!for ims = 1, 1 do begin
    
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
            
              
            !write(sonly_wind,'(I0)') only_wind
            
            !write(scst_diss,'(I0)') cst_diss
            
            !write(smass_planet,'(f0.2)') mass_planet(massp)
            
            !write(ssma_tabl,'(f0.2)') sma_tabl(sma)
            
            !write(smass_star,'(f0.2)') mass_star(ims)
            
            !filename1 = 'infotides_only_wind_'//trim(sonly_wind)//'_cstdiss_'// &
            !trim(scst_diss)//'_sma_' //trim(ssma_tabl)//'_Mp_'//trim(smass_planet) &
            !//'_Ms_'//trim(smass_star)//'_'//trim(diff)//'.dat'
            
            !filename2 = 'datatides_only_wind_'//trim(sonly_wind)//'_cstdiss_'// &
            !trim(scst_diss)//'_sma_' //trim(ssma_tabl)//'_Mp_'//trim(smass_planet) &
            !//'_Ms_'//trim(smass_star)//'_'//trim(diff)//'.dat'
            
            !filename3 = 'wconvevol_'//trim(sonly_wind)//'_cstdiss_'// &
            !trim(scst_diss)//'_sma_' //trim(ssma_tabl)//'_Mp_'//trim(smass_planet) &
            !//'_Ms_'//trim(smass_star)//'_'//trim(diff)//'.dat'
            
            
            !filename1 = trim(filename1)
            !filename2 = trim(filename2)
            !filename3 = trim(filename3)
            
            
            filename1 = 'infotides.dat'
            filename2 =	'datatides.dat'
            filename3 = 'wconvevol.dat'
            filename4 = 'wconvevolfinal.dat'
            filename5 = 'torques.dat'
                        
            
        
            !***********************************************************************
            ! Two modes : 
            ! - or rotp is taken independently to be equal to
            !   pseudo-synchronisation value : babaorum = 1.
            ! - or rotp evolves according to its equation, in this
            !   case, it will evolve quickly to the pseudo-synchronisation
            !   value : babaorum = 2.
            
            

            babaorum = 1

            !***********************************************************************
            ! Properties of the considered bodies
            !*******************      STAR  PARAMETERS   ********************

            Ms = Ms_star(ims)
            write(6,*)'Ms = ',Ms/Msun,' Msun'

            t_init = t_init_tabl(ims)
                        
            
            write(6,*)'t_init',t_init/yr
            Kwind = Kwind_tabl(ims)
            spin_sat = spin_sat_tabl(ims)
            
            
            allocate (toto(sizetot(ims)))
            
            !toto = toto_tabl(ims, *)
            do i=1,sizetot(ims)
                toto(i) = toto_tabl(ims,i)
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
                !write(6,*) toto(i),radius(i),Qs_inv(i)
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
            
            
            allocate (Mass_star_evol(size(Mass_star_tabl,2)))
            
            !Mass = Mass_star_tabl(ims, *)
            do i=1,size(Mass_star_tabl,2)
                 Mass_star_evol(i) = Mass_star_tabl(ims,i)
            enddo
            
            
            
        
            ! Indice corresponding to t_init
            do i = 1,size(toto)
               if (toto(i) .ge. (t_init/yr-1d0)) then 
                  indice0 = max(1, i-1)
                  exit
               endif   
            enddo
            ! Indice corresponding to Time (end simu)
            do i = 1,size(toto)
               if (toto(i) .ge. (Time/yr-1d0)) then
                  indicend = i
                  exit
               endif   
            enddo
            
            
                
            

            write(6,*)"Verification time zero", toto(indice0), toto(1)
            write(6,*)"Verification end time", toto(indicend) 
            Rs     =  radius(indice0) * Rsun    !m
            write(6,*)'Rs =',Rs,' m  = ',Rs/Rsun,' Rsun = ',Rs/AU,' AU'
            
            
        

            ! Initialization of spin star
            

            !0.56
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

 
open(unit=23,file=filename4,status='new')
close(23)  

open(unit=23,file=filename5,status='new')
close(23)    
         
if (ori .eq. 0) then
	
	deltat = toto(2)*yr
	numbertoto = 2
		
endif


do while (deltat .lt. t_init)
		
        Lum2(1) = Lum(numbertoto-1)
        Lum2(2) = Lum(numbertoto)

        Teff2(1) = Teff(numbertoto-1)
        Teff2(2) = Teff(numbertoto)

        tstar2(1) = toto(numbertoto-1)
        tstar2(2) = toto(numbertoto)

        Rstar2(1) = radius(numbertoto-1)
        Rstar2(2) = radius(numbertoto)

        r2conv2(1) = r2conv(numbertoto-1)
        r2conv2(2) = r2conv(numbertoto)

        r2rad2(1) = r2rad(numbertoto-1)
        r2rad2(2) = r2rad(numbertoto)

        Mrad2(1) = Mrad(numbertoto-1)
        Mrad2(2) = Mrad(numbertoto)

        Rrad2(1) = Rrad(numbertoto-1)
        Rrad2(2) = Rrad(numbertoto)
        
        Mstarin2(1) = Mass_star_evol(numbertoto-1)
        Mstarin2(2) = Mass_star_evol(numbertoto)
        
        
        
		
        !call rotevol with initial velocity wconvi and wradi
        !calculates new rotation rate rotstmp and wradf
        
            
        call rotevol(taudecinit,tdiskparam,wconvi,rotstmp,wradi,wradf,Lum2,Teff2,tstar2,Rstar2,r2conv2,r2rad2,Mrad2,Rrad2 &
            ,Wuprad,juprad,deltaJ,timedec,DAMbrake,DAMtide,Wtide,Wupconv,Wdownconv,Wjuprad,deltaWc,Mstarin2)

            
        !open(unit=23,file=filename3,access = 'append',status='old')
        !write(23,*) (deltat-t_init)/yr/1.e6, rotstmp/2.87e-6,wradf/2.87e-6,0.00,DAMbrake, DAMtide, &
        !Wtide,Wupconv,Wdownconv,Wjuprad,deltaWc
        
        grav = log10(Mass_star_evol(numbertoto)*Msun*1.e3*G*1.e3/((radius(numbertoto)*Rsun*100)**2.))
        Jrad = wradf * r2rad(numbertoto) *Mass_star_evol(numbertoto)*Msun*1000*(radius(numbertoto)*Rsun*100)**2
    	Jconv = rotstmp * r2conv(numbertoto)*Mass_star_evol(numbertoto)*Msun*1000*(radius(numbertoto)*Rsun*100)**2
                     
        open(unit=23,file=filename3,access = 'append',status='old')
		write(23,*) (deltat-t_init)/yr/1.e6, rotstmp/2.87e-6,wradf/2.87e-6,Mass_star_evol(numbertoto),Teff(numbertoto),grav &
			,Mrad(numbertoto),Rrad(numbertoto),radius(numbertoto)*Rsun,r2conv(numbertoto),r2rad(numbertoto)  &
			, Jconv,Jrad, indice_BGM 
        
        close(23)
        
		Irad = r2rad(numbertoto)*Mass_star_evol(numbertoto)*Msun*1000*(radius(numbertoto)*Rsun*100.0)**2.0
		Iconv = r2conv(numbertoto)*Mass_star_evol(numbertoto)*Msun*1000*(radius(numbertoto)*Rsun*100.0)**2.0
    	
    		
		DeltaJ_t = (Iconv*Jrad-Irad*Jconv) / (Iconv+Irad)
		term1 = DeltaJ_t / (Iconv*taudecinit*yr)
		dt_t = (toto(numbertoto) - toto(numbertoto-1))*yr
		dMrad = Mrad(numbertoto) - Mrad(numbertoto-1)
		dIconv = Iconv - r2conv(numbertoto-1)*Mass_star_evol(numbertoto-1)*Msun*1000*(radius(numbertoto-1)*Rsun*100.0)**2.0            
		term2 = (2.0 * (radius(numbertoto)*Rsun*100.0)**2.0 * rotstmp * dMrad*Msun*1000)/(3.0*Iconv*dt_t)
		term3 = dIconv*rotstmp/(dt_t*Iconv)
		DAMbrake = abs(Wdownconv) * Iconv / dt_t
		term4 = DAMbrake / Iconv
		open(unit=23,file=filename5,access = 'append',status='old')
		write(23,*) Wdownconv,(deltat-t_init)/yr/1.e6,DeltaJ &
		,taudecinit,dt_t,dMrad*Msun*1000/dt_t,dIconv/dt_t,wdownconv,term1,term2,term3,term4
		close(23)
        
        numbertoto = numbertoto +1
        deltat = toto(numbertoto)*yr
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

            
            

        
        	! Allocate array hamlet of 10 elements
        	
        	
            allocate (hamlet(size(findgen2)))
            hamlet  =  findgen2*1+t_init/yr
            

            
            do kk = 1,239 
                allocate (hamlettmp(size(hamlet)))
                hamlettmp = hamlet 
                
                if (kk .lt. 6 ) then
                
                    deallocate (hamlet)
                    allocate (hamlet(size(findgen)+size(hamlettmp)))
                    hamlet = (/hamlettmp,t_init/yr+findgen*10**(float(kk))/)  
                endif
                 
                if ( (kk .ge. 6) .and. (kk .le. 8)) then
                
                    deallocate (hamlet)
                    allocate (hamlet(size(findgen2)+size(hamlettmp)))
                    hamlet = (/hamlettmp,t_init/yr+findgen2*10**(float(kk))/)  
                endif    
                if ( (kk .gt. 8) .and. (kk .le. 189) ) then
                    deallocate (hamlet)
                    allocate (hamlet(size(hamlettmp)+1))
                	hamlet = (/hamlettmp,t_init/yr+1d9+(kk-9)*5.d7/) 
                endif 
                if ( kk .gt. 189) then
                    deallocate (hamlet)
                    allocate (hamlet(size(hamlettmp)+1))
                	hamlet = (/hamlettmp,t_init/yr+1d10+(kk-189)*1.d8/) 
                endif 
                                
                deallocate (hamlettmp)
            enddo
            !write(6,*)hamlet-t_init/yr
            !write(6,*) "size=", size(hamlet)
            !stop


        
            allocate( horatio(3,size(hamlet)))

            ! Second derivative in Rsun/yr**2
            
            
            allocate(totosmall(indicend-(indice0)+1))
            allocate(radiussmall(indicend-(indice0)+1))
            allocate(Qs_invsmall(indicend-(indice0)+1))
            allocate(Lumsmall(indicend-(indice0)+1))
            allocate(Teffsmall(indicend-(indice0)+1))
            allocate(r2convsmall(indicend-(indice0)+1))
            allocate(r2radsmall(indicend-(indice0)+1))
            allocate(Mradsmall(indicend-(indice0)+1))
            allocate(Rradsmall(indicend-(indice0)+1))
            allocate(Massstarsmall(indicend-(indice0)+1))
            
            allocate(d2RSdt(indicend-(indice0)+1))
            allocate (d2QSdt(indicend-(indice0)+1) )
            allocate (d2Ldt(indicend-(indice0)+1) )
            allocate (d2Teffdt(indicend-(indice0)+1) )
            allocate (d2r2convdt(indicend-(indice0)+1) )
            allocate (d2r2raddt(indicend-(indice0)+1) )
            allocate (d2Mraddt(indicend-(indice0)+1) )
            allocate (d2Rraddt(indicend-(indice0)+1) )
            allocate (d2Massstardt(indicend-(indice0)+1) )

            
            
            do k=1,indicend-indice0+1
                totosmall(k) = toto(k+indice0-1)
                radiussmall(k) = radius(k+indice0-1)
                Qs_invsmall(k) = Qs_inv(k+indice0-1)
                Lumsmall(k) = Lum(k+indice0-1)
                Teffsmall(k) = Teff(k+indice0-1)
                r2convsmall(k) = r2conv(k+indice0-1)
                r2radsmall(k) = r2rad(k+indice0-1)
                Mradsmall(k) = Mrad(k+indice0-1)
                Rradsmall(k) = Rrad(k+indice0-1)  
                Massstarsmall(k) = Mass_star_evol(k+indice0-1)    
            enddo
            
                         
            sizesmall = int(indicend-indice0)
            
            
            call splineme(totosmall,radiussmall,indicend-indice0+1,d2RSdt)
            
            
            do i = indice0, indicend-1
                if (toto(i+1).eq.toto(i)) then 
                    write(6,*) 'two times are identical, spline is not going to work'
                    write(6,*) i, toto(i), toto(i+1)
                    stop
                endif
            enddo

            ! Second derivative in Rsun/yr**2
            
            call splineme(totosmall,Qs_invsmall,indicend-indice0+1,d2QSdt)    
            call splineme(totosmall,Lumsmall,indicend-indice0+1,d2Ldt)
            call splineme(totosmall,Teffsmall,indicend-indice0+1,d2Teffdt)
            call splineme(totosmall,r2convsmall,indicend-indice0+1,d2r2convdt)
            call splineme(totosmall,r2radsmall,indicend-indice0+1,d2r2raddt)
            call splineme(totosmall,Mradsmall,indicend-indice0+1,d2Mraddt)
            call splineme(totosmall,Rradsmall,indicend-indice0+1,d2Rraddt)
            call splineme(totosmall,Massstarsmall,indicend-indice0+1,d2Massstardt)
            
                        
            
            
            ! Call function that calculates integration
!            if (only_wind.eq.1) then 
!                !??? Call?
!                call TIDESTEST_NO_EVOLUTION(hamlet,horatio,t_init,ap0,ep0,obls0,oblp0,wconvi,rotp0 &
!                    ,Mp,Ms,Rp,Rs,d2RSdt,d2QSdt,epsilon_squared,sigmap,sigmas,dttry,rg2s,rg2p &
!                    ,filename1,filename2,Time,indice0,indicend &
!                    ,toto,radius,Qs_inv,a_roche,babaorum,spin_sat,Kwind,cst_diss,sigmas_eq, only_wind)
!                    
!            endif   

     
             !sizearray = size(Lum)
             !sizederiv = size(d2RSdt)
             
             
             write(6,*) "dt try = " ,dttry, yr
             
                                 
     
            if (only_wind.eq.0) then 
                call TIDESTEST(hamlet,horatio,t_init,ap0,ep0,obls0,oblp0,wconvi,wradi,rotp0 &
                    ,Mp,Ms,Rp,Rs,d2RSdt,d2QSdt,epsilon_squared,sigmap,sigmas,dttry,rg2s,rg2p &
                    ,filename1,filename2,filename3,filename4,filename5,Time,indice0,indicend &
                    ,totosmall,radiussmall,Qs_invsmall,a_roche,babaorum,spin_sat,Kwind,cst_diss,sigmas_eq,Sun_dyn_freq, &
                    Lumsmall,d2Ldt,Teffsmall,d2Teffdt,r2convsmall,d2r2convdt,r2radsmall,d2r2raddt,Mradsmall,d2Mraddt, &
                    Rradsmall,d2Rraddt,Massstarsmall,d2Massstardt,indicend-(indice0)+1, &
                    indicend-(indice0)+1,ori,timedec,tideop,indice_BGM)
                    
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
