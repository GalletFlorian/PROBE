module tide

use constants
use fonctions
use rotevoldec

implicit none

contains

subroutine TIDESTEST(hamlet,horatio,t_init,ap0,ep0,obls0,oblp0,wconvi,wradi,rotp0 &
    ,Mp,Ms,Rp,Rs,d2RSdt,d2QSdt,epsilon_squared,sigmap,sigmas,dttry,rg2s,rg2p &
    ,filename1,filename2,filename3,Time,indice0,indicend &
    ,toto,radius,Qs_inv,a_roche,babaorum,spin_sat,Kwind,cst_diss,sigmas_eq,Sun_dyn_freq, &
                    Lum,d2Ldt,Teff,d2Teffdt,r2conv,d2r2convdt,r2rad,d2r2raddt,Mrad,d2Mraddt,Rrad,d2Rraddt,sizev, &
                    sizeder,ori,timedec,tideop)

use constants

implicit none

! Cash-Karp parameters for embedded Runga-Kutta Method 
! (Numerical recipes in fortran, the Art of Scientific
! Computing. W.H. Press, S.A. Teukolsky, W.T. Vetterling,
! B.P. Flannery)


!input

integer :: sizev,sizeder,ori,flagcrash,flagcrashout,out

real(dp) :: timedec

!Declare these param
real(dp) :: t_init,ap0,ep0,obls0,oblp0,wconvi,rotp0 &
    ,Mp,Ms,Rp,epsilon_squared,sigmap,sigmas,dttry,rg2s,rg2p &
    ,Time,a_roche,spin_sat,Kwind,sigmas_eq,Sun_dyn_freq
    
    
integer :: indice0,indicend,babaorum,cst_diss,tideop


real(dp), dimension(sizeder) :: d2RSdt,d2QSdt,d2Ldt,d2Teffdt,d2r2convdt,d2r2raddt,d2Mraddt,d2Rraddt
real(dp), dimension(sizev) :: Lum,Teff, tstar,r2conv,r2rad,Mrad,Rrad,toto,radius,Qs_inv


real(dp) :: dt, dtmin,dtmax,eps,safety,pgrow,pshrink,errcon
real(dp) :: Rs0,rg2s0,frotstmp,alphac,Q_str_inv


real(dp) :: Ls0,Teffs0, r2convs0,r2rads0,Rrads0,Mrads0
real(dp) :: Ls,Teffs, r2convs,r2rads,Rrads,Mrads
    
real(dp) :: hamlet(474),Rstar(10000),horatio(3,474)
real(dp) :: Rs

real(dp) :: tmp,delta_t,dttemp,errmax,fwindtmp,infinity,dtnext
real(dp) :: ka0,ka1,ka2,ka3,ka4,ka5,kas,ka
real(dp) :: ke0,ke1,ke2,ke3,ke4,ke5,kes,ke
real(dp) :: kop0,kop1,kop2,kop3,kop4,kop5,kops,kop
real(dp) :: kos0,kos1,kos2,kos3,kos4,kos5,koss,kos
real(dp) :: kp0,kp1,kp2,kp3,kp4,kp5,kps,kp
real(dp) :: ks0,ks1,ks2,ks3,ks4,ks5,kss,ks
real(dp) :: kw0,kw1,kw2,kw3,kw4,kw5,kws,kw
real(dp) :: rotstmp1,rotstmp2,rotstmp3,rotstmp4,rotstmp5


character :: filename1*100,filename2*100,filename3*100

real(dp) :: lag_angle,lop0,n_orb
real(dp) :: rg2s1,rg2s2,rg2s3,rg2s4,rg2s5,rg2sdot

real(dp) :: rotpnew,Rso
real(dp) :: yscal(7),error(7)

character :: stime*10
real(dp) :: enerdotp1,enerdots1

real(dp) :: alphaout,norbout

integer :: debug,j,i,int,k,kkk,numbertoto,comp,flagtinit

real(dp),dimension(6) :: aa,cc,ccs,sigmas_tabl,Rs_tabl
real(dp),dimension(6) :: Ls_tabl,Teffs_tabl,r2convs_tabl,r2rads_tabl,Mrads_tabl,Rrads_tabl,rotstmp_tabl
real(dp) :: bb(6,6)
real(dp), dimension(:), allocatable :: t,a,e,oblp,obls,rots,rotp,aequ,enerdotp,enerdots,sigmast,tenerdots

integer :: ndata_line

real(dp) :: ttmp,atmp,etmp,oblptmp,oblstmp,rotstmp,rotptmp,aequtmp,sigmastmp


real(dp) :: wradi,wradf,Wuprad,juprad,deltaJ,Iconvtmp
real(dp) :: DAMbrake, DAMtide,AMtottmp,AMtotinit,DAMtot,DAMwind
real(dp) :: Wtide,Wupconv,Wdownconv,Wjuprad,deltaWc


real(dp), dimension(2) :: Lum2,Teff2,tstar2,Rstar2,r2conv2,r2rad2,Mrad2,Rrad2




! Initial values 

out = 0

comp = 1
flagtinit = 0.0
flagcrashout = 0

!Used for rotevol

Wuprad = 0.0
juprad = 0.0
DeltaJ = 0.0

! Runge -Kutta parameters

aa(1)    =  0.0
aa(2)    =  1./5.
aa(3)    =  3./10.
aa(4)    =  3./5.
aa(5)    =  1.
aa(6)    =  7./8.

bb(1,2)  =  1./5.
bb(1,3)  =  3./40.
bb(1,4)  =  3./10.
bb(1,5)  =  -11./54.
bb(1,6)  =  1631./55296.

bb(2,3)  =  9./40.
bb(2,4)  =  -9./10.
bb(2,5)  =  5./2.
bb(2,6)  =  175./512.

bb(3,4)  =  6./5.
bb(3,5)  =  -70./27.
bb(3,6)  =  575./13824.

bb(4,5)  =  35./27.
bb(4,6)  =  44275./110592.

bb(5,6)  =  253./4096.

cc(1)    =  37./378.
cc(2)    =  0.0
cc(3)    =  250./621.
cc(4)    =  125./594.
cc(5)    =  0.0
cc(6)    =  512./1771.

ccs(1)   =  2825./27648.
ccs(2)   =  0.0
ccs(3)   =  18575./48384.
ccs(4)   =  13525./55296.
ccs(5)   =  277./14336.
ccs(6)   =  1./4.

!write(6,*) "Qs_inv=", Qs_inv

!Flag to know if the planet crash on the star 
!Initialization at 0
flagcrash = 0 

!*******************************************************************************
!*******************************************************************************
!*******************************************************************************

!!write(6,*)'Begin tables creation'
! Creation of tables of a, e, rot
ndata_line  =  size(hamlet)

allocate( t(ndata_line) )
allocate( a(ndata_line) )
allocate( e(ndata_line) )
allocate( oblp(ndata_line) )
allocate( obls(ndata_line) )
allocate( rots(ndata_line) )
allocate( rotp(ndata_line) )
allocate( aequ(ndata_line) )
allocate( enerdotp(ndata_line) )
allocate( enerdots(ndata_line) )
allocate( sigmast(ndata_line) )

if (ori .eq. 1) then
	DAMbrake = 0.0
	DAMtide = 0.0
endif




!Initial conditions (bibi)
ttmp          =  t_init
atmp          =  ap0
etmp          =  ep0
oblptmp       =  oblp0
oblstmp       =  obls0
rotstmp       =  wconvi
rotptmp       =  rotp0
call aeq(ep0,rotp0,wconvi,G,Mp,Ms,Rp,Rs,sigmap,sigmas,aequtmp)
sigmastmp     =  sigmas
frotstmp      = 0.
fwindtmp      = 0.

t(1)         =  ttmp
a(1)         =  atmp
e(1)         =  etmp
oblp(1)      =  oblptmp
obls(1)      =  oblstmp
rots(1)      =  rotstmp
rotp(1)      =  rotptmp
aequ(1)      =  aequtmp


call enerdot(ap0,ep0,rotp0,G,Mp,Ms,Rp,sigmap,enerdotp1) 
call enerdot(ap0,ep0,wconvi,G,Ms,Mp,Rs,sigmas,enerdots1) 

enerdotp(1) = enerdotp1
enerdots(1) = enerdots1
sigmast(1)   =  sigmastmp

        



!*******************************************************************************
!*******************************************************************************
!*******************************************************************************

! Stuff for Runge Kutta integration

dt       =  dttry               ! initial value of timestep
dtmin    =  1.0d-5 * yr         ! minimal value of timestep
dtmax    =  1.d4 * yr

if (tideop .eq. 0) then

dtmin = 0.5e4 * yr
endif

write(6,*) "Tideop =", tideop

eps      =  1d-3                !1.0d-2              ! 4d-3!

safety   =  0.9
pgrow    =  -0.2
pshrink  =  -0.25

errcon   =  (5/safety)**(1/pgrow)


! Interpolate Radius

call splint(toto,radius,d2RSdt,indicend-(indice0)+1,t_init/yr,Rs0)
Rs0 =Rs0 * Rsun
Rs = Rs0

! Interpolate of stellar quantities for rotevol

call splint(toto,Lum,d2Ldt,indicend-(indice0)+1,t_init/yr,Ls0)
call splint(toto,Teff,d2Teffdt,indicend-(indice0)+1,t_init/yr,Teffs0)
call splint(toto,r2conv,d2r2convdt,indicend-(indice0)+1,t_init/yr,r2convs0)
call splint(toto,r2rad,d2r2raddt,indicend-(indice0)+1,t_init/yr,r2rads0)
call splint(toto,Mrad,d2Mraddt,indicend-(indice0)+1,t_init/yr,Mrads0)
call splint(toto,Rrad,d2Rraddt,indicend-(indice0)+1,t_init/yr,Rrads0)



if (ori .eq. 1) rg2s0 = rg2s
if (ori .eq. 0) then
    rg2s0 = r2convs0
    rg2s  = r2convs0
endif

if (ori .eq. 1) then
    ! Call the routine that calculates 1/Omega*dOmega/dt (rg2s = cst, Rs = cst)
    ! due to tides
    call frots_dot(ap0,ep0,obls0,wconvi,G,Mp,Ms,Rs0,sigmas,rg2s0,ks)
else
    ! Call the routine that calculates dJ/dt due to tides
    call frots_dot2(ap0,ep0,obls0,wconvi,G,Mp,Ms,Rs0,sigmas,rg2s0,ks)
endif    
ks  =  dt * ks

call adot(ap0,ep0,oblp0,obls0,rotp0,wconvi &
        ,G,Mp,Ms,Rp,Rs0,sigmap,sigmas,ka)        
ka = dt * ka


!if (ori .eq. 1) then
    !! Normally we don't use that if we use frots_dot2
    !call frots_dot(ap0,ep0,obls0,wconvi,G,Mp,Ms,Rs0,sigmas,rg2s,ks)
!endif

! Different wind parametrization for slow vs fast rotating
if ( wconvi .le. spin_sat)  alphac = 3.d0
if ( wconvi .gt. spin_sat)  alphac = 1.d0


if (ori .eq. 1) then
    call fwind(wconvi,spin_sat,Kwind,Ms,Msun,Rs0,Rsun,rg2s,alphac,kw)
    kw  =  dt * kw
endif



!****************************** WRITE FILES ***********************************

close(1)
close(2)
!openw,1,filename1


open(unit=1,file=filename1,status='new')

write(1,*)'Information file on the simulation'
write(1,*)''
write(1,*)''

write(1,*)'Bodies properties :'
write(1,*)''
write(1,*)'Mp      = ',Mp/Mearth,'   M_earth'
write(1,*)'Rp      = ',Rp/Rearth,'   R_earth'
write(1,*)'sigmap  = ',sigmap,'   kg**-1.m**-2.s**-1'
write(1,*)'rg2p    = ',rg2p
write(1,*)''
write(1,*)'Ms      = ',Ms/Msun,'   M_sun'
write(1,*)'Rs      = ',Rs0/Rsun,'   R_sun','  = ',Rs0/AU,' AU' 
write(1,*)'sigmas  = ',sigmas,'   kg**-1.m**-2.s**-1'
write(1,*)'rg2s    = ',rg2s
write(1,*)''
write(1,*)'a_roche = ',a_roche/AU,' AU' 
write(1,*)''

call alpha(ep0,alphaout)
call norb(G,Mp,Ms,norbout)

write(1,*)'Initial conditions : '
write(1,*)''
write(1,*)'t_init        = ',t_init/yr,' yr'
write(1,*)'a_initial     = ',ap0/AU,' AU'
write(1,*)'a(da/dt=0)    = ',aequ(1)/AU,' AU'
write(1,*)'e_initial     = ',ep0
write(1,*)'oblp_initial  = ',oblp0*180.d0/Pi,' deg'
write(1,*)'obls_initial  = ',obls0*180.d0/Pi,' deg'
write(1,*)'rotp_initial  = ',rotp0,' s-1'
write(1,*)'pseudo_rot    = ',alphaout*norbout*ap0**(-1.5),' s-1'
write(1,*)'Tp_init       = ',2*Pi/(rotp0*hr),' hr'
write(1,*)'rots_initial  = ',wconvi,' s-1'
write(1,*)'Ts_init       = ',2*Pi/(wconvi*hr),' hr'
write(1,*)''
write(1,*)''

write(1,*)'final conditions : '
write(1,*)''
close(1)

open(unit=2,file=filename2,status='new')
write(2,*) ' ','t (yrs)',' ','a (AU)',' ','e',' ','oblp (deg)',' ','obls (deg)',' ' &
    ,'rotp (s-1)',' ','rots (s-1)',' ','Rp (Rearth)',' ','Rs (AU)',' ' &
    ,'sigmast',' ','1/tau_s (s-1)',' ','tau_s (yr)',' ' &
    ,'1/tau_a (s-1)',' ','tau_a (yr)',' ','d(rots)/dt (s-2)',' ','dt (yr)' 
    

    
!write(2,111) (ttmp)/yr,ap0/AU,ep0,oblp0*180.d0/Pi,obls0*180.d0/Pi &
!    ,rotp0,wconvi,Rp/Rearth,Rs0/AU,sigmas &
!    ,(ks+kw)/(rotstmp*dt),rotstmp*dt/((ks+kw)*yr) &
!    ,ka/(atmp*dt),(atmp*dt)/(ka*yr) &
!    ,(ks+kw)/dt ,dt/yr 
    
    
 111  format(10(1x,ES15.7),5(1x,ES16.7),1x(ES15.7))

close(2)


!stop

debug = 10
j  =  0

write(6,*)'j =',j
write(6,*)'t = ',(ttmp-t_init)/yr,'yrs'
write(6,*)'dt =',dt/yr

!*******************************************************************************
! Begin integration

AMtotinit = rg2s0 * mstar * msun*1.e3 * (Rs0*1.e2)**2*wconvi

ttmp = t_init


do while (ttmp .le. Time)


	if (ori .eq. 0) then
		frotstmp = 0.0
	endif 

	!classical calculation
	
	
    if (ori .eq. 1) then
        if ( rotstmp .le. spin_sat) then 
            alphac = 3.d0
            if ( j .eq. 0) write(6,*)'w<wsat, alpha =',alphac
        endif
        if ( rotstmp .gt. spin_sat) then 
            alphac = 1.d0
            if ( j .eq. 0) write(6,*)'w>wsat, alpha =',alphac
        endif
    endif

    if ( mod(j,100000) .eq. 0)  write(6,*)'t = ',(ttmp-t_init)/yr   !100000

    n_orb = sqrt(G*(Ms+Mp)/(atmp*atmp*atmp))

    ! taking into account the change in stellar radius with time
    ! Interpolation radius
    
    
    do int = 1, 6 
    
	
        call splint(toto,radius,d2RSdt,indicend-(indice0)+1,(ttmp+aa(int)*dt)/yr,Rs)
        Rs = Rs * Rsun
        Rs_tabl(int) = Rs
        

        if (ori .eq. 0) then
        
        	
            !Change of stellar parameters
            call splint(toto,Lum,d2Ldt,indicend-(indice0)+1,(ttmp+aa(int)*dt)/yr,Ls)
            Ls_tabl(int) = Ls
            
            call splint(toto,Teff,d2Teffdt,indicend-(indice0)+1,(ttmp+aa(int)*dt)/yr,Teffs)
            Teffs_tabl(int) = Teffs

            call splint(toto,r2conv,d2r2convdt,indicend-(indice0)+1,(ttmp+aa(int)*dt)/yr,r2convs)
            r2convs_tabl(int) = r2convs

            call splint(toto,r2rad,d2r2raddt,indicend-(indice0)+1,(ttmp+aa(int)*dt)/yr,r2rads)
            r2rads_tabl(int) = r2rads

            call splint(toto,Mrad,d2Mraddt,indicend-(indice0)+1,(ttmp+aa(int)*dt)/yr,Mrads)
            Mrads_tabl(int) = Mrads

            call splint(toto,Rrad,d2Rraddt,indicend-(indice0)+1,(ttmp+aa(int)*dt)/yr,Rrads)
            Rrads_tabl(int) = Rrads
        endif
        
                

        if ( cst_diss .eq. 0) then 
            if ( n_orb .le. 2.d0*rotstmp) then 

                ! Interpolation of the dissipation Q
                call splint(toto,Qs_inv,d2QSdt,indicend-(indice0),(ttmp+aa(int)*dt)/yr,Q_str_inv)

                epsilon_squared = rotstmp*rotstmp/Sun_dyn_freq
                lag_angle = 3.d0*epsilon_squared*Q_str_inv/4.d0
                ! 2*lag_angle = freq * Delta t, freq = 2*(n-Omega)
                delta_t = lag_angle/max(1.d-5,abs(n_orb-rotstmp))
                ! Dynamical dissipation + equilibrium dissipation
                sigmas_tabl(int) = 2.d0*G* delta_t /(3.d0*Rs*Rs*Rs*Rs*Rs) + sigmas_eq
            endif
            if ( n_orb .gt. 2.d0*rotstmp) sigmas_tabl(int) = sigmas_eq
        endif
        if ( cst_diss .eq. 1) sigmas_tabl(int) = sigmas_eq

    enddo
    

			if (Ls_tabl(5) .lt. 0.0) then

				write(6,*) "Ls_tabl(5) lt 0"
				write(6,*) tstar2,aa(5)*dt/yr
				write(6,*) Teffs0,Rs0,r2convs0,r2rads0,Mrads0,Rrads0
				do i=1,indicend-(indice0)+1
				 write(6,*) toto(i), Lum(i) , d2Ldt(i)
				enddo
				stop
			endif    
    
    
    ! Runge Kutta stuff necessary to compute newt timestep 
    ! If Ori = 1 need to use rotstmp not wconvi as in that case wconvi = cst = initial rotation
    ! If Ori = 0, wconvi evolves. 1st time rotstmp = wconvi, then wconvi updated at the end of loop and is
    ! equal to rotstmp. So before the call of rotevol, rotstmp (step i-1) is equal to the new wconvi (step i) value. 
    
    call adot(atmp,etmp,oblptmp,oblstmp,rotptmp,rotstmp &
        ,G,Mp,Ms,Rp,Rs_tabl(1),sigmap,sigmas_tabl(1),tmp)
    yscal(1) = tmp

    call edot(atmp,etmp,oblptmp,oblstmp,rotptmp,rotstmp &
        ,G,Mp,Ms,Rp,Rs_tabl(1),sigmap,sigmas_tabl(1),tmp)
    yscal(2) = tmp

    call obldot(atmp,etmp,oblptmp,rotptmp &
        ,G,Mp,Ms,Rp,sigmap,rg2p,tmp)
    yscal(3) = tmp

    call obldot(atmp,etmp,oblstmp,rotstmp &
        ,G,Ms,Mp,Rs_tabl(1),sigmas_tabl(1),rg2s,tmp)
    yscal(4) = tmp

    call rotp_dot(atmp,etmp,oblptmp,rotptmp &
        ,G,Mp,Ms,Rp,sigmap,rg2p,tmp)
    yscal(5) = tmp

    if (ori .eq. 1) then
        call frots_dot(atmp,etmp,oblstmp,rotstmp &
            ,G,Mp,Ms,Rs_tabl(1),sigmas_tabl(1),rg2s,tmp)
        yscal(6) = tmp

        call fwind(rotstmp,spin_sat,Kwind,Ms,Msun &
            ,Rs_tabl(1),Rsun,rg2s,alphac,tmp)
        yscal(7) = tmp
    else
        call frots_dot2(atmp,etmp,oblstmp,rotstmp &
            ,G,Mp,Ms,Rs_tabl(1),sigmas_tabl(1),rg2s,tmp)
        yscal(6) = tmp
        yscal(7) = 0.
    endif

    do i=1,7
        yscal(i) = yscal(i) * eps * dt
    enddo
        
    
    !*******************************************************************************
    ! Runge Kutta factors 1st
    !*******************************************************************************
    
    call adot(atmp,etmp,oblptmp,oblstmp,rotptmp,rotstmp &
        ,G,Mp,Ms,Rp,Rs_tabl(1),sigmap,sigmas_tabl(1),ka0)
    ka0 = ka0 * dt

    if ( etmp .gt. 0.0d0) then 
        call edot(atmp,etmp,oblptmp,oblstmp,rotptmp,rotstmp &
            ,G,Mp,Ms,Rp,Rs_tabl(1),sigmap,sigmas_tabl(1),ke0)
        ke0   = ke0 * dt 
    else 
        ke0   =  0.0d0
    endif

    if ( oblptmp .gt. 0.0d0) then 
        call obldot(atmp,etmp,oblptmp,rotptmp,G,Mp,Ms,Rp,sigmap,rg2p,kop0)
        kop0 = kop0 * dt
    else 
        kop0  =  0.0d0
    endif

    if ( oblstmp .gt. 0.0d0) then 
        call obldot(atmp,etmp,oblstmp,rotstmp,G,Ms,Mp,Rs_tabl(1),sigmas_tabl(1),rg2s,kos0)
        kos0  = kos0 * dt
    else
        kos0  =  0.0d0
    endif 


    if (ori .eq. 1) then
        call frots_dot(atmp,etmp,oblstmp,rotstmp,G,Mp,Ms,Rs_tabl(1),sigmas_tabl(1),rg2s,ks0)
        ks0  =  ks0 * dt
        call fwind(rotstmp,spin_sat,Kwind,Ms,Msun,Rs_tabl(1),Rsun,rg2s,alphac,kw0)
        kw0  =  kw0 * dt 
    else
        call frots_dot2(atmp,etmp,oblstmp,rotstmp,G,Mp,Ms,Rs_tabl(1),sigmas_tabl(1),rg2s,ks0)
        ks0  =  ks0 * dt
    endif

    !write(6,*) "kw0", rotstmp,rotstmp,spin_sat,Kwind,Ms,Msun,Rs_tabl(1),Rsun,rg2s,alphac,kw0
    !write(6,*) "ks0", atmp,etmp,oblstmp,rotstmp,G,Mp,Ms,Rs_tabl(1),sigmas_tabl(1),rg2s,ks0
    
    
    if ( babaorum .eq. 1)  kp0  =  0.0
    if ( babaorum .eq. 2) then 
        call rotp_dot(atmp,etmp,oblptmp,rotptmp,G,Mp,Ms,Rp,sigmap,rg2p,kp0)
        kp0  =  kp0 * dt
    endif   
    
    if (flagcrash .eq. 1 ) then
    	ka0 = 0.0
    	ke0 = 0.0
    	kop0 = 0.0
    	kos0 = 0.0
    	ks0 = 0.0  
    	kp0 = 0.0
    	frotstmp = 0.0  	
    endif

    !*******************************************************************************
    
	!*******************************************************************************
    ! 							HERE call rotevol !!!!!
    !		Rotevol calculates evolution of rots between two timestep
    !*******************************************************************************

    if (ori .eq. 0) then
        Lum2(1) = Ls0
        Lum2(2) = Ls_tabl(5)

        Teff2(1) = Teffs0
        Teff2(2) = Teffs_tabl(5)

        tstar2(1) = ttmp/yr
        tstar2(2) = (ttmp+aa(5)*dt)/yr

        Rstar2(1) = Rs0/Rsun
        Rstar2(2) = Rs_tabl(5)/Rsun

        r2conv2(1) = r2convs0
        r2conv2(2) = r2convs_tabl(5)

        r2rad2(1) = r2rads0
        r2rad2(2) = r2rads_tabl(5)

        Mrad2(1) = Mrads0
        Mrad2(2) = Mrads_tabl(5)

        Rrad2(1) = Rrads0
        Rrad2(2) = Rrads_tabl(5)

        !call rotevol with initial velocity wconvi and wradi
        !calculates new rotation rate rotstmp and wradf
        
        call rotevol(taudecinit,tdiskparam,wconvi,rotstmp,wradi,wradf,Lum2,Teff2,tstar2,Rstar2,r2conv2,r2rad2,Mrad2,Rrad2 &
            ,Wuprad,juprad,deltaJ,timedec,DAMbrake, DAMtide,Wtide,Wupconv,Wdownconv,Wjuprad,deltaWc)
            

			if (Ls0 .lt. 0.0) then

				write(6,*) "Ls0 lt 0"
				write(6,*) tstar2
				write(6,*) indicend,indice0,indicend-(indice0)+1,Teffs0,Rs0,r2convs0,r2rads0,Mrads0,Rrads0

				stop

			endif
            
        if (rotstmp .ne. rotstmp) then
        	write(6,*) "NaN? ", rotstmp
        	write(6,*) "Control donnee"
        	write(6,*) Lum2,Teff2,tstar2,Rstar2,r2conv2,r2rad2,Mrad2,Rrad2 &
            ,Wuprad,juprad,deltaJ,timedec
        
        	stop
        endif

        do int = 1, 6 
            ! Linear inertpolation of spin between input and output of rotevol
            ! (wconvi and rotstmp)
            rotstmp_tabl(int) = wconvi + (rotstmp-wconvi)*((ttmp+aa(int)*dt)/yr-tstar2(1))/(tstar2(2)-tstar2(1))
        enddo
    endif

    sigmas = sigmas_tabl(1)

    if (ori .eq. 1) then
        ! Here rg2s = cst
        rg2s  = rg2s
        rg2s1 = rg2s
        rg2s2 = rg2s
        rg2s3 = rg2s
        rg2s4 = rg2s
        rg2s5 = rg2s
        rg2sdot = 0.0
    else
        rg2s  = r2convs_tabl(1)
        rg2s1 = r2convs_tabl(2)
        rg2s2 = r2convs_tabl(3)
        rg2s3 = r2convs_tabl(4)
        rg2s4 = r2convs_tabl(5)
        rg2s5 = r2convs_tabl(6)
        rg2sdot = 0.0
    endif



    !*******************************************************************************
    ! Runge Kutta factors bis
    !*******************************************************************************

    
	!*******************************************************************************


    if (ori .eq. 1) then
        rotstmp1 = wconvi*rg2s0/rg2s1*(Rs0/Rs_tabl(2))**2 &
            *exp(frotstmp+bb(1,2)*ks0)*exp(fwindtmp+ bb(1,2)*kw0)
    else
        !interpolation + integration?
        Iconvtmp = r2convs_tabl(2) * mstar * msun * Rs_tabl(2)**2
        if (tideop .eq. 1) then
        	rotstmp1 = rotstmp_tabl(2) + (frotstmp + bb(1,2)*ks0)/Iconvtmp
        else 
        	rotstmp1 = rotstmp1
        endif		
    endif

    call adot(atmp + bb(1,2)*ka0 &
        ,etmp    + bb(1,2)*ke0 &
        ,oblptmp + bb(1,2)*kop0 &
        ,oblstmp + bb(1,2)*kos0 &
        ,rotptmp + bb(1,2)*kp0 &
        ,rotstmp1 &
        ,G,Mp,Ms,Rp,Rs_tabl(2),sigmap,sigmas_tabl(2),ka1)

    ka1 = ka1 * dt    

    if ( etmp .gt. 0.0d0) then
        call edot(atmp + bb(1,2)*ka0 &
            ,etmp    + bb(1,2)*ke0 &
            ,oblptmp + bb(1,2)*kop0 &
            ,oblstmp + bb(1,2)*kos0 &
            ,rotptmp + bb(1,2)*kp0 &
            ,rotstmp1 &
            ,G,Mp,Ms,Rp,Rs_tabl(2),sigmap,sigmas_tabl(2),ke1)
        ke1  =  dt * ke1
    else
        ke1  =  0.0d0
    endif 

    if ( oblptmp .gt. 0.0d0) then
        call obldot(atmp + bb(1,2)*ka0 &
            ,etmp    + bb(1,2)*ke0 &
            ,oblptmp + bb(1,2)*kop0 &
            ,rotptmp + bb(1,2)*kp0 &
            ,G,Mp,Ms,Rp,sigmap,rg2p,kop1) 
        kop1 = dt * kop1
    else
        kop1  =  0.0d0
    endif 

    if ( oblstmp .gt. 0.0d0) then 
        call obldot(atmp+ bb(1,2)*ka0 &
            ,etmp    + bb(1,2)*ke0 &
            ,oblstmp + bb(1,2)*kos0 &
            ,rotstmp1 &
            ,G,Ms,Mp,Rs_tabl(2),sigmas_tabl(2),rg2s1,kos1)
        kos1 = dt * kos1
    else
        kos1  =  0.0d0
    endif

    if (ori .eq. 1) then
        call frots_dot(atmp + bb(1,2)*ka0 &
            ,etmp    + bb(1,2)*ke0 &
            ,oblstmp + bb(1,2)*kos0 &
            ,rotstmp1 &
            ,G,Mp,Ms,Rs_tabl(2),sigmas_tabl(2),rg2s1,ks1) 
        ks1  =  dt *   ks1

        call fwind(rotstmp1 &
            ,spin_sat,Kwind,Ms,Msun &
            ,Rs_tabl(2),Rsun,rg2s1,alphac,kw1)
        kw1  =  dt * kw1
    else
        call frots_dot2(atmp + bb(1,2)*ka0 &
            ,etmp    + bb(1,2)*ke0 &
            ,oblstmp + bb(1,2)*kos0 &
            ,rotstmp1 &
            ,G,Mp,Ms,Rs_tabl(2),sigmas_tabl(2),rg2s1,ks1) 
        ks1  =  dt *   ks1
    endif


    if ( babaorum .eq. 1) kp1  =  0.0
    if ( babaorum .eq. 2) then 
        call rotp_dot(atmp + bb(1,2)*ka0 &
            ,oblptmp + bb(1,2)*kop0 &
            ,etmp + bb(1,2)*ke0 &
            ,rotptmp + bb(1,2)*kp0 &
            ,G,Mp,Ms,Rp,sigmap,rg2p,kp1)
        kp1  =  dt * kp1
    endif   
    
    
    if (flagcrash .eq. 1 .or. tideop .eq. 0 ) then
    	ka1 = 0.0
    	ke1 = 0.0
    	kop1 = 0.0
    	kos1 = 0.0
    	ks1 = 0.0  
    	kp1 = 0.0
    	frotstmp = 0.0  	
    endif 

    !*******************************************************************************

    if (ori .eq. 1) then
        rotstmp2 = wconvi*rg2s0/rg2s2*(Rs0/Rs_tabl(3))**2 &
            *exp(frotstmp + bb(1,3)*ks0 + bb(2,3)*ks1) &
            *exp(fwindtmp + bb(1,3)*kw0 + bb(2,3)*kw1)
    else
        !interpolation + integration?
        Iconvtmp = r2convs_tabl(3) * mstar * msun * Rs_tabl(3)**2
        if (tideop .eq. 1) then
        	rotstmp2 = rotstmp_tabl(3) &
            	+ (frotstmp + bb(1,3)*ks0 + bb(2,3)*ks1)/Iconvtmp
        else
        	rotstmp2 = rotstmp2
        endif    	
    endif

    call adot(atmp + bb(1,3)*ka0 + bb(2,3)*ka1 &
        ,etmp    + bb(1,3)*ke0  + bb(2,3)*ke1 &
        ,oblptmp + bb(1,3)*kop0 + bb(2,3)*kop1 &
        ,oblstmp + bb(1,3)*kos0 + bb(2,3)*kos1 &
        ,rotptmp + bb(1,3)*kp0  + bb(2,3)*kp1 &
        ,rotstmp2 &
        ,G,Mp,Ms,Rp,Rs_tabl(3),sigmap,sigmas_tabl(3),ka2)
    ka2  =  dt * ka2
    !write(6,*)'detail ka2'
    !write(6,*)atmp + bb(0,2)*ka0 + bb(1,2)*ka1
    !write(6,*)rotptmp + bb(0,2)*kp0  + bb(1,2)*kp1
    !write(6,*)rotstmp + bb(0,2)*ks0  + bb(1,2)*ks1
    !write(6,*)G,Mp,Ms,Rp,Rs_tabl(2),sigmap,sigmas_tabl(2)
    !write(6,*)'end detail ka2'

    if ( etmp .gt. 0.0d0) then 
        call edot(atmp + bb(1,3)*ka0 + bb(2,3)*ka1 &
            ,etmp    + bb(1,3)*ke0  + bb(2,3)*ke1 &
            ,oblptmp + bb(1,3)*kop0 + bb(2,3)*kop1 &
            ,oblstmp + bb(1,3)*kos0 + bb(2,3)*kos1 &
            ,rotptmp + bb(1,3)*kp0  + bb(2,3)*kp1 &
            ,rotstmp2 &
            ,G,Mp,Ms,Rp,Rs_tabl(3),sigmap,sigmas_tabl(3),ke2)
        ke2  =  dt * ke2
    else
        ke2  =  0.0d0
    endif 

    if ( oblptmp .gt. 0.0d0) then 
        call obldot(atmp + bb(1,3)*ka0 + bb(2,3)*ka1 &
            ,etmp    + bb(1,3)*ke0  + bb(2,3)*ke1 &
            ,oblptmp + bb(1,3)*kop0 + bb(2,3)*kop1 &
            ,rotptmp + bb(1,3)*kp0  + bb(2,3)*kp1 &
            ,G,Mp,Ms,Rp,sigmap,rg2p,kop2)
        kop2 =  dt * kop2
    else
        kop2  =  0.0d0
    endif 

    if ( oblstmp .gt. 0.0d0) then
        call obldot(atmp + bb(1,3)*ka0 + bb(2,3)*ka1 &
            ,etmp    + bb(1,3)*ke0  + bb(2,3)*ke1 &
            ,oblstmp + bb(1,3)*kos0 + bb(2,3)*kos1 &
            ,rotstmp2 &
            ,G,Ms,Mp,Rs_tabl(3),sigmas_tabl(3),rg2s2,kos2) 
        kos2 =  dt * kos2
    else
        kos2  =  0.0d0
    endif 

    if (ori .eq. 1) then
        call frots_dot(atmp + bb(1,3)*ka0 + bb(2,3)*ka1 &
            ,etmp    + bb(1,3)*ke0  + bb(2,3)*ke1 &
            ,oblstmp + bb(1,3)*kos0 + bb(2,3)*kos1 &
            ,rotstmp2 &
            ,G,Mp,Ms,Rs_tabl(3),sigmas_tabl(3),rg2s2,ks2)
        ks2  =  dt * ks2

        call fwind(rotstmp2 &
            ,spin_sat,Kwind,Ms,Msun &
            ,Rs_tabl(3),Rsun,rg2s2,alphac,kw2)
        kw2  =  dt * kw2
    else
        call frots_dot2(atmp + bb(1,3)*ka0 + bb(2,3)*ka1 &
            ,etmp    + bb(1,3)*ke0  + bb(2,3)*ke1 &
            ,oblstmp + bb(1,3)*kos0 + bb(2,3)*kos1 &
            ,rotstmp2 &
            ,G,Mp,Ms,Rs_tabl(3),sigmas_tabl(3),rg2s2,ks2)
        ks2  =  dt * ks2
    endif

    if ( babaorum .eq. 1) kp2  =  0.0
    if ( babaorum .eq. 2) then 
        call rotp_dot(atmp + bb(1,3)*ka0 + bb(2,3)*ka1 &
            ,etmp    + bb(1,3)*ke0  + bb(2,3)*ke1 &
            ,oblptmp + bb(1,3)*kop0 + bb(2,3)*kop1 &
            ,rotptmp + bb(1,3)*kp0  + bb(2,3)*kp1 &
            ,G,Mp,Ms,Rp,sigmap,rg2p,kp2)
        kp2  =  dt * kp2
    endif
    
    if (flagcrash .eq. 1 .or. tideop .eq. 0 ) then
    	ka2 = 0.0
    	ke2 = 0.0
    	kop2 = 0.0
    	kos2 = 0.0
    	ks2 = 0.0  
    	kp2 = 0.0
    	frotstmp = 0.0  	
    endif

    !*******************************************************************************

    if (ori .eq. 1) then
        rotstmp3 = wconvi*rg2s0/rg2s3*(Rs0/Rs_tabl(4))**2 &
            *exp(frotstmp + bb(1,4)*ks0 + bb(2,4)*ks1 + bb(3,4)*ks2) &
            *exp(fwindtmp + bb(1,4)*kw0 + bb(2,4)*kw1 + bb(3,4)*kw2)
    else
        Iconvtmp = r2convs_tabl(4) * mstar * msun * Rs_tabl(4)**2
        if (tideop .eq. 1) then
        	rotstmp3 = rotstmp_tabl(4) &
            	+ (frotstmp + bb(1,4)*ks0 + bb(2,4)*ks1 + bb(3,4)*ks2)/Iconvtmp
        else
        	rotstmp3 = rotstmp3
        endif    	
    endif

    call adot(atmp + bb(1,4)*ka0 + bb(2,4)*ka1 + bb(3,4)*ka2 &
        ,etmp    + bb(1,4)*ke0  + bb(2,4)*ke1  + bb(3,4)*ke2 &
        ,oblptmp + bb(1,4)*kop0 + bb(2,4)*kop1 + bb(3,4)*kop2 &
        ,oblstmp + bb(1,4)*kos0 + bb(2,4)*kos1 + bb(3,4)*kos2 &
        ,rotptmp + bb(1,4)*kp0  + bb(2,4)*kp1  + bb(3,4)*kp2 &
        ,rotstmp3 &
        ,G,Mp,Ms,Rp,Rs_tabl(4),sigmap,sigmas_tabl(4),ka3)
    ka3   =  dt * ka3
    
    !write(6,*)'detail ka3'
    !write(6,*)atmp + bb(0,3)*ka0 + bb(2,4)*ka1 + bb(3,4)*ka2
    !write(6,*)rotptmp + bb(0,3)*kp0  + bb(2,4)*kp1  + bb(3,4)*kp2
    !write(6,*)rotstmp + bb(0,3)*ks0  + bb(2,4)*ks1  + bb(3,4)*ks2
    !write(6,*)G,Mp,Ms,Rp,Rs_tabl(4),sigmap,sigmas_tabl(4)
    !write(6,*)'end detail ka3'

    if ( etmp .gt. 0.0d0) then 
        call edot(atmp + bb(1,4)*ka0 + bb(2,4)*ka1 + bb(3,4)*ka2 &
            ,etmp    + bb(1,4)*ke0  + bb(2,4)*ke1  + bb(3,4)*ke2 &
            ,oblptmp + bb(1,4)*kop0 + bb(2,4)*kop1 + bb(3,4)*kop2 &
            ,oblstmp + bb(1,4)*kos0 + bb(2,4)*kos1 + bb(3,4)*kos2 &
            ,rotptmp + bb(1,4)*kp0  + bb(2,4)*kp1  + bb(3,4)*kp2 &
            ,rotstmp3 &
            ,G,Mp,Ms,Rp,Rs_tabl(4),sigmap,sigmas_tabl(4),ke3)
        ke3   =  dt * ke3
    else
        ke3  =  0.0d0
    endif 

    if ( oblptmp .gt. 0.0d0) then
        call obldot(atmp + bb(1,4)*ka0 + bb(2,4)*ka1 + bb(3,4)*ka2 &
            ,etmp    + bb(1,4)*ke0  + bb(2,4)*ke1  + bb(3,4)*ke2 &
            ,oblptmp + bb(1,4)*kop0 + bb(2,4)*kop1 + bb(3,4)*kop2 &
            ,rotptmp + bb(1,4)*kp0  + bb(2,4)*kp1  + bb(3,4)*kp2 &
            ,G,Mp,Ms,Rp,sigmap,rg2p,kop3) 
        kop3  =  dt * kop3
    else
        kop3  =  0.0d0
    endif 

    if ( oblstmp .gt. 0.0d0) then 
        call obldot(atmp + bb(1,4)*ka0 + bb(2,4)*ka1 + bb(3,4)*ka2 &
            ,etmp    + bb(1,4)*ke0  + bb(2,4)*ke1  + bb(3,4)*ke2 &
            ,oblstmp + bb(1,4)*kos0 + bb(2,4)*kos1 + bb(3,4)*kos2 &
            ,rotstmp3 &
            ,G,Ms,Mp,Rs_tabl(4),sigmas_tabl(4),rg2s3,kos3)
        kos3  =  dt * kos3
    else
        kos3  =  0.0d0
    endif 

    if (ori .eq. 1) then
        call frots_dot(atmp + bb(1,4)*ka0 + bb(2,4)*ka1 + bb(3,4)*ka2 &
            ,etmp    + bb(1,4)*ke0  + bb(2,4)*ke1  + bb(3,4)*ke2 &
            ,oblstmp + bb(1,4)*kos0 + bb(2,4)*kos1 + bb(3,4)*kos2 &
            ,rotstmp3 &
            ,G,Mp,Ms,Rs_tabl(4),sigmas_tabl(4),rg2s3,ks3)
        ks3   =  dt * ks3

        call fwind(rotstmp3 &
            ,spin_sat,Kwind,Ms,Msun &
            ,Rs_tabl(4),Rsun,rg2s3,alphac,kw3)
        kw3  =  dt * kw3
    else
        call frots_dot2(atmp + bb(1,4)*ka0 + bb(2,4)*ka1 + bb(3,4)*ka2 &
            ,etmp    + bb(1,4)*ke0  + bb(2,4)*ke1  + bb(3,4)*ke2 &
            ,oblstmp + bb(1,4)*kos0 + bb(2,4)*kos1 + bb(3,4)*kos2 &
            ,rotstmp3 &
            ,G,Mp,Ms,Rs_tabl(4),sigmas_tabl(4),rg2s3,ks3)
        ks3   =  dt * ks3
    endif

    !write(6,*)'ka3  = ',ka3
    !write(6,*)'kos3 = ',kos3
    !write(6,*)'ks3 = ',ks3
    !write(6,*)'kw3 = ',kw3

    if ( babaorum .eq. 1) kp3  =  0.0
    if ( babaorum .eq. 2) then 
        call rotp_dot(atmp + bb(1,4)*ka0 + bb(2,4)*ka1 + bb(3,4)*ka2 &
            ,etmp    + bb(1,4)*ke0  + bb(2,4)*ke1  + bb(3,4)*ke2 &
            ,oblptmp + bb(1,4)*kop0 + bb(2,4)*kop1 + bb(3,4)*kop2 &
            ,rotptmp + bb(1,4)*kp0  + bb(2,4)*kp1  + bb(3,4)*kp2 &
            ,G,Mp,Ms,Rp,sigmap,rg2p,kp3)
        kp3  =  dt *  kp3
    endif
    
    if (flagcrash .eq. 1 .or. tideop .eq. 0 ) then
    	ka3 = 0.0
    	ke3 = 0.0
    	kop3 = 0.0
    	kos3 = 0.0
    	ks3 = 0.0  
    	kp3 = 0.0
    	frotstmp = 0.0  	
    endif

    !*******************************************************************************

    if (ori .eq. 1) then
        rotstmp4 = wconvi*rg2s0/rg2s4*(Rs0/Rs_tabl(5))**2 &
            *exp(frotstmp + bb(1,5)*ks0 + bb(2,5)*ks1 + bb(3,5)*ks2 + bb(4,5)*ks3) &
            *exp(fwindtmp + bb(1,5)*kw0 + bb(2,5)*kw1 + bb(3,5)*kw2 + bb(4,5)*kw3)
    else
        Iconvtmp = r2convs_tabl(5) * mstar * msun * Rs_tabl(5)**2
        if (tideop .eq. 1) then
        	rotstmp4 = rotstmp_tabl(5) &
            	+ (frotstmp + bb(1,5)*ks0 + bb(2,5)*ks1 + bb(3,5)*ks2 + bb(4,5)*ks3)/Iconvtmp
        else
        	rotstmp4 = rotstmp4        
        endif
    endif

    call adot(atmp   + bb(1,5)*ka0  + bb(2,5)*ka1  + bb(3,5)*ka2  + bb(4,5)*ka3 &
            ,etmp    + bb(1,5)*ke0  + bb(2,5)*ke1  + bb(3,5)*ke2  + bb(4,5)*ke3  &
            ,oblptmp + bb(1,5)*kop0 + bb(2,5)*kop1 + bb(3,5)*kop2 + bb(4,5)*kop3 &
            ,oblstmp + bb(1,5)*kos0 + bb(2,5)*kos1 + bb(3,5)*kos2 + bb(4,5)*kos3 &
            ,rotptmp + bb(1,5)*kp0  + bb(2,5)*kp1  + bb(3,5)*kp2  + bb(4,5)*kp3 &
            ,rotstmp4 &
            ,G,Mp,Ms,Rp,Rs_tabl(5),sigmap,sigmas_tabl(5),ka4)
    ka4  =  dt * ka4
        
    !write(6,*)'detail ka4'
    !write(6,*)- 1.d0/Temp(Mp,Ms,Rp,sigmap)*atmp**(-7.d0) 
    !write(6,*)(Na1(etmp)-cos(oblptmp)*(rotptmp/norb(G,Mp,Ms))*atmp**(1.5d0)*Na2(etmp))
    !write(6,*)- 1.d0/Temp(Ms,Mp,Rs_tabl(5),sigmas_tabl(5))*atmp**(-7.d0)
    !write(6,*)(Na1(etmp)-cos(oblstmp)*(rotstmp/norb(G,Mp,Ms))*atmp**(1.5d0)*Na2(etmp))
    !write(6,*)atmp    + bb(1,5)*ka0  + bb(2,5)*ka1  + bb(3,5)*ka2  + bb(4,5)*ka3
    !write(6,*)etmp    + bb(1,5)*ke0  + bb(2,5)*ke1  + bb(3,5)*ke2  + bb(4,5)*ke3
    !write(6,*)oblptmp + bb(1,5)*kop0 + bb(2,5)*kop1 + bb(3,5)*kop2 + bb(4,5)*kop3
    !write(6,*)oblstmp + bb(1,5)*kos0 + bb(2,5)*kos1 + bb(3,5)*kos2 + bb(4,5)*kos3
    !write(6,*)rotptmp + bb(1,5)*kp0  + bb(2,5)*kp1  + bb(3,5)*kp2  + bb(4,5)*kp3
    !write(6,*)rotstmp + bb(1,5)*ks0  + bb(2,5)*ks1  + bb(3,5)*ks2  + bb(4,5)*ks3
    !write(6,*)G,Mp,Ms,Rp,Rs_tabl(5),sigmap,sigmas_tabl(5)
    !write(6,*)'end detail ka4'

    if ( etmp .gt. 0.0d0) then 
        call edot(atmp    + bb(1,5)*ka0  + bb(2,5)*ka1  + bb(3,5)*ka2  + bb(4,5)*ka3 &
            ,etmp    + bb(1,5)*ke0  + bb(2,5)*ke1  + bb(3,5)*ke2  + bb(4,5)*ke3 &
            ,oblptmp + bb(1,5)*kop0 + bb(2,5)*kop1 + bb(3,5)*kop2 + bb(4,5)*kop3 &
            ,oblstmp + bb(1,5)*kos0 + bb(2,5)*kos1 + bb(3,5)*kos2 + bb(4,5)*kos3 &
            ,rotptmp + bb(1,5)*kp0  + bb(2,5)*kp1  + bb(3,5)*kp2  + bb(4,5)*kp3 &
            ,rotstmp4 &
            ,G,Mp,Ms,Rp,Rs_tabl(5),sigmap,sigmas_tabl(5),ke4)
        ke4  =  dt * ke4
    else
        ke4  =  0.0d0
    endif 

    if ( oblptmp .gt. 0.0d0) then 
        call obldot(atmp    + bb(1,5)*ka0  + bb(2,5)*ka1  + bb(3,5)*ka2  + bb(4,5)*ka3 &
            ,etmp    + bb(1,5)*ke0  + bb(2,5)*ke1  + bb(3,5)*ke2  + bb(4,5)*ke3 &
            ,oblptmp + bb(1,5)*kop0 + bb(2,5)*kop1 + bb(3,5)*kop2 + bb(4,5)*kop3 &
            ,rotptmp + bb(1,5)*kp0  + bb(2,5)*kp1  + bb(3,5)*kp2  + bb(4,5)*kp3 &
            ,G,Mp,Ms,Rp,sigmap,rg2p,kop4)
        kop4  =  dt * kop4
    else
        kop4  =  0.0d0
    endif 

    if ( oblstmp .gt. 0.0d0) then
        call obldot(atmp    + bb(1,5)*ka0  + bb(2,5)*ka1  + bb(3,5)*ka2  + bb(4,5)*ka3 &
            ,etmp    + bb(1,5)*ke0  + bb(2,5)*ke1  + bb(3,5)*ke2  + bb(4,5)*ke3 &
            ,oblstmp + bb(1,5)*kos0 + bb(2,5)*kos1 + bb(3,5)*kos2 + bb(4,5)*kos3 &
            ,rotstmp4 &
            ,G,Ms,Mp,Rs_tabl(5),sigmas_tabl(5),rg2s4,kos4) 
        kos4  =  dt * kos4
    else
        kos4  =  0.0d0
    endif 

    if (ori .eq. 1) then
        call frots_dot(atmp    + bb(1,5)*ka0  + bb(2,5)*ka1  + bb(3,5)*ka2  + bb(4,5)*ka3 &
            ,etmp    + bb(1,5)*ke0  + bb(2,5)*ke1  + bb(3,5)*ke2  + bb(4,5)*ke3 &
            ,oblstmp + bb(1,5)*kos0 + bb(2,5)*kos1 + bb(3,5)*kos2 + bb(4,5)*kos3 &
            ,rotstmp4 &
            ,G,Mp,Ms,Rs_tabl(5),sigmas_tabl(5),rg2s4,ks4)
        ks4  =  dt * ks4

        call fwind(rotstmp4 &
                 ,spin_sat,Kwind,Ms,Msun &
                 ,Rs_tabl(5),Rsun,rg2s4,alphac,kw4)
        kw4  =  dt * kw4
    else
        call frots_dot2(atmp    + bb(1,5)*ka0  + bb(2,5)*ka1  + bb(3,5)*ka2  + bb(4,5)*ka3 &
            ,etmp    + bb(1,5)*ke0  + bb(2,5)*ke1  + bb(3,5)*ke2  + bb(4,5)*ke3 &
            ,oblstmp + bb(1,5)*kos0 + bb(2,5)*kos1 + bb(3,5)*kos2 + bb(4,5)*kos3 &
            ,rotstmp4 &
            ,G,Mp,Ms,Rs_tabl(5),sigmas_tabl(5),rg2s4,ks4)
        ks4  =  dt * ks4
    endif

    if ( babaorum .eq. 1) kp4  =  0.0
    if ( babaorum .eq. 2) then
        call rotp_dot(atmp    + bb(1,5)*ka0  + bb(2,5)*ka1  + bb(3,5)*ka2  + bb(4,5)*ka3 &
            ,etmp    + bb(1,5)*ke0  + bb(2,5)*ke1  + bb(3,5)*ke2  + bb(4,5)*ke3 &
            ,oblptmp + bb(1,5)*kop0 + bb(2,5)*kop1 + bb(3,5)*kop2 + bb(4,5)*kop3 &
            ,rotptmp + bb(1,5)*kp0  + bb(2,5)*kp1  + bb(3,5)*kp2  + bb(4,5)*kp3 &
            ,G,Mp,Ms,Rp,sigmap,rg2p,kp4)  
        kp4  =  dt * kp4
    endif
    
    
    if (flagcrash .eq. 1 .or. tideop .eq. 0 ) then
    	ka4 = 0.0
    	ke4 = 0.0
    	kop4 = 0.0
    	kos4 = 0.0
    	ks4 = 0.0  
    	kp4 = 0.0
    	frotstmp = 0.0  	
    endif

    !*******************************************************************************

    if (ori .eq. 1) then
        rotstmp5 = wconvi*rg2s0/rg2s5*(Rs0/Rs_tabl(6))**2 &
            *exp(frotstmp + bb(1,6)*ks0 + bb(2,6)*ks1 + bb(3,6)*ks2 + bb(4,6)*ks3 + bb(5,6)*ks4) &
            *exp(fwindtmp + bb(1,6)*kw0 + bb(2,6)*kw1 + bb(3,6)*kw2 + bb(4,6)*kw3 + bb(5,6)*kw4)
    else
        Iconvtmp = r2convs_tabl(6) * mstar * msun * Rs_tabl(6)**2
        if (tideop .eq. 1) then
        	rotstmp5 = rotstmp_tabl(6) &
            	+ (frotstmp + bb(1,6)*ks0 + bb(2,6)*ks1 + bb(3,6)*ks2 + bb(4,6)*ks3 + bb(5,6)*ks4)/Iconvtmp
        else
        	rotstmp5 =rotstmp5        
        endif    	
    endif

    call adot(atmp    + bb(1,6)*ka0  + bb(2,6)*ka1  + bb(3,6)*ka2  + bb(4,6)*ka3  + bb(5,6)*ka4 &
        ,etmp    + bb(1,6)*ke0  + bb(2,6)*ke1  + bb(3,6)*ke2  + bb(4,6)*ke3  + bb(5,6)*ke4 &
        ,oblptmp + bb(1,6)*kop0 + bb(2,6)*kop1 + bb(3,6)*kop2 + bb(4,6)*kop3 + bb(5,6)*kop4 &
        ,oblstmp + bb(1,6)*kos0 + bb(2,6)*kos1 + bb(3,6)*kos2 + bb(4,6)*kos3 + bb(5,6)*kos4 &
        ,rotptmp + bb(1,6)*kp0  + bb(2,6)*kp1  + bb(3,6)*kp2  + bb(4,6)*kp3  + bb(5,6)*kp4 &
        ,rotstmp5 &
        ,G,Mp,Ms,Rp,Rs_tabl(6),sigmap,sigmas_tabl(6),ka5)
    ka5  =  dt * ka5

    if ( etmp .gt. 0.0d0) then
        call edot(atmp    + bb(1,6)*ka0  + bb(2,6)*ka1  + bb(3,6)*ka2  + bb(4,6)*ka3  + bb(5,6)*ka4 &
            ,etmp    + bb(1,6)*ke0  + bb(2,6)*ke1  + bb(3,6)*ke2  + bb(4,6)*ke3  + bb(5,6)*ke4 &
            ,oblptmp + bb(1,6)*kop0 + bb(2,6)*kop1 + bb(3,6)*kop2 + bb(4,6)*kop3 + bb(5,6)*kop4 &
            ,oblstmp + bb(1,6)*kos0 + bb(2,6)*kos1 + bb(3,6)*kos2 + bb(4,6)*kos3 + bb(5,6)*kos4 &
            ,rotptmp + bb(1,6)*kp0  + bb(2,6)*kp1  + bb(3,6)*kp2  + bb(4,6)*kp3  + bb(5,6)*kp4 &
            ,rotstmp5 &
            ,G,Mp,Ms,Rp,Rs_tabl(6),sigmap,sigmas_tabl(6),ke5) 
        ke5  =  dt * ke5
    else
        ke5  =  0.0d0
    endif 

    if ( oblptmp .gt. 0.0d0) then 
        call obldot(atmp    + bb(1,6)*ka0  + bb(2,6)*ka1  + bb(3,6)*ka2  + bb(4,6)*ka3  + bb(5,6)*ka4 &
            ,etmp    + bb(1,6)*ke0  + bb(2,6)*ke1  + bb(3,6)*ke2  + bb(4,6)*ke3  + bb(5,6)*ke4 &
            ,oblptmp + bb(1,6)*kop0 + bb(2,6)*kop1 + bb(3,6)*kop2 + bb(4,6)*kop3 + bb(5,6)*kop4 &
            ,rotptmp + bb(1,6)*kp0  + bb(2,6)*kp1  + bb(3,6)*kp2  + bb(4,6)*kp3  + bb(5,6)*kp4 &
            ,G,Mp,Ms,Rp,sigmap,rg2p,kop5)
        kop5  =  dt * kop5
    else
        kop5  =  0.0d0
    endif 

    if ( oblstmp .gt. 0.0d0) then
        call obldot(atmp    + bb(1,6)*ka0  + bb(2,6)*ka1  + bb(3,6)*ka2  + bb(4,6)*ka3  + bb(5,6)*ka4 &
            ,etmp    + bb(1,6)*ke0  + bb(2,6)*ke1  + bb(3,6)*ke2  + bb(4,6)*ke3  + bb(5,6)*ke4 &
            ,oblstmp + bb(1,6)*kos0 + bb(2,6)*kos1 + bb(3,6)*kos2 + bb(4,6)*kos3 + bb(5,6)*kos4 &
            ,rotstmp5 &
            ,G,Ms,Mp,Rs_tabl(6),sigmas_tabl(6),rg2s5,kos5) 
        kos5  =  dt * kos5
    else
        kos5  =  0.0d0
    endif 

    if (ori .eq. 1) then
        call frots_dot(atmp    + bb(1,6)*ka0  + bb(2,6)*ka1  + bb(3,6)*ka2  + bb(4,6)*ka3  + bb(5,6)*ka4 &
            ,etmp    + bb(1,6)*ke0  + bb(2,6)*ke1  + bb(3,6)*ke2  + bb(4,6)*ke3  + bb(5,6)*ke4 &
            ,oblstmp + bb(1,6)*kos0 + bb(2,6)*kos1 + bb(3,6)*kos2 + bb(4,6)*kos3 + bb(5,6)*kos4 &
            ,rotstmp5 &
            ,G,Mp,Ms,Rs_tabl(6),sigmas_tabl(6),rg2s5,ks5)
        ks5  =  dt * ks5

        call fwind(rotstmp5 &
            ,spin_sat,Kwind,Ms,Msun &
            ,Rs_tabl(6),Rsun,rg2s5,alphac,kw5)
        kw5  =  dt * kw5
    else
        call frots_dot2(atmp    + bb(1,6)*ka0  + bb(2,6)*ka1  + bb(3,6)*ka2  + bb(4,6)*ka3  + bb(5,6)*ka4 &
            ,etmp    + bb(1,6)*ke0  + bb(2,6)*ke1  + bb(3,6)*ke2  + bb(4,6)*ke3  + bb(5,6)*ke4 &
            ,oblstmp + bb(1,6)*kos0 + bb(2,6)*kos1 + bb(3,6)*kos2 + bb(4,6)*kos3 + bb(5,6)*kos4 &
            ,rotstmp5 &
            ,G,Mp,Ms,Rs_tabl(6),sigmas_tabl(6),rg2s5,ks5)
        ks5  =  dt * ks5
    endif

    !write(6,*)'ka5  = ',ka5
    !write(6,*)'kos5 = ',kos5
    !write(6,*)'ks5 = ',ks5
    !write(6,*)'kw5 = ',kw5

    if ( babaorum .eq. 1) kp5  =  0.0
    if ( babaorum .eq. 2) then 
        call rotp_dot(atmp    + bb(1,6)*ka0  + bb(2,6)*ka1  + bb(3,6)*ka2  + bb(4,6)*ka3  + bb(5,6)*ka4 &
            ,etmp    + bb(1,6)*ke0  + bb(2,6)*ke1  + bb(3,6)*ke2  + bb(4,6)*ke3  + bb(5,6)*ke4 &
            ,oblptmp + bb(1,6)*kop0 + bb(2,6)*kop1 + bb(3,6)*kop2 + bb(4,6)*kop3 + bb(5,6)*kop4 &
            ,rotptmp + bb(1,6)*kp0  + bb(2,6)*kp1  + bb(3,6)*kp2  + bb(4,6)*kp3  + bb(5,6)*kp4 &
            ,G,Mp,Ms,Rp,sigmap,rg2p,kp5)
        kp5  =  dt * kp5
    endif
    
    
    if (flagcrash .eq. 1 .or. tideop .eq. 0 ) then
    	ka5 = 0.0
    	ke5 = 0.0
    	kop5 = 0.0
    	kos5 = 0.0
    	ks5 = 0.0  
    	kp5 = 0.0
    	frotstmp = 0.0  	
    endif

    if ( etmp .lt. 1d-6) then 
        ke = 0. 
        kes = 0.
    endif
    if ( oblptmp .lt. 1d-4) then 
        kop = 0. 
        kops = 0.
    endif    
    if ( oblstmp .lt. 1d-4) then 
        kos = 0. 
        koss = 0.
    endif

    ka   =  cc(1)*ka0 + cc(2)*ka1 + cc(3)*ka2 + cc(4)*ka3 + cc(5)*ka4 + cc(6)*ka5
    ke   =  cc(1)*ke0 + cc(2)*ke1 + cc(3)*ke2 + cc(4)*ke3 + cc(5)*ke4 + cc(6)*ke5
    kop  =  cc(1)*kop0 + cc(2)*kop1 + cc(3)*kop2 + cc(4)*kop3 + cc(5)*kop4 + cc(6)*kop5
    kos  =  cc(1)*kos0 + cc(2)*kos1 + cc(3)*kos2 + cc(4)*kos3 + cc(5)*kos4 + cc(6)*kos5
    ks   =  cc(1)*ks0 + cc(2)*ks1 + cc(3)*ks2 + cc(4)*ks3 + cc(5)*ks4 + cc(6)*ks5
    kp   =  cc(1)*kp0 + cc(2)*kp1 + cc(3)*kp2 + cc(4)*kp3 + cc(5)*kp4 + cc(6)*kp5
    if (ori .eq. 1) kw   =  cc(1)*kw0 + cc(2)*kw1 + cc(3)*kw2 + cc(4)*kw3 + cc(5)*kw4 + cc(6)*kw5

    kas   =  ccs(1)*ka0 + ccs(2)*ka1 + ccs(3)*ka2 + ccs(4)*ka3 + ccs(5)*ka4 + ccs(6)*ka5
    kes   =  ccs(1)*ke0 + ccs(2)*ke1 + ccs(3)*ke2 + ccs(4)*ke3 + ccs(5)*ke4 + ccs(6)*ke5
    kops  =  ccs(1)*kop0 + ccs(2)*kop1 + ccs(3)*kop2 + ccs(4)*kop3 + ccs(5)*kop4 + ccs(6)*kop5
    koss  =  ccs(1)*kos0 + ccs(2)*kos1 + ccs(3)*kos2 + ccs(4)*kos3 + ccs(5)*kos4 + ccs(6)*kos5
    kss   =  ccs(1)*ks0 + ccs(2)*ks1 + ccs(3)*ks2 + ccs(4)*ks3 + ccs(5)*ks4 + ccs(6)*ks5
    kps   =  ccs(1)*kp0 + ccs(2)*kp1 + ccs(3)*kp2 + ccs(4)*kp3 + ccs(5)*kp4 + ccs(6)*kp5
    if (ori .eq. 1) kws   =  ccs(1)*kw0 + ccs(2)*kw1 + ccs(3)*kw2 + ccs(4)*kw3 + ccs(5)*kw4 + ccs(6)*kw5
    

    !*******************************************************************************

    ! Error estimate :
    
    !write(6,*) cc 
    
    error(1)  =  ka - kas
    error(2)  =  ke - kes
    error(3)  =  kop - kops
    error(4)  =  kos - koss
    error(5)  =  kp - kps
    error(6)  =  ks - kss
    if (ori .eq. 1) error(7)  =  kw - kws
    if (ori .eq. 0) error(7)  =  0.
    

    ! Runge Kutta integration. Next time step : 

    if ( ka .eq. ka .and. tideop .eq. 1) then 
        atmp     =  atmp + ka
    else
        atmp     =  atmp
    endif 

    if ( ke .eq. ke .and. tideop .eq. 1) then 
        etmp     =  etmp + ke
    else
        etmp     =  etmp
    endif 

    if ( kop .eq. kop .and. tideop .eq. 1) then 
        oblptmp  =  oblptmp + kop
    else
        oblptmp  =  oblptmp
    endif 

    if ( kos .eq. kos .and. tideop .eq. 1) then 
        oblstmp  =  oblstmp + kos
    else
        oblstmp  =  oblstmp
    endif 

    if ( ks .eq. ks .and. tideop .eq. 1) then 
        frotstmp =  frotstmp + ks
    else
        frotstmp =  frotstmp
    endif 

    if (ori .eq. 1) then
        fwindtmp =  fwindtmp + kw
        rotstmp  =  wconvi*rg2s0/rg2s*(Rs0/Rs_tabl(1))**2*exp(frotstmp)*exp(fwindtmp)
        
        Iconvtmp = rg2s * Ms*1e3 * (Rs_tabl(1)*100)**2 
        call fwind(rotstmp,spin_sat,Kwind,Ms*1e3,Msun*1e3,Rs_tabl(1)*1e3,Rsun*1e3,rg2s,alphac,DAMwind)
        AMtottmp = Iconvtmp*rotstmp
        DAMbrake  = abs(fwindtmp*Iconvtmp*rotstmp*1.e-3*1.e-4*1e7)
        DAMtot = abs(AMtottmp - AMtotinit)
        DAMtide = abs(DAMtot-DAMbrake)
        DAMtide = abs(frotstmp*Iconvtmp*rotstmp)
        
        
        !write(6,*) "COUCOU TEST CACA", AMtottmp,AMtotinit,abs(AMtottmp - AMtotinit),DAMbrake,DAMtide,rotstmp

        
    else
        !Need to express frotstmp as a variation of rotation rate =>  frotstmp = frotstmp / Iconv(j)
        
        !Moment of inertia of the star. frotstmp from frots_dot2 gives dJ/dt due to tide.
        !We need to translate this AM evolution into Omega evolution. 
        ! I = radius of gyration M R**2
        ! dOmega /dt approx 1/Iconv *  dJ/dt
        Iconvtmp = r2convs_tabl(5) * mstar * msun * Rs_tabl(5)**2 
    
        frotstmp = frotstmp / Iconvtmp

        !include tides?    
        !dOmega approx 1/Iconv *  dJ/dt * dt
        ! TO VERIFY, but if I'm not mistaken, frotstmp already was multiplied by
        ! dt, no no need to do it again?
        
        if (tideop .eq. 1) then
        
        	rotstmp = rotstmp + frotstmp
        else
        	rotstmp = rotstmp
        endif	
        
        
        if (flagcrash .eq. 1) then
        	rotstmp = rotstmp
        endif

    endif    
        

      if ( babaorum .eq. 1) then 
        call pseudorot(etmp,G,Mp,Ms,rotptmp)
        rotptmp  = rotptmp *atmp**(-1.5d0)
    endif
        
    if ( babaorum .eq. 2) rotptmp  =  rotptmp +  kp
    

    call pseudorot(etmp,G,Mp,Ms,rotpnew)
    rotpnew = rotpnew*atmp**(-1.5d0)

    ! if ( oblp = 0 .and. synchronization, don't calculate spin anymore
    if ( (oblptmp .eq. 0.0d0) .and. (abs(rotptmp-rotpnew) .le. 1d-9)) then
        babaorum = 1
    endif

    call aeq(etmp,rotptmp,rotstmp,G,Mp,Ms,Rp,Rs_tabl(1),sigmap,sigmas,aequtmp)

    if ( mod(j,10000) .eq. 0) then 
        write(6,*)''
        write(6,*)'time =',(ttmp-t_init)/yr
        call obldot(atmp,etmp,oblptmp,rotptmp,G,Mp,Ms,Rp,sigmap,rg2p,tmp)
        write(6,*) "oblptmp",oblptmp, tmp
        write(6,*)'tau_oblp =',abs(oblptmp/(tmp*yr)),' yr'
        call obldot(atmp,etmp,oblstmp,rotstmp,G,Ms,Mp,Rs_tabl(1),sigmas,rg2s,tmp)
        write(6,*)'tau_obls =',abs(oblstmp/(tmp*yr)),' yr'
        call adot(atmp,etmp,oblptmp,oblstmp &
            ,rotptmp,rotstmp,G &
            ,Mp,Ms,Rp,Rs_tabl(1),sigmap,sigmas,tmp)
        write(6,*)'tau_a =',abs(atmp/(tmp*yr)),' yr'
        call edot(atmp,etmp,oblptmp,oblstmp,rotptmp,rotstmp &
            ,G,Mp,Ms,Rp,Rs_tabl(1),sigmap,sigmas,tmp)
        write(6,*)'tau_e =',abs(etmp/(tmp*yr)),' yr'
        call rotp_dot(atmp,etmp,oblptmp,rotptmp &
            ,G,Mp,Ms,Rp,sigmap,rg2p,tmp)
        write(6,*)'tau_rotp =',abs(rotptmp/(tmp*yr)),' yr'
        !write(6,*)'rotpdot =',abs(rotp_dot(atmp,etmp,oblptmp,rotptmp &
            !,G,Mp,Ms,Rp,sigmap,rg2p))
        call frots_dot(atmp,etmp,oblstmp &
            ,rotstmp,G,Mp,Ms &
            ,Rs_tabl(1),sigmas,rg2s,tmp)
        write(6,*)'tau_rots_tides =',abs(1./(tmp*yr)), frotstmp*dt/(ks*yr),' yr'
        write(6,*)'tau_rots_wind =',fwindtmp*dt/(kw*yr),' yr'
        write(6,*)'tau_rots_both =',rotstmp*dt/((ks+kw)*yr),' yr'
        write(6,*)'error',error
        write(6,*)'dt',dt/yr
        write(6,*)''
    endif

    ! Write in table : ! superzoyotte
    
   
    
    !if ( ((ttmp .gt. (tdiskparam-0.05*tdiskparam)*yr) .and. (ttmp .lt. (tdiskparam+0.05*tdiskparam)*yr))  &
    !.or. ((ttmp .gt. (timedec-0.05*timedec)*yr) .and. (ttmp .lt. (timedec+0.05*timedec)*yr)) &
    !.or. ((ttmp .gt. 7.e6*yr) .and. (ttmp .lt. 30.e6*yr)) ) then
    
    if ( ((ttmp .gt. (tdiskparam-0.05*tdiskparam)*yr) .and. (ttmp .lt. (tdiskparam+0.05*tdiskparam)*yr))  &
    .or. ((ttmp .gt. (timedec-0.05*timedec)*yr) .and. (ttmp .lt. (timedec+0.05*timedec)*yr)) .or. ((out .eq. 0) &
     .and. (flagcrash .eq. 1)) .or. (ttmp .eq. t(1))  ) then
    
    if (flagcrash .eq. 1  ) then
    	out = 1
            	open(unit=2,file=filename2,access = 'append',status='old')
	            write(2,111) (ttmp-t(1))/yr &
    	            ,atmp/AU  &
        	        ,etmp  &
            	    ,oblptmp*180.d0/Pi  &
                	,oblstmp*180.d0/Pi  &
	                ,rotptmp ,rotstmp  &
    	            ,Rp/Rearth ,Rs_tabl(1)/AU  &
        	        ,sigmas  &
            	    !,frotstmp ,1.d0/(frotstmp*yr)  &
                	,(ks+kw)/(rotstmp*dt) ,rotstmp*dt/((ks+kw)*yr)  &
	                ,ka/(atmp*dt) ,(atmp*dt)/(ka*yr)  &
    	            !,frotstmp*rotstmp  &
        	        ,(ks+kw)/dt  &
            	    ,dt/yr 

            	close(2)
    endif	
     !.or. (ttmp .gt. 1.e10*yr) 
    
    	
		!write(6,*) "Super Connard 1",  ttmp   , tdiskparam*yr

        if ( etmp .gt. 1.0d-98) then 
        
        	if (flagcrash .eq. 0) then
            	open(unit=2,file=filename2,access = 'append',status='old')
	            write(2,111) (ttmp-t(1))/yr &
    	            ,atmp/AU  &
        	        ,etmp  &
            	    ,oblptmp*180.d0/Pi  &
                	,oblstmp*180.d0/Pi  &
	                ,rotptmp ,rotstmp  &
    	            ,Rp/Rearth ,Rs_tabl(1)/AU  &
        	        ,sigmas  &
            	    !,frotstmp ,1.d0/(frotstmp*yr)  &
                	,(ks+kw)/(rotstmp*dt) ,rotstmp*dt/((ks+kw)*yr)  &
	                ,ka/(atmp*dt) ,(atmp*dt)/(ka*yr)  &
    	            !,frotstmp*rotstmp  &
        	        ,(ks+kw)/dt  &
            	    ,dt/yr 

            	close(2)
            endif
          
             ! open(unit=23,file=filename3,access = 'append',status='old')
             ! write(23,*) (ttmp-t(1))/yr/1.e6,rotstmp/2.87e-6,wradf/2.87e-6,frotstmp,DAMbrake, DAMtide, &
             ! frotstmp,Wupconv,Wdownconv,Wjuprad,deltaWc
             
              open(unit=23,file=filename3,access = 'append',status='old')
              write(23,*) (ttmp-t(1))/yr/1.e6,rotstmp/2.87e-6,wradf/2.87e-6,frotstmp,DAMbrake, DAMtide, &
            frotstmp,Wupconv,Wdownconv,Wjuprad,deltaWc
              close(23)

        else
        	
        	if (flagcrash .eq. 0) then
	            open(unit=2,file=filename2,access = 'append',status='old')
    	        write(2,111) (ttmp-t(1))/yr &
        	        ,atmp/AU  &
            	    ,0.0d0  &
                	,oblptmp*180.d0/Pi  &
	                ,oblstmp*180.d0/Pi  &
    	            ,rotptmp ,rotstmp  &
        	        ,Rp/Rearth ,Rs_tabl(1)/AU  &
            	    ,sigmas  &
                	!,frotstmp ,1.d0/(frotstmp*yr)  &
	                ,(ks+kw)/(rotstmp*dt) ,rotstmp*dt/((ks+kw)*yr)  &
    	            ,ka/(atmp*dt) ,(atmp*dt)/(ka*yr)  &
        	        !,frotstmp*rotstmp  &
            	    ,(ks+kw)/dt  &
                	,dt/yr 
            	close(2)
            endif
            
            open(unit=23,file=filename3,access = 'append',status='old')
            write(23,*) (ttmp-t(1))/yr/1.e6, rotstmp/2.87e-6,wradf/2.87e-6,frotstmp,DAMbrake, DAMtide, &
            frotstmp,Wupconv,Wdownconv,Wjuprad,deltaWc
            close(23)
            
        endif 
    endif    
    
    
    if ( ttmp .eq. t_init) kkk = 2
    
    if (ori .eq. 0) then
    	 !if ( (ttmp .ge. t_init) .and. (flagtinit .eq. 0)) then
    	 !	kkk = 2
    	 !	flagtinit = 1
    	 !endif	
    endif
    
    if ( ttmp .ge. hamlet(kkk)*yr ) then 
    
    
    
        t(kkk+1)     =  ttmp
        a(kkk+1)     =  atmp
        e(kkk+1)     =  etmp
        oblp(kkk+1)  =  oblptmp
        obls(kkk+1)  =  oblstmp
        rotp(kkk+1)  =  rotptmp
        rots(kkk+1)  =  rotstmp
        aequ(kkk+1)  =  aequtmp
        sigmast(kkk+1)  =  sigmas
        


        !write(6,*) "rotstmp =" ,rotstmp, "rots(kkk+1)",rots(kkk+1)

        if ( e(kkk+1) .gt. 1.0d-98) then 
        
        	if (flagcrash .eq. 0) then
            	open(unit=2,file=filename2,access = 'append',status='old')
	            write(2,111) (t(kkk+1)-t(1))/yr &
    	            ,a(kkk+1)/AU  &
        	        ,e(kkk+1)  &
            	    ,oblp(kkk+1)*180.d0/Pi  &
                	,obls(kkk+1)*180.d0/Pi  &
	                ,rotp(kkk+1) ,rots(kkk+1)  &
    	            ,Rp/Rearth ,Rs_tabl(1)/AU  &
        	        ,sigmast(kkk+1)  &
            	    !,frotstmp ,1.d0/(frotstmp*yr)  &
                	,(ks+kw)/(rotstmp*dt) ,rotstmp*dt/((ks+kw)*yr)  &
	                ,ka/(atmp*dt) ,(atmp*dt)/(ka*yr)  &
    	            !,frotstmp*rotstmp  &
        	        ,(ks+kw)/dt  &
            	    ,dt/yr 

            	close(2)
            endif
          
              open(unit=23,file=filename3,access = 'append',status='old')
              write(23,*) (t(kkk+1)-t(1))/yr/1.e6,rotstmp/2.87e-6,wradf/2.87e-6,frotstmp,DAMbrake, DAMtide, &
              frotstmp,Wupconv,Wdownconv,Wjuprad,deltaWc
              close(23)

			
        else
        	if (flagcrash .eq. 0) then
	            open(unit=2,file=filename2,access = 'append',status='old')
    	        write(2,111) (t(kkk+1)-t(1))/yr &
        	        ,a(kkk+1)/AU  &
            	    ,0.0d0  &
                	,oblp(kkk+1)*180.d0/Pi  &
	                ,obls(kkk+1)*180.d0/Pi  &
    	            ,rotp(kkk+1) ,rots(kkk+1)  &
        	        ,Rp/Rearth ,Rs_tabl(1)/AU  &
            	    ,sigmast(kkk+1)  &
                	!,frotstmp ,1.d0/(frotstmp*yr)  &
	                ,(ks+kw)/(rotstmp*dt) ,rotstmp*dt/((ks+kw)*yr)  &
    	            ,ka/(atmp*dt) ,(atmp*dt)/(ka*yr)  &
        	        !,frotstmp*rotstmp  &
            	    ,(ks+kw)/dt  &
                	,dt/yr 
            	close(2)
            endif
            
            open(unit=23,file=filename3,access = 'append',status='old')
            write(23,*) (t(kkk+1)-t(1))/yr/1.e6, rotstmp/2.87e-6,wradf/2.87e-6,frotstmp,DAMbrake, DAMtide, &
            frotstmp,Wupconv,Wdownconv,Wjuprad,deltaWc
            close(23)
            
        endif 
        kkk = kkk+1
        !break
        
        if (kkk .ge. size(hamlet)) then
        	write(6,*) "Outside bound, stop"
        	stop
        endif
        
    endif
    

    ttmp = ttmp + dt

    ! Calculation of value of next timestep
    errmax   =  0. 
    do k = 1,size(error)
        if ( error(k) .eq. 0.) errmax = max(errmax,0.)
        if ( error(k) .ne.0.) errmax = max(errmax,abs(error(k)/yscal(k)))
    enddo
    errmax = errmax/eps

    if ( mod(j,10000) .eq. 0) then 

        do k=1,7
            write(6,*)'yscal(k) =',yscal(k)
        enddo
        
    !write(6,*) "ka's= " ,ka0,ka1,ka2,ka3,ka4,ka5
    !write(6,*) "ke's= " ,ke0,ke1,ke2,ke3,ke4,ke5
    !write(6,*) "kop's= " ,kop0,kop1,kop2,kop3,kop4,kop5
    !write(6,*) "kos's= " ,kos0,kos1,kos2,kos3,kos4,kos5
    !write(6,*) "ks's= " ,ks0,ks1,ks2,ks3,ks4,ks5
    !write(6,*) "kw's= " ,kw0,kw1,kw2,kw3,kw4,kw5
    !write(6,*) "kp's= " ,kp0,kp1,kp2,kp3,kp4,kp5
    
    !write(6,*) "ka, kas =", ka, kas
	!write(6,*) "kp, kps =", kp, kps
        
    !write(6,*) "error's = ", error
    
       
       write(6,*)'abs(error(k)/yscal(k)) =',abs(error((1))/yscal(1)) &
             ,abs(error(2)/yscal(2)),abs(error(3)/yscal(3)) &
             ,abs(error(4)/yscal(4)),abs(error(5)/yscal(5)) &
             ,abs(error(6)/yscal(6)),abs(error(7)/yscal(7))
       write(6,*)'Calculation of value of next timestep'
       write(6,*)'errmax = ',errmax
       write(6,*)'errcon = ',errcon
       
    endif
    

    if ( errmax .gt. 1.) then 

        dttemp  =  safety*dt*(errmax**(pshrink))
        dt  =  max(abs(dttemp),0.1*abs(dt),dtmin) 
        if ( mod(j,10000) .eq. 0) then 
            write(6,*)'errmax > 1, dt must be decreased'
            write(6,*)'dttemp = ',dttemp/yr
            write(6,*)'dt = ',dt/yr
        endif

        if ( dt .eq. 0) then 
            write(6,*)'stepsize underflow'
            stop
        endif
    else  
       !! if ( j mod 100000 .eq. 0) then 
       !!    write(6,*)'errmax < 1, dt must be increased'
       !! endif
        if ( errmax .gt. errcon) then 
            dtnext  = min(safety*dt*(errmax**(pgrow)), dtmax)
            !! if ( j mod 100000 .eq. 0) then 
            !!    write(6,*)'errmax > errcon, dt must be increased'
            !!    write(6,*)'dt = dtnext = ',dtnext/yr
            !! endif
        else
            dtnext  = min(5.*dt, dtmax)
            !! if ( j mod 100000 .eq. 0) then 
            !!    write(6,*)'errmax < errcon, dt must be increased'
            !!    write(6,*)'dt = dtnext = ',dtnext/yr
            !! endif
        endif 

        !if ((ttmp-dt)/yr/1.e6 .gt. 10.0) stop

		
        dt      =  dtnext
        
    endif 

    !*******************************************************************************
    !*******************************************************************************
    !*******************************************************************************

    ! ERROR/INFORMATION MESSAGES, 
    
    if (flagcrash .eq. 0) then

    if ( atmp*(1-etmp) .le. a_roche) then  
        write(6,*)'ROCHE LIMIT REACHED -----> PLANET DESTROYED !!!'
        write(6,*)'time of impact = ',ttmp/yr,'  yrs'
        write(6,*)'atmp =',atmp/AU,' AU'
        write(6,*)'a_roche = ',a_roche/AU,' AU'
        write(6,*)'j = ',j
        flagcrash = 1
    endif

    if ( etmp .lt. 0.0) then 
        write(6,*)'WARNING : e<0, you need to shorten the time step !'
        write(6,*)'OR, this could mean : CRASH !!!! PLANET DESTROYED !!'
        write(6,*)'time of impact = ',ttmp/yr,'  yrs'
        write(6,*)'j = ',j
        flagcrash = 1
    endif

    if ( atmp*(1-etmp) .le. Rs) then 
        write(6,*)'Planet reaches surface star !!!!!!!!'
        write(6,*)'time of impact = ',ttmp/yr,'  yrs'
        write(6,*)'atmp =',atmp/AU,' AU'
        write(6,*)'Rs = ',Rs_tabl(1)/AU,' AU'
        write(6,*)'a_roche = ',a_roche/AU,' AU'
        flagcrash = 1
    endif

    infinity = 5d8!5d7 !1d9!

    if ( abs(ka) .ne. abs(ka) ) then 
        write(6,*)'WARNING : derivative of a = infinity' 
        write(6,*)'This could mean : CRASH !!!! PLANET DESTROYED !!'
        write(6,*)'time of impact = ',ttmp/yr,'  yrs'
        write(6,*)'j = ',j
        write(6,*) "ka = ",ka
        flagcrash = 1
    endif 
    

    if ( abs(ke) .ne. abs(ke)  ) then 
        write(6,*)'WARNING : derivative of e = infinity' 
        write(6,*)'This could mean : CRASH !!!! PLANET DESTROYED !!'
        write(6,*)'time of impact = ',ttmp/yr,'  yrs'
        write(6,*)'j = ',j
        flagcrash = 1
    endif 

    if ( abs(kp) .ne. abs(kp) ) then 
        write(6,*)'WARNING : derivative of rotp = infinity' 
        write(6,*)'This could mean : CRASH !!!! PLANET DESTROYED !!'
        write(6,*)'time of impact = ',ttmp/yr,'  yrs'
        write(6,*)'j = ',j
        flagcrash = 1
    endif 

    if ( abs(ks) .ne. abs(ks) ) then 
        write(6,*)'WARNING : derivative of rots = infinity' 
        write(6,*)'This could mean : CRASH !!!! PLANET DESTROYED !!'
        write(6,*)'time of impact = ',ttmp/yr,'  yrs'
        write(6,*)'j = ',j
        flagcrash = 1
    endif 

    if ( rotptmp .lt. 0.0) then 
        write(6,*)'WARNING : rotp<0 !'
        write(6,*)'this could mean : CRASH !!!! PLANET DESTROYED !!'
        write(6,*)'time of impact',ttmp/yr,'  yrs'
        write(6,*)'j = ',j
        flagcrash = 1
    endif 
    endif 
    
    if (flagcrash .eq. 1) then
    	tideop = 0
    	flagcrashout = 1
    endif 

    !************************** END MESSAGES ******************************

    j = float(j) + 1

    !update new initial values

    if (ori .eq. 0) then
        
        Rs0      = Rs_tabl(5)
        Teffs0   = Teffs_tabl(5)
        Ls0      = Ls_tabl(5)
        r2convs0 = r2convs_tabl(5)
        r2rads0  = r2rads_tabl(5)
        Mrads0   = Mrads_tabl(5)
        Rrads0   = Rrads_tabl(5)
        wconvi   = rotstmp
        wradi    = wradf
    endif
    
    
    if (ori .eq. 1) then
    
    	AMtotinit = AMtottmp
    endif
        
        !write(6,*) "wconvi final= ", wconvi,rotstmp
        
    if(comp .eq. 3) then
    	!stop
    endif
    
    comp = comp +1
    

enddo


12 write(6,*) "On bypasse le calcule"

	write(6,*) "Pas de pb à la fin du program"


!! write(6,*)'j+1    = ',j+1
!! write(6,*)'tf     = ',ttmp/yr,'  yrs'
!! write(6,*)'af     = ',atmp/AU,'  AU'
!! write(6,*)'ef     = ',etmp
!! write(6,*)'rotsf  = ',rotstmp,'  s-1'
!! write(6,*)'rotpf  = ',rotptmp,'  s-1'

!****************************** WRITE FILES ***********************************

open(unit=23,file=filename3,access = 'append',status='old')
write(23,*) (ttmp-t(1))/yr/1.e6, rotstmp/2.87e-6,wradf/2.87e-6,frotstmp,DAMbrake, DAMtide, &
frotstmp,Wupconv,Wdownconv,Wjuprad,deltaWc
close(23)


open(unit=1,file=filename1,access = 'append',status='old')
write(1,*)'t_final     = ',ttmp/yr,' yrs'
write(1,*)'a_final     = ',atmp/AU,' AU'
write(1,*)'e_final     = ',etmp
write(1,*)'oblp_final  = ',oblptmp*180.d0/Pi,' deg'
write(1,*)'obls_final  = ',oblstmp*180.d0/Pi,' deg'
write(1,*)'rotp_final  = ',rotptmp,' s-1'
write(1,*)'rots_final  = ',rotstmp,' s-1'
write(1,*)'Rs_final      = ',Rs_tabl(1)/Rsun, ' Rsun' &
       ,'  = ',Rs_tabl(1)/AU,' AU'
write(1,*)'sigmas_final      = ',sigmas
close(1)

if (flagcrash .eq. 0) then
	open(unit=2,file=filename2,access = 'append',status='old')
	write(2,111) (ttmp-t(1))/yr,atmp/AU,etmp,oblptmp*180.d0/Pi &
 	   	,oblstmp*180.d0/Pi,rotptmp,rotstmp,Rp/Rearth,Rs_tabl(1)/AU &
		,sigmas &
		!,frotstmp,1.d0/(frotstmp*yr) &
  	  	,(ks+kw)/(rotstmp*dt),rotstmp*dt/((ks+kw)*yr) &
    	,ka/(atmp*dt),(atmp*dt)/(ka*yr) &
    	!,frotstmp*rotstmp &
    	,(ks+kw)/dt &
    	,dt/yr 
	close(2)
endif


end subroutine TIDESTEST

end module tide
