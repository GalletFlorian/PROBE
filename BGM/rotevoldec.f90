module rotevoldec

use constants
	  
	  
	  contains

	  subroutine rotevol(taudecinittmp,tdiskparam,wconvi,wconvf,wradi,wradf,Lum2,Teff2,tstar2,Rstar2,k2conv2,k2rad2,Mrad2,Rrad2,&
	  		Wuprad,juprad,deltaJ,timedec,DAMbrake,DAMtide,Wtide,Wupconv,Wdownconv,Wjuprad,deltaWc,mstarin)


      implicit none

      common/const/mstar2,msol,rsol,Wsol,taudecinit,mdotsun,alphaS,conS
      common/const2/numtest
      common/var/rstar,tstar,rrad,mrad,Irad,Iconv
      common/var/k2rad,k2conv,Lum,Teff
      common/freinage/ksk,kmm,ksc,kmp,ksat,wcrit,wsat
      common/freinage/G,mdotstar,flag2,Bstar,Ro,ff,mdotstarC,mdotstarJ
      common/indexs/indexms,indexdec,j,nn,n,indexsun
      common/coeff/coeffri,coeffkir,coeffkic,coeffrri,coeffmri


      INTEGER NSTEP,NTRACK,NUMTEST,NN,NBTRACK

      PARAMETER(NTRACK=6) !Number of disk lifetime

! -----------------------------------------------------------
! |if you change Ntrack, you also have to change format #110|
! -----------------------------------------------------------
!
      PARAMETER(NSTEP= 2)
      PARAMETER(NBTRACK=1) !To select one of the Ntrack disk lifetime
      						! NBTRACK=1 => fast
      						! NBTRACK=3 => slow/median  

! added 12.07.07 Solar angular momentum 
	  REAL*8 test
       real(dp) :: wconvi,wconvf,wradi,wradf,timedisk
       real(dp) :: timedec
      
      
!       real(dp) :: taudecinit2,  tdiskparam2

      real(dp), dimension(2) :: Lum2,Teff2,tstar2,Rstar2
      real(dp),dimension(2) :: k2conv2,k2rad2,Mrad2,Rrad2
      real(dp),dimension(2) :: mstarin
      
      real(dp), dimension(2) :: Lum,Teff,tstar,Rstar,mstar2
      real(dp),dimension(2) :: k2conv,k2rad,Mrad,Rrad
	  
	  
      REAL*8 AMSOL,ISOL
      REAL*8 flag2
      REAL*8 mdotsun,wsun,mdotstar,Period,Bstar,Ro,ff
      REAL*8 G
      REAL*8 DAMSOL,lsol
      REAL*8 wbreakup(NSTEP),vbreakup(NSTEP)
      REAL*8 K2star(NSTEP)
      REAL* 8 ICONV(NSTEP),IRAD(NSTEP)
      REAL*8 TDISK,DT,DAM(NSTEP)
      REAL*8 wconviPAR,JUPRAD,WUPRAD,WUPCONV,WDOWNCONV,Wtide,Wjuprad
      REAL*8 WRAD(NSTEP),WCONV(NSTEP),WCRIT,WSAT,AMconv(nstep), DJ(nstep),AMrad(nstep)
      REAL*8 RINIT,DELTAWR,DELTAWC, deltaJ
!      real*8 TAUDEC
      real*8 taudec(NSTEP),breaking(NSTEP)
      REAL*8 DPI, KMS, MSOL, YEAR
      REAL*8 KSK,KMM,KSC,KMP,KSAT,a,b
      REAL*8 DUMM, RSOL, day2, WINTERFIN
      INTEGER IDUMM, INDEX, N, ITRACK,INDEXDEC,NDEBUT, INDEXMS, FLAG
      integer indexsun,index1,j,index2,i,l
      real*8 agesun
      real*8 pinit,winterfinr,wsol,taudecfin
      CHARACTER MODEL*30, HEAD*72
      CHARACTER BRAKING_LAW*2,IDENT*19
      real*8 coeffri(nstep),coeffkir(nstep),coeffkic(nstep) !interpolation spline coefficients
      real*8 coeffrri(nstep),coeffmri(nstep)                !interpolation spline coefficients
      REAL*8 q
      REAL*8 mdotcra,omegacra,mdotyr,mdotstarC,mdotstarJ
      REAL*8 acran,bcran,dt1,dt2
      INTEGER cramin , cramax , bo
      REAL*8 alphaS,conS
      CHARACTER*100 command
      
      REAL*8 Itot, AMtot,DAMtot,DAMbrake, DAMtide
  	  real(dp) ::  taudecinit,tdiskparam,taudecinittmp
  	  
      real*8 Wtest
      integer ninter,n1,n2,n3
  	  
  	 ! REAL*8 K,K1MP,K2MP,mmp,K3,K4,mvic
  	  
  	  !integer :: brklaw
  	  !real :: K,mmp,K1MP,K2MP,K3,K4,mvic
  	  
      NN=NTRACK
      
      
  	  

! ---------------------------------------------
! UNITS CGS: MSOL in G, RSOL in CM, WSOL in S-1
! ---------------------------------------------

		!write(6,*) "Dans rotevol",taudecinittmp,tdiskparam,mvic,mmp,brklaw


	  
	  wconvi = wconvi !/ (3600.0*24.0)
	  wradi = wradi !/ (3600.0*24.0)	  
	  
	  
	  !write(6,*) "Wconv = ", wconvi,wradi
	  
	
! added 12.07.07 Solar angular momentum
      AMSOL = 1.63D48
      !AMSOL = 1.84e48 !Pinto et al. 2010
      ISOL = 6.411d53
      DAMSOL = 7.169D30
      !DAMSOL = 3.67785E+30 !Calcul perso en prenant Matt et al. 2012 pour RA et Bsun = 2 Gauss
      MSOL = 1.989D33
      RSOL = 6.9599D+10
      WSOL=2.87D-6
      lsol = 3.826d+33
      kms = 1.1d5 
      day2 = 2.4d1*36.e2
      year = 365.25*day2
      numtest=1
      !pi = 3.14159
      dpi = 2.0*pi
      agesun=4.5e9
      mdotyr = 6.30276d+25 !mdotyr in Msol/yr
      alphaS = 0.0 !Classic
!      alphaS = 0.076 !Index Spada
      conS = 0.2 !Spada
!      alphaS = 0.2
      conS = 10
	 

! Added by F.G
!-----------------------------------------------------------------------
      G = 6.67300d-8
      mdotsun = 1.31d12  !Comes from Holzwarth & Jardine 2007 units : g/s
!-----------------------------------------------------------------------



 99   format(a72)
 98   format(a19,d10.4)
 200  format(a19,d10.6)
 97   format(a19,a30)
 96   format(a19,I1)

		ksk = 0.0
		kmm = 0.0
		ksc = 0.0
		kmp = 1.0
	  
      
	  
!_______Control braking law_______ 
!	  brklaw = 3   !brklaw = braking law si magnbrak.f brklaw=3 = Matt 2015
!      K=1.0
!      m=0.22
!      K1MP = 1.7
!      K2MP = 0.0506d0
!      K3 = 2.8
!      K4 = 0.2
!      mvic = 0.11
!_________________________________ 
      
	itrack = 0
      
	  do j=1,2
	  
	  	 Lum(j) = Lum2(j)
	  	 Teff(j) = Teff2(j)
	  	 tstar(j) = tstar2(j)
	  	 rstar(j)= Rstar2(j)
	  	 k2conv(j) = k2conv2(j)
	  	 k2rad(j) = k2rad2(j)
	  	 Mrad(j) = Mrad2(j)
	  	 Rrad(j) = Rrad2(j)
	  	 taudecinit = taudecinittmp
	  	 !tdiskparam = tdiskparam2
	  	 mstar2(j) = mstarin(j)
	  	 
            K2star(j) = k2conv(j) + k2rad(j)
            Rstar(j)=Rstar(j)*rsol
            Rrad(j)=Rrad(j)*Rstar(j)
            Irad(j)=k2rad(j)*mstar2(j)*msol*(Rstar(j))**2
            Iconv(j)=k2conv(j)*mstar2(j)*msol*(Rstar(j))**2
            
            
            wbreakup(j)=(2./3.)**(3./2.)*(G*mstar2(j)*msol)**(1./2.) /(rstar(j)**(3./2.))
            vbreakup(j)= (2.*G*mstar2(j)*msol/(3.*rstar(j)))**(1./2.)
            if (tstar(j).le.tdisk) index=j
            if (tstar(j).le. agesun) indexsun=j
            if (mstar2(j).eq. 1. .or. mstar2(j).eq. 0.5) then
               if (tstar(j).lt. 4.8e5) ndebut=j
            else
               if (tstar(j).le. 3.e5) ndebut=j
            endif
            if (tstar(j).lt. 6.3e7) indexms=j
            if ((k2rad(j).gt. 1.d-9).and.(k2rad(j-1).lt. 1.d-9)) indexdec=j 
            tstar(j)=tstar(j)*year
	  enddo

      
      if ((ksk.eq. 0.) .and. (kmm.eq. 0.) .and. (kmp.ne. 0.) .and. (brklaw .eq. 0)) then 
      !write(6,*) '******************'
      !write(6,*) 'Matt et al. (2012)'
      !write(6,*) '******************'
      endif
      
      
      if ((brklaw .eq. 1) .or. (brklaw .eq. 2) ) then
      !write(6,*) '***********************'
      !write(6,"(a24,I1)") 'Reville et al. (2014) v.',brklaw
      !write(6,*) '***********************'
      endif
      
      if (brklaw .eq. 3) then
      !!write(6,*) '******************'
      !!write(6,*) 'Matt et al. (2015)'
      !!write(6,*) '******************'
      endif
 


      ksat = 0.d0
      wcrit=0.d0


 1001 format(f4.1,2x,d10.4,2x,d10.4,2x,d10.4)
 1002 format(d10.4)


      !!write (*,*)'braking law constant :'
      !!write(6,95) K1MP,K2MP,K3
 95   format(' K1=',f5.2,' K2=',f10.6, ' K3=',f5.2)

      !!write (*,1111)taudecinit
      
! F.G 6/10/2014
! Taudec Spada
      taudecinit = taudecinit*year
      taudec (1) = taudecinit
      

!      taudec=taudec*year !taudec is now in second

 1111 format('coupling time=',e7.2,' years')


! Open the braking law backup file 


 512    format(1x,f10.4,1x,f12.6,1x,ES14.7E3,1x,f10.4)
 513    format(1x,f10.4,1x,f12.6,1x,f10.4,1x,f10.4)
 514    format(1x,f10.4,1x,f10.7,1x,f10.4,1x,f12.6)
 515    format(1x,f10.4,1x,ES14.7E3)


!-----------------------------------------------------------------------
! input: disk lifetime distribution : t_disk(i),i=1,ntrack
! (ntrack = number of track that will be calculated)
      
      tdisk= tdiskparam
!-----------------------------------------------------------------------



         index = 0 

!------------------Reiners braking law test------------------------
! If flag2 = 1. the model use the Reiners & Mohanty 2012 braking law
!
      flag2 = 0.
!------------------------------------------------------------------

		n = 2
 
 		do j = 1,n
 				WRAD(j) = 0.0
 				WCONV(j) = 0.0 
 				AMconv(j) = 0.0
 				DJ(j) = 0.0
 				AMrad(j) = 0.0
 		enddo
 

! the initial conditions refer to the moment at which the disk disappear
! because before that the star and the disk are couple 
! P_star = P_init => Wconv=wconvi ; Pconv=pinit

         rinit = rstar(1)
!         wconvipar = dpi/(pinit*day2)


!---------------------------------------------
! \m/ |ROTATIONAL EVOLUTION  computation| \m/
!---------------------------------------------

         if (indexdec.ge.index) then
! case A: the disk disappear before the core/envelope decoupling
            index1=index
            index2=indexdec-1
         else
! case B: the disk disappear after the core/envelope decoupling 
            index1=indexdec-1
            index2=index
         endif

		 wconv(1) = wconvi
		 wrad(1) = wradi
         AMconv(1)=iconv(1)*Wconv(1)
         AMrad(1)=Irad(1)*Wrad(1)
		 
    do j=2,2
         
         dt1= tstar(j)/year

		 !faire un truc avec tdisk		 
		 !Genre si t* < tdisk ET Irad = 0.0

		 
		
		 !if (dt1 .le. tdiskparam .and. Irad(j) .eq. 0.0) then
		 if (dt1 .le. tdiskparam .and. dt1 .lt. timedec) then
		 
		 !if (j .le. index1) then
            Wconv(j)=Wconv(j-1)
            Wrad(j)=0.
            dj(j)=0.
            AMconv(j)=iconv(j)*Wconv(j)
            AMrad(j)=Irad(j)*Wrad(j)
            taudec(j) = taudecinit
            
            
            
         !else if (dt1 .le. tdiskparam .and. Irad(j) .ne. 0.0) then 
         else if (dt1 .le. tdiskparam .and. dt1 .ge. timedec) then 
         
         	
! II B. the disk disappear after the decoupling
!       ***************************************
!      Core-enveloppe decoupling
!      -------------------------           
            !wdownconv=0.d0
            !Wupconv=0.d0
            !deltaWc=0.d0
            !deltaJ = 0.0
            
! F.G 6/10/2014
! Taudec Spada
     
                    
          taudec(j) = taudecinit
                    


!    2.B.1) disk-star coupling 
!           ------------------
               Wconv(j)=wconvi
               

!    2.B.2) instantaneous spin-up from contraction
!           --------------------------------------
               Wuprad=Wrad(j-1)*(Irad(j-1)/Irad(j) - 1.)
               

              
!    2.B.3) core developps
!           --------------
               Juprad=2./3.*Rrad(j)**2.*Wconv(j-1)*(Mrad(j)-Mrad(j-1))*msol
                              
        	if (Irad(j) .ne. 0.0) then
            	deltaWr=deltaJ*(tstar(j)-tstar(j-1))/(taudec(j)*Irad(j))
            else
            	deltaWr= 0.0
            endif

               flag= 0
               if (tstar(j) .gt. timedec*year) then !call interpolation routine if velocity change .gt. 10% of initial velocity
              ! write(6,*) "Avant interpol" ,Juprad/Irad(j),deltaJ,deltaWr,Wuprad 
               	call interpolB(Juprad,deltaJ,deltaWr,Wuprad,Wrad(j-1),Wconv(j-1),Winterfinr,itrack,flag,taudecfin)
               endif
!               flag = 0 !Remove interpolation
               if (flag .eq. 1) then
                  Wrad(j) = Winterfinr
                  taudec(j) = taudecfin
                  write(6,*) "InteprolB",tstar(j)/year/1.e6, Winterfinr
               else
                  Wrad(j) = Wrad(j-1)+Juprad/Irad(j)-deltaWr+Wuprad  
                  write(6,*) "Pas InteprolB",tstar(j)/year/1.e6, Wrad(j),Juprad/Irad(j),deltaWr,Wuprad                            
               endif
                 
            
               
               !write(6,*) "Rotevol",Wrad(j),Wrad(j-1),Juprad/Irad(j),deltaWr,Wuprad,flag  
               


!----------------------------------------------------------------------------------------------
!| We impose to the core velocity to be equal to the envelope velocity when the core/envelope |
!| decoupling starts.                                                                         |
!----------------------------------------------------------------------------------------------
               
               !if (k2rad(j) .ne. 0.0 .and. k2rad(j-1) .eq. 0.0) then             ! added JB 16.06.08
               !   Wrad(j) = Wconv(j) ! force the inital core velocity to be the convective velocity
               !   !write(6,*) "Condition", Wrad(j), Wconv(j)
               !end if
               write(6,*) "Time caca= ",tstar(1), tstar(2),timedec*year
               
           	   if ( (tstar(j) .ge. timedec*year)  .and. (  tstar(j-1) .lt. timedec*year ) ) then
           	   		Wrad(j) = Wconv(j)
           		write(6,*) "Condition, wrad 2",tstar(j)/year, wrad(j)/wsol, Wuprad,Juprad/Irad(j),deltaWr
           		!stop
           	   endif
           	              	            
           	   if(tstar(j) .ge. timedec*year) then
           	   		if ( Irad(j)/Irad(j-1) .lt. 1e-4) then
           				Wrad(j) = Wconv(j)
           			endif
           	   endif
               
                        
               dj(j)=deltaJ
               AMrad(j)=Irad(j)*Wrad(j)
               AMconv(j)=Iconv(j)*Wconv(j)

 109           format(f7.4,2x,2(F8.4,2x),7(d10.4,2x))

 
!     2.B.4) core-enveloppe angular momentum exchange
!            ----------------------------------------
               dt1= tstar(j)/(year*1.0e6)
               dt2= tstar(j-1)/(year*1.0e6)
               
               if (Wconv(j) .lt. Wrad(j)) then
                  deltaJ = Irad(j)*Iconv(j)/(Irad(j)+Iconv(j))*(Wrad(j)-Wconv(j))
               else if (Wconv(j) .eq. Wrad(j)) then
                  deltaJ = 0.0
               else if (Wconv(j) .gt. Wrad(j)) then
                  deltaJ= -Irad(j)*Iconv(j)/(Irad(j)+Iconv(j))* (Wconv(j)-Wrad(j)) !If the core rotates slower than the envelope
               	  !deltaJ = 0.0
               endif
               
               
               

		!else if (dt1 .gt. tdiskparam .and. Irad(j) .eq. 0.0) then
		else if (dt1 .gt. tdiskparam .and. dt1 .le. timedec ) then

		!Ensuite si t* > tdisk ET Irad = 0.0
		!copier
! II. A. The disk disappear before the decoupling (indexdec > index)
!     **************************************************************
!      Solid-body rotation
!      Disk-star uncoupled: spin-up contraction + magnetic braking
!      -----------------------------------------------------------
         !deltaj=0.d0
         !Wuprad=0.d0
         !juprad=0.d0
         
!    2.A.1) dJ/dt = 0, instantaneous spin-up from contraction 
!           -------------------------------------------------
               Wupconv = Wconv(j-1) * (Iconv(j-1)/Iconv(j)- 1.)
                

!    2.A.2) Instantaneous magnetic braking 
!           ------------------------------         
			
               call magnbrak(Wdownconv,Wconv(j-1),rstar(j-1),Lum(j-1),tstar(j),tstar(j-1),Iconv(j-1),itrack,Teff(j-1))
               
!  call interpolation routine if velocity change .gt. 10% of initial velocity            
               flag = 0

               call interpolA(Wdownconv,Wupconv,Wconv(j-1),Winterfin,itrack,braking_law,flag)
               
               if (flag .eq. 1) then
                  Wconv(j) = Winterfin
                  Wrad(j)=0.
                  !write(6,*) "Entrée dans interpolA"
               else
                  Wconv(j) = Wconv(j-1)+Wupconv-Wdownconv
                  Wrad(j)=0.
                  !write(6,*) "Pas entrée dans interpolA"
               endif

               dt1= tstar(j)/(year*1.0e6)
               dt2= tstar(j-1)/(year*1.0e6)

               dt=(tstar(j)-tstar(j-1))/year
               dj(j)=deltaJ 
               

     
               AMconv(j)=iconv(j)*Wconv(j)
               AMrad(j)=Irad(j)*Wrad(j)

 108           format(f7.4,2x,2(f8.3,2x),7(d10.4,2x))
 111           format(f7.4,2x,2(f7.3,2x),6(d10.4,2x))
 
 				


		else 
!   III Core-envelope decoupling & star-disk uncoupled
!       **********************************************


!     3.1) instantaneous spin-up from contraction 
!          --------------------------------------
            Wupconv = Wconv(j-1) * (Iconv(j-1)/Iconv(j)-1.)
            if (Irad(j) .ne. 0.0) then
              Wuprad = Wrad(j-1) * (Irad(j-1)/Irad(j)-1.)
            else
              Wuprad = 0.0
            endif
              
!     3.2) Instantaneous magnetic braking 
!          ------------------------------
            call magnbrak(Wdownconv,Wconv(j-1),rstar(j-1),Lum(j-1),tstar(j),tstar(j-1),Iconv(j-1),itrack,Teff(j-1))
				

           dt1= tstar(j)/(year*1.0e6)
               
! F.G 6/10/2014
! Taudec Spada


                 
         taudec(j) = taudecinit


     
               
!     3.3) Core develops
!          -------------
            Juprad=2./3.*Rrad(j)**2.*Wconv(j-1)*(Mrad(j)-Mrad(j-1))*msol  
            
                        
            if (Irad(j) .ne. 0.0) then
            	deltaWr=deltaJ*(tstar(j)-tstar(j-1))/(taudec(j)*Irad(j))
            else
            	deltaWr= 0.0
            endif
            deltaWc=deltaJ*(tstar(j)-tstar(j-1))/(taudec(j)*Iconv(j))


               
            flag=0
            !if (k2rad(j) .ne. 0.0 .and. k2rad(j-1) .ge. 0.0 ) then
            
            if (tstar(j) .gt. timedec*year) then
            
    	        ninter = 0
    	        n1 = 0
    	        n2 = 0
    	        n3 = 0
      			Wtest=max(dabs(Wdownconv),dabs(Juprad/Iconv(j)),dabs(deltaWc))
      				if (Wtest.gt.(0.05*Wconv(j-1))) then
         				if ((Wconv(j-1)/Wsol) .gt. 1.d-4) then
            				n1=int(Wtest/(0.05d0*Wconv(j-1)))
       	  				else
        	    			n1=10
         				endif
       				endif

  	     		Wtest=max(abs(Juprad/Irad(j)),abs(deltaWr))
   	    			if(Wtest.gt.(0.05d0*Wconv(j-1))) then
   	       				if ((Wrad(j-1)/Wsol) .gt. 1e-4) then
   	        	  			n2=int(Wtest/(0.05*Wrad(j-1)))
   	    	   			else
   		          			n2=10
 	         			endif
       				endif
 	      		if ((tstar(j)-tstar(j-1)) .ge. 0.5e14) then
   		       		n3=int((tstar(j)-tstar(j-1))/0.5e14)
       			endif

       			ninter=max(n1,n2,n3)
       		  
       			if (ninter.gt.999) then
          			ninter=999
    	   		endif
    	   		
    	
       	     	if (ninter .ne. 0) then 
       	     		if (ninter .lt. 400) then
       	     			ninter = 400
       	     		endif      	
    		
       		     	call interpol3(Wdownconv,Juprad,Wconv(j-1),Wrad(j-1),deltaWr,deltaWc,deltaJ, &
            		Winterfin,Winterfinr,Wuprad,Wupconv,itrack,braking_law,flag,taudecfin,ninter)
            		
            
            	
            	!else
            	
            	!	ninter = 10
       		    ! 	call interpol3(Wdownconv,Juprad,Wconv(j-1),Wrad(j-1),deltaWr,deltaWc,deltaJ, &
            	!	Winterfin,Winterfinr,Wuprad,Wupconv,itrack,braking_law,flag,taudecfin,ninter)
            		
           	 	endif
            endif

! We pass here only once, when j = index2+1 = indexdec 
            if (flag .eq. 0 ) then
               Wconv(j)=Wconv(j-1)+Wupconv-Wdownconv-Juprad/Iconv(j)+deltaWc
                              
               Juprad=2./3.*Rrad(j)**2.*Wconv(j)*(Mrad(j)-Mrad(j-1))*msol
     
    			if (Irad(j) .ne. 0.0) then
               	Wrad(j) = Wrad(j-1)+Wuprad+Juprad/Irad(j)-deltaWr  
     			else
               	Wrad(j) = 0.0
            	endif
               
            else
               Wconv(j)=Winterfin
               Wrad(j)=Winterfinr
               taudec(j) = taudecfin


            endif
                        
           ! force the inital core velocity to be the convective velocity using timedec
           if (   (tstar(j) .ge. timedec*year)  .and. (  tstar(j-1) .lt. timedec*year ) ) then
           		Wrad(j) = Wconv(j)
           		write(6,*) "Condition, wrad", tstar(j)/year/1.e6, wrad(j)/wsol, Wuprad,Juprad/Irad(j),deltaWr
           		!stop	
           endif 
           
           
           if(tstar(j) .gt. timedec*year) then
           		if ( Irad(j)/Irad(j-1) .lt. 1e-3) then
           			Wrad(j) = Wconv(j)
           		endif
           
           endif
           

           
            !if ( (Irad(j-1)/Iconv(j-1) .le. 1.e-3) .and. (Irad(j)/Iconv(j) .gt. 1.e-3)) then    
     			!Wrad(j) = Wconv(j)
    		!endif
           

            dt1= tstar(j)/(year*1.0e6)
            dt2= tstar(j-1)/(year*1.0e6)
            
            
           !if (Irad(j)/Iconv(j) .le. 1e-3) then
           		!Wrad(j) = Wconv(j)
           		!write(6,*) "Condition ok, time =", dt1 
           !endif
            	     
            dj(j)=deltaJ
            AMconv(j)=iconv(j)*Wconv(j)
            AMrad(j)=Irad(j)*Wrad(j)
            
            

 112        format(5(f8.5,2x))

! 107        format(f7.4,2x,2(F8.4,2x),9(d10.4,2x))

!     3.4) core-enveloppe angular momentum exchange
!          ----------------------------------------
            if (Wconv(j) .lt. Wrad(j)) then
            	deltaJ = Irad(j)*Iconv(j)/(Irad(j)+Iconv(j))*(Wrad(j)-Wconv(j))
            else if (Wconv(j) .eq. Wrad(j)) then
                 deltaJ = 0.0
            else if (Wconv(j) .gt. Wrad(j)) then
               deltaJ= -Irad(j)*Iconv(j)/(Irad(j)+Iconv(j))*(Wconv(j)-Wrad(j))            
            endif
            

         endif 
         
         

    enddo
         


!      end do !continue with next evolutionary track

	  wconvf = Wconv(2) 
	  wradf = Wrad(2) 
	  DAMtot = abs(AMconv(2) + Amrad(2) - AMconv(1) - Amrad(1))
	  DAMbrake	= abs(Wdownconv * Iconv(j-1))
	  DAMtide = abs(DAMtot - DAMbrake)
	  
	  DAMtot = (AMconv(2) + Amrad(2) - AMconv(1) - Amrad(1))/ (tstar(2)-tstar(1))
	  DAMbrake	= (Wdownconv * Iconv(1) /  (tstar(2)-tstar(1)))
	  DAMtide = (DAMtot - DAMbrake)

	  Wtide = DAMtide *  (tstar(2)-tstar(1))/ Iconv(1)
	  Wjuprad = Juprad/Iconv(2)
	  
	  !Wupconv-Wdownconv-Juprad/Iconv(j)+deltaWc
	  
	  
	  DAMbrake = (AMconv(2) + Amrad(2) - AMconv(1) - Amrad(1))/ ((tstar(2)-tstar(1))*Iconv(2))
	  DAMtide = wconv(2)* (Iconv(2)-Iconv(1))/ ((tstar(2)-tstar(1))*Iconv(2))

 144  format(6(1x,d12.6))

! [ 9(1x,f6.2) stands for ntrack(1x,f6.2) ]
 110  format(1x,f14.8,6(1x,d12.6))
 172  format(2(1x,f10.4),6(1x,d12.6))
 143  format(1x,f10.4,1x,d12.6)
      


      !!write (*,'(1x,a)') char(7)
      return
      end subroutine rotevol
      
      
! linear interpolation routine for rotevoldec.f
!***********************************************************************

      subroutine interpolA(Wdownconv,Wup,W0, & 
      Winterfin,itrack,braking_law,flag)

      common/const/mstar2,msol,rsol,Wsol,taudecinit,mdotsun,alphaS,conS
      common/const2/numtest
      common/var/rstar,tstar,rrad,mrad,Irad,Iconv
      common/var/k2rad,k2conv,Lum,Teff
      common/indexs/indexms,indexdec,j,nn,n,indexsun


      integer numtest,indexdec,indexms,j,n,nn,indexsun,flag
      real*8 mstar2(2),msol,rsol,Wsol,mdotsun,conS
      real*8 rstar(2),tstar(2),rrad(2),mrad(2), & 
     k2rad(2),k2conv(2),Irad(2),Iconv(2)
      real*8 W0,Wtest,Wdownconv,Winterfin,Wup,Wi(20000),kic(20000), & 
     ti(20000),kir(20000),ri(20000),ici(20000),iri(20000)
      real*8 taudeci(20000),taudecfin
      integer k,itrack,ninter
      character*2 braking_law
      real*8 xx,yy,zz,tt,uu,sum1,sum2,year,alphaS
      REAL*8 Lum(2), Li(20000),Teff(2),Teffi(20000)
      
      real(dp) ::  taudecinit

      
      Wtest=max(dabs(Wdownconv),dabs(Wup))

      if (Wtest.gt.(0.1*W0)) then
         flag = 1
         
         !!write(6,666) itrack,n, W0/Wsol,Wup/Wsol, Wdownconv/Wsol
 666     format(1x,'interpolA: track=',i2,2x,i4, ' Wo=',d10.4,' Wup=' & 
        ,d10.4,' Wdown=',d10.4)
!         ninter=10
         if ((W0/Wsol) .gt. 1d-4) then
             ninter=int(Wtest/(0.1*W0))
         else
            ninter=10
         endif
         if (ninter.gt. 1000) then 
            ninter=1000
         endif
         
       !Nécessaire sinon pas de temps trop long par rapport à couplage
       if (ninter .lt. dabs((tstar(j)-tstar(j-1))/taudec)) then
       		ninter = 4*int((tstar(j)-tstar(j-1))/taudec)
       		write(6,*) "InterpolA, ninter optimum = ", ninter
       		if (ninter .gt. 5000) then
       			ninter = 5000
      		endif
       endif 
                  
         !print *,tstar(j),ninter
         do k=1,ninter+1
            ti(k)=tstar(j-1)+(k-1.)*(tstar(j)-tstar(j-1))/ninter
            ri(k)=rstar(j-1)+(k-1.)*(rstar(j)-rstar(j-1))/ninter
            Li(k)=Lum(j-1)+(k-1)*(Lum(j)-Lum(j-1))/ninter
            Teffi(k)=Teff(j-1)+(k-1)*(Teff(j)-Teff(j-1))/ninter
            kir(k)=k2rad(j-1)+(k-1.)*(k2rad(j)-k2rad(j-1))/ninter
            kic(k)=k2conv(j-1)+(k-1.)*(k2conv(j)-k2conv(j-1))/ninter
            Iri(k)=kir(k)*mstar2(j-1)*msol*(ri(k))**2
            Ici(k)=kic(k)*mstar2(j-1)*msol*(ri(k))**2
         enddo
      
         Wi(1)=W0
         sum1=0
         sum2=0
         xx=0
         yy=0
         zz=0
         tt=0
         uu=0
         year=3600*24*365
         do k=2,ninter+1
            Wup = Wi(k-1) * (Ici(k-1)/Ici(k) - 1.)
            call magnbrak(Wdownconv,Wi(k-1),ri(k-1),Li(k-1), & 
           ti(k),ti(k-1),ici(k-1),itrack,Teffi(k-1))
            Wi(k)=Wi(k-1)+Wup-Wdownconv
            sum1=sum1+Wdownconv
            sum2=sum2+Wup
         enddo
         Winterfin=Wi(ninter+1)
         wdownconv = sum1
         
		write(6,*) "InterpolA condition OK", W0/Wsol, Winterfin/Wsol, tstar(j-1)/yr/1.0e6


 108        format(f7.4,2x,2(F8.3,2x),7(d10.4,2x))

      end if
      
      return
      end subroutine interpolA
!***********************************************************************

      subroutine interpol3(Wdownconv,Juprad,Wc0,Wr0,dWr, & 
     dWc,dJ,Winterfin,Winterfinr,Wur,Wuc,itrack, & 
     braking_law,flag,taudecfin,ninter)

      common/const/mstar2,msol,rsol,Wsol,taudecinit,mdotsun,alphaS,conS
      common/const2/numtest
      common/var/rstar,tstar,rrad,mrad,Irad,Iconv
      common/var/k2rad,k2conv,Lum,Teff
      common/indexs/indexms,indexdec,j,nn,n,indexsun

      integer numtest,indexms,indexdec,nn,n,indexsun,flag
      real*8 mstar2(2),msol,rsol,Wsol,mdotsun
      real*8 Rstar(2),tstar(2),Rrad(2),Mrad(2), & 
     k2rad(2),k2conv(2),Irad(2),Iconv(2)
      real*8 Wdownconv,Wc0,Wr0,dWr,dWc,dJ,Juprad,mdotsat(20000)
      real*8 Winterfin,Winterfinr,Wur,Wuc,kic(20000), & 
     ti(20000),kir(20000),ri(20000),rri(20000),mri(20000), & 
     Ici(20000),Iri(20000),Wtest,Wci(20000),Wri(20000)
      character*2 braking_law
      integer itrack,j,k,ninter,n1,n2,n3
      real*8 sum1,sum2,sum3,sum4,sum5,sum6,year
      REAL*8 Lum(2), Li(20000),Teff(2),Teffi(20000)
      real*8 taudeci(20000),taudecfin,AMconv0,DAMconv,alphaS,conS

	  real(dp) ::  taudecinit

      REAL*8 grav,G,M

      G = 6.6732d-8
      M = mstar2(j-1)*msol

      


!      Wtest=max(dabs(Wdownconv),dabs(Juprad/Iconv(j)),dabs(dWc))
!      if (Wtest.gt.(0.05*Wc0)) then
!         if ((Wc0/Wsol) .gt. 1.d-4) then
!            n1=int(Wtest/(0.05d0*Wc0))
!         else
!            n1=10
!         endif
!       endif

!       Wtest=max(abs(Juprad/Irad(j)),abs(dWr))
!       if(Wtest.gt.(0.05d0*Wc0)) then
!          if ((Wr0/Wsol) .gt. 1e-4) then
!             n2=int(Wtest/(0.05*Wr0))
!          else
!             n2=10
!          endif
!       endif
!       if ((tstar(j)-tstar(j-1)) .ge. 0.5e14) then
!          n3=int((tstar(j)-tstar(j-1))/0.5e14)
!       endif

!       ninter=max(n1,n2,n3)
              
!       if (ninter.gt.999) then
!          ninter=999
!       endif
              
       !Nécessaire sinon pas de temps trop long par rapport au couplage
       if (ninter .lt. dabs((tstar(j)-tstar(j-1))/taudecinit)) then
       ninter = 4*int((tstar(j)-tstar(j-1))/taudecinit)
       		if (ninter .gt. 3000) then
       		   ninter = 3000
             endif
       endif 
       
	   ninter = 5000
       
       if (ninter .ne. 0) then
          flag=1
          

          do k=1,ninter+1
             ti(k)=tstar(j-1)+(k-1)*(tstar(j)-tstar(j-1))/ninter
             ri(k)=rstar(j-1)+(k-1)*(rstar(j)-rstar(j-1))/ninter
             Li(k)=Lum(j-1)+(k-1)*(Lum(j)-Lum(j-1))/ninter
             Teffi(k)=Teff(j-1)+(k-1)*(Teff(j)-Teff(j-1))/ninter
             grav = M*G/(ri(k)**2.)
             kir(k)=k2rad(j-1)+(k-1.)*(k2rad(j)-k2rad(j-1))/ninter
             kic(k)=k2conv(j-1)+(k-1.)*(k2conv(j)-k2conv(j-1))/ninter
             rri(k)=rrad(j-1)+(k-1)*(rrad(j)-rrad(j-1))/ninter
             mri(k)=mrad(j-1)+(k-1)*(mrad(j)-mrad(j-1))/ninter
             Iri(k)=kir(k)*mstar2(j-1)*msol*(ri(k))**2.
             Ici(k)=kic(k)*mstar2(j-1)*msol*(ri(k))**2.
          enddo
 13       continue


 667      format(1x,'interpol3: j=',i3,i4 ,' Wco=',e10.4,' Wro=',e10.4, & 
        ' Wup=',e10.4,' Wdown=',e10.4,' dWc=',e10.4,' dWr=',e10.4, & 
        ' Wuprad=',e10.4,' Wupconv=',e10.4)

          Wci(1)=Wc0
          Wri(1)=Wr0
          sum1=0
          sum2=0
          sum3=0
          sum4=0
          sum5=0
          sum6=0
         year=3600*24*365
         

          do k=2,ninter+1
             wuc=wci(k-1)*(Ici(k-1)/Ici(k)-1)
        if (iri(k) .ne. 0.0) then
             wur=wri(k-1)*(Iri(k-1)/Iri(k)-1)
        else
        	wur = 0.0
        endif
        
             call magnbrak(Wdownconv, & 
         Wci(k-1),ri(k-1),Li(k-1),ti(k),ti(k-1),ici(k-1),itrack, & 
         Teffi(k-1))
 
		

! F.G 6/10/2014
! Taudec Spada
         if (wri(k-1) .eq. wci(k-1)) then
            taudeci(k) = taudecinit*10000. !si vitesse egale : couplage faible?
         else
            taudeci(k)=taudecinit* & 
      ((conS*wsol)/abs(wri(k-1)-wci(k-1)))**alphaS
        endif
        
          if ( (taudeci(k)/year) .gt. 5000e6) then
          
           taudeci(k) = 5000e6 * year
          
          endif
          

        
            taudeci(k) = taudecinit
          

             Juprad=2./3.*rri(k)**2.*Wci(k-1)*(mri(k)-mri(k-1))*msol
             dWc=dJ*(ti(k)-ti(k-1))/(taudeci(k)*ici(k))
        if (iri(k) .ne. 0.0) then
        	dWr=dJ*(ti(k)-ti(k-1))/(taudeci(k)*iri(k))
        else
        	dWr = 0.0
        endif	

             Wci(k)=Wci(k-1)-Wdownconv+Wuc-Juprad/Ici(k)+dWc
 116         format(6(f8.5,2x))
             Juprad=2./3.*rri(k)**2.*Wci(k)*(mri(k)-mri(k-1))*msol
    		if (Irad(j) .ne. 0.0) then             
             Wri(k)=Wri(k-1)+Wur+Juprad/Iri(k)-dWr
            else
             Wri(k)= 0.0
            endif
             if (Wci(k) .lt. Wri(k)) then
             	dJ=iri(k)*ici(k)/(iri(k)+ici(k))*(Wri(k)-Wci(k))
             
             else if (Wci(k) .eq. Wri(k)) then
                dJ=0.d0
             else if (Wci(k) .gt. Wri(k)) then
                dJ=-iri(k)*ici(k)/(iri(k)+ici(k))*(Wci(k)-Wri(k))
             endif
             sum1=sum1+wur
             sum2=sum2+wuc
             sum3=sum3+wdownconv
             sum4=sum4+juprad/iri(k)
             sum5=sum5+dWr
             sum6=sum6+dWc
             
             
          enddo
 114      format(5(f8.5,2x))         
          Winterfin=Wci(ninter+1)
          Winterfinr=Wri(ninter+1)
          taudecfin = taudeci(ninter+1)
          wdownconv = sum3

          
 108         format(f7.4,2x,2(F8.3,2x),7(d10.4,2x))

       endif
       
      return
      end subroutine interpol3

!***********************************************************************

      subroutine interpolB(Juprad,dJ,dWr,Wur,Wr0,Wc0,Winterfinr, & 
     itrack,flag,taudecfin)

      common/const/mstar2,msol,rsol,Wsol,taudecinit,mdotsun,alphaS,conS
      common/const2/numtest
      common/var/rstar,tstar,rrad,mrad,Irad,Iconv
      common/var/k2rad,k2conv,Lum,Teff
      common/indexs/indexms,indexdec,j,nn,n,indexsun

      integer indexms,indexdec,nn,n,indexsun,flag,numtest
      real*8 mstar2(2),msol,rsol,Wsol,mdotsun
      real*8 rstar(2),tstar(2),rrad(2),mrad(2), & 
     k2rad(2),k2conv(2),Irad(2),Iconv(2)
      real*8 Wc0,Wr0,dWr,dJ,Juprad,Wur
      real*8 Winterfinr,kir(20000),kic(20000)
      real*8 ti(20000),rri(20000),mri(20000),ri(20000)
      real*8 Ici(20000),Iri(20000),Wtest,Wri(20000),Lum(2)
      real*8 Teff(2),taudeci(20000),taudecfin
      integer itrack,j,k,ninter
      real*8 xx,yy,zz,sum1,sum2,sum3,year,alphaS,conS
      
      real(dp) ::  taudecinit


      Wtest=max(abs(Juprad/Irad(j)),abs(dWr),abs(Wur))
      
      Wtest = 0.5*Wr0
      
      if (Wtest.ge.(0.1*Wr0)) then
         flag=1
         if ((Wr0/Wsol) .gt. 1e-4) then
            ninter=int(Wtest/(0.1*Wr0))
         else
            ninter=10
         endif
         if (ninter.gt. 999) then
            ninter=999
         endif
         
         
        ninter = 200 

 668     format(1x,'interpolB: j=',i3,2x, ' Wro=',e10.4,' Wup=',e10.4, & 
        'Juprad=',e10.4,' dWr=',e10.4)

         do k=1,ninter+1
            ti(k)=tstar(j-1)+(k-1)*(tstar(j)-tstar(j-1))/ninter 
            ri(k)=rstar(j-1)+(k-1)*(rstar(j)-rstar(j-1))/ninter
            kir(k)=k2rad(j-1)+(k-1)*(k2rad(j)-k2rad(j-1))/ninter
            kic(k)=k2conv(j-1)+(k-1)*(k2conv(j)-k2conv(j-1))/ninter
            rri(k)=rrad(j-1)+(k-1)*(rrad(j)-rrad(j-1))/ninter
            mri(k)=mrad(j-1)+(k-1)*(mrad(j)-mrad(j-1))/ninter
            Iri(k)=kir(k)*mstar2(j-1)*msol*(ri(k))**2.
            Ici(k)=kic(k)*mstar2(j-1)*msol*(ri(k))**2.
         enddo

         Wri(1)=Wr0
         sum1=0
         sum2=0
         xx=0
         yy=0
         zz=0
         sum3=0
         year=3600.*24.*365.
         do k=2,ninter+1
         
! F.G 6/10/2014
! Taudec Spada
            if (Wri(k-1) .eq. Wc0) then
            	taudeci(k) = taudecinit*10000. !si vitesse egale : couplage faible?
            else
            	taudeci(k)=taudecinit*((conS*wsol)/abs(Wri(k-1)-Wc0))**alphaS
            endif
            
          
          	if ( (taudeci(k)/year) .gt. 5000e6) then
          
           		taudeci(k) = 5000e6 * year
          
          	endif
                    
            taudeci(k) = taudecinit  
      
            wur=wri(k-1)*(Iri(k-1)/Iri(k)-1.)
            Juprad=2./3.*rri(k)**2.*Wc0*(mri(k)-mri(k-1))*msol
            dWr=dJ*(ti(k)-ti(k-1))/(taudeci(k)*iri(k))
            Wri(k)=Wri(k-1)+Juprad/Iri(k)-dWr+Wur
            
            if (Wc0 .lt. Wri(k)) then
            	dJ=iri(k)*ici(k)/(iri(k)+ici(k))*(Wri(k)-Wc0)
            else if (Wc0 .eq. Wri(k)) then
            	dJ=0.d0
            else
               dJ=-iri(k)*ici(k)/(iri(k)+ici(k))*(Wc0-Wri(k))
            endif
            sum1=sum1+juprad/iri(k)
            sum2=sum2+dWr
            sum3=sum3+wur
         enddo

         Winterfinr=Wri(ninter+1)
         taudecfin = taudeci(ninter+1)
!         if (itrack.eq.numtest) then
!            !write(20,108)log10(tstar(j)/year),Wc0/Wsol,
!     *           Winterfinr/Wsol,sum3/Wsol
!     *           ,xx,yy,sum1/Wsol,sum2/Wsol,zz
!         endif
 108        format(f7.4,2x,2(F8.3,2x),7(d10.4,2x))

      endif
      return
      end subroutine interpolB
      
      
!***********************************************************************

      subroutine interpolAM(t1,t2,AMconv1,AMconv2,AMrad1,AMrad2)

      common/const/mstar2,msol,rsol,Wsol,taudecinit,mdotsun,alphaS,conS
      common/const2/numtest
      common/var/rstar,tstar,rrad,mrad,Irad,Iconv
      common/var/k2rad,k2conv,Lum,Teff
      common/indexs/indexms,indexdec,j,nn,n,indexsun

      integer numtest,indexms,indexdec,nn,n,indexsun,flag
      real*8 mstar2(2),msol,rsol,Wsol,mdotsun
      real*8 rstar(2),tstar(2),rrad(2),mrad(2), & 
     k2rad(2),k2conv(2),Irad(2),Iconv(2)
      real*8 Wdownconv,Wc0,Wr0,dWr,dWc,dJ,Juprad,mdotsat(20000)
      real*8 Winterfin,Winterfinr,Wur,Wuc,kic(20000), & 
     ti(20000),kir(20000),ri(20000),rri(20000),mri(20000), & 
     Ici(20000),Iri(20000),Wtest,Wci(20000),Wri(20000)
      character*2 braking_law
      integer itrack,j,k,ninter,n1,n2,n3
      real*8 sum1,sum2,sum3,sum4,sum5,sum6,year
      REAL*8 Lum(2), Li(20000),Teff(2),Teffi(20000)
      real*8 taudeci(20000),taudecfin
      REAL*8 alphaS,conS
      
      
      REAL*8 critere,t1,t2,AMconv1,AMconv2,AMrad1,AMrad2,dt
      REAL*8 AMci(20000),AMri(20000)
      
      real(dp) ::  taudecinit
      
      !???
      REAL*8 DAMconv(20000)


      REAL*8 grav,G,M

		critere = 0.05
		
		
		ninter = ((t2-t1)/t1) * ((1-critere)/critere)
		
		ninter = 2000		

          do k=1,ninter+1
             ti(k)=tstar(j-1)+(k-1)*(tstar(j)-tstar(j-1))/ninter
             AMci(k) = AMconv1+(k-1)*(AMconv2-AMconv1)/ninter
             AMri(k) = AMrad1+(k-1)*(AMrad2-AMrad1)/ninter 
          enddo


         year=3600*24*365

          do k=2,ninter+1

            DAMconv(k) = abs(AMci(k)+AMri(k)-AMci(k-1)-AMri(k-1))  /     & 
 (ti(k)-ti(k-1))
     
     		dt = ti(k) /(year*1.0e6) 
     		!write(26,*) dt,DAMconv(k),"inter"
               
             
          enddo

      return
      end subroutine interpolAM
      
      
! Magnetic braking calculation subroutine
!****************************************

      subroutine magnbrak(Wdownconv,W0,r,L,t1,t2,Ic,itrack,Teff)
      
      implicit none

      common/const/mstar2,msol,rsol,Wsol,taudecinit,mdotsun,alphaS,conS
      common/const2/numtest
      common/freinage/ksk,kmm,ksc,kmp,ksat,wcrit,wsat
      common/freinage/G,mdotstar,flag2,Bstar,Ro,ff,mdotstarC,mdotstarJ
      common/indexs/indexms,indexdec,j,nn,n,indexsun


      integer numtest,indexms,indexdec,j,nn,n,itrack,indexsun
      !integer brklaw
      
      real(dp) ::  taudecinit
      !real*8 K3,K4,mvic,,K,K1MP,K2MP,m,pi
      real*8 ksk,kmm,ksc,kmp,ksat,a,b
      real*8 mstar2(2),wcrit,wsat,msol,rsol,Wsol,taudec
      real*8 t1,t2,r,Wdownconv,W0,Ic,G,mdotsun,mdotstar,Bstar,L,ff
      real qq
      REAL*8 Wcr,C,flag2,Teff
      REAL*8 omegacra, mdotcra , dumm , year , dt, Period , mdotyr
      REAL*8 Ro
      CHARACTER*100 command
      CHARACTER HEAD*72
      REAL*8 mdotstarC,mdotstarJ,conS,Wsun,alphaS
      REAL*8 omegasatJ,WdownconvJ,Wdownconv1,Wdownconv2,Wdownconv3
      REAL*8 tcz,Rosun,tczsun,T0,factor,gamma,ratio
      
      
      
      !real(dp) :: chi,p 

      !pi = 3.14159
      !year = 365.0*2.4d1*36.e2
      mdotyr = 6.30276d+25


         mdotstar = 0.0
         bstar = 0.0
        call mdot(W0,r,L,mdotstar,bstar,Ro,ff,Teff)
        
        
        if (brklaw .eq. 0) then
        
        
        mdotstarC = mdotstar
        


!------------------------Mdot/Braking Johnstone-------------------------

       omegasatJ = 15 * wsol * (mstar2(j-1)/1)**2.3


       if (W0 .ge. omegasatJ) then
        mdotstarJ = (mdotsun * (r/rsol)**2. * (omegasatJ/wsol)**1.33) / (mstar2(j-1))**3.36
     
        WdownconvJ = 7.15e30*(15)**1.89*(W0/wsol)*(mstar2(j-1))**4.42 * ((t1-t2)/Ic) * K1MP
     
     
        else 
         mdotstarJ = (mdotsun * (r/rsol)**2. * (W0/wsol)**1.33) / (mstar2(j-1))**3.36
        WdownconvJ = 7.15e30*(W0/wsol)**2.89*((t1-t2)/Ic)*K1MP
       endif
        

        
        
        mdotstar = mdotstarC


!------------Mass loss Johnstone------------
!
!          mdotstar = mdotstarJ
!
!-------------------------------------------




!-----------------------------------------------------------------------


        
!----Test pour voir quel Mdot il faut si K_1 = 1.7 pour faible masse----        
!        mdotstar = mdotstar * 8.
!-----------------------------------------------------------------------

!------------Utilisation of Mdot and Bstar from Cranmer 2011------------

          Wdownconv = K1MP**2.*kmp*Bstar**(4.*mmp)*  & 
          r**(5.*mmp+2.)*w0**(1.)*(mdotstar)**(1.-2.*mmp) * & 
          (t1-t2)/(Ic *(K2MP**2. *2.*G*mstar2(j-1)*msol + Kbr*w0**(2.) * & 
          r**(3.))**(mmp))
     
! Avec facteur 2/3 pour la vitesse de liberation! Facteur 3**m 
!          Wdownconv = K1MP**2.*kmp*Bstar**(4.*m)* 
!     *        r**(5.*m+2.)*w0**(1.)*(mdotstar)**(1.-2.*m) *
!     *        (t1-t2)/(Ic *(K2MP**2. *2./3.*G*mstar2(j-1)*msol + K*w0**(2.) *
!     *        r**(3.))**(m))
     

     
!------------Freinage Johnstone------------
!          
!          Wdownconv = WdownconvJ
!
!------------------------------------------


	   else if (brklaw .eq. 1) then
!     Reville et al. (2015)	   
	      Wdownconv = K3**2.*Bstar**(4.*mvic) *  & 
       r**(5.*mvic+2.)*w0**(1.)*(mdotstar)**(1.-2.*mvic) * & 
          (t1-t2)/ (Ic *(2.*G*mstar2(j-1)*msol + (2/K4**2.)*w0**(2.) * & 
          r**(3.))**(mvic)) * & 
        (4.*pi)**(4.*mvic)
     
         
	   else if (brklaw .eq. 2) then
	   
	   
!	        if (t1/year .lt. 4.e7) then
	        	K3 = 2.  
	            mvic = 0.235
	      Wdownconv1 = K3**2.*Bstar**(4.*mvic)  & 
       *r**(5.*mvic+2.)*w0**(1.)*(mdotstar)**(1.-2.*mvic) * & 
          (t1-t2)/ (Ic *(2.*G*mstar2(j-1)*msol + (2/K4**2.)*w0**(2.) * & 
          r**(3.))**(mvic))
	            
	            
	            
!	        else if (t1/year .lt. 2.e8) then
	            K3 = 1.7
	            mvic = 0.15
	      Wdownconv2 = K3**2.*Bstar**(4.*mvic)  & 
       *r**(5.*mvic+2.)*w0**(1.)*(mdotstar)**(1.-2.*mvic) * & 
          (t1-t2)/ (Ic *(2.*G*mstar2(j-1)*msol + (2/K4**2.)*w0**(2.) * & 
          r**(3.))**(mvic))
     	        
!	        else 
	            K3 = 1.7
	            mvic = 0.11
	      Wdownconv3 = K3**2.*Bstar**(4.*mvic)  & 
       *r**(5.*mvic+2.)*w0**(1.)*(mdotstar)**(1.-2.*mvic) * & 
          (t1-t2)/ (Ic *(2.*G*mstar2(j-1)*msol + (2/K4**2.)*w0**(2.) * & 
          r**(3.))**(mvic))	            
!	        endif            
	   
!	      Wdownconv = K3**2.*Bstar**(4.*mvic) 
!     *     *r**(5.*mvic+2.)*w0**(1.)*(mdotstar)**(1.-2.*mvic) *
!     *        (t1-t2)/ (Ic *(2.*G*mstar2(j-1)*msol + (2/K4**2.)*w0**(2.) *
!     *        r**(3.))**(mvic))

          Wdownconv =Wdownconv1*2.0 + Wdownconv2*2.0 + Wdownconv3*1.
     
        
        else if (brklaw .eq. 3) then
!      Matt et al. 2015

            tcz = 314.24*exp(-(Teff/1952.5)-(Teff/6250.)**18.)+ 0.002
            
            Rosun =  1.96
!            chi = 10.0
            Wsun = 2.87e-6
            tczsun = 12.9
!            factor = 0.16 !slow
            factor = 0.16
                        
            ratio = W0*(r)**(3/2)*(G*mstar2(j-1)*msol)**(-0.5)
            gamma = (1+(ratio/0.072)**2)**0.5 
        	T0 = factor* 5.0e31*(r/rsol)**3.1 * mstar2(j-1)**0.5*gamma**(2*mmp)
!        	p = 2.55 !slow
!            p = 2.
            	if (Ro .le. (Rosun/chi) ) then
            		
            		Wdownconv = (T0 * chi**p * (W0/Wsun))*((t1-t2)/Ic)
            	else 
            		Wdownconv = (T0*(tcz/tczsun)**p * (W0/Wsun)**(p+1)) & 
         			*((t1-t2)/Ic)
            	endif	
       
       
       


!-----------------------------------------------------------------------
!
!
!------------------Reiners & Mohanty 2012 braking law-------------------

!        if (w0.lt.Wsat) then
!          Wdownconv = K1MP**2.*kmp*3**(4.*m)*Wsol**(-b*4.*m)* 
!     *        r**(5.*m+2.)*w0**(1.+b*4.*m) * (mdotstar)**(1.-2.*m) *
!     *        (t1-t2)/(Ic *(K2MP**2. *2.*G*mstar2(j-1)*msol + K*w0**(2.) *
!     *        r**(3.))**(m))
!        else
!          Wdownconv = K1MP**(2.)*3**(4.*m)*ksat* Wsol**(-b*4.*m)*
!     *        r**(5.*m+2.)*w0 * (mdotstar)**(1.-2.*m) * 
!     *        (t1-t2) / (Ic * (K2MP**2. *2.*G*mstar2(j-1)*msol + 
!     *        K*wsat**(2.)*r**(3.))**(m))
!        end if   

!-----------------------------------------------------------------------


! ---------------------------------------------
! Add of the Reiners & Mohanty 2012 braking law
      else if (brklaw .eq. 4) then

        if (flag2 .eq. 0.) then
          !write (6,*) 'Reiners!'
          flag2 = 1.
        endif
       
        Wcr = 8.56d-6  ! Wcrit from Reiners & Mohanty  
        C = 2.66d3
        if (w0.lt.Wcr) then
          Wdownconv = C*((w0/Wcr)**4*w0*(r**16/(mstar2(j-1)*msol)**2)**(1/3))* & 
          ((t1-t2)/Ic)   
        else
          Wdownconv = C*(w0*(r**16/(mstar2(j-1)*msol)**2)**(1/3))*((t1-t2)/Ic)    
        end if
!----------------------------------------------------------

      else 
         !!write(6,1005) itrack,j
 1005    format('unknown braking-law for itrack =',I2,' and step=',I3)
         !print *,'Wconv=',W0/wsol,'  Wcrit=',Wcrit/wsol
         close(21)
         stop
      end if

!-------------No braking case-------------     
!          Wdownconv = 0.0
!_________________________________________
     


      qq=j/4.
      if ((itrack.eq.1.or.itrack.eq.nn).and.(j.gt.indexms.and. & 
       j.lt.indexsun).and.(qq.eq.aint(qq).or.abs(w0/wcrit).le.0.05 & 
       .or.abs(w0/wsat).le.0.05)) then
         !!write (21,*)w0/wsol,log10(wdownconv*Ic/(t1-t2))
      endif
      return
      end subroutine magnbrak
      
      
! Provides the mass loss rate and the mean magnetic field as a function
! of the stellar parameters

      subroutine mdot(wstar,rstar,lstar,mdotstar,bstarf,Ro,fstar,Teff)
      
      
      implicit none

      common/const/mstar2,msol,rsol,Wsol,taudecinit,mdotsun,alphaS,conS

      common/const2/numtest
!      common/freinage/G,mdotstar,flag2,Bstar

      REAL*8 taudec,mdotsun,flag2,Bstarf
      REAL*8 pi,rstar,mstar2(2),msol,rsol,fstar,period,grav,lgrav,Teff,G,M,L
      REAL*8 wsol,vesc,vesc2,day2,year,tauc,Ro,Rosol,x,kb,stefan
      REAL*8 fstarTR,lstar,lsol,F0,T0,eps,FA,theta,wstar
      real(dp) ::  taudecinit
      REAL*8 mdotstar,mh,mu,lammax,lammaxsol,ratioZ,rhostar,Qstar,QTR
      REAL*8 lamperstar,lampersol,Hsol,alpha,alphaTR,ratioalpha,Bstar
      REAL*8 VA,vper,FluxhTR,Fluxcond,h,crad,PTR,rhoTR,sqbrac,alphaS
      REAL*8 Cons
      REAL*8 rho(13000),Temp(13000),gravit(13000),a,b
      REAL*8 rhotemp, rhogravit,test1,test2,gravsol,gratio,test3,test4
      REAL*8 mdotstarhot, mdotstarcold
      REAL*8 rcrit , ucrit , rhocrit , Bcrit , VAcrit , MAcrit , macfac
      REAL*8 vpercrit,FATR, BTR , uinf , TTR, Reflcoef , Reflcoefnew
      REAL*8 VATR , fmin , fmax , fmod, wsatnorm,mdotsat
      REAL*8 C0,C1,xteff !use to calculate the density
      REAL*8 uTR_hot,MATR

      INTEGER j,n,numtest,flag
      


      pi = 3.14159
      G = 6.6732d-8
      day2 = 2.4d1*36.e2
      year = 365.0*day2
      Rosol = 1.96
      kb  = 1.380622d-16
      stefan = 5.66961d-5
      lsol = 3.826d+33
      mh = 1.67333d-24
      lampersol = 3.d7
      lammaxsol = 7.4d-23 + 4.2d-22
      Hsol = 1.39d7
      h = 0.5  ! h = [0.5-1.5]


      M = mstar2(1)*msol
      L = lstar*lsol
      grav = M*G/(rstar**2.)
      gravsol = msol*G/(rsol**2.)
      lgrav = log10(grav)
      !Teff = (L/(4.*pi*rstar**2.*stefan))**0.25
      theta = 1./3.
      ratioZ = 1.
      alpha = 0.5
      ratioalpha = 0.9997794  !comes from Boreas
      alphaTR = ratioalpha * alpha
      vesc = (2.*G*M/rstar)**0.5
      vesc2 = (2.*G*M/rstar)

!     OK for   3300 < Teff < 7000 K
      tauc = 314.24*exp(-(Teff/1952.5)-(Teff/6250.)**18.)+ 0.002
      gratio = gravsol/grav
      if (gratio .gt. 1.) then
        tauc = tauc * gratio**0.18
      endif

!      tauc = -250.* log10(Teff) + 955.

      Period = 2.*pi/(wstar*day2) 
      Ro = Period/tauc
      x = Ro/Rosol
      fmin = 0.5 /(1.+(x/0.16)**2.6)**1.3
      fmax =  1. / (1. + (x/0.31)**2.5)

!     -----1 Msol----
      wsatnorm = 0.16
      fmod = 0.55 /(1.+(x/wsatnorm)**2.3)**1.22
!     ---------------

!     -----0.8 Msol----
!      wsatnorm = 0.25 
!      fmod = 0.55 /(1.+(x/wsatnorm)**2.15)**1.
!     ---------------

!     -----1.3 Msol----
!      wsatnorm = 0.3
!      fmod = 0.55 /(1.+(x/wsatnorm)**2.9)**1.5
!     ---------------


      fstar = fmod
      fstarTR = fstar**theta


!-----------------------------------------------------------------------
!     ********************************************************
!     * Calculation of the density by using an approximation *
!     *                  provided by Cranmer                 *
!     ********************************************************
       xteff = Teff / 1000.
       C0=4.0872049-48.979961*xteff+47.135345*(xteff**2.)- & 
        20.204336*(xteff**3.)+4.4110828*(xteff**4.)- & 
        0.48112223*(xteff**5.)+0.020825121*(xteff**6.)
       C1=24.562656-25.713078*xteff+9.8731818*(xteff**2.)- & 
        1.3825986*(xteff**3.)-0.055190235*(xteff**4.) & 
        + 0.031616983*(xteff**5.)-0.0021585992*(xteff**6.)

      rhostar =10.**(C0+C1*lgrav)
!-----------------------------------------------------------------------



!
!------------------------Calculation of Mdot_hot------------------------
!                        ***********************

      F0 = 5.724*exp(-lgrav/11.48)*1.d9
      T0 = 1000.*(5.624+0.6002*lgrav)
      eps = 6.774+0.5057*lgrav
      FA = F0*(Teff/T0)**eps *exp(-(Teff/T0)**25.) 
!-----with alpha = 1.5 and Bstar/Beq = 1.13 => FA is about 5 times------ 
!------------------------greater than it should-------------------------
!----------------------------FA = FA / 4.9566---------------------------
!      FA = FA / 4.9566
      FA = FA / 2.5   ! Best fit for 1Msol

      lammax = 1.d-23 * (7.4+42.*ratioZ**1.13)

      mu = 7./4. + 0.5*tanh((3500.-Teff)/600.)

      lamperstar = (lampersol/Hsol) * kb*Teff/(mu*mh*grav)
      Bstar = 1.13 * sqrt(8.*pi*rhostar*kb*Teff/(mu*mh))
      Bstarf = Bstar*fstar
      VA = Bstar / (4.*pi*rhostar)**0.5
      vper = sqrt(FA/(rhostar*VA))

      Qstar = alpha * rhostar * vper**3. / lamperstar

      BTR = Bstar * fstar / fstarTR
      uinf = vesc
      TTR = 2.0d5
      Reflcoef = 0.5

      do j=1,50
      ratioalpha = ReflCoef*(1.+ReflCoef)/(1.+ReflCoef**2)**1.5*sqrt(2.) 
      sqbrac = ratioalpha*Qstar*mh**2./ (rhostar**0.25 * lammax)
      rhoTR = sqbrac**(4./7.) * fstar**(2.*(1.-theta)/7.)
      QTR = Qstar*ratioalpha*((rhoTR/rhostar)**0.25)*sqrt(BTR/Bstar)
      VATR = BTR / sqrt(4.*pi*rhoTR)
      PTR = rhoTR * kb * TTR / (mu*mh)
      Reflcoefnew = abs((VATR-uinf)/(VATR+uinf))      
      Reflcoef = sqrt(Reflcoefnew*Reflcoef) 
      enddo
      
      FluxhTR = QTR*rstar*h
      FATR = FA * fstar/fstarTR

      if (FluxhTR .gt. FATR) then
        FluxhTR = FATR
      endif

      crad = 1.4d6 * sqrt(lammax/lammaxsol)
      Fluxcond = crad * PTR

	  if (Fluxcond .gt. 0.9) Fluxcond = 0.9
	  
      mdotstarhot = (4.*pi*rstar**2.*fstarTR/vesc2)*(FluxhTR - Fluxcond)
    
	  uTR_hot = mdotstarhot / (rhoTR*4.*pi*rstar**2.*fstarTR)
	  
      MATR  = uTR_hot / VATR	
 	
!
!------------------------Calculation of Mdot_cold-----------------------
!                        ************************

      rcrit = rstar * 7. / (4.*(1. + (vper/vesc2)**2.))
      ucrit = (G*M/(2.*rcrit))**0.5
 
      Bcrit = Bstar * (rstar/rcrit)**2. * fstar  

      vpercrit = 2.*ucrit 
      rhocrit = 4.*pi*(rhostar*vper**2.*VA*fstar*4.*pi*rstar**2./ & 
                            (vpercrit**2.*Bcrit*4.*pi*rcrit**2.))**2.

      do j=1,50
        VAcrit = Bcrit / (4.*pi*rhocrit)**0.5
        MAcrit = ucrit / VAcrit
        macfac = (1. + 3.*MAcrit)/(1. + MAcrit)
        vpercrit = 2.*ucrit / sqrt(macfac)
        rhocrit = rhostar*vper**2.*VA*fstar*rstar**2./(vpercrit**2. & 
                                     *VAcrit*rcrit**2.*(1.+MAcrit)*2.)
      enddo


      mdotstarcold = 4.*pi*rcrit**2.*ucrit*rhocrit

!      mdotstar = mdotstarhot + mdotstarcold
      
      mdotstar = mdotstarhot*exp(-4.*MATR**2) + mdotstarcold
      

 
  
      close(10)
      return
      end subroutine mdot
      
      
      SUBROUTINE SPLINE2(X,Y,N,Y2)

      PARAMETER (NMAX=100)

      real*8 X(N),Y(N),Y2(N),U(NMAX),sig,p,qn,un
      integer n,i,KS

      Y2(1)=0.
      U(1)=0.
      DO  I=2,N-1
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
        P=SIG*Y2(I-1)+2.
        Y2(I)=(SIG-1.)/P
        U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1)) / & 
        	(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
      enddo

      QN=0.
      UN=0.

      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.)
      DO KS=N-1,1,-1
        Y2(KS)=Y2(KS)*Y2(KS+1)+U(KS)
      enddo
      
      RETURN
      END subroutine SPLINE2 


end module rotevoldec  
