! Hansen 2010

module fonctions

use constants

implicit none

contains 

! The 2 eccentricity dependant factors in the equation in a 
subroutine Na1 (e,val)
    real(dp) :: e,val
    val = (1.d0+31.d0/2.d0*e**2+255.d0/8.d0*e**4 &
        +185.d0/16.d0*e**6+25.d0/64.d0*e**8)/(1.d0-e**2)**(15.d0/2.d0)     
    return
end subroutine Na1

subroutine Na2(e,val)
    real(dp) :: e,val
    val = (1d0+15/2.d0*e**2+45/8.d0*e**4+5/16.d0*e**6)/(1d0-e**2)**6
    return
end subroutine Na2


! The 2 eccentricity dependant factors in the equation in e
subroutine Ne1 (e,val)
    real(dp) :: e,val
    val = (1d0+15d0/4.d0*e**2+15d0/8.d0*e**4+5d0/64.d0*e**6)/(1-e**2)**(13d0/2.d0)
    return
end subroutine Ne1

subroutine Ne2(e,val)
    real(dp) :: e,val
    val =  (1d0+3d0/2.d0*e**2+1d0/8.d0*e**4)/(1d0-e**2)**5
    return
end subroutine Ne2


! The 2 eccentricity dependant factors in the equation in Ohmega and
! in i
subroutine No1(e,val)
    real(dp) :: e,val
    val =   (1d0+15d0/2.d0*e**2+45d0/8.d0*e**4+5d0/16.d0*e**6)/(1-e**2)**(13d0/2.d0)
    return
end subroutine No1

subroutine No2(e,val)
    real(dp) :: e,val
    val =    (1d0+3d0*e**2+3d0/8.d0*e**4)/(1d0-e**2)**5
    return
end subroutine No2


! Desperate attempt to fix code : 
subroutine Nob(e,val)
    real(dp) :: e,val
    val =    (1.d0+3.d0*e**2+3.d0/8.d0*e**4)/(1.d0-e**2)**(9.d0/2.d0)
    return
end subroutine Nob


! Mean orbital angular velocity without the a dependance         (m**3/2.s-1)
subroutine norb(G,Mp,Ms,val)
    real(dp) :: G,Mp,Ms,val
    val = sqrt(G)*sqrt(Mp+Ms)
    return
end subroutine norb



! Timescales without the a dependance (exchange p <-> s)         (s.m**-8)
subroutine Temp(Ms,Mp,Rs,sigmas,val)
    real(dp) :: Ms,Mp,Rs,sigmas,val
    val = Ms/(9.d0*Mp*(Mp+Ms)*Rs**10*sigmas)
    
    !write(6,*) "Temp", Ms,Mp,Rs,sigmas,val
    
    return
end subroutine Temp

! Orbital angular momentum
! without the a dependance
subroutine horb(e,Mp,Ms,G,val)
    real(dp) :: e,Mp,Ms,G,val
    val = sqrt(G*(1d0-e**2)/(Ms+Mp))*Ms*Mp
    return
end subroutine horb


! Ratio of orbital angular momentum on spin angular momentum
! (m**-1/2)
! without the a dependance, I = rg2*(M*R**2)
subroutine gamma(e,rotp,Mp,Ms,Rp,G,rg2p,val)
    real(dp) :: e,rotp,Mp,Ms,Rp,G,rg2p,val
    val = sqrt(G*(1d0-e**2)/(Ms+Mp))*Ms*Mp &
        /(rg2p*Mp*Rp**2*rotp)
    return
end subroutine gamma


! ratio between orbital angular momentum and spin angular momentum for
! pseudo synchronisation
subroutine alpha(e,val)
    real(dp) :: e,val
    val = (1d0+15/2.d0*e**2+45/8.d0*e**4+5/16.d0*e**6)* &
        1.d0/(1d0+3d0*e**2+3/8.d0*e**4)*1./(1-e**2)**1.5d0
    return
end subroutine alpha


! pseudo synchronisation rotation without the a dependance:
subroutine pseudorot(e,G,Mp,Ms,val)
    real(dp) :: e,G,Mp,Ms,val
    real(dp) :: tmp,tmp2

    call alpha(e,tmp)
    call norb(G,Mp,Ms,tmp2)

    val = tmp*tmp2
    return
end subroutine pseudorot


! Derivative of a
subroutine adot(a,e,oblp,obls,rotp,rots,G,Mp,Ms,Rp,Rs,sigmap,sigmas,val)

    real(dp) :: a,e,oblp,obls,rotp,rots,G,Mp,Ms,Rp,Rs,sigmap,sigmas,val
    real(dp) :: tmp,tmp2,tmp3,tmp4,tmp5
   
    call norb(G,Mp,Ms,tmp)
    call Na1(e,tmp2)
    call Na2(e,tmp3)
    call Temp(Ms,Mp,Rs,sigmas,tmp4)
    call Temp(Mp,Ms,Rp,sigmap,tmp5)

    val = - (1.d0/tmp5*a**(-7.d0) &
        *(tmp2-cos(oblp)*(rotp/tmp)*a**(1.5d0)*tmp3)) &
        - (1.d0/tmp4*a**(-7.d0) &
        *(tmp2-cos(obls)*(rots/tmp)*a**(1.5d0)*tmp3)) 
        
        
     !write(6,*) "a...=",a,cos(obls),rots,rotp
        
    return
end subroutine adot


!Derivative of e
subroutine edot(a,e,oblp,obls,rotp,rots,G,Mp,Ms,Rp,Rs,sigmap,sigmas,val)
    real(dp) :: a,e,oblp,obls,rotp,rots,G,Mp,Ms,Rp,Rs,sigmap,sigmas,val
    real(dp) :: tmp,tmp2,tmp3,tmp4,tmp5

    call Temp(Mp,Ms,Rp,sigmap,tmp)
    call Temp(Ms,Mp,Rs,sigmas,tmp5)
    call Ne1(e,tmp2)
    call norb(G,Mp,Ms,tmp3)
    call Ne2(e,tmp4)

    val = - (9.d0/2.d0)*1.d0/tmp*e*a**(-8.d0) &
        *(tmp2-(11.d0/18.d0) &
        *cos(oblp)*(rotp/tmp3)*a**(1.5d0)*tmp4) &
        - (9.d0/2.d0)*1.d0/tmp5*e*a**(-8.d0) &
        *(tmp2-(11.d0/18.d0) &
        *cos(obls)*(rots/tmp3)*a**(1.5d0)*tmp4)
    return
end subroutine edot


!Derivative of obl
subroutine obldot(a,e,oblp,rotp,G,Mp,Ms,Rp,sigmap,rg2p,val)
    real(dp) :: a,e,oblp,rotp,G,Mp,Ms,Rp,sigmap,rg2p,val
      real(dp) :: tmp,tmp2,tmp3,tmp4,tmp5

    
      call Temp(Mp,Ms,Rp,sigmap,tmp)
      call Nob(e,tmp2)
      call norb(G,Mp,Ms,tmp3)
      call Na2(e,tmp4)
      call gamma(e,rotp,Mp,Ms,Rp,G,rg2p,tmp5)
    
    
  val = 1.d0/(4.d0*tmp)*sin(oblp)*a**(-15.d0/2.d0) &
          *G*Mp*Ms/(Mp*rg2p*rotp*tmp3*Rp**2) &
          *((cos(oblp)-a**(-1.d0/2.d0)/tmp5)*tmp2 &
            *(rotp/tmp3)*a**(1.5d0)-2.d0*tmp4)
    return
end subroutine obldot



!Derivative of rot_star: 1/Omega*dOmega/dt (rg2s = cst, Rs = cst)
subroutine frots_dot(a,e,obls,rots,G,Mp,Ms,Rs,sigmas,rg2s,val)
    real(dp) ::  a,e,obls,rots,G,Mp,Ms,Rs,sigmas,rg2s,val  
    real(dp) :: tmp,tmp2,tmp3,tmp4,tmp5

    !call horb(e,Ms,Mp,G,tmp)
    call gamma(e,rots,Ms,Mp,Rs,G,rg2s,tmp)
    
    call Temp(Ms,Mp,Rs,sigmas,tmp2)
    call No1(e,tmp3)
    call norb(G,Mp,Ms,tmp4)
    call No2(e,tmp5)

    val = tmp/(4.d0*tmp2) &
        *a**(-15.d0/2.d0)*(2.d0*cos(obls)*tmp3 &
        -(1.d0+cos(obls)*cos(obls)) &
        *(rots/tmp4)*a**(1.5d0)*tmp5)
 return
end subroutine frots_dot


!Derivative of J: dJ/dt (rg2s = cst, Rs = cst)
subroutine frots_dot2(a,e,obls,rots,G,Mp,Ms,Rs,sigmas,rg2s,val)
    real(dp) ::  a,e,obls,rots,G,Mp,Ms,Rs,sigmas,rg2s,val  
    real(dp) :: tmp,tmp2,tmp3,tmp4,tmp5

    call horb(e,Ms,Mp,G,tmp)
    !call gamma(e,rots,Ms,Mp,Rs,G,rg2s,tmp)
    
    call Temp(Ms,Mp,Rs,sigmas,tmp2)
    call No1(e,tmp3)
    call norb(G,Mp,Ms,tmp4)
    call No2(e,tmp5)

    val = tmp/(4.d0*tmp2) &
        *a**(-15.d0/2.d0)*(2.d0*cos(obls)*tmp3 &
        -(1.d0+cos(obls)*cos(obls)) &
        *(rots/tmp4)*a**(1.5d0)*tmp5)
 return
end subroutine frots_dot2



!Thing to take into account wind braking (Bouvier 97)
subroutine fwind(rots,spin_sat,Kwind,Ms,Msun,Rs,Rsun,rg2s,alpha,val)
    real(dp) :: rots,spin_sat,Kwind,Ms,Msun,Rs,Rsun,rg2s,alpha,val

    val =-1.d0/(Ms*rg2s*Rs**(2.))*Kwind*rots**(alpha-1.d0) &
        *spin_sat**(3.d0-alpha)*(Rs/Rsun)**(0.5)*(Ms/Msun)**(-0.5)
    return
end subroutine fwind

!Derivative of rot_plan
subroutine rotp_dot(a,e,oblp,rotp,G,Mp,Ms,Rp,sigmap,rg2p,val)
    real(dp) :: a,e,oblp,rotp,G,Mp,Ms,Rp,sigmap,rg2p,val
    real(dp) :: tmp,tmp2,tmp3,tmp4,tmp5

    call gamma(e,rotp,Mp,Ms,Rp,G,rg2p,tmp)
    call Temp(Mp,Ms,Rp,sigmap,tmp2)
    call No1(e,tmp3)
    call norb(G,Mp,Ms,tmp4)
    call No2(e,tmp5)

    val = tmp/(4.d0*tmp2)*rotp &
        *a**(-15.d0/2.d0)*(2.d0*cos(oblp)*tmp3 &
        -(1.d0+cos(oblp)*cos(oblp)) &
        *(rotp/tmp4)*a**(1.5d0)*tmp5)
    return
end subroutine rotp_dot


!a for which da/dt =0
subroutine aeq(e,rotp,rots,G,Mp,Ms,Rp,Rs,sigmap,sigmas,val)
    real(dp) :: e,rotp,rots,G,Mp,Ms,Rp,Rs,sigmap,sigmas,val
    real(dp) :: tmp,tmp2,tmp3,tmp4,tmp5

    call Na1(e,tmp)
    call norb(G,Mp,Ms,tmp2)
    call Temp(Mp,Ms,Rp,sigmap,tmp3)
    call Temp(Ms,Mp,Rs,sigmas,tmp4)
    call Na2(e,tmp5)

    val = (tmp*tmp2* (tmp3+tmp4)/ (tmp5*(tmp3*rotp + tmp4*rots)))**(2/3.)

    return                 
end subroutine aeq



!! !Energy due to tidal dissipation in bodies 
subroutine enerdot(a,e,rotp,G,Mp,Ms,Rp,sigmap,val)
    real(dp) :: a,e,rotp,G,Mp,Ms,Rp,sigmap,val
    real(dp) :: tmp,tmp2,tmp3,tmp4,tmp5

    call Temp(Mp,Ms,Rp,sigmap,tmp)
    call Na1(e,tmp2)
    call Na2(e,tmp3)
    call norb(G,Mp,Ms,tmp4)
    call No2(e,tmp5)
    
    val = 1.d0/(2.d0*tmp)*G*Mp*Ms*a**(-9.0) &
        * (tmp2-2.d0*tmp3*(rotp/tmp4)*a**(3d0/2.d0) &
        + tmp5*sqrt(1-e**2.)*(rotp/tmp4)**2.0*a**3.0)
    return
end subroutine enerdot

!Jeremy's Ki factor : 
subroutine Kplan(k2deltat_plan,G,Mp,Ms,Rp,a,val)
    real(dp) :: k2deltat_plan,G,Mp,Ms,Rp,a,val
    real(dp) :: tmp

    call norb(G,Mp,Ms,tmp)

    val = 3.d0/2.d0*k2deltat_plan*(G*Mp**2/Rp)*(Ms/Mp)**2 &
        *(Rp/a)**6*(tmp*a**(-1.5d0))**2.
    return
end subroutine Kplan


!Energy due to tidal dissipation in bodies Jeremy's exact formula
!subroutine enerdot,a,e,rotp,oblp,G,Mp,Ms,Rp,k2deltat_plan
  !val = 2.d0* Kplan(k2deltat_plan,G,Mp,Ms,Rp,a) &
          !* (Na1(e)-2.d0*Na2(e)*cos(oblp)*(rotp/(norb(G,Mp,Ms)*a**(-1.5d0))) &
             !+ (1.0d0+cos(oblp)**2)/2.d0*No2(e)*sqrt(1.d0-e**2)*(rotp/(norb(G,Mp,Ms)*a**(-1.5d0)))**2)
!end

! ********** Second derivative calculation *************

subroutine splineme(x,y,n,val)
 
    integer :: n
    integer :: Nmax,i
    real(dp), dimension(n) :: y2,x,y,val
    real(dp), dimension(n) :: u
    real(dp) :: sig,p,qn,un
  
    Nmax = n
    !Nmax = n
    !u  = dblarr(Nmax)
    !y2 = dblarr(n)
    
    y2(1) = 0.0
    u(1)  = 0.0

    do i = 2,n-1 
       sig   = (x(i)-x(i-1))/(x(i+1)-x(i-1))
       p     = sig * y2(i-1)+2
       y2(i) = (sig-1.)/p
       u(i)  = (6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) &
                    /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
    enddo

    qn = 0.0
    un = 0.0
    y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.)

    do i = 1,n-1 
       y2(n-1-i+1) = y2(n-1-i+1)*y2(n-1-i+2)+u(n-1-i+1)
    enddo

    val = y2

end subroutine splineme

! ************* Interpolation *************

subroutine splint(xa,ya,y2a,n,x,val)

    integer :: klo,khi,n,k
    real(dp), dimension(n) :: xa,ya,y2a
    real(dp) :: val,x,h,a,b,y
  
    klo = 1
    khi = n
    
    do while ((khi-klo) .gt. 1)
       k = int((khi+klo)/2.)
       if (xa(k) .gt. x) then 
          khi = k
       else 
          klo = k
       endif
    enddo
    h = xa(khi)-xa(klo)
    if (h .eq. 0.) write(6,*) 'Bad xa input in splint'
    a = (xa(khi)-x)/h
    b = (x-xa(klo))/h
    y = a*ya(klo)+b*ya(khi) &
        + ((a**3.-a)*y2a(klo)+(b**3.-b)*y2a(khi))*(h**2)/6.
        
    if ( (ya(klo) .eq. 0.0) .and. (ya(khi) .eq. 0.0)) then
        val = 0.0
    else   
        val = y
    endif    
end subroutine splint

! ********** First derivative calculation *************

subroutine derivee(xa,ya,y2a,n,x ,val)

    integer :: klo,khi,n,k,i
    real(dp), dimension(n) :: xa,ya,y2a
    real(dp) :: val,x,h,dydx1,dydx2,dydx3

  do i = 1,n-1
     klo = 0
     khi = n
     do while ((khi-klo) .gt. 1)
        k = int((khi+klo)/2.)
        if (xa(k) .gt. x) then 
           khi = k
        else 
           klo = k
        endif
     enddo    
  enddo
  h = xa(khi)-xa(klo)
  dydx1 = (ya(khi)-ya(klo))/(xa(khi)-xa(klo))
  dydx2 = -(3*((xa(khi)-x)/(xa(khi)-xa(klo)))**2-1)/6.*(xa(khi)-xa(klo))*y2a(klo)
  dydx3 = (3*((x-xa(klo))/(xa(khi)-xa(klo)))**2-1)/6.*(xa(khi)-xa(klo))*y2a(khi) 
  val = dydx1+dydx2+dydx3
end subroutine derivee

end module fonctions
