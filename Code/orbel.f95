!**********************************************************************
!                    ORBEL_FGET.F
!**********************************************************************
!     PURPOSE:  Solves Kepler's eqn. for hyperbola using hybrid approach.  
!
!             Input:
!                           e ==> eccentricity anomaly. (real scalar)
!                        capn ==> hyperbola mean anomaly. (real scalar)
!             Returns:
!                  orbel_fget ==>  eccentric anomaly. (real scalar)
!
!     ALGORITHM: Based on pp. 70-72 of Fitzpatrick's book "Principles of
!           Cel. Mech. ".  Quartic convergence from Danby's book.
!     REMARKS: 
!     AUTHOR: M. Duncan 
!     DATE WRITTEN: May 11, 1992.
!     REVISIONS: 2/26/93 hfl
!     Modified by JEC
!**********************************************************************

    real(8) function orbel_fget(e,capn)

      include 'swift.inc'

!...  Inputs Only: 
    real(8) e,capn

!...  Internals:
    integer i,IMAX
    real(8) tmp,x,shx,chx
    real(8) esh,ech,f,fp,fpp,fppp,dx
    PARAMETER (IMAX = 10)

!----
!...  Executable code 

! Function to solve "Kepler's eqn" for F (here called
! x) for given e and CAPN. 

!  begin with a guess proposed by Danby    
    if( capn .lt. 0.d0) then
       tmp = -2.d0*capn/e + 1.8d0
       x = -log(tmp)
    else
       tmp = +2.d0*capn/e + 1.8d0
       x = log( tmp)
    endif

    orbel_fget = x

    do i = 1,IMAX
          call mco_sinh (x,shx,chx)
      esh = e*shx
      ech = e*chx
      f = esh - x - capn
!      write(6,*) 'i,x,f : ',i,x,f
      fp = ech - 1.d0  
      fpp = esh 
      fppp = ech 
      dx = -f/fp
      dx = -f/(fp + dx*fpp/2.d0)
      dx = -f/(fp + dx*fpp/2.d0 + dx*dx*fppp/6.d0)
      orbel_fget = x + dx
!   If we have converged here there's no point in going on
      if(abs(dx) .le. TINY) RETURN
      x = orbel_fget
    enddo    

    write(6,*) 'FGET : RETURNING WITHOUT COMPLETE CONVERGENCE' 
    return
    end   ! orbel_fget
!------------------------------------------------------------------
!
!**********************************************************************
!                    ORBEL_FHYBRID.F
!**********************************************************************
!     PURPOSE:  Solves Kepler's eqn. for hyperbola using hybrid approach.  
!
!             Input:
!                           e ==> eccentricity anomaly. (real scalar)
!                           n ==> hyperbola mean anomaly. (real scalar)
!             Returns:
!               orbel_fhybrid ==>  eccentric anomaly. (real scalar)
!
!     ALGORITHM: For abs(N) < 0.636*ecc -0.6 , use FLON 
!             For larger N, uses FGET
!     REMARKS: 
!     AUTHOR: M. Duncan 
!     DATE WRITTEN: May 26,1992.
!     REVISIONS: 
!     REVISIONS: 2/26/93 hfl
!**********************************************************************

    real(8) function orbel_fhybrid(e,n)

      include 'swift.inc'

!...  Inputs Only: 
    real(8) e,n

!...  Internals:
    real(8) abn
        real(8) orbel_flon,orbel_fget

!----
!...  Executable code 

    abn = n
    if(n.lt.0.d0) abn = -abn

    if(abn .lt. 0.636d0*e -0.6d0) then
      orbel_fhybrid = orbel_flon(e,n)
    else 
      orbel_fhybrid = orbel_fget(e,n)
    endif   

    return
    end  ! orbel_fhybrid
!-------------------------------------------------------------------
!
!**********************************************************************
!                    ORBEL_FLON.F
!**********************************************************************
!     PURPOSE:  Solves Kepler's eqn. for hyperbola using hybrid approach.  
!
!             Input:
!                           e ==> eccentricity anomaly. (real scalar)
!                        capn ==> hyperbola mean anomaly. (real scalar)
!             Returns:
!                  orbel_flon ==>  eccentric anomaly. (real scalar)
!
!     ALGORITHM: Uses power series for N in terms of F and Newton,s method
!     REMARKS: ONLY GOOD FOR LOW VALUES OF N (N < 0.636*e -0.6)
!     AUTHOR: M. Duncan 
!     DATE WRITTEN: May 26, 1992.
!     REVISIONS: 
!**********************************************************************

    real(8) function orbel_flon(e,capn)

      include 'swift.inc'

!...  Inputs Only: 
    real(8) e,capn

!...  Internals:
    integer iflag,i,IMAX
    real(8) a,b,sq,biga,bigb
    real(8) x,x2
    real(8) f,fp,dx
    real(8) diff
    real(8) a0,a1,a3,a5,a7,a9,a11
    real(8) b1,b3,b5,b7,b9,b11
    PARAMETER (IMAX = 10)
    PARAMETER (a11 = 156.d0,a9 = 17160.d0,a7 = 1235520.d0)
    PARAMETER (a5 = 51891840.d0,a3 = 1037836800.d0)
    PARAMETER (b11 = 11.d0*a11,b9 = 9.d0*a9,b7 = 7.d0*a7)
    PARAMETER (b5 = 5.d0*a5, b3 = 3.d0*a3)

!----
!...  Executable code 


! Function to solve "Kepler's eqn" for F (here called
! x) for given e and CAPN. Only good for smallish CAPN 

    iflag = 0
    if( capn .lt. 0.d0) then
       iflag = 1
       capn = -capn
    endif

    a1 = 6227020800.d0 * (1.d0 - 1.d0/e)
    a0 = -6227020800.d0*capn/e
    b1 = a1

!  Set iflag nonzero if capn < 0., in which case solve for -capn
!  and change the sign of the final answer for F.
!  Begin with a reasonable guess based on solving the cubic for small F    


    a = 6.d0*(e-1.d0)/e
    b = -6.d0*capn/e
    sq = sqrt(0.25*b*b +a*a*a/27.d0)
    biga = (-0.5*b + sq)**0.3333333333333333d0
    bigb = -(+0.5*b + sq)**0.3333333333333333d0
    x = biga + bigb
!    write(6,*) 'cubic = ',x**3 +a*x +b
    orbel_flon = x
! If capn is tiny (or zero) no need to go further than cubic even for
! e =1.
    if( capn .lt. TINY) go to 100

    do i = 1,IMAX
      x2 = x*x
      f = a0 +x*(a1+x2*(a3+x2*(a5+x2*(a7+x2*(a9+x2*(a11+x2))))))
      fp = b1 +x2*(b3+x2*(b5+x2*(b7+x2*(b9+x2*(b11 + 13.d0*x2)))))   
      dx = -f/fp
!      write(6,*) 'i,dx,x,f : '
!      write(6,432) i,dx,x,f
432      format(1x,i3,3(2x,1p1e22.15))
      orbel_flon = x + dx
!   If we have converged here there's no point in going on
      if(abs(dx) .le. TINY) go to 100
      x = orbel_flon
    enddo    

! Abnormal return here - we've gone thru the loop 
! IMAX times without convergence
    if(iflag .eq. 1) then
       orbel_flon = -orbel_flon
       capn = -capn
    endif
    write(6,*) 'FLON : RETURNING WITHOUT COMPLETE CONVERGENCE' 
      diff = e*sinh(orbel_flon) - orbel_flon - capn
      write(6,*) 'N, F, ecc*sinh(F) - F - N : '
      write(6,*) capn,orbel_flon,diff
    return

!  Normal return here, but check if capn was originally negative
100    if(iflag .eq. 1) then
       orbel_flon = -orbel_flon
       capn = -capn
    endif

    return
    end     ! orbel_flon
!
!***********************************************************************
!                    ORBEL_ZGET.F
!***********************************************************************
!     PURPOSE:  Solves the equivalent of Kepler's eqn. for a parabola 
!          given Q (Fitz. notation.)
!
!             Input:
!                           q ==>  parabola mean anomaly. (real scalar)
!             Returns:
!                  orbel_zget ==>  eccentric anomaly. (real scalar)
!
!     ALGORITHM: p. 70-72 of Fitzpatrick's book "Princ. of Cel. Mech."
!     REMARKS: For a parabola we can solve analytically.
!     AUTHOR: M. Duncan 
!     DATE WRITTEN: May 11, 1992.
!     REVISIONS: May 27 - corrected it for negative Q and use power
!          series for small Q.
!**********************************************************************

    real(8) function orbel_zget(q)

      include 'swift.inc'

!...  Inputs Only: 
    real(8) q

!...  Internals:
    integer iflag
    real(8) x,tmp

!----
!...  Executable code 

    iflag = 0
    if(q.lt.0.d0) then
      iflag = 1
      q = -q
    endif

    if (q.lt.1.d-3) then
       orbel_zget = q*(1.d0 - (q*q/3.d0)*(1.d0 -q*q))
    else
       x = 0.5d0*(3.d0*q + sqrt(9.d0*(q**2) +4.d0))
       tmp = x**(1.d0/3.d0)
       orbel_zget = tmp - 1.d0/tmp
    endif

    if(iflag .eq.1) then
           orbel_zget = -orbel_zget
       q = -q
    endif
    
    return
    end    ! orbel_zget
!----------------------------------------------------------------------
