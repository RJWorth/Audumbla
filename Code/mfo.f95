!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MFO_USER.FOR    (ErikSoft   2 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Applies an arbitrary force, defined by the user.
!
! If using with the symplectic algorithm MAL_MVS, the force should be
! small compared with the force from the central object.
! If using with the conservative Bulirsch-Stoer algorithm MAL_BS2, the
! force should not be a function of the velocities.
!
! N.B. All coordinates and velocities must be with respect to central body
! ===
!------------------------------------------------------------------------------
!
! Galactic tides and dispersing gas cloud.
!

      subroutine mfo_user (time,jcen,nbod,nbig,m,x,v,a)
!
      implicit none
      include 'mercury95.inc'
!
! Input/Output
      integer nbod, nbig
      real*8 time,jcen(3),m(nbod),x(3,nbod),v(3,nbod),a(3,nbod)
!  **************   DGV Dimitri Veras added these parameters ************
      real*8 Omegag,rhog,delta
      parameter (Omegag=7.66963e-11,rhog=1.05629e-17,delta=-0.0573185)
!
! From Veras & Evans (2013) MNRAS
!
! R(kpc)      Omega         Delta          Rho
!  0.1    4.37096*10^-9    0.377521   9.06232*10^-15
!  0.2    2.75482*10^-9    0.283677   3.26618*10^-15
!  0.5    1.32134*10^-9    0.105695   6.48106*10^-16
!  1.0    6.78356*10^-10  -0.0214944  1.94041*10^-16
!  2.0    3.27146*10^-10  -0.0641641  8.42905*10^-17
!  3.0    2.13134*10^-10  -0.0469292  5.53658*10^-17
!  5.0    1.25414*10^-10  -0.0362197  2.77995*10^-17
!  8.0    7.66963*10^-11  -0.0573185  1.05629*10^-17
! 12.0    4.99957*10^-11  -0.0431656  3.12614*10^-18
! 20.0    2.99106*10^-11   0.0313698  4.43947*10^-19
! 35.0    1.74928*10^-11   0.0346311  8.98535*10^-20
!
! Local
      integer j
!  **************   #rjw# Rachel Worth added these parameters ************
      real*8 aa(3,nbod),ab(3,nbod),cm(3)
      real*8 mtot,mgas0,mgas,rgas,gaslife,gasfrac
!
!  **************   DGV Dimitri Veras added this subroutine ************
! #rjw# changed a to aa
!------------------------------------------------------------------------------
!
      do j = 2, nbod
        aa(1,j)= (Omegag**2.0) &
     &         *((1.0-delta)*cos(2.0*Omegag*time) - delta) &
     &         *x(1,j) + &
     &          (Omegag**2.0)*(1.0-delta)*sin(2.0*Omegag*time) &
               *x(2,j)

        aa(2,j)= (Omegag**2.0)*(1.0-delta)*sin(2.0*Omegag*time)*x(1,j) - &
     &          (Omegag**2.0)*((1.0-delta)*cos(2.0*Omegag*time) + delta) &
     &         *x(2,j)

        aa(3,j)= -(4.0*3.1415926535897932*K2*rhog - &
     &           2.0*(Omegag**2.0)*delta)*x(3,j) 
      end do
!
!------------------------------------------------------------------------------
!  **************   DGV Dimitri Veras end ************
!
!  **************   #rjw# Rachel Worth added this subroutine ************
!------------------------------------------------------------------------------
! Initial mass of cloud, solar masses*K2
      mgas0 = 10*K2
! Lifetime of gas cloud (this and 'time' are in days)
      gaslife = 440000*365.25
! Mass over time, as cloud dissipates (...What about star accretion?)
      gasfrac= max(0.0, (gaslife-time)/gaslife)
      mgas = mgas0*gasfrac
! Radius of cloud core, in AU
      rgas = 7500

! Calculate total mass
      do j = 1, nbod
        mtot= mtot+m(j)
      end do

! Calculate location of center of mass
      do j = 2, nbod
        cm(1)= cm(1)+(1/mtot)*(m(j)*x(1,j))
!        write(*,*) cm(1)
        cm(2)= cm(2)+(1/mtot)*(m(j)*x(2,j))

        cm(3)= cm(3)+(1/mtot)*(m(j)*x(3,j))
      end do
      
! Calculate force from remnant gas cloud
      do j = 2, nbod
        ab(1,j)= -mgas*x(1,j) / (x(1,j)**2.0 + rgas**2.0)**1.5
!         ab(1,j)= 0.0
        ab(2,j)= -mgas*x(2,j) / (x(2,j)**2.0 + rgas**2.0)**1.5
!         ab(2,j)= 0.0
        ab(3,j)= -mgas*x(3,j) / (x(3,j)**2.0 + rgas**2.0)**1.5
!         ab(3,j)= 0.0
      end do
    
! Combine galactic tide and gas cloud effects
      do j = 2, nbod
        a(1,j)= aa(1,j) + ab(1,j)

        a(2,j)= aa(2,j) + ab(2,j)

        a(3,j)= aa(3,j) + ab(3,j)
      end do


!------------------------------------------------------------------------------
!  **************   #rjw# Rachel Worth end ************
!
      return
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MFO_ALL.FOR    (ErikSoft   2 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Calculates accelerations on a set of NBOD bodies (of which NBIG are Big)
! due to Newtonian gravitational perturbations, post-Newtonian
! corrections (if required), cometary non-gravitational forces (if required)
! and user-defined forces (if required).
!
! N.B. Input/output must be in coordinates with respect to the central body.
! ===
!
!------------------------------------------------------------------------------
!
      subroutine mfo_all (time,jcen,nbod,nbig,m,x,v,s,rcrit,a,stat,ngf, &
       ngflag,opt,nce,ice,jce)
!
      implicit none
      include 'mercury95.inc'
!
! Input/Output
      integer nbod,nbig,ngflag,stat(nbod),opt(8),nce,ice(nce),jce(nce)
      real(8) time,jcen(3),m(nbod),x(3,nbod),v(3,nbod),s(3,nbod)
      real(8) a(3,nbod),ngf(4,nbod),rcrit(nbod)
!
! Local
      integer j
      real(8) acor(3,NMAX),acen(3)
!
!------------------------------------------------------------------------------
!
! Newtonian gravitational forces
      call mfo_grav (nbod,nbig,m,x,v,a,stat)
!
! Correct for oblateness of the central body
      if (jcen(1).ne.0.or.jcen(2).ne.0.or.jcen(3).ne.0) then
        call mfo_obl (jcen,nbod,m,x,acor,acen)
        do j = 2, nbod
          a(1,j) = a(1,j) + (acor(1,j) - acen(1))
          a(2,j) = a(2,j) + (acor(2,j) - acen(2))
          a(3,j) = a(3,j) + (acor(3,j) - acen(3))
        end do
      end if
!
! Include non-gravitational (cometary jet) accelerations if necessary
      if (ngflag.eq.1.or.ngflag.eq.3) then
        call mfo_ngf (nbod,x,v,acor,ngf)
        do j = 2, nbod
          a(1,j) = a(1,j) + acor(1,j)
          a(2,j) = a(2,j) + acor(2,j)
          a(3,j) = a(3,j) + acor(3,j)
        end do
      end if
!
! Include radiation pressure/Poynting-Robertson drag if necessary
      if (ngflag.eq.2.or.ngflag.eq.3) then
        call mfo_pr (nbod,nbig,m,x,v,acor,ngf)
        do j = 2, nbod
          a(1,j) = a(1,j) + acor(1,j)
          a(2,j) = a(2,j) + acor(2,j)
          a(3,j) = a(3,j) + acor(3,j)
        end do
      end if
!
! Include post-Newtonian corrections if required
      if (opt(7).eq.1) then
        call mfo_pn (nbod,nbig,m,x,v,acor)
        do j = 2, nbod
          a(1,j) = a(1,j) + acor(1,j)
          a(2,j) = a(2,j) + acor(2,j)
          a(3,j) = a(3,j) + acor(3,j)
        end do
      end if
!
! Include user-defined accelerations if required
      if (opt(8).eq.1) then
        call mfo_user (time,jcen,nbod,nbig,m,x,v,acor)
        do j = 2, nbod
          a(1,j) = a(1,j) + acor(1,j)
          a(2,j) = a(2,j) + acor(2,j)
          a(3,j) = a(3,j) + acor(3,j)
        end do
      end if
!
!------------------------------------------------------------------------------
!
      return
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MFO_GRAV.FOR    (ErikSoft   3 October 2000)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Calculates accelerations on a set of NBOD bodies (NBIG of which are Big)
! due to gravitational perturbations by all the other bodies, except that
! Small bodies do not interact with one another.
!
! The positions and velocities are stored in arrays X, V with the format
! (x,y,z) and (vx,vy,vz) for each object in succession. The accelerations 
! are stored in the array A (ax,ay,az).
!
! N.B. All coordinates and velocities must be with respect to central body!!!!
! ===
!------------------------------------------------------------------------------
!
      subroutine mfo_grav (nbod,nbig,m,x,v,a,stat)
!
      implicit none
      include 'mercury95.inc'
!
! Input/Output
      integer nbod, nbig, stat(nbod)
      real(8) m(nbod), x(3,nbod), v(3,nbod), a(3,nbod)
!
! Local
      integer i, j
      real(8) sx, sy, sz, dx, dy, dz, tmp1, tmp2, s_1, s2, s_3, r3(NMAX)
!
!------------------------------------------------------------------------------
!
      sx = 0.d0
      sy = 0.d0
      sz = 0.d0
      do i = 2, nbod
        a(1,i) = 0.d0
        a(2,i) = 0.d0
        a(3,i) = 0.d0
        s2 = x(1,i)*x(1,i) + x(2,i)*x(2,i) + x(3,i)*x(3,i)
        s_1  = 1.d0 / sqrt(s2)
        r3(i) = s_1 * s_1 * s_1
      end do
!
      do i = 2, nbod
        tmp1 = m(i) * r3(i)
        sx = sx  -  tmp1 * x(1,i)
        sy = sy  -  tmp1 * x(2,i)
        sz = sz  -  tmp1 * x(3,i)
      end do
!
! Direct terms
      do i = 2, nbig
        do j = i + 1, nbod
          dx = x(1,j) - x(1,i)
          dy = x(2,j) - x(2,i)
          dz = x(3,j) - x(3,i)
          s2 = dx*dx + dy*dy + dz*dz
          s_1 = 1.d0 / sqrt(s2)
          s_3 = s_1 * s_1 * s_1
          tmp1 = s_3 * m(i)
          tmp2 = s_3 * m(j)
          a(1,j) = a(1,j)  -  tmp1 * dx
          a(2,j) = a(2,j)  -  tmp1 * dy
          a(3,j) = a(3,j)  -  tmp1 * dz
          a(1,i) = a(1,i)  +  tmp2 * dx
          a(2,i) = a(2,i)  +  tmp2 * dy
          a(3,i) = a(3,i)  +  tmp2 * dz
        end do
      end do
!
! Indirect terms (add these on last to reduce roundoff error)
      do i = 2, nbod
        tmp1 = m(1) * r3(i)
        a(1,i) = a(1,i)  +  sx  -  tmp1 * x(1,i)
        a(2,i) = a(2,i)  +  sy  -  tmp1 * x(2,i)
        a(3,i) = a(3,i)  +  sz  -  tmp1 * x(3,i)
      end do
!
!------------------------------------------------------------------------------
!
      return
      end
!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!     MFO_DRCT.FOR    (ErikSoft   27 February 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Calculates direct accelerations between bodies in the interaction part
! of the Hamiltonian of a symplectic integrator that partitions close
! encounter terms (e.g. hybrid symplectic algorithms or SyMBA).
! The routine calculates accelerations between all pairs of bodies with
! indices I >= I0.
!
!------------------------------------------------------------------------------
!
      subroutine mfo_drct (i0,nbod,nbig,m,x,rcrit,a,stat)
!
      implicit none
      include 'mercury95.inc'
!
! Input/Output
      integer i0, nbod, nbig, stat(nbod)
      real(8) m(nbod), x(3,nbod), a(3,nbod), rcrit(nbod)
!
! Local
      integer i,j
      real(8) dx,dy,dz,s,s_1,s2,s_3,rc,rc2,q,q2,q3,q4,q5,tmp2,faci,facj
!
!------------------------------------------------------------------------------
!
      if (i0.le.0) i0 = 2
!
      do i = i0, nbig
        do j = i + 1, nbod
          dx = x(1,j) - x(1,i)
          dy = x(2,j) - x(2,i)
          dz = x(3,j) - x(3,i)
          s2 = dx * dx  +  dy * dy  +  dz * dz
          rc = max(rcrit(i), rcrit(j))
          rc2 = rc * rc
!
          if (s2.ge.rc2) then
            s_1 = 1.d0 / sqrt(s2)
            tmp2 = s_1 * s_1 * s_1
          else if (s2.le.0.01*rc2) then
            tmp2 = 0.d0
          else
            s_1 = 1.d0 / sqrt(s2)
            s   = 1.d0 / s_1
            s_3 = s_1 * s_1 * s_1
            q = (s - 0.1d0*rc) / (0.9d0 * rc)
            q2 = q  * q
            q3 = q  * q2
            q4 = q2 * q2
            q5 = q2 * q3
            tmp2 = (10.d0*q3 - 15.d0*q4 + 6.d0*q5) * s_3
          end if
!
          faci = tmp2 * m(i)
          facj = tmp2 * m(j)
          a(1,j) = a(1,j)  -  faci * dx
          a(2,j) = a(2,j)  -  faci * dy
          a(3,j) = a(3,j)  -  faci * dz
          a(1,i) = a(1,i)  +  facj * dx
          a(2,i) = a(2,i)  +  facj * dy
          a(3,i) = a(3,i)  +  facj * dz
        end do
      end do
!
!------------------------------------------------------------------------------
!
      return
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!     MFO_HY.FOR    (ErikSoft   2 October 2000)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Calculates accelerations due to the Interaction part of the Hamiltonian 
! of a hybrid symplectic integrator for a set of NBOD bodies (NBIG of which 
! are Big), where Small bodies do not interact with one another.
!
!------------------------------------------------------------------------------
!
      subroutine mfo_hy (jcen,nbod,nbig,m,x,rcrit,a,stat)
!
      implicit none
      include 'mercury95.inc'
!
! Input/Output
      integer nbod, nbig, stat(nbod)
      real(8) jcen(3), m(nbod), x(3,nbod), a(3,nbod), rcrit(nbod)
!
! Local
      integer k
      real(8) aobl(3,NMAX),acen(3)
!
!------------------------------------------------------------------------------
!
! Initialize accelerations to zero
      do k = 1, nbod
        a(1,k) = 0.d0
        a(2,k) = 0.d0
        a(3,k) = 0.d0
      end do
!
! Calculate direct terms
      call mfo_drct (2,nbod,nbig,m,x,rcrit,a,stat)
!
! Add accelerations due to oblateness of the central body
      if (jcen(1).ne.0.or.jcen(2).ne.0.or.jcen(3).ne.0) then
        call mfo_obl (jcen,nbod,m,x,aobl,acen)
        do k = 2, nbod
          a(1,k) = a(1,k) + aobl(1,k) - acen(1)
          a(2,k) = a(2,k) + aobl(2,k) - acen(2)
          a(3,k) = a(3,k) + aobl(3,k) - acen(3)
        end do
      end if
!
!------------------------------------------------------------------------------
!
      return
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!     MFO_HKCE.FOR    (ErikSoft   27 February 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Calculates accelerations due to the Keplerian part of the Hamiltonian 
! of a hybrid symplectic integrator, when close encounters are taking place,
! for a set of NBOD bodies (NBIG of which are Big). Note that Small bodies
! do not interact with one another.
!
!------------------------------------------------------------------------------
!
      subroutine mfo_hkce (time,jcen,nbod,nbig,m,x,v,spin,rcrit,a,stat, &
       ngf,ngflag,opt,nce,ice,jce)
!
      implicit none
      include 'mercury95.inc'
!
! Input/Output
      integer nbod,nbig,stat(nbod),ngflag,opt(8),nce,ice(nce),jce(nce)
      real(8) time,jcen(3),rcrit(nbod),ngf(4,nbod),m(nbod)
      real(8) x(3,nbod),v(3,nbod),a(3,nbod),spin(3,nbod)
!
! Local
      integer i, j, k
      real(8) tmp2,dx,dy,dz,s,s_1,s2,s_3,faci,facj,rc,rc2,q,q2,q3,q4,q5
!
!------------------------------------------------------------------------------
!
! Initialize accelerations
      do j = 1, nbod
        a(1,j) = 0.d0
        a(2,j) = 0.d0
        a(3,j) = 0.d0
      end do
!
! Direct terms
      do k = 1, nce
        i = ice(k)
        j = jce(k)
        dx = x(1,j) - x(1,i)
        dy = x(2,j) - x(2,i)
        dz = x(3,j) - x(3,i)
        s2 = dx * dx  +  dy * dy  +  dz * dz
        rc = max (rcrit(i), rcrit(j))
        rc2 = rc * rc
!
        if (s2.lt.rc2) then
          s_1 = 1.d0 / sqrt(s2)
          s_3 = s_1 * s_1 * s_1
          if (s2.le.0.01*rc2) then
            tmp2 = s_3
          else
            s = 1.d0 / s_1
            q = (s - 0.1d0*rc) / (0.9d0 * rc)
            q2 = q * q
            q3 = q * q2
            q4 = q2 * q2
            q5 = q2 * q3
            tmp2 = (1.d0 - 10.d0*q3 + 15.d0*q4 - 6.d0*q5) * s_3
          end if
!
          faci = tmp2 * m(i)
          facj = tmp2 * m(j)
          a(1,j) = a(1,j)  -  faci * dx
          a(2,j) = a(2,j)  -  faci * dy
          a(3,j) = a(3,j)  -  faci * dz
          a(1,i) = a(1,i)  +  facj * dx
          a(2,i) = a(2,i)  +  facj * dy
          a(3,i) = a(3,i)  +  facj * dz
        end if
      end do
!
! Solar terms
      do i = 2, nbod
        s2 = x(1,i)*x(1,i) + x(2,i)*x(2,i) + x(3,i)*x(3,i)
        s_1 = 1.d0 / sqrt(s2)
        tmp2 = m(1) * s_1 * s_1 * s_1
        a(1,i) = a(1,i)  -  tmp2 * x(1,i)
        a(2,i) = a(2,i)  -  tmp2 * x(2,i)
        a(3,i) = a(3,i)  -  tmp2 * x(3,i)
      end do
!
!------------------------------------------------------------------------------
!
      return
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!     MFO_MVS.FOR    (ErikSoft   2 October 2000)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Calculates accelerations on a set of NBOD bodies (of which NBIG are Big)
! due to gravitational perturbations by all the other bodies.
! This routine is designed for use with a mixed-variable symplectic
! integrator using Jacobi coordinates.
!
! Based upon routines from Levison and Duncan's SWIFT integrator.
!
!------------------------------------------------------------------------------
!
      subroutine mfo_mvs (jcen,nbod,nbig,m,x,xj,a,stat)
!
      implicit none
      include 'mercury95.inc'
!
! Input/Output
      integer nbod, nbig, stat(nbod)
      real(8) jcen(3), m(nbod), x(3,nbod), xj(3,nbod), a(3,nbod)
!
! Local
      integer i,j,k,k1
      real(8) fac0,fac1,fac12,fac2,minside,dx,dy,dz,s_1,s2,s_3,faci,facj
      real(8) a0(3),a0tp(3),a1(3,NMAX),a2(3,NMAX),a3(3,NMAX),aobl(3,NMAX)
      real(8) r,r2,r3,rj,rj2,rj3,q,q2,q3,q4,q5,q6,q7,acen(3)
!
!------------------------------------------------------------------------------
!
! Initialize variables
      a0(1) = 0.d0
      a0(2) = 0.d0
      a0(3) = 0.d0
      a1(1,2) = 0.d0
      a1(2,2) = 0.d0
      a1(3,2) = 0.d0
      a2(1,2) = 0.d0
      a2(2,2) = 0.d0
      a2(3,2) = 0.d0
      minside = 0.d0
!
! Calculate acceleration terms
      do k = 3, nbig
        k1 = k - 1
        minside = minside + m(k1)
        r2   = x(1,k)  * x(1,k)  +  x(2,k) * x(2,k)  +  x(3,k) * x(3,k)
        rj2  = xj(1,k) * xj(1,k) + xj(2,k) * xj(2,k) + xj(3,k) * xj(3,k)
        r  = 1.d0 / sqrt(r2)
        rj = 1.d0 / sqrt(rj2)
        r3  = r  * r  * r
        rj3 = rj * rj * rj
!
        fac0 = m(k) * r3
        fac12 = m(1) * rj3
        fac2 = m(k) * fac12 / (minside + m(1))
        q = (r2 - rj2) * .5d0 / rj2
        q2 = q  * q
        q3 = q  * q2
        q4 = q2 * q2
        q5 = q2 * q3
        q6 = q3 * q3
        q7 = q3 * q4
        fac1 = 402.1875d0*q7 - 187.6875d0*q6 + 86.625d0*q5 &
            - 39.375d0*q4 + 17.5d0*q3 - 7.5d0*q2 + 3.d0*q - 1.d0
!
! Add to A0 term
        a0(1) = a0(1)  -  fac0 * x(1,k)
        a0(2) = a0(2)  -  fac0 * x(2,k)
        a0(3) = a0(3)  -  fac0 * x(3,k)
!
! Calculate A1 for this body
        a1(1,k) = fac12 * (xj(1,k) + fac1*x(1,k))
        a1(2,k) = fac12 * (xj(2,k) + fac1*x(2,k))
        a1(3,k) = fac12 * (xj(3,k) + fac1*x(3,k))
!
! Calculate A2 for this body
        a2(1,k) = a2(1,k1)  +  fac2 * xj(1,k)
        a2(2,k) = a2(2,k1)  +  fac2 * xj(2,k)
        a2(3,k) = a2(3,k1)  +  fac2 * xj(3,k)
      end do
!
      r2   = x(1,2)  * x(1,2)  +  x(2,2) * x(2,2)  +  x(3,2) * x(3,2)
      r  = 1.d0 / sqrt(r2)
      r3  = r  * r  * r
      fac0 = m(2) * r3
      a0tp(1) = a0(1)  -  fac0 * x(1,2)
      a0tp(2) = a0(2)  -  fac0 * x(2,2)
      a0tp(3) = a0(3)  -  fac0 * x(3,2)
!
! Calculate A3 (direct terms)
      do k = 2, nbod
        a3(1,k) = 0.d0
        a3(2,k) = 0.d0
        a3(3,k) = 0.d0
      end do
      do i = 2, nbig
        do j = i + 1, nbig
          dx = x(1,j) - x(1,i)
          dy = x(2,j) - x(2,i)
          dz = x(3,j) - x(3,i)
          s2 = dx*dx + dy*dy + dz*dz
          s_1 = 1.d0 / sqrt(s2)
          s_3 = s_1 * s_1 * s_1
          faci = m(i) * s_3
          facj = m(j) * s_3
          a3(1,j) = a3(1,j)  -  faci * dx
          a3(2,j) = a3(2,j)  -  faci * dy
          a3(3,j) = a3(3,j)  -  faci * dz
          a3(1,i) = a3(1,i)  +  facj * dx
          a3(2,i) = a3(2,i)  +  facj * dy
          a3(3,i) = a3(3,i)  +  facj * dz
        end do
!
        do j = nbig + 1, nbod
          dx = x(1,j) - x(1,i)
          dy = x(2,j) - x(2,i)
          dz = x(3,j) - x(3,i)
          s2 = dx*dx + dy*dy + dz*dz
          s_1 = 1.d0 / sqrt(s2)
          s_3 = s_1 * s_1 * s_1
          faci = m(i) * s_3
          a3(1,j) = a3(1,j)  -  faci * dx
          a3(2,j) = a3(2,j)  -  faci * dy
          a3(3,j) = a3(3,j)  -  faci * dz
        end do
      end do
!
! Big-body accelerations
      do k = 2, nbig
        a(1,k) = a0(1) + a1(1,k) + a2(1,k) + a3(1,k)
        a(2,k) = a0(2) + a1(2,k) + a2(2,k) + a3(2,k)
        a(3,k) = a0(3) + a1(3,k) + a2(3,k) + a3(3,k)
      end do
!
! Small-body accelerations
      do k = nbig + 1, nbod
        a(1,k) = a0tp(1) + a3(1,k)
        a(2,k) = a0tp(2) + a3(2,k)
        a(3,k) = a0tp(3) + a3(3,k)
      end do
!
! Correct for oblateness of the central body
      if (jcen(1).ne.0.or.jcen(2).ne.0.or.jcen(3).ne.0) then
        call mfo_obl (jcen,nbod,m,x,aobl,acen)
        do k = 2, nbod
          a(1,k) = a(1,k) + (aobl(1,k) - acen(1))
          a(2,k) = a(2,k) + (aobl(2,k) - acen(2))
          a(3,k) = a(3,k) + (aobl(3,k) - acen(3))
        end do
      end if
!
!------------------------------------------------------------------------------
!
      return
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MFO_NGF.FOR    (ErikSoft  29 November 1999)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Calculates accelerations on a set of NBOD bodies due to cometary
! non-gravitational jet forces. The positions and velocities are stored in
! arrays X, V with the format (x,y,z) and (vx,vy,vz) for each object in
! succession. The accelerations are stored in the array A (ax,ay,az). The
! non-gravitational accelerations follow a force law described by Marsden
! et al. (1973) Astron. J. 211-225, with magnitude determined by the
! parameters NGF(1,2,3) for each object.
!
! N.B. All coordinates and velocities must be with respect to central body!!!!
! ===
!------------------------------------------------------------------------------
!
      subroutine mfo_ngf (nbod,x,v,a,ngf)
!
      implicit none
      include 'mercury95.inc'
!
! Input/Output
      integer nbod
      real(8) x(3,nbod), v(3,nbod), a(3,nbod), ngf(4,nbod)
!
! Local
      integer j
      real(8) r2,r,rv,q,g,tx,ty,tz,nx,ny,nz,a1,a2,a3
!
!------------------------------------------------------------------------------
!
      do j = 2, nbod
        r2 = x(1,j)*x(1,j) + x(2,j)*x(2,j) +x(3,j)*x(3,j)
!
! Only calculate accelerations if body is close to the Sun (R < 9.36 AU), 
! or if the non-gravitational force parameters are exceptionally large.
        if (r2.lt.88.d0.or.abs(ngf(1,j)).gt.1d-7 &
         .or.abs(ngf(2,j)).gt.1d-7.or.abs(ngf(3,j)).gt.1d-7) then
          r = sqrt(r2)
          rv = x(1,j)*v(1,j) + x(2,j)*v(2,j) + x(3,j)*v(3,j)
!
! Calculate Q = R / R0, where R0 = 2.808 AU
          q = r * .3561253561253561d0
          g = .111262d0 * q**(-2.15d0) * (1.d0+q**5.093d0)**(-4.6142d0)
!
! Within-orbital-plane transverse vector components
          tx = r2*v(1,j) - rv*x(1,j)
          ty = r2*v(2,j) - rv*x(2,j)
          tz = r2*v(3,j) - rv*x(3,j)
!
! Orbit-normal vector components
          nx = x(2,j)*v(3,j) - x(3,j)*v(2,j)
          ny = x(3,j)*v(1,j) - x(1,j)*v(3,j)
          nz = x(1,j)*v(2,j) - x(2,j)*v(1,j)
!
! Multiplication factors
          a1 = ngf(1,j) * g / r
          a2 = ngf(2,j) * g / sqrt(tx*tx + ty*ty + tz*tz)
          a3 = ngf(3,j) * g / sqrt(nx*nx + ny*ny + nz*nz)
!
! X,Y and Z components of non-gravitational acceleration
          a(1,j) = a1*x(1,j) + a2*tx + a3*nx
          a(2,j) = a1*x(2,j) + a2*ty + a3*ny
          a(3,j) = a1*x(3,j) + a2*tz + a3*nz
        else
          a(1,j) = 0.d0
          a(2,j) = 0.d0
          a(3,j) = 0.d0
        end if
      end do
!
!------------------------------------------------------------------------------
!
      return
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!     MFO_OBL.FOR    (ErikSoft   2 October 2000)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Calculates barycentric accelerations of NBOD bodies due to oblateness of
! the central body. Also returns the corresponding barycentric acceleration
! of the central body.
!
! N.B. All coordinates must be with respect to the central body!!!!
! ===
!------------------------------------------------------------------------------
!
      subroutine mfo_obl (jcen,nbod,m,x,a,acen)
!
      implicit none
      include 'mercury95.inc'
!
! Input/Output
      integer nbod
      real(8) jcen(3), m(nbod), x(3,nbod), a(3,nbod), acen(3)
!
! Local
      integer i
      real(8) jr2,jr4,jr6,r2,r_1,r_2,r_3,u2,u4,u6,tmp1,tmp2,tmp3,tmp4
!
!------------------------------------------------------------------------------
!
      acen(1) = 0.d0
      acen(2) = 0.d0
      acen(3) = 0.d0
!
      do i = 2, nbod
!
! Calculate barycentric accelerations on the objects
        r2 = x(1,i)*x(1,i) + x(2,i)*x(2,i) + x(3,i)*x(3,i)
        r_1 = 1.d0 / sqrt(r2)
        r_2 = r_1 * r_1
        r_3 = r_2 * r_1
        jr2 = jcen(1) * r_2
        jr4 = jcen(2) * r_2 * r_2
        jr6 = jcen(3) * r_2 * r_2 * r_2
        u2 = x(3,i) * x(3,i) * r_2
        u4 = u2 * u2
        u6 = u4 * u2
!
        tmp1 = m(1) * r_3
        tmp2 =jr2*(7.5d0*u2 - 1.5d0) &
            +jr4*(39.375d0*u4 - 26.25d0*u2 + 1.875d0) &
            +jr6*(187.6875d0*u6 -216.5625d0*u4 +59.0625d0*u2 -2.1875d0)
        tmp3 = jr2*3.d0 + jr4*(17.5d0*u2 - 7.5d0) &
            + jr6*(86.625d0*u4 - 78.75d0*u2 + 13.125d0)
!
        a(1,i) = x(1,i) * tmp1 * tmp2
        a(2,i) = x(2,i) * tmp1 * tmp2
        a(3,i) = x(3,i) * tmp1 * (tmp2 - tmp3)
!
! Calculate barycentric accelerations on the central body
        tmp4 = m(i) / m(1)
        acen(1) = acen(1)  -  tmp4 * a(1,i)
        acen(2) = acen(2)  -  tmp4 * a(2,i)
        acen(3) = acen(3)  -  tmp4 * a(3,i)
      end do
!
!------------------------------------------------------------------------------
!
      return
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MFO_PN.FOR    (ErikSoft   3 October 2000)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! ****** To be completed at a later date ******
!
! Calculates post-Newtonian relativistic corrective accelerations for a set
! of NBOD bodies (NBIG of which are Big).
!
! This routine should not be called from the symplectic algorithm MAL_MVS
! or the conservative Bulirsch-Stoer algorithm MAL_BS2.
!
! N.B. All coordinates and velocities must be with respect to central body!!!!
! ===
!------------------------------------------------------------------------------
!
      subroutine mfo_pn (nbod,nbig,m,x,v,a)
!
      implicit none
      include 'mercury95.inc'
!
! Input/Output
      integer nbod, nbig
      real(8) m(nbod), x(3,nbod), v(3,nbod), a(3,nbod)
!
! Local
      integer j
!
!------------------------------------------------------------------------------
!
      do j = 1, nbod
        a(1,j) = 0.d0
        a(2,j) = 0.d0
        a(3,j) = 0.d0
      end do
!
!------------------------------------------------------------------------------
!
      return
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MFO_PR.FOR    (ErikSoft   3 October 2000)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! ****** To be completed at a later date ******
!
! Calculates radiation pressure and Poynting-Robertson drag for a set
! of NBOD bodies (NBIG of which are Big).
!
! This routine should not be called from the symplectic algorithm MAL_MVS
! or the conservative Bulirsch-Stoer algorithm MAL_BS2.
!
! N.B. All coordinates and velocities must be with respect to central body!!!!
! ===
!------------------------------------------------------------------------------
!
      subroutine mfo_pr (nbod,nbig,m,x,v,a,ngf)
!
      implicit none
      include 'mercury95.inc'
!
! Input/Output
      integer nbod, nbig
      real(8) m(nbod), x(3,nbod), v(3,nbod), a(3,nbod), ngf(4,nbod)
!
! Local
      integer j
!
!------------------------------------------------------------------------------
!
      do j = 1, nbod
        a(1,j) = 0.d0
        a(2,j) = 0.d0
        a(3,j) = 0.d0
      end do
!
!------------------------------------------------------------------------------
!
      return
      end
!

