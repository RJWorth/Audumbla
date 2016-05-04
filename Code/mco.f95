!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MCO_ACSH.FOR    (ErikSoft  2 March 1999)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Calculates inverse hyperbolic cosine of an angle X (in radians).
!
!------------------------------------------------------------------------------
!
      function mco_acsh (x)
!
      implicit none
!
! Input/Output
      real(8) x,mco_acsh
!
!------------------------------------------------------------------------------
!
      if (x.ge.1.d0) then
        mco_acsh = log (x + sqrt(x*x - 1.d0))
      else
        mco_acsh = 0.d0
      end if
!
!------------------------------------------------------------------------------
!
      return
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MCO_B2H.FOR    (ErikSoft   2 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Converts barycentric coordinates to coordinates with respect to the central
! body.
!
!------------------------------------------------------------------------------
!
      subroutine mco_b2h (time,jcen,nbod,nbig,h,m,x,v,xh,vh,ngf,ngflag, &
       opt)
!
      implicit none
!
! Input/Output
      integer nbod,nbig,ngflag,opt(8)
      real(8) time,jcen(3),h,m(nbod),x(3,nbod),v(3,nbod),xh(3,nbod)
      real(8) vh(3,nbod),ngf(4,nbod)
!
! Local
      integer j
!
!------------------------------------------------------------------------------
!
      do j = 2, nbod
        xh(1,j) = x(1,j) - x(1,1)
        xh(2,j) = x(2,j) - x(2,1)
        xh(3,j) = x(3,j) - x(3,1)
        vh(1,j) = v(1,j) - v(1,1)
        vh(2,j) = v(2,j) - v(2,1)
        vh(3,j) = v(3,j) - v(3,1)
      enddo
!
!------------------------------------------------------------------------------
!
      return
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MCO_DH2H.FOR    (ErikSoft   2 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Converts democratic heliocentric coordinates to coordinates with respect
! to the central body.
!
!------------------------------------------------------------------------------
!
      subroutine mco_dh2h (time,jcen,nbod,nbig,h,m,x,v,xh,vh,ngf,ngflag, &
       opt)
!
      implicit none
!
! Input/Output
      integer nbod,nbig,ngflag,opt(8)
      real(8) time,jcen(3),h,m(nbod),x(3,nbod),v(3,nbod),xh(3,nbod)
      real(8) vh(3,nbod),ngf(4,nbod)
!
! Local
      integer j
      real(8) mvsum(3),temp
!
!------------------------------------------------------------------------------
!
      mvsum(1) = 0.d0
      mvsum(2) = 0.d0
      mvsum(3) = 0.d0
!
      do j = 2, nbod
        xh(1,j) = x(1,j)
        xh(2,j) = x(2,j)
        xh(3,j) = x(3,j)
        mvsum(1) = mvsum(1)  +  m(j) * v(1,j)
        mvsum(2) = mvsum(2)  +  m(j) * v(2,j)
        mvsum(3) = mvsum(3)  +  m(j) * v(3,j)
      end do
!
      temp = 1.d0 / m(1)
      mvsum(1) = temp * mvsum(1)
      mvsum(2) = temp * mvsum(2)
      mvsum(3) = temp * mvsum(3)
!
      do j = 2, nbod
        vh(1,j) = v(1,j) + mvsum(1)
        vh(2,j) = v(2,j) + mvsum(2)
        vh(3,j) = v(3,j) + mvsum(3)
      end do
!
!------------------------------------------------------------------------------
!
      return
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MCO_MVS2H.FOR    (ErikSoft   28 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Applies a symplectic corrector, which converts coordinates for a second-
! order mixed-variable symplectic integrator to ones with respect to the
! central body.
!
!------------------------------------------------------------------------------
!
      subroutine mco_mvs2h (time,jcen,nbod,nbig,h,m,x,v,xh,vh,ngf, &
       ngflag,opt)
!
      implicit none
      include 'mercury95.inc'
!
! Input/Output
      integer nbod,nbig,ngflag,opt(8)
      real(8) time,jcen(3),h,m(nbod),x(3,nbod),v(3,nbod),xh(3,nbod)
      real(8) vh(3,nbod),ngf(4,nbod)
!
! Local
      integer j,k,iflag,stat(NMAX)
      real(8) minside,msofar,gm(NMAX),a(3,NMAX),xj(3,NMAX),vj(3,NMAX)
      real(8) ha(2),hb(2),rt10,angf(3,NMAX),ausr(3,NMAX)
!
!------------------------------------------------------------------------------
!
      rt10 = sqrt(10.d0)
      ha(1) =  h * rt10 * 3.d0 / 10.d0
      hb(1) = -h * rt10 / 72.d0
      ha(2) =  h * rt10 / 5.d0
      hb(2) =  h * rt10 / 24.d0
      do j = 2, nbod
        angf(1,j) = 0.d0
        angf(2,j) = 0.d0
        angf(3,j) = 0.d0
        ausr(1,j) = 0.d0
        ausr(2,j) = 0.d0
        ausr(3,j) = 0.d0
      end do
      call mco_iden (jcen,nbod,nbig,h,m,x,v,xh,vh)
!
! Calculate effective central masses for Kepler drifts
      minside = m(1)
      do j = 2, nbig
        msofar = minside + m(j)
        gm(j) = m(1) * msofar / minside
        minside = msofar
      end do
!
! Two step corrector
      do k = 1, 2
!
! Advance Keplerian Hamiltonian (Jacobi/helio coords for Big/Small bodies)
        call mco_h2j(jcen,nbig,nbig,h,m,xh,vh,xj,vj)
        do j = 2, nbig
          iflag = 0
          call drift_one (gm(j),xj(1,j),xj(2,j),xj(3,j),vj(1,j), &
           vj(2,j),vj(3,j),ha(k),iflag)
        end do
        do j = nbig + 1, nbod
          iflag = 0
          call drift_one (m(1),xh(1,j),xh(2,j),xh(3,j),vh(1,j),vh(2,j), &
           vh(3,j),ha(k),iflag)
        end do
!
! Advance Interaction Hamiltonian
        call mco_j2h(time,jcen,nbig,nbig,h,m,xj,vj,xh,vh,ngf,ngflag,opt)
        call mfo_mvs (jcen,nbod,nbig,m,xh,xj,a,stat)
!
! If required, apply non-gravitational and user-defined forces
        if (opt(8).eq.1) call mfo_user (time,jcen,nbod,nbig,m,xh,vh, &
         ausr)
        if (ngflag.eq.1.or.ngflag.eq.3) call mfo_ngf (nbod,xh,vh,angf, &
         ngf)
!
        do j = 2, nbod
          vh(1,j) = vh(1,j)  -  hb(k) * (angf(1,j) + ausr(1,j) + a(1,j))
          vh(2,j) = vh(2,j)  -  hb(k) * (angf(2,j) + ausr(2,j) + a(2,j))
          vh(3,j) = vh(3,j)  -  hb(k) * (angf(3,j) + ausr(3,j) + a(3,j))
        end do
!
! Advance Keplerian Hamiltonian (Jacobi/helio coords for Big/Small bodies)
        call mco_h2j(jcen,nbig,nbig,h,m,xh,vh,xj,vj)
        do j = 2, nbig
          iflag = 0
          call drift_one (gm(j),xj(1,j),xj(2,j),xj(3,j),vj(1,j), &
           vj(2,j),vj(3,j),-2.d0*ha(k),iflag)
        end do
        do j = nbig + 1, nbod
          iflag = 0
         call drift_one (m(1),xh(1,j),xh(2,j),xh(3,j),vh(1,j),vh(2,j), &
          vh(3,j),-2.d0*ha(k),iflag)
        end do
!
! Advance Interaction Hamiltonian
        call mco_j2h(time,jcen,nbig,nbig,h,m,xj,vj,xh,vh,ngf,ngflag,opt)
        call mfo_mvs (jcen,nbod,nbig,m,xh,xj,a,stat)
!
! If required, apply non-gravitational and user-defined forces
        if (opt(8).eq.1) call mfo_user (time,jcen,nbod,nbig,m,xh,vh, &
         ausr)
        if (ngflag.eq.1.or.ngflag.eq.3) call mfo_ngf (nbod,xh,vh,angf, &
         ngf)
!
        do j = 2, nbod
          vh(1,j) = vh(1,j)  +  hb(k) * (angf(1,j) + ausr(1,j) + a(1,j))
          vh(2,j) = vh(2,j)  +  hb(k) * (angf(2,j) + ausr(2,j) + a(2,j))
          vh(3,j) = vh(3,j)  +  hb(k) * (angf(3,j) + ausr(3,j) + a(3,j))
        end do
!
! Advance Keplerian Hamiltonian (Jacobi/helio coords for Big/Small bodies)
        call mco_h2j(jcen,nbig,nbig,h,m,xh,vh,xj,vj)
        do j = 2, nbig
          iflag = 0
          call drift_one (gm(j),xj(1,j),xj(2,j),xj(3,j),vj(1,j), &
           vj(2,j),vj(3,j),ha(k),iflag)
        end do
        do j = nbig + 1, nbod
          iflag = 0
          call drift_one (m(1),xh(1,j),xh(2,j),xh(3,j),vh(1,j),vh(2,j), &
           vh(3,j),ha(k),iflag)
        end do
        call mco_j2h(time,jcen,nbig,nbig,h,m,xj,vj,xh,vh,ngf,ngflag,opt)
      end do
!
!------------------------------------------------------------------------------
!
      return
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MCO_H2MVS.FOR    (ErikSoft   28 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Applies an inverse symplectic corrector, which converts coordinates with
! respect to the central body to integrator coordinates for a second-order
! mixed-variable symplectic integrator.
!
!------------------------------------------------------------------------------
!
      subroutine mco_h2mvs (time,jcen,nbod,nbig,h,m,xh,vh,x,v,ngf, &
       ngflag,opt)
!
      implicit none
      include 'mercury95.inc'
!
! Input/Output
      integer nbod,nbig,ngflag,opt(8)
      real(8) time,jcen(3),h,m(nbod),xh(3,nbod),vh(3,nbod),x(3,nbod)
      real(8) v(3,nbod),ngf(4,nbod)
!
! Local
      integer j,k,iflag,stat(NMAX)
      real(8) minside,msofar,gm(NMAX),a(3,NMAX),xj(3,NMAX),vj(3,NMAX)
      real(8) ha(2),hb(2),rt10,angf(3,NMAX),ausr(3,NMAX)
!
!------------------------------------------------------------------------------
!
      rt10 = sqrt(10.d0)
      ha(1) = -h * rt10 / 5.d0
      hb(1) = -h * rt10 / 24.d0
      ha(2) = -h * rt10 * 3.d0 / 10.d0
      hb(2) =  h * rt10 / 72.d0
      do j = 2, nbod
        angf(1,j) = 0.d0
        angf(2,j) = 0.d0
        angf(3,j) = 0.d0
        ausr(1,j) = 0.d0
        ausr(2,j) = 0.d0
        ausr(3,j) = 0.d0
      end do
      call mco_iden (jcen,nbod,nbig,h,m,xh,vh,x,v)
!
! Calculate effective central masses for Kepler drifts
      minside = m(1)
      do j = 2, nbig
        msofar = minside + m(j)
        gm(j) = m(1) * msofar / minside
        minside = msofar
      end do
!
      do k = 1, 2
!
! Advance Keplerian Hamiltonian (Jacobi/helio coords for Big/Small bodies)
        call mco_h2j (jcen,nbig,nbig,h,m,x,v,xj,vj)
        do j = 2, nbig
          iflag = 0
          call drift_one (gm(j),xj(1,j),xj(2,j),xj(3,j),vj(1,j), &
           vj(2,j),vj(3,j),ha(k),iflag)
        end do
        do j = nbig + 1, nbod
          iflag = 0
          call drift_one (m(1),x(1,j),x(2,j),x(3,j),v(1,j),v(2,j), &
           v(3,j),ha(k),iflag)
        end do
!
! Advance Interaction Hamiltonian
        call mco_j2h (time,jcen,nbig,nbig,h,m,xj,vj,x,v,ngf,ngflag,opt)
        call mfo_mvs (jcen,nbod,nbig,m,x,xj,a,stat)
!
! If required, apply non-gravitational and user-defined forces
        if (opt(8).eq.1) call mfo_user (time,jcen,nbod,nbig,m,x,v,ausr)
        if (ngflag.eq.1.or.ngflag.eq.3) call mfo_ngf (nbod,x,v,angf,ngf)
!
        do j = 2, nbod
          v(1,j) = v(1,j)  +  hb(k) * (angf(1,j) + ausr(1,j) + a(1,j))
          v(2,j) = v(2,j)  +  hb(k) * (angf(2,j) + ausr(2,j) + a(2,j))
          v(3,j) = v(3,j)  +  hb(k) * (angf(3,j) + ausr(3,j) + a(3,j))
        end do
!
! Advance Keplerian Hamiltonian (Jacobi/helio coords for Big/Small bodies)
        call mco_h2j (jcen,nbig,nbig,h,m,x,v,xj,vj)
        do j = 2, nbig
          iflag = 0
          call drift_one (gm(j),xj(1,j),xj(2,j),xj(3,j),vj(1,j), &
           vj(2,j),vj(3,j),-2.d0*ha(k),iflag)
        end do
        do j = nbig + 1, nbod
          iflag = 0
          call drift_one (m(1),x(1,j),x(2,j),x(3,j),v(1,j),v(2,j), &
           v(3,j),-2.d0*ha(k),iflag)
        end do
!
! Advance Interaction Hamiltonian
        call mco_j2h (time,jcen,nbig,nbig,h,m,xj,vj,x,v,ngf,ngflag,opt)
        call mfo_mvs (jcen,nbod,nbig,m,x,xj,a,stat)
!
! If required, apply non-gravitational and user-defined forces
        if (opt(8).eq.1) call mfo_user (time,jcen,nbod,nbig,m,x,v,ausr)
        if (ngflag.eq.1.or.ngflag.eq.3) call mfo_ngf (nbod,x,v,angf,ngf)
!
        do j = 2, nbod
          v(1,j) = v(1,j)  -  hb(k) * (angf(1,j) + ausr(1,j) + a(1,j))
          v(2,j) = v(2,j)  -  hb(k) * (angf(2,j) + ausr(2,j) + a(2,j))
          v(3,j) = v(3,j)  -  hb(k) * (angf(3,j) + ausr(3,j) + a(3,j))
        end do
!
! Advance Keplerian Hamiltonian (Jacobi/helio coords for Big/Small bodies)
        call mco_h2j (jcen,nbig,nbig,h,m,x,v,xj,vj)
        do j = 2, nbig
          iflag = 0
          call drift_one (gm(j),xj(1,j),xj(2,j),xj(3,j),vj(1,j), &
           vj(2,j),vj(3,j),ha(k),iflag)
        end do
        do j = nbig + 1, nbod
          iflag = 0
          call drift_one (m(1),x(1,j),x(2,j),x(3,j),v(1,j),v(2,j), &
           v(3,j),ha(k),iflag)
        end do
        call mco_j2h (time,jcen,nbig,nbig,h,m,xj,vj,x,v,ngf,ngflag,opt)
      end do
!
!------------------------------------------------------------------------------
!
      return
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MCO_H2DH.FOR    (ErikSoft   2 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Convert coordinates with respect to the central body to democratic
! heliocentric coordinates.
!
!------------------------------------------------------------------------------
!
      subroutine mco_h2dh (time,jcen,nbod,nbig,h,m,xh,vh,x,v,ngf,ngflag, &
       opt)
!
      implicit none
!
! Input/Output
      integer nbod,nbig,ngflag,opt(8)
      real(8) time,jcen(3),h,m(nbod),xh(3,nbod),vh(3,nbod),x(3,nbod)
      real(8) v(3,nbod),ngf(4,nbod)
!
! Local
      integer j
      real(8) mtot,temp,mvsum(3)
!
!------------------------------------------------------------------------------
!
      mtot = 0.d0
      mvsum(1) = 0.d0
      mvsum(2) = 0.d0
      mvsum(3) = 0.d0
!
      do j = 2, nbod
        x(1,j) = xh(1,j)
        x(2,j) = xh(2,j)
        x(3,j) = xh(3,j)
        mtot = mtot + m(j)
        mvsum(1) = mvsum(1)  +  m(j) * vh(1,j)
        mvsum(2) = mvsum(2)  +  m(j) * vh(2,j)
        mvsum(3) = mvsum(3)  +  m(j) * vh(3,j)
      end do
!
      temp = 1.d0 / (m(1) + mtot)
!
      mvsum(1) = temp * mvsum(1)
      mvsum(2) = temp * mvsum(2)
      mvsum(3) = temp * mvsum(3)
!
      do j = 2, nbod
        v(1,j) = vh(1,j) - mvsum(1)
        v(2,j) = vh(2,j) - mvsum(2)
        v(3,j) = vh(3,j) - mvsum(3)
      end do
!
!------------------------------------------------------------------------------
!
      return
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MCO_J2H.FOR    (ErikSoft   2 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Converts Jacobi coordinates to coordinates with respect to the central
! body.
!
! N.B. The Jacobi coordinates of the small bodies are assumed to be equal
! ===  to their coordinates with respect to the central body.
!
!------------------------------------------------------------------------------
!
      subroutine mco_j2h (time,jcen,nbod,nbig,h,m,x,v,xh,vh,ngf,ngflag, &
       opt)
!
      implicit none
!
! Input/Output
      integer nbod,nbig,ngflag,opt(8)
      real(8) time,jcen(3),h,m(nbod),x(3,nbod),v(3,nbod),xh(3,nbod)
      real(8) vh(3,nbod),ngf(4,nbod)
!
! Local
      integer j
      real(8) mtot, mx, my, mz, mu, mv, mw, temp
!
!------------------------------------------------------------------------------
!
      xh(1,2) = x(1,2)
      xh(2,2) = x(2,2)
      xh(3,2) = x(3,2)
      vh(1,2) = v(1,2)
      vh(2,2) = v(2,2)
      vh(3,2) = v(3,2)
      mtot = m(2)
      temp = m(2) / (mtot + m(1))
      mx = temp * x(1,2)
      my = temp * x(2,2)
      mz = temp * x(3,2)
      mu = temp * v(1,2)
      mv = temp * v(2,2)
      mw = temp * v(3,2)
!
      do j = 3, nbig - 1
        xh(1,j) = x(1,j) + mx
        xh(2,j) = x(2,j) + my
        xh(3,j) = x(3,j) + mz
        vh(1,j) = v(1,j) + mu
        vh(2,j) = v(2,j) + mv
        vh(3,j) = v(3,j) + mw
        mtot = mtot + m(j)
        temp = m(j) / (mtot + m(1))
        mx = mx  +  temp * x(1,j)
        my = my  +  temp * x(2,j)
        mz = mz  +  temp * x(3,j)
        mu = mu  +  temp * v(1,j)
        mv = mv  +  temp * v(2,j)
        mw = mw  +  temp * v(3,j)
      enddo
!
      if (nbig.gt.2) then
        xh(1,nbig) = x(1,nbig) + mx
        xh(2,nbig) = x(2,nbig) + my
        xh(3,nbig) = x(3,nbig) + mz
        vh(1,nbig) = v(1,nbig) + mu
        vh(2,nbig) = v(2,nbig) + mv
        vh(3,nbig) = v(3,nbig) + mw
      end if
!
      do j = nbig + 1, nbod
        xh(1,j) = x(1,j)
        xh(2,j) = x(2,j)
        xh(3,j) = x(3,j)
        vh(1,j) = v(1,j)
        vh(2,j) = v(2,j)
        vh(3,j) = v(3,j)
      end do
!
!------------------------------------------------------------------------------
!
      return
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MCO_X2A.FOR    (ErikSoft   4 October 2000)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Calculates an object's orbital semi-major axis given its Cartesian coords.
!
!------------------------------------------------------------------------------
!
      subroutine mco_x2a (gm,x,y,z,u,v,w,a,r,v2)
!
      implicit none
!
! Input/Output
      real(8) gm,x,y,z,u,v,w,a,r,v2
!
!------------------------------------------------------------------------------
!
      r  = sqrt(x * x  +  y * y  +  z * z)
      v2 =      u * u  +  v * v  +  w * w
      a  = gm * r / (2.d0 * gm  -  r * v2)
!
!------------------------------------------------------------------------------
!
      return
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MCO_X2OV.FOR    (ErikSoft   20 February 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Calculates output variables for an object given its coordinates and
! velocities. The output variables are:
!  r = the radial distance
!  theta = polar angle
!  phi = azimuthal angle
!  fv = 1 / [1 + 2(ke/be)^2], where be and ke are the object's binding and
!                             kinetic energies. (Note that 0 < fv < 1).
!  vtheta = polar angle of velocity vector
!  vphi = azimuthal angle of the velocity vector
!
!------------------------------------------------------------------------------
!
      subroutine mco_x2ov (rcen,rmax,mcen,m,x,y,z,u,v,w,fr,theta,phi,fv, &
       vtheta,vphi)
!
      implicit none
      include 'mercury95.inc'
!
! Input/Output
      real(8) rcen,rmax,mcen,m,x,y,z,u,v,w,fr,theta,phi,fv,vtheta,vphi
!
! Local
      real(8) r,v2,v1,be,ke,temp
!
!------------------------------------------------------------------------------
!
        r = sqrt(x*x + y*y + z*z)
        v2 =     u*u + v*v + w*w
        v1 = sqrt(v2)
        be = (mcen + m) / r
        ke = .5d0 * v2
!
        fr = log10 (min(max(r, rcen), rmax) / rcen)
        temp = ke / be
        fv = 1.d0 / (1.d0 + 2.d0*temp*temp)
!
        theta  = mod (acos (z / r) + TWOPI, TWOPI)
        vtheta = mod (acos (w / v1) + TWOPI, TWOPI)
        phi  = mod (atan2 (y, x) + TWOPI, TWOPI)
        vphi = mod (atan2 (v, u) + TWOPI, TWOPI)
!
!------------------------------------------------------------------------------
!
      return
      end
!

