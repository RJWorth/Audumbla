!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MDT_BS1.FOR    (ErikSoft   2 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Integrates NBOD bodies (of which NBIG are Big) for one timestep H0
! using the Bulirsch-Stoer method. The accelerations are calculated using the 
! subroutine FORCE. The accuracy of the step is approximately determined 
! by the tolerance parameter TOL.
!
! N.B. Input/output must be in coordinates with respect to the central body.
! ===
!
!------------------------------------------------------------------------------
!
      subroutine mdt_bs1 (time,h0,hdid,tol,jcen,nbod,nbig,mass,x0,v0,s, &
       rphys,rcrit,ngf,stat,dtflag,ngflag,opt,nce,ice,jce,force)
!
      implicit none
      include 'mercury95.inc'
!
      real(8) SHRINK,GROW
      parameter (SHRINK=.55d0,GROW=1.3d0)
!
! Input/Output
      integer nbod, nbig, opt(8), stat(nbod), dtflag, ngflag
      integer nce, ice(nce), jce(nce)
      real(8) time,h0,hdid,tol,jcen(3),mass(nbod),x0(3,nbod),v0(3,nbod)
      real(8) s(3,nbod),ngf(4,nbod),rphys(nbod),rcrit(nbod)
      external force
!
! Local
      integer j, j1, k, n
      real(8) tmp0,tmp1,tmp2,errmax,tol2,h,hx2,h2(8)
      real(8) x(3,NMAX),v(3,NMAX),xend(3,NMAX),vend(3,NMAX)
      real(8) a(3,NMAX),a0(3,NMAX),d(6,NMAX,8),xscal(NMAX),vscal(NMAX)
!
!------------------------------------------------------------------------------
!
      tol2 = tol * tol
!
! Calculate arrays used to scale the relative error (R^2 for position and
! V^2 for velocity).
      do k = 2, nbod
        tmp1 = x0(1,k)*x0(1,k) + x0(2,k)*x0(2,k) + x0(3,k)*x0(3,k)
        tmp2 = v0(1,k)*v0(1,k) + v0(2,k)*v0(2,k) + v0(3,k)*v0(3,k)
        xscal(k) = 1.d0 / tmp1
        vscal(k) = 1.d0 / tmp2
      end do
!
! Calculate accelerations at the start of the step
      call force (time,jcen,nbod,nbig,mass,x0,v0,s,rcrit,a0,stat,ngf, &
       ngflag,opt,nce,ice,jce)
!
 100  continue
!
! For each value of N, do a modified-midpoint integration with 2N substeps
      do n = 1, 8
        h = h0 / (2.d0 * float(n))
        h2(n) = .25d0 / (n*n)
        hx2 = h * 2.d0
!
        do k = 2, nbod
          x(1,k) = x0(1,k) + h*v0(1,k)
          x(2,k) = x0(2,k) + h*v0(2,k)
          x(3,k) = x0(3,k) + h*v0(3,k)
          v(1,k) = v0(1,k) + h*a0(1,k)
          v(2,k) = v0(2,k) + h*a0(2,k)
          v(3,k) = v0(3,k) + h*a0(3,k)
        end do
        call force (time,jcen,nbod,nbig,mass,x,v,s,rcrit,a,stat,ngf, &
         ngflag,opt,nce,ice,jce)
        do k = 2, nbod
          xend(1,k) = x0(1,k) + hx2*v(1,k)
          xend(2,k) = x0(2,k) + hx2*v(2,k)
          xend(3,k) = x0(3,k) + hx2*v(3,k)
          vend(1,k) = v0(1,k) + hx2*a(1,k)
          vend(2,k) = v0(2,k) + hx2*a(2,k)
          vend(3,k) = v0(3,k) + hx2*a(3,k)
        end do
!
        do j = 2, n
          call force (time,jcen,nbod,nbig,mass,xend,vend,s,rcrit,a,stat, &
           ngf,ngflag,opt,nce,ice,jce)
          do k = 2, nbod
            x(1,k) = x(1,k) + hx2*vend(1,k)
            x(2,k) = x(2,k) + hx2*vend(2,k)
            x(3,k) = x(3,k) + hx2*vend(3,k)
            v(1,k) = v(1,k) + hx2*a(1,k)
            v(2,k) = v(2,k) + hx2*a(2,k)
            v(3,k) = v(3,k) + hx2*a(3,k)
          end do
          call force (time,jcen,nbod,nbig,mass,x,v,s,rcrit,a,stat,ngf, &
           ngflag,opt,nce,ice,jce)
          do k = 2, nbod
            xend(1,k) = xend(1,k) + hx2*v(1,k)
            xend(2,k) = xend(2,k) + hx2*v(2,k)
            xend(3,k) = xend(3,k) + hx2*v(3,k)
            vend(1,k) = vend(1,k) + hx2*a(1,k)
            vend(2,k) = vend(2,k) + hx2*a(2,k)
            vend(3,k) = vend(3,k) + hx2*a(3,k)
          end do
        end do
!
        call force (time,jcen,nbod,nbig,mass,xend,vend,s,rcrit,a,stat, &
         ngf,ngflag,opt,nce,ice,jce)
!
        do k = 2, nbod
          d(1,k,n) = .5d0*(xend(1,k) + x(1,k) + h*vend(1,k))
          d(2,k,n) = .5d0*(xend(2,k) + x(2,k) + h*vend(2,k))
          d(3,k,n) = .5d0*(xend(3,k) + x(3,k) + h*vend(3,k))
          d(4,k,n) = .5d0*(vend(1,k) + v(1,k) + h*a(1,k))
          d(5,k,n) = .5d0*(vend(2,k) + v(2,k) + h*a(2,k))
          d(6,k,n) = .5d0*(vend(3,k) + v(3,k) + h*a(3,k))
        end do
!
! Update the D array, used for polynomial extrapolation
        do j = n - 1, 1, -1
          j1 = j + 1
          tmp0 = 1.d0 / (h2(j) - h2(n))
          tmp1 = tmp0 * h2(j1)
          tmp2 = tmp0 * h2(n)
          do k = 2, nbod
            d(1,k,j) = tmp1 * d(1,k,j1)  -  tmp2 * d(1,k,j)
            d(2,k,j) = tmp1 * d(2,k,j1)  -  tmp2 * d(2,k,j)
            d(3,k,j) = tmp1 * d(3,k,j1)  -  tmp2 * d(3,k,j)
            d(4,k,j) = tmp1 * d(4,k,j1)  -  tmp2 * d(4,k,j)
            d(5,k,j) = tmp1 * d(5,k,j1)  -  tmp2 * d(5,k,j)
            d(6,k,j) = tmp1 * d(6,k,j1)  -  tmp2 * d(6,k,j)
          end do
        end do
!
! After several integrations, test the relative error on extrapolated values
        if (n.gt.3) then
          errmax = 0.d0
!
! Maximum relative position and velocity errors (last D term added)
          do k = 2, nbod
            tmp1 = max( d(1,k,1)*d(1,k,1), d(2,k,1)*d(2,k,1), &
                       d(3,k,1)*d(3,k,1) )
            tmp2 = max( d(4,k,1)*d(4,k,1), d(5,k,1)*d(5,k,1), &
                       d(6,k,1)*d(6,k,1) )
            errmax = max(errmax, tmp1*xscal(k), tmp2*vscal(k))
          end do
!
! If error is smaller than TOL, update position and velocity arrays, and exit
          if (errmax.le.tol2) then
            do k = 2, nbod
              x0(1,k) = d(1,k,1)
              x0(2,k) = d(2,k,1)
              x0(3,k) = d(3,k,1)
              v0(1,k) = d(4,k,1)
              v0(2,k) = d(5,k,1)
              v0(3,k) = d(6,k,1)
            end do
!
            do j = 2, n
              do k = 2, nbod
                x0(1,k) = x0(1,k) + d(1,k,j)
                x0(2,k) = x0(2,k) + d(2,k,j)
                x0(3,k) = x0(3,k) + d(3,k,j)
                v0(1,k) = v0(1,k) + d(4,k,j)
                v0(2,k) = v0(2,k) + d(5,k,j)
                v0(3,k) = v0(3,k) + d(6,k,j)
              end do
            end do
!
! Save the actual stepsize used
            hdid = h0
!
! Recommend a new stepsize for the next call to this subroutine
            if (n.eq.8) h0 = h0 * SHRINK
            if (n.lt.7) h0 = h0 * GROW
            return
          end if
        end if
!
      end do
!
! If errors were too large, redo the step with half the previous step size.
      h0 = h0 * .5d0
      goto 100
!
!------------------------------------------------------------------------------
!
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MDT_BS2.FOR    (ErikSoft   2 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Integrates NBOD bodies (of which NBIG are Big) for one timestep H0
! using the Bulirsch-Stoer method. The accelerations are calculated using the 
! subroutine FORCE. The accuracy of the step is approximately determined 
! by the tolerance parameter TOL.
!
! N.B. This version only works for conservative systems (i.e. force is a
! ===  function of position only) !!!! Hence, non-gravitational forces
!      and post-Newtonian corrections cannot be used.
!
! N.B. Input/output must be in coordinates with respect to the central body.
! ===
!
!------------------------------------------------------------------------------
!
      subroutine mdt_bs2 (time,h0,hdid,tol,jcen,nbod,nbig,mass,x0,v0,s, &
       rphys,rcrit,ngf,stat,dtflag,ngflag,opt,nce,ice,jce,force)
!
      implicit none
      include 'mercury95.inc'
!
      real(8) SHRINK,GROW
      parameter (SHRINK=.55d0,GROW=1.3d0)
!
! Input/Output
      integer nbod, nbig, opt(8), stat(nbod), dtflag, ngflag
      real(8) time,h0,hdid,tol,jcen(3),mass(nbod),x0(3,nbod),v0(3,nbod)
      real(8) s(3,nbod),ngf(4,nbod),rphys(nbod),rcrit(nbod)
      integer nce,ice(nce),jce(nce)
      external force
!
! Local
      integer j,j1,k,n
      real(8) tmp0,tmp1,tmp2,errmax,tol2,h,h2(12),hby2,h2by2
      real(8) xend(3,NMAX),b(3,NMAX),c(3,NMAX)
      real(8) a(3,NMAX),a0(3,NMAX),d(6,NMAX,12),xscal(NMAX),vscal(NMAX)
!
!------------------------------------------------------------------------------
!
      tol2 = tol * tol
!
! Calculate arrays used to scale the relative error (R^2 for position and
! V^2 for velocity).
      do k = 2, nbod
        tmp1 = x0(1,k)*x0(1,k) + x0(2,k)*x0(2,k) + x0(3,k)*x0(3,k)
        tmp2 = v0(1,k)*v0(1,k) + v0(2,k)*v0(2,k) + v0(3,k)*v0(3,k)
        xscal(k) = 1.d0 / tmp1
        vscal(k) = 1.d0 / tmp2
      end do
!
! Calculate accelerations at the start of the step
      call force (time,jcen,nbod,nbig,mass,x0,v0,s,rcrit,a0,stat,ngf, &
       ngflag,opt,nce,ice,jce)
!
 100  continue
!
! For each value of N, do a modified-midpoint integration with N substeps
      do n = 1, 12
        h = h0 / (dble(n))
        hby2  = .5d0 * h
        h2(n) = h * h
        h2by2 = .5d0 * h2(n)
!
        do k = 2, nbod
          b(1,k) = .5d0*a0(1,k)
          b(2,k) = .5d0*a0(2,k)
          b(3,k) = .5d0*a0(3,k)
          c(1,k) = 0.d0
          c(2,k) = 0.d0
          c(3,k) = 0.d0
          xend(1,k) = h2by2 * a0(1,k)  +  h * v0(1,k)  +  x0(1,k)
          xend(2,k) = h2by2 * a0(2,k)  +  h * v0(2,k)  +  x0(2,k)
          xend(3,k) = h2by2 * a0(3,k)  +  h * v0(3,k)  +  x0(3,k)
        end do
!
        do j = 2, n
          call force (time,jcen,nbod,nbig,mass,xend,v0,s,rcrit,a,stat, &
           ngf,ngflag,opt,nce,ice,jce)
          tmp0 = h * dble(j)
          do k = 2, nbod
            b(1,k) = b(1,k) + a(1,k)
            b(2,k) = b(2,k) + a(2,k)
            b(3,k) = b(3,k) + a(3,k)
            c(1,k) = c(1,k) + b(1,k)
            c(2,k) = c(2,k) + b(2,k)
            c(3,k) = c(3,k) + b(3,k)
            xend(1,k) = h2(n)*c(1,k) + h2by2*a0(1,k) + tmp0*v0(1,k) &
                     + x0(1,k)
            xend(2,k) = h2(n)*c(2,k) + h2by2*a0(2,k) + tmp0*v0(2,k) &
                     + x0(2,k)
            xend(3,k) = h2(n)*c(3,k) + h2by2*a0(3,k) + tmp0*v0(3,k) &
                     + x0(3,k)
          end do
        end do
!
        call force (time,jcen,nbod,nbig,mass,xend,v0,s,rcrit,a,stat,ngf, &
         ngflag,opt,nce,ice,jce)
!
        do k = 2, nbod
          d(1,k,n) = xend(1,k)
          d(2,k,n) = xend(2,k)
          d(3,k,n) = xend(3,k)
          d(4,k,n) = h*b(1,k) + hby2*a(1,k) + v0(1,k)
          d(5,k,n) = h*b(2,k) + hby2*a(2,k) + v0(2,k)
          d(6,k,n) = h*b(3,k) + hby2*a(3,k) + v0(3,k)
        end do
!
! Update the D array, used for polynomial extrapolation
        do j = n - 1, 1, -1
          j1 = j + 1
          tmp0 = 1.d0 / (h2(j) - h2(n))
          tmp1 = tmp0 * h2(j1)
          tmp2 = tmp0 * h2(n)
          do k = 2, nbod
            d(1,k,j) = tmp1 * d(1,k,j1)  -  tmp2 * d(1,k,j)
            d(2,k,j) = tmp1 * d(2,k,j1)  -  tmp2 * d(2,k,j)
            d(3,k,j) = tmp1 * d(3,k,j1)  -  tmp2 * d(3,k,j)
            d(4,k,j) = tmp1 * d(4,k,j1)  -  tmp2 * d(4,k,j)
            d(5,k,j) = tmp1 * d(5,k,j1)  -  tmp2 * d(5,k,j)
            d(6,k,j) = tmp1 * d(6,k,j1)  -  tmp2 * d(6,k,j)
          end do
        end do
!
! After several integrations, test the relative error on extrapolated values
        if (n.gt.3) then
          errmax = 0.d0
!
! Maximum relative position and velocity errors (last D term added)
          do k = 2, nbod
            tmp1 = max( d(1,k,1)*d(1,k,1), d(2,k,1)*d(2,k,1), &
                       d(3,k,1)*d(3,k,1) )
            tmp2 = max( d(4,k,1)*d(4,k,1), d(5,k,1)*d(2,k,1), &
                       d(6,k,1)*d(6,k,1) )
            errmax = max( errmax, tmp1*xscal(k), tmp2*vscal(k) )
          end do
!
! If error is smaller than TOL, update position and velocity arrays and exit
          if (errmax.le.tol2) then
            do k = 2, nbod
              x0(1,k) = d(1,k,1)
              x0(2,k) = d(2,k,1)
              x0(3,k) = d(3,k,1)
              v0(1,k) = d(4,k,1)
              v0(2,k) = d(5,k,1)
              v0(3,k) = d(6,k,1)
            end do
!
            do j = 2, n
              do k = 2, nbod
                x0(1,k) = x0(1,k) + d(1,k,j)
                x0(2,k) = x0(2,k) + d(2,k,j)
                x0(3,k) = x0(3,k) + d(3,k,j)
                v0(1,k) = v0(1,k) + d(4,k,j)
                v0(2,k) = v0(2,k) + d(5,k,j)
                v0(3,k) = v0(3,k) + d(6,k,j)
              end do
            end do
!
! Save the actual stepsize used
            hdid = h0
!
! Recommend a new stepsize for the next call to this subroutine
            if (n.ge.8) h0 = h0 * SHRINK
            if (n.lt.7) h0 = h0 * GROW
            return
          end if
        end if
!
      end do
!
! If errors were too large, redo the step with half the previous step size.
      h0 = h0 * .5d0
      goto 100
!
!------------------------------------------------------------------------------
!
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MDT_HY.FOR    (ErikSoft   2 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Integrates NBOD bodies (of which NBIG are Big) for one timestep H
! using a second-order hybrid-symplectic integrator algorithm
!
! DTFLAG = 0 implies first ever call to this subroutine, 
!        = 1 implies first call since number/masses of objects changed.
!        = 2 normal call
!
! N.B. Input/output must be in democratic heliocentric coordinates.
! ===
!
!------------------------------------------------------------------------------
!
      subroutine mdt_hy (time,tstart,h0,tol,rmax,en,am,jcen,rcen,nbod, &
       nbig,m,x,v,s,rphys,rcrit,rce,stat,id,ngf,algor,opt,dtflag, &
       ngflag,opflag,colflag,nclo,iclo,jclo,dclo,tclo,ixvclo,jxvclo, &
       outfile,mem,lmem)
!
      implicit none
      include 'mercury95.inc'
!
! Input/Output
      integer nbod,nbig,stat(nbod),algor,opt(8),dtflag,ngflag,opflag
      integer colflag,lmem(NMESS),nclo,iclo(CMAX),jclo(CMAX)
      real(8) time,tstart,h0,tol,rmax,en(3),am(3),jcen(3),rcen
      real(8) m(nbod),x(3,nbod),v(3,nbod),s(3,nbod),rphys(nbod)
      real(8) rce(nbod),rcrit(nbod),ngf(4,nbod),tclo(CMAX),dclo(CMAX)
      real(8) ixvclo(6,CMAX),jxvclo(6,CMAX)
      character(80) outfile(3),mem(NMESS)
      character(8) id(nbod)
!
! Local
      integer j,nce,ice(NMAX),jce(NMAX),ce(NMAX),iflag
      real(8) a(3,NMAX),hby2,hrec,x0(3,NMAX),v0(3,NMAX),mvsum(3),temp
      real(8) angf(3,NMAX),ausr(3,NMAX)
      external mfo_hkce
!
!------------------------------------------------------------------------------
!
      save a, hrec, angf, ausr
      hby2 = h0 * .5d0
      nclo = 0
      colflag = 0
!
! If accelerations from previous call are not valid, calculate them now
      if (dtflag.ne.2) then
        if (dtflag.eq.0) hrec = h0
        call mfo_hy (jcen,nbod,nbig,m,x,rcrit,a,stat)
        dtflag = 2
        do j = 2, nbod
          angf(1,j) = 0.d0
          angf(2,j) = 0.d0
          angf(3,j) = 0.d0
          ausr(1,j) = 0.d0
          ausr(2,j) = 0.d0
          ausr(3,j) = 0.d0
        end do
! If required, apply non-gravitational and user-defined forces
        if (opt(8).eq.1) call mfo_user (time,jcen,nbod,nbig,m,x,v,ausr)
        if (ngflag.eq.1.or.ngflag.eq.3) call mfo_ngf (nbod,x,v,angf,ngf)
      end if
!
! Advance interaction Hamiltonian for H/2
      do j = 2, nbod
        v(1,j) = v(1,j)  +  hby2 * (angf(1,j) + ausr(1,j) + a(1,j))
        v(2,j) = v(2,j)  +  hby2 * (angf(2,j) + ausr(2,j) + a(2,j))
        v(3,j) = v(3,j)  +  hby2 * (angf(3,j) + ausr(3,j) + a(3,j))
      end do
!
! Advance solar Hamiltonian for H/2
      mvsum(1) = 0.d0
      mvsum(2) = 0.d0
      mvsum(3) = 0.d0
      do j = 2, nbod
        mvsum(1) = mvsum(1)  +  m(j) * v(1,j)
        mvsum(2) = mvsum(2)  +  m(j) * v(2,j)
        mvsum(3) = mvsum(3)  +  m(j) * v(3,j)
      end do
!
      temp = hby2 / m(1)
      mvsum(1) = temp * mvsum(1)
      mvsum(2) = temp * mvsum(2)
      mvsum(3) = temp * mvsum(3)
      do j = 2, nbod
        x(1,j) = x(1,j)  +  mvsum(1)
        x(2,j) = x(2,j)  +  mvsum(2)
        x(3,j) = x(3,j)  +  mvsum(3)
      end do
!
! Save the current coordinates and velocities
      call mco_iden (jcen,nbod,nbig,h0,m,x,v,x0,v0)
!
! Advance H_K for H
      do j = 2, nbod
        iflag = 0
        call drift_one (m(1),x(1,j),x(2,j),x(3,j),v(1,j),v(2,j), &
         v(3,j),h0,iflag)
      end do
!
! Check whether any object separations were < R_CRIT whilst advancing H_K
      call mce_snif (h0,2,nbod,nbig,x0,v0,x,v,rcrit,ce,nce,ice,jce)
!
! If objects had close encounters, advance H_K using Bulirsch-Stoer instead
      if (nce.gt.0) then
        do j = 2, nbod
          if (ce(j).ne.0) then
            x(1,j) = x0(1,j)
            x(2,j) = x0(2,j)
            x(3,j) = x0(3,j)
            v(1,j) = v0(1,j)
            v(2,j) = v0(2,j)
            v(3,j) = v0(3,j)
          end if
        end do
        call mdt_hkce (time,tstart,h0,hrec,tol,rmax,en(3),jcen,rcen, &
         nbod,nbig,m,x,v,s,rphys,rcrit,rce,stat,id,ngf,algor,opt, &
         ngflag,colflag,ce,nce,ice,jce,nclo,iclo,jclo,dclo,tclo,ixvclo, &
         jxvclo,outfile,mem,lmem,mfo_hkce)
      end if
!
! Advance solar Hamiltonian for H/2
      mvsum(1) = 0.d0
      mvsum(2) = 0.d0
      mvsum(3) = 0.d0
      do j = 2, nbod
        mvsum(1) = mvsum(1)  +  m(j) * v(1,j)
        mvsum(2) = mvsum(2)  +  m(j) * v(2,j)
        mvsum(3) = mvsum(3)  +  m(j) * v(3,j)
      end do
!
      temp = hby2 / m(1)
      mvsum(1) = temp * mvsum(1)
      mvsum(2) = temp * mvsum(2)
      mvsum(3) = temp * mvsum(3)
      do j = 2, nbod
        x(1,j) = x(1,j)  +  mvsum(1)
        x(2,j) = x(2,j)  +  mvsum(2)
        x(3,j) = x(3,j)  +  mvsum(3)
      end do
!
! Advance interaction Hamiltonian for H/2
      call mfo_hy (jcen,nbod,nbig,m,x,rcrit,a,stat)
      if (opt(8).eq.1) call mfo_user (time,jcen,nbod,nbig,m,x,v,ausr)
      if (ngflag.eq.1.or.ngflag.eq.3) call mfo_ngf (nbod,x,v,angf,ngf)
!
      do j = 2, nbod
        v(1,j) = v(1,j)  +  hby2 * (angf(1,j) + ausr(1,j) + a(1,j))
        v(2,j) = v(2,j)  +  hby2 * (angf(2,j) + ausr(2,j) + a(2,j))
        v(3,j) = v(3,j)  +  hby2 * (angf(3,j) + ausr(3,j) + a(3,j))
      end do
!
!------------------------------------------------------------------------------
!
      return
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MDT_HKCE.FOR    (ErikSoft   1 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Integrates NBOD bodies (of which NBIG are Big) for one timestep H under
! the Hamiltonian H_K, including close-encounter terms.
!
!------------------------------------------------------------------------------
!
      subroutine mdt_hkce (time,tstart,h0,hrec,tol,rmax,elost,jcen, &
       rcen,nbod,nbig,m,x,v,s,rphy,rcrit,rce,stat,id,ngf,algor,opt, &
       ngflag,colflag,ce,nce,ice,jce,nclo,iclo,jclo,dclo,tclo,ixvclo, &
       jxvclo,outfile,mem,lmem,force)
!
      implicit none
      include 'mercury95.inc'
!
! Input/Output
      integer nbod,nbig,nce,ice(nce),jce(nce),stat(nbod),ngflag,ce(nbod)
      integer algor,opt(8),colflag,lmem(NMESS),nclo,iclo(CMAX)
      integer jclo(CMAX)
      real(8) time,tstart,h0,hrec,tol,rmax,elost,jcen(3),rcen
      real(8) m(nbod),x(3,nbod),v(3,nbod),s(3,nbod)
      real(8) rce(nbod),rphy(nbod),rcrit(nbod),ngf(4,nbod)
      real(8) tclo(CMAX),dclo(CMAX),ixvclo(6,CMAX),jxvclo(6,CMAX)
      character(80) outfile(3),mem(NMESS)
      character(8) id(nbod)
      external force
!
! Local
      integer iback(NMAX),index(NMAX),ibs(NMAX),jbs(NMAX),nclo_old
      integer i,j,k,nbs,nbsbig,statbs(NMAX)
      integer nhit,ihit(CMAX),jhit(CMAX),chit(CMAX),nowflag,dtflag
      real(8) tlocal,hlocal,hdid,tmp0
      real(8) mbs(NMAX),xbs(3,NMAX),vbs(3,NMAX),sbs(3,NMAX)
      real(8) rcritbs(NMAX),rcebs(NMAX),rphybs(NMAX)
      real(8) ngfbs(4,NMAX),x0(3,NMAX),v0(3,NMAX)
      real(8) thit(CMAX),dhit(CMAX),thit1,temp
      character(8) idbs(NMAX)
!
!------------------------------------------------------------------------------
!
! N.B. Don't set nclo to zero!!
      nbs = 1
      nbsbig = 0
      mbs(1) = m(1)
      if (algor.eq.11) mbs(1) = m(1) + m(2)
      sbs(1,1) = s(1,1)
      sbs(2,1) = s(2,1)
      sbs(3,1) = s(3,1)
!
! Put data for close-encounter bodies into local arrays for use with BS routine
      do j = 2, nbod
        if (ce(j).ne.0) then
          nbs = nbs + 1
          if (j.le.nbig) nbsbig = nbs
          mbs(nbs)   = m(j)
          xbs(1,nbs) = x(1,j)
          xbs(2,nbs) = x(2,j)
          xbs(3,nbs) = x(3,j)
          vbs(1,nbs) = v(1,j)
          vbs(2,nbs) = v(2,j)
          vbs(3,nbs) = v(3,j)
          sbs(1,nbs) = s(1,j)
          sbs(2,nbs) = s(2,j)
          sbs(3,nbs) = s(3,j)
          rcebs(nbs) = rce(j)
          rphybs(nbs) = rphy(j)
          statbs(nbs) = stat(j)
          rcritbs(nbs) = rcrit(j)
          idbs(nbs) = id(j)
          index(nbs) = j
          iback(j) = nbs
        end if
      end do
!
      do k = 1, nce
        ibs(k) = iback(ice(k))
        jbs(k) = iback(jce(k))
      end do
!
      tlocal = 0.d0
      hlocal = sign(hrec,h0)
!
! Begin the Bulirsch-Stoer integration
  50  continue
        tmp0 = abs(h0) - abs(tlocal)
        hrec = hlocal
        if (abs(hlocal).gt.tmp0) hlocal = sign (tmp0, h0)
!
! Save old coordinates and integrate
        call mco_iden (jcen,nbs,0,h0,mbs,xbs,vbs,x0,v0)
        call mdt_bs2 (time,hlocal,hdid,tol,jcen,nbs,nbsbig,mbs,xbs,vbs, &
         sbs,rphybs,rcritbs,ngfbs,statbs,dtflag,ngflag,opt,nce, &
         ibs,jbs,force)
        tlocal = tlocal + hdid
!
! Check for close-encounter minima
        nclo_old = nclo
        temp = time + tlocal
        call mce_stat (temp,hdid,rcen,nbs,nbsbig,mbs,x0,v0,xbs,vbs, &
         rcebs,rphybs,nclo,iclo,jclo,dclo,tclo,ixvclo,jxvclo,nhit,ihit, &
         jhit,chit,dhit,thit,thit1,nowflag,statbs,outfile(3),mem,lmem)
!
! If collisions occurred, resolve the collision and return a flag
        if (nhit.gt.0.and.opt(2).ne.0) then
          do k = 1, nhit
            if (chit(k).eq.1) then
              i = ihit(k)
              j = jhit(k)
              call mce_coll (thit(k),tstart,elost,jcen,i,j,nbs,nbsbig, &
               mbs,xbs,vbs,sbs,rphybs,statbs,idbs,opt,mem,lmem, &
               outfile(3))
              colflag = colflag + 1
            end if
          end do
        end if
!
! If necessary, continue integrating objects undergoing close encounters
      if ((tlocal - h0)*h0.lt.0) goto 50
!
! Return data for the close-encounter objects to global arrays
      do k = 2, nbs
        j = index(k)
        m(j)   = mbs(k)
        x(1,j) = xbs(1,k)
        x(2,j) = xbs(2,k)
        x(3,j) = xbs(3,k)
        v(1,j) = vbs(1,k)
        v(2,j) = vbs(2,k)
        v(3,j) = vbs(3,k)
        s(1,j) = sbs(1,k)
        s(2,j) = sbs(2,k)
        s(3,j) = sbs(3,k)
        stat(j) = statbs(k)
      end do
      do k = 1, nclo
        iclo(k) = index(iclo(k))
        jclo(k) = index(jclo(k))
      end do
!
!------------------------------------------------------------------------------
!
      return
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MDT_MVS.FOR    (ErikSoft   28 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Integrates NBOD bodies (of which NBIG are Big) for one timestep H
! using a second-order mixed-variable symplectic integrator.
!
! DTFLAG = 0 implies first ever call to this subroutine, 
!        = 1 implies first call since number/masses of objects changed.
!        = 2 normal call
!
! N.B. Input/output must be in coordinates with respect to the central body.
! ===
!
!------------------------------------------------------------------------------
!
      subroutine mdt_mvs (time,tstart,h0,tol,rmax,en,am,jcen,rcen,nbod, &
       nbig,m,x,v,s,rphys,rcrit,rce,stat,id,ngf,algor,opt,dtflag, &
       ngflag,opflag,colflag,nclo,iclo,jclo,dclo,tclo,ixvclo,jxvclo, &
       outfile,mem,lmem)
!
      implicit none
      include 'mercury95.inc'
!
! Input/Output
      integer nbod,nbig,stat(nbod),algor,opt(8),dtflag,ngflag,opflag
      integer colflag,lmem(NMESS),nclo,iclo(CMAX),jclo(CMAX)
      real(8) time,tstart,h0,tol,rmax,en(3),am(3),jcen(3),rcen
      real(8) m(nbod),x(3,nbod),v(3,nbod),s(3,nbod),rphys(nbod)
      real(8) rce(nbod),rcrit(nbod),ngf(4,nbod),tclo(CMAX),dclo(CMAX)
      real(8) ixvclo(6,CMAX),jxvclo(6,CMAX)
      character(80) outfile(3),mem(NMESS)
      character(8) id(nbod)
!
! Local
      integer j,iflag,nhit,ihit(CMAX),jhit(CMAX),chit(CMAX),nowflag
      real(8) xj(3,NMAX),vj(3,NMAX),a(3,NMAX),gm(NMAX),hby2,thit1,temp
      real(8) msofar,minside,x0(3,NMAX),v0(3,NMAX),dhit(CMAX),thit(CMAX)
      real(8) angf(3,NMAX),ausr(3,NMAX)
!
!------------------------------------------------------------------------------
!
      save a, xj, gm, angf, ausr
      hby2 = .5d0 * h0
      nclo = 0
!
! If accelerations from previous call are not valid, calculate them now,
! and also the Jacobi coordinates XJ, and effective central masses GM.
      if (dtflag.ne.2) then
        dtflag = 2
        call mco_h2j (jcen,nbig,nbig,h0,m,x,v,xj,vj)
        call mfo_mvs (jcen,nbod,nbig,m,x,xj,a,stat)
!
        minside = m(1)
        do j = 2, nbig
          msofar = minside + m(j)
          gm(j) = m(1) * msofar / minside
          minside = msofar
          angf(1,j) = 0.d0
          angf(2,j) = 0.d0
          angf(3,j) = 0.d0
          ausr(1,j) = 0.d0
          ausr(2,j) = 0.d0
          ausr(3,j) = 0.d0
        end do
! If required, apply non-gravitational and user-defined forces
        if (opt(8).eq.1) call mfo_user (time,jcen,nbod,nbig,m,x,v,ausr)
        if (ngflag.eq.1.or.ngflag.eq.3) call mfo_ngf (nbod,x,v,angf,ngf)
      end if
!
! Advance interaction Hamiltonian for H/2
      do j = 2, nbod
        v(1,j) = v(1,j)  +  hby2 * (angf(1,j) + ausr(1,j) + a(1,j))
        v(2,j) = v(2,j)  +  hby2 * (angf(2,j) + ausr(2,j) + a(2,j))
        v(3,j) = v(3,j)  +  hby2 * (angf(3,j) + ausr(3,j) + a(3,j))
      end do
!
! Save current coordinates and velocities
      call mco_iden (jcen,nbod,nbig,h0,m,x,v,x0,v0)
!
! Advance Keplerian Hamiltonian (Jacobi/helio coords for Big/Small bodies)
      call mco_h2j (jcen,nbig,nbig,h0,m,x,v,xj,vj)
      do j = 2, nbig
        iflag = 0
        call drift_one (gm(j),xj(1,j),xj(2,j),xj(3,j),vj(1,j), &
         vj(2,j),vj(3,j),h0,iflag)
      end do
      do j = nbig + 1, nbod
        iflag = 0
        call drift_one (m(1),x(1,j),x(2,j),x(3,j),v(1,j),v(2,j), &
         v(3,j),h0,iflag)
      end do
      call mco_j2h (time,jcen,nbig,nbig,h0,m,xj,vj,x,v,ngf,ngflag,opt)
!
! Check for close-encounter minima during drift step
      temp = time + h0
      call mce_stat (temp,h0,rcen,nbod,nbig,m,x0,v0,x,v,rce,rphys,nclo, &
       iclo,jclo,dclo,tclo,ixvclo,jxvclo,nhit,ihit,jhit,chit,dhit,thit, &
       thit1,nowflag,stat,outfile(3),mem,lmem)
!
! Advance interaction Hamiltonian for H/2
      call mfo_mvs (jcen,nbod,nbig,m,x,xj,a,stat)
      if (opt(8).eq.1) call mfo_user (time,jcen,nbod,nbig,m,x,v,ausr)
      if (ngflag.eq.1.or.ngflag.eq.3) call mfo_ngf (nbod,x,v,angf,ngf)
!
      do j = 2, nbod
        v(1,j) = v(1,j)  +  hby2 * (angf(1,j) + ausr(1,j) + a(1,j))
        v(2,j) = v(2,j)  +  hby2 * (angf(2,j) + ausr(2,j) + a(2,j))
        v(3,j) = v(3,j)  +  hby2 * (angf(3,j) + ausr(3,j) + a(3,j))
      end do
!
!------------------------------------------------------------------------------
!
      return
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MDT_RA15.FOR    (ErikSoft   2 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Integrates NBOD bodies (of which NBIG are Big) for one timestep H0 using
! Everhart's RA15 integrator algorithm. The accelerations are calculated
! using the subroutine FORCE. The accuracy of the step is approximately 
! determined by the tolerance parameter TOL.
!
! Based on RADAU by E. Everhart, Physics Department, University of Denver.
! Comments giving equation numbers refer to Everhart (1985) ``An Efficient
! Integrator that Uses Gauss-Radau Spacings'', in The Dynamics of Comets:
! Their Origin and Evolution, p185-202, eds. A. Carusi & G. B. Valsecchi,
! pub Reidel. (A listing of the original subroutine is also given in this 
! paper.)
!
! DTFLAG = 0 implies first ever call to this subroutine, 
!        = 1 implies first call since number/masses of objects changed.
!        = 2 normal call
!
! N.B. Input/output must be in coordinates with respect to the central body.
! ===
!
!------------------------------------------------------------------------------
!
      subroutine mdt_ra15 (time,t,tdid,tol,jcen,nbod,nbig,mass,x1,v1, &
       spin,rphys,rcrit,ngf,stat,dtflag,ngflag,opt,nce,ice,jce,force)
!
      implicit none
      include 'mercury95.inc'
!
! Input/Output
      integer nbod,nbig,dtflag,ngflag,opt(8),stat(nbod)
      integer nce,ice(nce),jce(nce)
      real(8) time,t,tdid,tol,jcen(3),mass(nbod)
      real(8) x1(3*nbod),v1(3*nbod),spin(3*nbod)
      real(8) ngf(4,nbod),rphys(nbod),rcrit(nbod)
      external force
!
! Local
      integer nv,niter,j,k,n
      real(8) x(3*NMAX),v(3*NMAX),a(3*NMAX),a1(3*NMAX)
      real(8) g(7,3*NMAX),b(7,3*NMAX),e(7,3*NMAX)
      real(8) h(8),xc(8),vc(7),c(21),d(21),r(28),s(9)
      real(8) q,q2,q3,q4,q5,q6,q7,temp,gk
!
!------------------------------------------------------------------------------
!
      save h,xc,vc,c,d,r,b,e
!
! Gauss-Radau spacings for substeps within a sequence, for the 15th order 
! integrator. The sum of the H values should be 3.733333333333333
!
      data h/          0.d0,.0562625605369221d0,.1802406917368924d0, &
       .3526247171131696d0,.5471536263305554d0,.7342101772154105d0, &
       .8853209468390958d0,.9775206135612875d0/
!
! Constant coefficients used in series expansions for X and V
!  XC: 1/2,  1/6,  1/12, 1/20, 1/30, 1/42, 1/56, 1/72
!  VC: 1/2,  1/3,  1/4,  1/5,  1/6,  1/7,  1/8
      data xc/.5d0,.1666666666666667d0,.08333333333333333d0,.05d0, &
       .03333333333333333d0,.02380952380952381d0,.01785714285714286d0, &
       .01388888888888889d0/
      data vc/.5d0,.3333333333333333d0,.25d0,.2d0, &
       .1666666666666667d0,.1428571428571429d0,.125d0/
!
! If this is first call to the subroutine, set values of the constant arrays
! (R = R21, R31, R32, R41, R42, R43 in Everhart's paper.)
      if (dtflag.eq.0) then
        n = 0
        do j = 2, 8
          do k = 1, j - 1
            n = n + 1
            r(n) = 1.d0 / (h(j) - h(k))
          end do
        end do
!
! Constants to convert between B and G arrays (C = C21, C31, C32, C41, C42...)
        c(1) = - h(2)
        d(1) =   h(2)
        n = 1
        do j = 3, 7
          n = n + 1
          c(n) = -h(j) * c(n-j+2)
          d(n) =  h(2) * d(n-j+2)
          do k = 3, j - 1
            n = n + 1
            c(n) = c(n-j+1)  -  h(j) * c(n-j+2)
            d(n) = d(n-j+1)  +  h(k) * d(n-j+2)
          end do
          n = n + 1
          c(n) = c(n-j+1) - h(j)
          d(n) = d(n-j+1) + h(j)
        end do
!
        dtflag = 1
      end if
!
      nv = 3 * nbod
  100 continue
!
! If this is first call to subroutine since number/masses of objects changed
! do 6 iterations and initialize B, E arrays, otherwise do 2 iterations.
      if (dtflag.eq.1) then
        niter = 6
        do j = 4, nv
          do k = 1, 7
            b (k,j) = 0.d0
            e (k,j) = 0.d0
          end do
        end do
      else
        niter = 2
      end if
!
! Calculate forces at the start of the sequence
      call force (time,jcen,nbod,nbig,mass,x1,v1,spin,rcrit,a1,stat,ngf, &
       ngflag,opt,nce,ice,jce)
!
! Find G values from B values predicted at the last call (Eqs. 7 of Everhart)
      do k = 4, nv
        g(1,k) = b(7,k)*d(16) + b(6,k)*d(11) + b(5,k)*d(7) &
              + b(4,k)*d(4)  + b(3,k)*d(2)  + b(2,k)*d(1)  + b(1,k)
        g(2,k) = b(7,k)*d(17) + b(6,k)*d(12) + b(5,k)*d(8) &
              + b(4,k)*d(5)  + b(3,k)*d(3)  + b(2,k)
        g(3,k) = b(7,k)*d(18) + b(6,k)*d(13) + b(5,k)*d(9) &
              + b(4,k)*d(6)  + b(3,k)
        g(4,k) = b(7,k)*d(19) + b(6,k)*d(14) + b(5,k)*d(10) + b(4,k)
        g(5,k) = b(7,k)*d(20) + b(6,k)*d(15) + b(5,k)
        g(6,k) = b(7,k)*d(21) + b(6,k)
        g(7,k) = b(7,k)
      end do
!
!------------------------------------------------------------------------------
!
!  MAIN  LOOP  STARTS  HERE
!
! For each iteration (six for first call to subroutine, two otherwise)...
      do n = 1, niter
!
! For each substep within a sequence...
        do j = 2, 8
!
! Calculate position predictors using Eqn. 9 of Everhart
          s(1) = t * h(j)
          s(2) = s(1) * s(1) * .5d0
          s(3) = s(2) * h(j) * .3333333333333333d0
          s(4) = s(3) * h(j) * .5d0
          s(5) = s(4) * h(j) * .6d0
          s(6) = s(5) * h(j) * .6666666666666667d0
          s(7) = s(6) * h(j) * .7142857142857143d0
          s(8) = s(7) * h(j) * .75d0
          s(9) = s(8) * h(j) * .7777777777777778d0
!
          do k = 4, nv
            x(k) = s(9)*b(7,k) + s(8)*b(6,k) + s(7)*b(5,k) &
                + s(6)*b(4,k) + s(5)*b(3,k) + s(4)*b(2,k) &
                + s(3)*b(1,k) + s(2)*a1(k)  + s(1)*v1(k) + x1(k)
          end do
!
! If necessary, calculate velocity predictors too, from Eqn. 10 of Everhart
          if (ngflag.ne.0) then
            s(1) = t * h(j)
            s(2) = s(1) * h(j) * .5d0
            s(3) = s(2) * h(j) * .6666666666666667d0
            s(4) = s(3) * h(j) * .75d0
            s(5) = s(4) * h(j) * .8d0
            s(6) = s(5) * h(j) * .8333333333333333d0
            s(7) = s(6) * h(j) * .8571428571428571d0
            s(8) = s(7) * h(j) * .875d0
!
            do k = 4, nv
              v(k) = s(8)*b(7,k) + s(7)*b(6,k) + s(6)*b(5,k) &
                  + s(5)*b(4,k) + s(4)*b(3,k) + s(3)*b(2,k) &
                  + s(2)*b(1,k) + s(1)*a1(k)  + v1(k)
            end do
          end if
!
! Calculate forces at the current substep
          call force (time,jcen,nbod,nbig,mass,x,v,spin,rcrit,a,stat, &
           ngf,ngflag,opt,nce,ice,jce)
!
! Update G values using Eqs. 4 of Everhart, and update B values using Eqs. 5
          if (j.eq.2) then
            do k = 4, nv
              temp = g(1,k)
              g(1,k) = (a(k) - a1(k)) * r(1)
              b(1,k) = b(1,k) + g(1,k) - temp
            end do
            goto 300
          end if
          if (j.eq.3) then
            do k = 4, nv
              temp = g(2,k)
              gk = a(k) - a1(k)
              g(2,k) = (gk*r(2) - g(1,k))*r(3)
              temp = g(2,k) - temp
              b(1,k) = b(1,k)  +  temp * c(1)
              b(2,k) = b(2,k)  +  temp
            end do
            goto 300
          end if
          if (j.eq.4) then
            do k = 4, nv
              temp = g(3,k)
              gk = a(k) - a1(k)
              g(3,k) = ((gk*r(4) - g(1,k))*r(5) - g(2,k))*r(6)
              temp = g(3,k) - temp
              b(1,k) = b(1,k)  +  temp * c(2)
              b(2,k) = b(2,k)  +  temp * c(3)
              b(3,k) = b(3,k)  +  temp
            end do
            goto 300
          end if
          if (j.eq.5) then
            do k = 4, nv
              temp = g(4,k)
              gk = a(k) - a1(k)
              g(4,k) = (((gk*r(7) - g(1,k))*r(8) - g(2,k))*r(9) &
                    - g(3,k))*r(10)
              temp = g(4,k) - temp
              b(1,k) = b(1,k)  +  temp * c(4)
              b(2,k) = b(2,k)  +  temp * c(5)
              b(3,k) = b(3,k)  +  temp * c(6)
              b(4,k) = b(4,k)  +  temp
            end do
            goto 300
          end if
          if (j.eq.6) then
            do k = 4, nv
              temp = g(5,k)
              gk = a(k) - a1(k)
              g(5,k) =  ((((gk*r(11) - g(1,k))*r(12) - g(2,k))*r(13) &
                    - g(3,k))*r(14) - g(4,k))*r(15)
              temp = g(5,k) - temp
              b(1,k) = b(1,k)  +  temp * c(7)
              b(2,k) = b(2,k)  +  temp * c(8)
              b(3,k) = b(3,k)  +  temp * c(9)
              b(4,k) = b(4,k)  +  temp * c(10)
              b(5,k) = b(5,k)  +  temp
            end do
            goto 300
          end if
          if (j.eq.7) then
            do k = 4, nv
              temp = g(6,k)
              gk = a(k) - a1(k)
              g(6,k) = (((((gk*r(16) - g(1,k))*r(17) - g(2,k))*r(18) &
                    - g(3,k))*r(19) - g(4,k))*r(20) - g(5,k))*r(21)
              temp = g(6,k) - temp
              b(1,k) = b(1,k)  +  temp * c(11)
              b(2,k) = b(2,k)  +  temp * c(12)
              b(3,k) = b(3,k)  +  temp * c(13)
              b(4,k) = b(4,k)  +  temp * c(14)
              b(5,k) = b(5,k)  +  temp * c(15)
              b(6,k) = b(6,k)  +  temp
            end do
            goto 300
          end if
          if (j.eq.8) then
            do k = 4, nv
              temp = g(7,k)
              gk = a(k) - a1(k)
              g(7,k) = ((((((gk*r(22) - g(1,k))*r(23) - g(2,k))*r(24) &
                     - g(3,k))*r(25) - g(4,k))*r(26) - g(5,k))*r(27) &
                     - g(6,k))*r(28)
              temp = g(7,k) - temp
              b(1,k) = b(1,k)  +  temp * c(16)
              b(2,k) = b(2,k)  +  temp * c(17)
              b(3,k) = b(3,k)  +  temp * c(18)
              b(4,k) = b(4,k)  +  temp * c(19)
              b(5,k) = b(5,k)  +  temp * c(20)
              b(6,k) = b(6,k)  +  temp * c(21)
              b(7,k) = b(7,k)  +  temp
            end do
          end if
 300      continue
        end do
      end do
!
!------------------------------------------------------------------------------
!
!  END  OF  MAIN  LOOP
!
! Estimate suitable sequence size for the next call to subroutine (Eqs. 15, 16)
      temp = 0.d0
      do k = 4, nv
        temp = max( temp, abs( b(7,k) ) )
      end do
      temp = temp / (72.d0 * abs(t)**7)
      tdid = t
      if (temp.eq.0) then
        t = tdid * 1.4d0
      else
        t = sign( (tol/temp)**(1.d0/9.d0), tdid )
      end if
!
! If sequence size for the first subroutine call is too big, go back and redo
! the sequence using a smaller size.
      if (dtflag.eq.1.and.abs(t/tdid).lt.1) then
        t = t * .8d0
        goto 100
      end if
!
! If new sequence size is much bigger than the current one, reduce it
      if (abs(t/tdid).gt.1.4d0) t = tdid * 1.4d0
!
! Find new position and velocity values at end of the sequence (Eqs. 11, 12)
      temp = tdid * tdid
      do k = 4 , nv
        x1(k) = (xc(8)*b(7,k) + xc(7)*b(6,k) + xc(6)*b(5,k) &
             +  xc(5)*b(4,k) + xc(4)*b(3,k) + xc(3)*b(2,k) &
             +  xc(2)*b(1,k) + xc(1)*a1(k))*temp + v1(k)*tdid + x1(k)
!
        v1(k) = (vc(7)*b(7,k) + vc(6)*b(6,k) + vc(5)*b(5,k) &
             +  vc(4)*b(4,k) + vc(3)*b(3,k) + vc(2)*b(2,k) &
             +  vc(1)*b(1,k) + a1(k))*tdid + v1(k)
      end do
!
! Predict new B values to use at the start of the next sequence. The predicted
! values from the last call are saved as E. The correction, BD, between the
! actual and predicted values of B is applied in advance as a correction.
      q = t / tdid
      q2 = q  * q
      q3 = q  * q2
      q4 = q2 * q2
      q5 = q2 * q3
      q6 = q3 * q3
      q7 = q3 * q4
!
      do k = 4, nv
        s(1) = b(1,k) - e(1,k)
        s(2) = b(2,k) - e(2,k)
        s(3) = b(3,k) - e(3,k)
        s(4) = b(4,k) - e(4,k)
        s(5) = b(5,k) - e(5,k)
        s(6) = b(6,k) - e(6,k)
        s(7) = b(7,k) - e(7,k)
!
! Estimate B values for the next sequence (Eqs. 13 of Everhart).
        e(1,k) = q* (b(7,k)* 7.d0 + b(6,k)* 6.d0 + b(5,k)* 5.d0 &
              +     b(4,k)* 4.d0 + b(3,k)* 3.d0 + b(2,k)*2.d0 + b(1,k))
        e(2,k) = q2*(b(7,k)*21.d0 + b(6,k)*15.d0 + b(5,k)*10.d0 &
              +     b(4,k)* 6.d0 + b(3,k)* 3.d0 + b(2,k))
        e(3,k) = q3*(b(7,k)*35.d0 + b(6,k)*20.d0 + b(5,k)*10.d0 &
              +     b(4,k)*4.d0  + b(3,k))
        e(4,k) = q4*(b(7,k)*35.d0 + b(6,k)*15.d0 + b(5,k)*5.d0 + b(4,k))
        e(5,k) = q5*(b(7,k)*21.d0 + b(6,k)*6.d0  + b(5,k))
        e(6,k) = q6*(b(7,k)*7.d0  + b(6,k))
        e(7,k) = q7* b(7,k)
!
        b(1,k) = e(1,k) + s(1)
        b(2,k) = e(2,k) + s(2)
        b(3,k) = e(3,k) + s(3)
        b(4,k) = e(4,k) + s(4)
        b(5,k) = e(5,k) + s(5)
        b(6,k) = e(6,k) + s(6)
        b(7,k) = e(7,k) + s(7)
      end do
      dtflag = 2
!
!------------------------------------------------------------------------------
!
      return
      end
!

