!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    MCE_BOX.FOR    (ErikSoft   30 September 2000)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Given initial and final coordinates and velocities, the routine returns
! the X and Y coordinates of a box bounding the motion in between the
! end points.
!
! If the X or Y velocity changes sign, the routine performs a quadratic
! interpolation to estimate the corresponding extreme value of X or Y.
!
!------------------------------------------------------------------------------
!
      subroutine mce_box (nbod,h,x0,v0,x1,v1,xmin,xmax,ymin,ymax)
!
      implicit none
      include 'mercury95.inc'
!
! Input/Output
      integer nbod
      real(8) h,x0(3,nbod), x1(3,nbod), v0(3,nbod),v1(3,nbod)
      real(8)   xmin(nbod), xmax(nbod), ymin(nbod),ymax(nbod)
!
! Local
      integer j
      real(8) temp
!
!------------------------------------------------------------------------------
!
      do j = 2, nbod
        xmin(j) = min (x0(1,j), x1(1,j))
        xmax(j) = max (x0(1,j), x1(1,j))
        ymin(j) = min (x0(2,j), x1(2,j))
        ymax(j) = max (x0(2,j), x1(2,j))
!
! If velocity changes sign, do an interpolation
        if ((v0(1,j).lt.0.and.v1(1,j).gt.0).or. &
           (v0(1,j).gt.0.and.v1(1,j).lt.0)) then
          temp = (v0(1,j)*x1(1,j) - v1(1,j)*x0(1,j)  &
                 - .5d0*h*v0(1,j)*v1(1,j)) / (v0(1,j) - v1(1,j))
          xmin(j) = min (xmin(j),temp)
          xmax(j) = max (xmax(j),temp)
        end if
!
        if ((v0(2,j).lt.0.and.v1(2,j).gt.0).or. &
           (v0(2,j).gt.0.and.v1(2,j).lt.0)) then
          temp = (v0(2,j)*x1(2,j) - v1(2,j)*x0(2,j)  &
                 - .5d0*h*v0(2,j)*v1(2,j)) / (v0(2,j) - v1(2,j))
          ymin(j) = min (ymin(j),temp)
          ymax(j) = max (ymax(j),temp)
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
!    MCE_CENT.FOR    (ErikSoft   4 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Checks all objects with index I >= I0, to see if they have had a collision
! with the central body in a time interval H, when given the initial and 
! final coordinates and velocities. The routine uses cubic interpolation
! to estimate the minimum separations.
!
! N.B. All coordinates & velocities must be with respect to the central body!!
! ===
!------------------------------------------------------------------------------
!
      subroutine mce_cent (time,h,rcen,jcen,i0,nbod,nbig,m,x0,v0,x1,v1, &
       nhit,jhit,thit,dhit,algor,ngf,ngflag)
!
      implicit none
      include 'mercury95.inc'
!
! Input/Output
      integer i0, nbod, nbig, nhit, jhit(CMAX), algor, ngflag
      real(8) time,h,rcen,jcen(3),m(nbod),x0(3,nbod),v0(3,nbod)
      real(8) x1(3,nbod),v1(3,nbod),thit(CMAX),dhit(CMAX),ngf(4,nbod)
!
! Local
      integer j
      real(8) rcen2,mco_acsh,a,q,u0,uhit,m0,mhit,mm,r0,mcen
      real(8) hx,hy,hz,h2,p,rr0,rr1,rv0,rv1,temp,e,v2
      real(8) xu0(3,NMAX),xu1(3,NMAX),vu0(3,NMAX),vu1(3,NMAX)
!
!------------------------------------------------------------------------------
!
      if (i0.le.0) i0 = 2
      nhit = 0
      rcen2 = rcen * rcen
      mcen = m(1)
!
! If using close-binary code, convert to coords with respect to the binary
!      if (algor.eq.11) then
!        mcen = m(1) + m(2)
!        call mco_h2ub (temp,jcen,nbod,nbig,h,m,x0,v0,xu0,vu0,ngf,ngflag)
!        call mco_h2ub (temp,jcen,nbod,nbig,h,m,x1,v1,xu1,vu1,ngf,ngflag)
!      end if
!
! Check for collisions with the central body
      do j = i0, nbod
        if (algor.eq.11) then
          rr0 = xu0(1,j)*xu0(1,j) + xu0(2,j)*xu0(2,j) +xu0(3,j)*xu0(3,j)
          rr1 = xu1(1,j)*xu1(1,j) + xu1(2,j)*xu1(2,j) +xu1(3,j)*xu1(3,j)
          rv0 = vu0(1,j)*xu0(1,j) + vu0(2,j)*xu0(2,j) +vu0(3,j)*xu0(3,j)
          rv1 = vu1(1,j)*xu1(1,j) + vu1(2,j)*xu1(2,j) +vu1(3,j)*xu1(3,j)
        else
          rr0 = x0(1,j)*x0(1,j) + x0(2,j)*x0(2,j) + x0(3,j)*x0(3,j)
          rr1 = x1(1,j)*x1(1,j) + x1(2,j)*x1(2,j) + x1(3,j)*x1(3,j)
          rv0 = v0(1,j)*x0(1,j) + v0(2,j)*x0(2,j) + v0(3,j)*x0(3,j)
          rv1 = v1(1,j)*x1(1,j) + v1(2,j)*x1(2,j) + v1(3,j)*x1(3,j)
        end if
!
! If inside the central body, or passing through pericentre, use 2-body approx.
        if ((rv0*h.le.0.and.rv1*h.ge.0).or.min(rr0,rr1).le.rcen2) then
          if (algor.eq.11) then
            hx = xu0(2,j) * vu0(3,j)  -  xu0(3,j) * vu0(2,j)
            hy = xu0(3,j) * vu0(1,j)  -  xu0(1,j) * vu0(3,j)
            hz = xu0(1,j) * vu0(2,j)  -  xu0(2,j) * vu0(1,j)
            v2 = vu0(1,j)*vu0(1,j) +vu0(2,j)*vu0(2,j) +vu0(3,j)*vu0(3,j)
          else
            hx = x0(2,j) * v0(3,j)  -  x0(3,j) * v0(2,j)
            hy = x0(3,j) * v0(1,j)  -  x0(1,j) * v0(3,j)
            hz = x0(1,j) * v0(2,j)  -  x0(2,j) * v0(1,j)
            v2 = v0(1,j)*v0(1,j) + v0(2,j)*v0(2,j) + v0(3,j)*v0(3,j)
          end if
          h2 = hx*hx + hy*hy + hz*hz
          p = h2 / (mcen + m(j))
          r0 = sqrt(rr0)
          temp = 1.d0 + p*(v2/(mcen + m(j)) - 2.d0/r0)
          e = sqrt( max(temp,0.d0) )
          q = p / (1.d0 + e)
!
! If the object hit the central body
! ##ajm## 28-05-12 edited this to check that collision
! occurs within the timestep
          if (q.le.rcen) then
!     Time of impact relative to the end of the timestep
            if (e.lt.1) then
              a = q / (1.d0 - e)
              uhit = sign (acos((1.d0 - rcen/a)/e), -h)
              u0   = sign (acos((1.d0 - r0/a  )/e), rv0)
              mhit = mod (uhit - e*sin(uhit) + PI, TWOPI) - PI
              m0   = mod (u0   - e*sin(u0)   + PI, TWOPI) - PI
            else
              a = q / (e - 1.d0)
              uhit = sign (mco_acsh((1.d0 - rcen/a)/e), -h)
              u0   = sign (mco_acsh((1.d0 - r0/a  )/e), rv0)
              mhit = mod (uhit - e*sinh(uhit) + PI, TWOPI) - PI
              m0   = mod (u0   - e*sinh(u0)   + PI, TWOPI) - PI
            end if
            mm = sqrt((mcen + m(j)) / (a*a*a))
            if ((mhit-m0)/mm.le.h) then
               nhit = nhit + 1
               jhit(nhit) = j
               dhit(nhit) = rcen
               thit(nhit) = (mhit - m0) / mm + time
            endif
         end if
! ##ajm##
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
!      MCE_COLL.FOR    (ErikSoft   2 October 2000)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Resolves a collision between two objects, using the collision model chosen
! by the user. Also writes a message to the information file, and updates the
! value of ELOST, the change in energy due to collisions and ejections.
!
! N.B. All coordinates and velocities must be with respect to central body.
! ===
!
!------------------------------------------------------------------------------
!
      subroutine mce_coll (time,tstart,elost,jcen,i,j,nbod,nbig,m,xh, &
       vh,s,rphys,stat,id,opt,mem,lmem,outfile)
!
      implicit none
      include 'mercury95.inc'
!
! Input/Output
      integer i,j,nbod,nbig,stat(nbod),opt(8),lmem(NMESS)
      real(8) time,tstart,elost,jcen(3)
      real(8) m(nbod),xh(3,nbod),vh(3,nbod),s(3,nbod),rphys(nbod)
      character(80) outfile,mem(NMESS)
      character(8) id(nbod)
!
! Local
      integer year,month,itmp
      real(8) t1
      character(38) flost,fcol
      character(6) tstring
!
!------------------------------------------------------------------------------
!
! If two bodies collided, check that the less massive one is removed
! (unless the more massive one is a Small body)
      if (i.ne.0) then
        if (m(j).gt.m(i).and.j.le.nbig) then
          itmp = i
          i = j
          j = itmp
        end if
      end if
!
! Write message to info file (I=0 implies collision with the central body)
  10  open (23, file=outfile, status='old', access='append', err=10)
!
      if (opt(3).eq.1) then
        call mio_jd2y (time,year,month,t1)
        if (i.eq.0) then
          flost = '(1x,a8,a,i10,1x,i2,1x,f8.5)'
          write (23,flost) id(j),mem(67)(1:lmem(67)),year,month,t1
        else
          fcol  = '(1x,a8,a,a8,a,i10,1x,i2,1x,f4.1)'
          write (23,fcol) id(i),mem(69)(1:lmem(69)),id(j), &
           mem(71)(1:lmem(71)),year,month,t1
        end if
      else
        if (opt(3).eq.3) then
          t1 = (time - tstart) / 365.25d0
          tstring = mem(2)
          flost = '(1x,a8,a,f18.7,a)'
          fcol  = '(1x,a8,a,a8,a,1x,f14.3,a)'
        else
          if (opt(3).eq.0) t1 = time
          if (opt(3).eq.2) t1 = time - tstart
          tstring = mem(1)
          flost = '(1x,a8,a,f18.5,a)'
          fcol  = '(1x,a8,a,a8,a,1x,f14.1,a)'
        end if
        if (i.eq.0.or.i.eq.1) then
          write (23,flost) id(j),mem(67)(1:lmem(67)),t1,tstring
        else
          write (23,fcol) id(i),mem(69)(1:lmem(69)),id(j), &
           mem(71)(1:lmem(71)),t1,tstring
        end if
      end if
      close (23)
!
! Do the collision (inelastic merger)
      call mce_merg (jcen,i,j,nbod,nbig,m,xh,vh,s,stat,elost)
!
!------------------------------------------------------------------------------
!
      return
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MCE_HILL.FOR    (ErikSoft   4 October 2000)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Calculates the Hill radii for all objects given their masses, M,
! coordinates, X, and velocities, V; plus the mass of the central body, M(1)
! Where HILL = a * (m/3*m(1))^(1/3)
!
! If the orbit is hyperbolic or parabolic, the Hill radius is calculated using:
!       HILL = r * (m/3*m(1))^(1/3)
! where R is the current distance from the central body.
!
! The routine also gives the semi-major axis, A, of each object's orbit.
!
! N.B. Designed to use heliocentric coordinates, but should be adequate using
! ===  barycentric coordinates.
!
!------------------------------------------------------------------------------
!
      subroutine mce_hill (nbod,m,x,v,hill,a)
!
      implicit none
      include 'mercury95.inc'
      real(8) THIRD
      parameter (THIRD = .3333333333333333d0)
!
! Input/Output
      integer nbod
      real(8) m(nbod),x(3,nbod),v(3,nbod),hill(nbod),a(nbod)
!
! Local
      integer j
      real(8) r, v2, gm
!
!------------------------------------------------------------------------------
!
      do j = 2, nbod
        gm = m(1) + m(j)
        call mco_x2a (gm,x(1,j),x(2,j),x(3,j),v(1,j),v(2,j),v(3,j),a(j), &
         r,v2)
! If orbit is hyperbolic, use the distance rather than the semi-major axis
        if (a(j).le.0) a(j) = r
        hill(j) = a(j) * (THIRD * m(j) / m(1))**THIRD
      end do
!
!------------------------------------------------------------------------------
!
      return
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MCE_INIT.FOR    (ErikSoft   28 February 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Calculates close-approach limits RCE (in AU) and physical radii RPHYS
! (in AU) for all objects, given their masses M, coordinates X, velocities
! V, densities RHO, and close-approach limits RCEH (in Hill radii).
!
! Also calculates the changeover distance RCRIT, used by the hybrid
! symplectic integrator. RCRIT is defined to be the larger of N1*HILL and
! N2*H*VMAX, where HILL is the Hill radius, H is the timestep, VMAX is the
! largest expected velocity of any body, and N1, N2 are parameters (see
! section 4.2 of Chambers 1999, Monthly Notices, vol 304, p793-799).
!
! N.B. Designed to use heliocentric coordinates, but should be adequate using
! ===  barycentric coordinates.
!
!------------------------------------------------------------------------------
!
      subroutine mce_init (tstart,algor,h,jcen,rcen,rmax,cefac,nbod, &
       nbig,m,x,v,s,rho,rceh,rphys,rce,rcrit,id,opt,outfile,rcritflag)
!
      implicit none
      include 'mercury95.inc'
!
      real(8) N2,THIRD
      parameter (N2=.4d0,THIRD=.3333333333333333d0)
!
! Input/Output
      integer nbod,nbig,algor,opt(8),rcritflag
      real(8) tstart,h,jcen(3),rcen,rmax,cefac,m(nbod),x(3,nbod)
      real(8) v(3,nbod),s(3,nbod),rho(nbod),rceh(nbod),rphys(nbod)
      real(8) rce(nbod),rcrit(nbod)
      character(8) id(nbod)
      character(80) outfile
!
! Local
      integer j
      real(8) a(NMAX),hill(NMAX),temp,amin,vmax,k_2,rhocgs,rcen_2
      character(80) header,c(NMAX)
      character(8) mio_re2c, mio_fl2c
!
!------------------------------------------------------------------------------
!
      rhocgs = AU * AU * AU * K2 / MSUN
      k_2 = 1.d0 / K2
      rcen_2 = 1.d0 / (rcen * rcen)
      amin = HUGE
!
! Calculate the Hill radii
      call mce_hill (nbod,m,x,v,hill,a)
!
! Determine the maximum close-encounter distances, and the physical radii
      temp = 2.25d0 * m(1) / PI
      do j = 2, nbod
        rce(j)   = hill(j) * rceh(j)
        rphys(j) = hill(j) / a(j) * (temp/rho(j))**THIRD
        amin = min (a(j), amin)
      end do
!
! If required, calculate the changeover distance used by hybrid algorithm
      if (rcritflag.eq.1) then
        vmax = sqrt (m(1) / amin)
        temp = N2 * h * vmax
        do j = 2, nbod
          rcrit(j) = max(hill(j) * cefac, temp)
        end do
      end if
!
! Write list of object's identities to close-encounter output file
      header(1:8)   = mio_fl2c (tstart)
      header(9:16)  = mio_re2c (dble(nbig - 1),   0.d0, 11239423.99d0)
      header(12:19) = mio_re2c (dble(nbod - nbig),0.d0, 11239423.99d0)
      header(15:22) = mio_fl2c (m(1) * k_2)
      header(23:30) = mio_fl2c (jcen(1) * rcen_2)
      header(31:38) = mio_fl2c (jcen(2) * rcen_2 * rcen_2)
      header(39:46) = mio_fl2c (jcen(3) * rcen_2 * rcen_2 * rcen_2)
      header(47:54) = mio_fl2c (rcen)
      header(55:62) = mio_fl2c (rmax)
!
      do j = 2, nbod
        c(j)(1:8) = mio_re2c (dble(j - 1), 0.d0, 11239423.99d0)
        c(j)(4:11) = id(j)
        c(j)(12:19) = mio_fl2c (m(j) * k_2)
        c(j)(20:27) = mio_fl2c (s(1,j) * k_2)
        c(j)(28:35) = mio_fl2c (s(2,j) * k_2)
        c(j)(36:43) = mio_fl2c (s(3,j) * k_2)
        c(j)(44:51) = mio_fl2c (rho(j) / rhocgs)
      end do
!
! Write compressed output to file
! AVIMANDELL
! Commenting collision output
!  50  open (22, file=outfile, status='old', access='append', err=50)
!      write (22,'(a1,a2,i2,a62,i1)') char(12),'6a',algor,header(1:62),
!     %  opt(4)
!      do j = 2, nbod
!        write (22,'(a51)') c(j)(1:51)
!      enddo
!      close (22)
!
!------------------------------------------------------------------------------
!
      return
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MCE_MERG.FOR    (ErikSoft   2 October 2000)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Author: John E. Chambers
!
! Merges objects I and J inelastically to produce a single new body by 
! conserving mass and linear momentum.
!   If J <= NBIG, then J is a Big body
!   If J >  NBIG, then J is a Small body
!   If I = 0, then I is the central body
!
! N.B. All coordinates and velocities must be with respect to central body.
! ===
!
!------------------------------------------------------------------------------
!
      subroutine mce_merg (jcen,i,j,nbod,nbig,m,xh,vh,s,stat,elost)
!
      implicit none
      include 'mercury95.inc'
!
! Input/Output
      integer i, j, nbod, nbig, stat(nbod)
      real(8) jcen(3),m(nbod),xh(3,nbod),vh(3,nbod),s(3,nbod),elost
!
! Local
      integer k
      real(8) tmp1, tmp2, dx, dy, dz, du, dv, dw, msum, mredu, msum_1
      real(8) e0, e1, l2
!
!------------------------------------------------------------------------------
!
! If a body hits the central body
      if (i.le.1) then
        call mxx_en (jcen,nbod,nbig,m,xh,vh,s,e0,l2)
!
! If a body hit the central body...
        msum   = m(1) + m(j)
        msum_1 = 1.d0 / msum
        mredu  = m(1) * m(j) * msum_1
        dx = xh(1,j)
        dy = xh(2,j)
        dz = xh(3,j)
        du = vh(1,j)
        dv = vh(2,j)
        dw = vh(3,j)
!
! Calculate new spin angular momentum of the central body
        s(1,1) = s(1,1)  +  s(1,j)  +  mredu * (dy * dw  -  dz * dv)
        s(2,1) = s(2,1)  +  s(2,j)  +  mredu * (dz * du  -  dx * dw)
        s(3,1) = s(3,1)  +  s(3,j)  +  mredu * (dx * dv  -  dy * du)
!
! Calculate shift in barycentric coords and velocities of central body
        tmp2 = m(j) * msum_1
        xh(1,1) = tmp2 * xh(1,j)
        xh(2,1) = tmp2 * xh(2,j)
        xh(3,1) = tmp2 * xh(3,j)
        vh(1,1) = tmp2 * vh(1,j)
        vh(2,1) = tmp2 * vh(2,j)
        vh(3,1) = tmp2 * vh(3,j)
        m(1) = msum
        m(j) = 0.d0
        s(1,j) = 0.d0
        s(2,j) = 0.d0
        s(3,j) = 0.d0
!
! Shift the heliocentric coordinates and velocities of all bodies
        do k = 2, nbod
          xh(1,k) = xh(1,k) - xh(1,1)
          xh(2,k) = xh(2,k) - xh(2,1)
          xh(3,k) = xh(3,k) - xh(3,1)
          vh(1,k) = vh(1,k) - vh(1,1)
          vh(2,k) = vh(2,k) - vh(2,1)
          vh(3,k) = vh(3,k) - vh(3,1)
        end do
!
! Calculate energy loss due to the collision
        call mxx_en (jcen,nbod,nbig,m,xh,vh,s,e1,l2)
        elost = elost + (e0 - e1)
      else
!
! If two bodies collided...
        msum   = m(i) + m(j)
        msum_1 = 1.d0 / msum
        mredu  = m(i) * m(j) * msum_1
        dx = xh(1,i) - xh(1,j)
        dy = xh(2,i) - xh(2,j)
        dz = xh(3,i) - xh(3,j)
        du = vh(1,i) - vh(1,j)
        dv = vh(2,i) - vh(2,j)
        dw = vh(3,i) - vh(3,j)
!
! Calculate energy loss due to the collision
        elost = elost  +  .5d0 * mredu * (du*du + dv*dv + dw*dw) &
             -  m(i) * m(j) / sqrt(dx*dx + dy*dy + dz*dz)
!
! Calculate spin angular momentum of the new body
        s(1,i) = s(1,i)  +  s(1,j)  +  mredu * (dy * dw  -  dz * dv)
        s(2,i) = s(2,i)  +  s(2,j)  +  mredu * (dz * du  -  dx * dw)
        s(3,i) = s(3,i)  +  s(3,j)  +  mredu * (dx * dv  -  dy * du)
!
! Calculate new coords and velocities by conserving centre of mass & momentum
        tmp1 = m(i) * msum_1
        tmp2 = m(j) * msum_1
        xh(1,i) = xh(1,i) * tmp1  +  xh(1,j) * tmp2
        xh(2,i) = xh(2,i) * tmp1  +  xh(2,j) * tmp2
        xh(3,i) = xh(3,i) * tmp1  +  xh(3,j) * tmp2
        vh(1,i) = vh(1,i) * tmp1  +  vh(1,j) * tmp2
        vh(2,i) = vh(2,i) * tmp1  +  vh(2,j) * tmp2
        vh(3,i) = vh(3,i) * tmp1  +  vh(3,j) * tmp2
        m(i) = msum
      end if
!
! Flag the lost body for removal, and move it away from the new body
      stat(j) = -2
      xh(1,j) = -xh(1,j)
      xh(2,j) = -xh(2,j)
      xh(3,j) = -xh(3,j)
      vh(1,j) = -vh(1,j)
      vh(2,j) = -vh(2,j)
      vh(3,j) = -vh(3,j)
      m(j)   = 0.d0
      s(1,j) = 0.d0
      s(2,j) = 0.d0
      s(3,j) = 0.d0
!
!------------------------------------------------------------------------------
!
      return
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MCE_MIN.FOR    (ErikSoft  1 December 1998)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Calculates minimum value of a quantity D, within an interval H, given initial
! and final values D0, D1, and their derivatives D0T, D1T, using third-order
! (i.e. cubic) interpolation.
!
! Also calculates the value of the independent variable T at which D is a
! minimum, with respect to the epoch of D1.
!
! N.B. The routine assumes that only one minimum is present in the interval H.
! ===
!------------------------------------------------------------------------------
!
      subroutine mce_min (d0,d1,d0t,d1t,h,d2min,tmin)
!
      implicit none
!
! Input/Output
      real(8) d0,d1,d0t,d1t,h,d2min,tmin
!
! Local
      real(8) a,b,c,temp,tau
!
!------------------------------------------------------------------------------
!
      if (d0t*h.gt.0.or.d1t*h.lt.0) then
        if (d0.le.d1) then
          d2min = d0
          tmin = -h
        else
          d2min = d1
          tmin = 0.d0
        end if
      else
        temp = 6.d0*(d0 - d1)
        a = temp + 3.d0*h*(d0t + d1t)
        b = temp + 2.d0*h*(d0t + 2.d0*d1t)
        c = h * d1t
!
        temp =-.5d0*(b + sign (sqrt(max(b*b - 4.d0*a*c,0.d0)), b) )
        if (temp.eq.0) then
          tau = 0.d0
        else
          tau = c / temp
        end if
!
! Make sure TAU falls in the interval -1 < TAU < 0
        tau = min(tau, 0.d0)
        tau = max(tau, -1.d0)
!
! Calculate TMIN and D2MIN
        tmin = tau * h
        temp = 1.d0 + tau
        d2min = tau*tau*((3.d0+2.d0*tau)*d0 + temp*h*d0t) &
             + temp*temp*((1.d0-2.d0*tau)*d1 + tau*h*d1t)
!
! Make sure D2MIN is not negative
        d2min = max(d2min, 0.d0)
      end if
!
!------------------------------------------------------------------------------
!
      return
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    MCE_SNIF.FOR    (ErikSoft   3 October 2000)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Given initial and final coordinates and velocities X and V, and a timestep 
! H, the routine estimates which objects were involved in a close encounter
! during the timestep. The routine examines all objects with indices I >= I0.
!
! Returns an array CE, which for each object is: 
!                           0 if it will undergo no encounters
!                           2 if it will pass within RCRIT of a Big body
!
! Also returns arrays ICE and JCE, containing the indices of each pair of
! objects estimated to have undergone an encounter.
!
! N.B. All coordinates must be with respect to the central body!!!!
! ===
!
!------------------------------------------------------------------------------
!
      subroutine mce_snif (h,i0,nbod,nbig,x0,v0,x1,v1,rcrit,ce,nce,ice, &
       jce)
!
      implicit none
      include 'mercury95.inc'
!
! Input/Output
      integer i0,nbod,nbig,ce(nbod),nce,ice(NMAX),jce(NMAX)
      real(8) x0(3,nbod),v0(3,nbod),x1(3,nbod),v1(3,nbod),h,rcrit(nbod)
!
! Local
      integer i,j
      real(8) d0,d1,d0t,d1t,d2min,temp,tmin,rc,rc2
      real(8) dx0,dy0,dz0,du0,dv0,dw0,dx1,dy1,dz1,du1,dv1,dw1
      real(8) xmin(NMAX),xmax(NMAX),ymin(NMAX),ymax(NMAX)
!
!------------------------------------------------------------------------------
!
      if (i0.le.0) i0 = 2
      nce = 0
      do j = 2, nbod
        ce(j) = 0
      end do
!
! Calculate maximum and minimum values of x and y coordinates
      call mce_box (nbod,h,x0,v0,x1,v1,xmin,xmax,ymin,ymax)
!
! Adjust values for the Big bodies by symplectic close-encounter distance
      do j = i0, nbig
        xmin(j) = xmin(j) - rcrit(j)
        xmax(j) = xmax(j) + rcrit(j)
        ymin(j) = ymin(j) - rcrit(j)
        ymax(j) = ymax(j) + rcrit(j)
      end do
!
! Identify pairs whose X-Y boxes overlap, and calculate minimum separation
      do i = i0, nbig
        do j = i + 1, nbod
          if (xmax(i).ge.xmin(j).and.xmax(j).ge.xmin(i) &
           .and.ymax(i).ge.ymin(j).and.ymax(j).ge.ymin(i)) then
!
! Determine the maximum separation that would qualify as an encounter
            rc = max(rcrit(i), rcrit(j))
            rc2 = rc * rc
!
! Calculate initial and final separations
            dx0 = x0(1,i) - x0(1,j)
            dy0 = x0(2,i) - x0(2,j)
            dz0 = x0(3,i) - x0(3,j)
            dx1 = x1(1,i) - x1(1,j)
            dy1 = x1(2,i) - x1(2,j)
            dz1 = x1(3,i) - x1(3,j)
            d0 = dx0*dx0 + dy0*dy0 + dz0*dz0
            d1 = dx1*dx1 + dy1*dy1 + dz1*dz1
!
! Check for a possible minimum in between
            du0 = v0(1,i) - v0(1,j)
            dv0 = v0(2,i) - v0(2,j)
            dw0 = v0(3,i) - v0(3,j)
            du1 = v1(1,i) - v1(1,j)
            dv1 = v1(2,i) - v1(2,j)
            dw1 = v1(3,i) - v1(3,j)
            d0t = (dx0*du0 + dy0*dv0 + dz0*dw0) * 2.d0
            d1t = (dx1*du1 + dy1*dv1 + dz1*dw1) * 2.d0
!
! If separation derivative changes sign, find the minimum separation
            d2min = HUGE
            if (d0t*h.le.0.and.d1t*h.ge.0) call mce_min (d0,d1,d0t,d1t, &
             h,d2min,tmin)
!
! If minimum separation is small enough, flag this as a possible encounter
            temp = min (d0,d1,d2min)
            if (temp.le.rc2) then
              ce(i) = 2
              ce(j) = 2
              nce = nce + 1
              ice(nce) = i
              jce(nce) = j
            end if
          end if
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
!    MCE_STAT.FOR    (ErikSoft   1 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Returns details of all close-encounter minima involving at least one Big
! body during a timestep. The routine estimates minima using the initial
! and final coordinates X(0),X(1) and velocities V(0),V(1) of the step, and
! the stepsize H.
!
!  ICLO, JCLO contain the indices of the objects
!  DCLO is their minimum separation
!  TCLO is the time of closest approach relative to current time
!
! The routine also checks for collisions/near misses given the physical radii 
! RPHYS, and returns the time THIT of the collision/near miss closest to the
! start of the timestep, and the identities IHIT and JHIT of the objects
! involved.
!
!  NHIT = +1 implies a collision
!         -1    "    a near miss
!
! N.B. All coordinates & velocities must be with respect to the central body!!
! ===
!------------------------------------------------------------------------------
!
      subroutine mce_stat (time,h,rcen,nbod,nbig,m,x0,v0,x1,v1,rce, &
       rphys,nclo,iclo,jclo,dclo,tclo,ixvclo,jxvclo,nhit,ihit,jhit, &
       chit,dhit,thit,thit1,nowflag,stat,outfile,mem,lmem)
!
      implicit none
      include 'mercury95.inc'
!
! Input/Output
      integer nbod,nbig,stat(nbod),nowflag
      integer nclo,iclo(CMAX),jclo(CMAX)
      integer nhit,ihit(CMAX),jhit(CMAX),chit(CMAX),lmem(NMESS)
      real(8) time,h,rcen,m(nbod),x0(3,nbod),v0(3,nbod)
      real(8) x1(3,nbod),v1(3,nbod),rce(nbod),rphys(nbod)
      real(8) dclo(CMAX),tclo(CMAX),thit(CMAX),dhit(CMAX),thit1
      real(8) ixvclo(6,CMAX),jxvclo(6,CMAX)
      character(80) outfile,mem(NMESS)
!
! Local
      integer i,j
      real(8) d0,d1,d0t,d1t,hm1,tmp0,tmp1
      real(8) dx0,dy0,dz0,du0,dv0,dw0,dx1,dy1,dz1,du1,dv1,dw1
      real(8) xmin(NMAX),xmax(NMAX),ymin(NMAX),ymax(NMAX)
      real(8) d2min,d2ce,d2near,d2hit,temp,tmin
!
!------------------------------------------------------------------------------
!
      nhit = 0
      thit1 = sign(HUGE, h)
      hm1 = 1.d0 / h
!
! Calculate maximum and minimum values of x and y coords for each object
      call mce_box (nbod,h,x0,v0,x1,v1,xmin,xmax,ymin,ymax)
!
! Adjust values by the maximum close-encounter radius plus a fudge factor
      do j = 2, nbod
        temp = rce(j) * 1.2d0
        xmin(j) = xmin(j)  -  temp
        xmax(j) = xmax(j)  +  temp
        ymin(j) = ymin(j)  -  temp
        ymax(j) = ymax(j)  +  temp
      end do
!
! Check for close encounters between each pair of objects
      do i = 2, nbig
        do j = i + 1, nbod
          if (   xmax(i).ge.xmin(j).and.xmax(j).ge.xmin(i) &
           .and.ymax(i).ge.ymin(j).and.ymax(j).ge.ymin(i) &
           .and.stat(i).ge.0.and.stat(j).ge.0) then
!
! If the X-Y boxes for this pair overlap, check circumstances more closely
            dx0 = x0(1,i) - x0(1,j)
            dy0 = x0(2,i) - x0(2,j)
            dz0 = x0(3,i) - x0(3,j)
            du0 = v0(1,i) - v0(1,j)
            dv0 = v0(2,i) - v0(2,j)
            dw0 = v0(3,i) - v0(3,j)
            d0t = (dx0*du0 + dy0*dv0 + dz0*dw0) * 2.d0
!
            dx1 = x1(1,i) - x1(1,j)
            dy1 = x1(2,i) - x1(2,j)
            dz1 = x1(3,i) - x1(3,j)
            du1 = v1(1,i) - v1(1,j)
            dv1 = v1(2,i) - v1(2,j)
            dw1 = v1(3,i) - v1(3,j)
            d1t = (dx1*du1 + dy1*dv1 + dz1*dw1) * 2.d0
!
! Estimate minimum separation during the time interval, using interpolation
            d0 = dx0*dx0 + dy0*dy0 + dz0*dz0
            d1 = dx1*dx1 + dy1*dy1 + dz1*dz1
            call mce_min (d0,d1,d0t,d1t,h,d2min,tmin)
            d2ce  = max (rce(i), rce(j))
            d2hit = rphys(i) + rphys(j)
            d2ce   = d2ce  * d2ce
            d2hit  = d2hit * d2hit
            d2near = d2hit * 4.d0
!
! If the minimum separation qualifies as an encounter or if a collision
! is in progress, store details
            if ((d2min.le.d2ce.and.d0t*h.le.0.and.d1t*h.ge.0) &
             .or.(d2min.le.d2hit)) then
              nclo = nclo + 1
              if (nclo.gt.CMAX) then
 230            open (23,file=outfile,status='old',access='append', &
                 err=230)
                write (23,'(/,2a,/,a)') mem(121)(1:lmem(121)), &
                 mem(132)(1:lmem(132)),mem(82)(1:lmem(82))
                close (23)
              else
                tclo(nclo) = tmin + time
                dclo(nclo) = sqrt (max(0.d0,d2min))
                iclo(nclo) = i
                jclo(nclo) = j
!
! Make sure the more massive body is listed first
                if (m(j).gt.m(i).and.j.le.nbig) then
                  iclo(nclo) = j
                  jclo(nclo) = i
                end if
!
! Make linear interpolation to get coordinates at time of closest approach
                tmp0 = 1.d0 + tmin*hm1
                tmp1 = -tmin*hm1
                ixvclo(1,nclo) = tmp0 * x0(1,i)  +  tmp1 * x1(1,i)
                ixvclo(2,nclo) = tmp0 * x0(2,i)  +  tmp1 * x1(2,i)
                ixvclo(3,nclo) = tmp0 * x0(3,i)  +  tmp1 * x1(3,i)
                ixvclo(4,nclo) = tmp0 * v0(1,i)  +  tmp1 * v1(1,i)
                ixvclo(5,nclo) = tmp0 * v0(2,i)  +  tmp1 * v1(2,i)
                ixvclo(6,nclo) = tmp0 * v0(3,i)  +  tmp1 * v1(3,i)
                jxvclo(1,nclo) = tmp0 * x0(1,j)  +  tmp1 * x1(1,j)
                jxvclo(2,nclo) = tmp0 * x0(2,j)  +  tmp1 * x1(2,j)
                jxvclo(3,nclo) = tmp0 * x0(3,j)  +  tmp1 * x1(3,j)
                jxvclo(4,nclo) = tmp0 * v0(1,j)  +  tmp1 * v1(1,j)
                jxvclo(5,nclo) = tmp0 * v0(2,j)  +  tmp1 * v1(2,j)
                jxvclo(6,nclo) = tmp0 * v0(3,j)  +  tmp1 * v1(3,j)
              end if
            end if
!
! Check for a near miss or collision
            if (d2min.le.d2near) then
              nhit = nhit + 1
              ihit(nhit) = i
              jhit(nhit) = j
              thit(nhit) = tmin + time
              dhit(nhit) = sqrt(d2min)
              chit(nhit) = -1
              if (d2min.le.d2hit) chit(nhit) = 1
!
! Make sure the more massive body is listed first
              if (m(jhit(nhit)).gt.m(ihit(nhit)).and.j.le.nbig) then
                ihit(nhit) = j
                jhit(nhit) = i
              end if
!
! Is this the collision closest to the start of the time step?
              if ((tmin-thit1)*h.lt.0) then
                thit1 = tmin
                nowflag = 0
                if (d1.le.d2hit) nowflag = 1
              end if
            end if
          end if
!
! Move on to the next pair of objects
        end do
      end do
!
!------------------------------------------------------------------------------
!
      return
      end
!

