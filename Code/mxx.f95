!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MXX_EJEC.FOR    (ErikSoft   2 November 2000)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Calculates the distance from the central body of each object with index
! I >= I0. If this distance exceeds RMAX, the object is flagged for ejection 
! (STAT set to -3). If any object is to be ejected, EJFLAG = 1 on exit,
! otherwise EJFLAG = 0.
!
! Also updates the values of EN(3) and AM(3)---the change in energy and
! angular momentum due to collisions and ejections.
!
!
! N.B. All coordinates must be with respect to the central body!!
! ===
!
!------------------------------------------------------------------------------
!
      subroutine mxx_ejec (time,tstart,rmax,en,am,jcen,i0,nbod,nbig,m,x, &
       v,s,stat,id,opt,ejflag,outfile,mem,lmem)
!
      implicit none
      include 'mercury95.inc'
!
! Input/Output
      integer i0, nbod, nbig, stat(nbod), opt(8), ejflag, lmem(NMESS)
      real(8) time, tstart, rmax, en(3), am(3), jcen(3)
      real(8) m(nbod), x(3,nbod), v(3,nbod), s(3,nbod)
      character(80) outfile, mem(NMESS)
      character(8) id(nbod)
!
! Local
      integer j, year, month
      real(8) r2,rmax2,t1,e,l
      character(38) flost
      character(6) tstring
!
!------------------------------------------------------------------------------
!
      if (i0.le.0) i0 = 2
      ejflag = 0
      rmax2 = rmax * rmax
!
! Calculate initial energy and angular momentum
      call mxx_en (jcen,nbod,nbig,m,x,v,s,e,l)
!
! Flag each object which is ejected, and set its mass to zero
      do j = i0, nbod
        r2 = x(1,j)*x(1,j) + x(2,j)*x(2,j) + x(3,j)*x(3,j)
        if (r2.gt.rmax2) then
          ejflag = 1
          stat(j) = -3
          m(j) = 0.d0
          s(1,j) = 0.d0
          s(2,j) = 0.d0
          s(3,j) = 0.d0
!
! Write message to information file
  20      open  (23,file=outfile,status='old',access='append',err=20)
          if (opt(3).eq.1) then
            call mio_jd2y (time,year,month,t1)
            flost = '(1x,a8,a,i10,1x,i2,1x,f8.5)'
            write (23,flost) id(j),mem(68)(1:lmem(68)),year,month,t1
          else
            if (opt(3).eq.3) then
              t1 = (time - tstart) / 365.25d0
              tstring = mem(2)
              flost = '(1x,a8,a,f18.7,a)'
            else
              if (opt(3).eq.0) t1 = time
              if (opt(3).eq.2) t1 = time - tstart
              tstring = mem(1)
              flost = '(1x,a8,a,f18.5,a)'
            end if
            write (23,flost) id(j),mem(68)(1:lmem(68)),t1,tstring
          end if
          close (23)
        end if
      end do
!
! If ejections occurred, update ELOST and LLOST
      if (ejflag.ne.0) then
        call mxx_en (jcen,nbod,nbig,m,x,v,s,en(2),am(2))
        en(3) = en(3) + (e - en(2))
        am(3) = am(3) + (l - am(2))
      end if
!
!------------------------------------------------------------------------------
!
      return
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MXX_ELIM.FOR    (ErikSoft   13 February 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Removes any objects with STAT < 0 (i.e. those that have been flagged for 
! removal) and reindexes all the appropriate arrays for the remaining objects.
!
!------------------------------------------------------------------------------
!
      subroutine mxx_elim (nbod,nbig,m,x,v,s,rho,rceh,rcrit,ngf,stat, &
       id,mem,lmem,outfile,nelim)
!
      implicit none
      include 'mercury95.inc'
!
! Input/Output
      integer nbod, nbig, nelim, stat(nbod), lmem(NMESS)
      real(8) m(nbod), x(3,nbod), v(3,nbod), s(3,nbod)
      real(8) rho(nbod), rceh(nbod), rcrit(nbod), ngf(4,nbod)
      character(8) id(nbod)
      character(80) outfile, mem(NMESS)
!
! Local
      integer j, k, l, nbigelim, elim(NMAX+1)
!
!------------------------------------------------------------------------------
!
! Find out how many objects are to be removed
      nelim = 0
      nbigelim = 0
      do j = 2, nbod
        if (stat(j).lt.0) then
          nelim = nelim + 1
          elim(nelim) = j
          if (j.le.nbig) nbigelim = nbigelim + 1
        end if
      end do
      elim(nelim+1) = nbod + 1
!
! Eliminate unwanted objects
      do k = 1, nelim
        do j = elim(k) - k + 1, elim(k+1) - k - 1
          l = j + k
          x(1,j) = x(1,l)
          x(2,j) = x(2,l)
          x(3,j) = x(3,l)
          v(1,j) = v(1,l)
          v(2,j) = v(2,l)
          v(3,j) = v(3,l)
          m(j)   = m(l)
          s(1,j) = s(1,l)
          s(2,j) = s(2,l)
          s(3,j) = s(3,l)
          rho(j) = rho(l)
          rceh(j) = rceh(l)
          stat(j) = stat(l)
          id(j) = id(l)
          ngf(1,j) = ngf(1,l)
          ngf(2,j) = ngf(2,l)
          ngf(3,j) = ngf(3,l)
          ngf(4,j) = ngf(4,l)
        end do
      end do
!
! Update total number of bodies and number of Big bodies
      nbod = nbod - nelim
      nbig = nbig - nbigelim
!
! If no massive bodies remain, stop the integration
      if (nbig.lt.1) then
  10    open (23,file=outfile,status='old',access='append',err=10)
        write (23,'(2a)') mem(81)(1:lmem(81)),mem(124)(1:lmem(124))
        close (23)
        stop
      end if
!
!------------------------------------------------------------------------------
!
      return
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!     MXX_EN.FOR    (ErikSoft   21 February 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Calculates the total energy and angular-momentum for a system of objects
! with masses M, coordinates X, velocities V and spin angular momenta S.
!
! N.B. All coordinates and velocities must be with respect to the central
! ===  body.
!
!------------------------------------------------------------------------------
!
      subroutine mxx_en  (jcen,nbod,nbig,m,xh,vh,s,e,l2)
!
      implicit none
      include 'mercury95.inc'
!
! Input/Output
      integer nbod,nbig
      real(8) jcen(3),m(nbod),xh(3,nbod),vh(3,nbod),s(3,nbod),e,l2
!
! Local
      integer j,k,iflag,itmp(8)
      real(8) x(3,NMAX),v(3,NMAX),temp,dx,dy,dz,r2,tmp,ke,pe,l(3)
      real(8) r_1,r_2,r_4,r_6,u2,u4,u6,tmp2(4,NMAX)
!
!------------------------------------------------------------------------------
!
      ke = 0.d0
      pe = 0.d0
      l(1) = 0.d0
      l(2) = 0.d0
      l(3) = 0.d0
!
! Convert to barycentric coordinates and velocities
      call mco_h2b(jcen,nbod,nbig,temp,m,xh,vh,x,v)
!
! Do the spin angular momenta first (probably the smallest terms)
      do j = 1, nbod
        l(1) = l(1) + s(1,j)
        l(2) = l(2) + s(2,j)
        l(3) = l(3) + s(3,j)
      end do
!
! Orbital angular momentum and kinetic energy terms
      do j = 1, nbod
        l(1) = l(1)  +  m(j)*(x(2,j) * v(3,j)  -  x(3,j) * v(2,j))
        l(2) = l(2)  +  m(j)*(x(3,j) * v(1,j)  -  x(1,j) * v(3,j))
        l(3) = l(3)  +  m(j)*(x(1,j) * v(2,j)  -  x(2,j) * v(1,j))
        ke = ke + m(j)*(v(1,j)*v(1,j)+v(2,j)*v(2,j)+v(3,j)*v(3,j))
      end do
!
! Potential energy terms due to pairs of bodies
      do j = 2, nbod
        tmp = 0.d0
        do k = j + 1, nbod
          dx = x(1,k) - x(1,j)
          dy = x(2,k) - x(2,j)
          dz = x(3,k) - x(3,j)
          r2 = dx*dx + dy*dy + dz*dz
          if (r2.ne.0) tmp = tmp + m(k) / sqrt(r2)
        end do
        pe = pe  -  tmp * m(j)
      end do
!
! Potential energy terms involving the central body
      do j = 2, nbod
        dx = x(1,j) - x(1,1)
        dy = x(2,j) - x(2,1)
        dz = x(3,j) - x(3,1)
        r2 = dx*dx + dy*dy + dz*dz
        if (r2.ne.0) pe = pe  -  m(1) * m(j) / sqrt(r2)
      end do
!
! Corrections for oblateness
      if (jcen(1).ne.0.or.jcen(2).ne.0.or.jcen(3).ne.0) then
        do j = 2, nbod
          r2 = xh(1,j)*xh(1,j) + xh(2,j)*xh(2,j) + xh(3,j)*xh(3,j)
          r_1 = 1.d0 / sqrt(r2)
          r_2 = r_1 * r_1
          r_4 = r_2 * r_2
          r_6 = r_4 * r_2
          u2 = xh(3,j) * xh(3,j) * r_2
          u4 = u2 * u2
          u6 = u4 * u2
          pe = pe + m(1) * m(j) * r_1 &
            * (jcen(1) * r_2 * (1.5d0*u2 - 0.5d0) &
            +  jcen(2) * r_4 * (4.375d0*u4 - 3.75d0*u2 + .375d0) &
            +  jcen(3) * r_6 &
            *(14.4375d0*u6 - 19.6875d0*u4 + 6.5625d0*u2 - .3125d0))
        end do
      end if
!
      e = .5d0 * ke  +  pe
      l2 = sqrt(l(1)*l(1) + l(2)*l(2) + l(3)*l(3))
!
!------------------------------------------------------------------------------
!
      return    
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!     MXX_JAC.FOR    (ErikSoft   2 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Calculates the Jacobi constant for massless particles. This assumes that
! there are only 2 massive bodies (including the central body) moving on
! circular orbits.
!
! N.B. All coordinates and velocities must be heliocentric!!
! ===
!
!------------------------------------------------------------------------------
!
      subroutine mxx_jac (jcen,nbod,nbig,m,xh,vh,jac)
!
      implicit none
      include 'mercury95.inc'
!
! Input/Output
      integer nbod,nbig
      real(8) jcen(3),m(nbod),xh(3,nbod),vh(3,nbod)
!
! Local
      integer j,itmp(8),iflag
      real(8) x(3,NMAX),v(3,NMAX),temp,dx,dy,dz,r,d,a2,n,jac(NMAX)
      real(8) tmp2(4,NMAX)
!
!------------------------------------------------------------------------------
!
      call mco_h2b(jcen,nbod,nbig,temp,m,xh,vh,x,v)
      dx = x(1,2) - x(1,1)
      dy = x(2,2) - x(2,1)
      dz = x(3,2) - x(3,1)
      a2 = dx*dx + dy*dy + dz*dz
      n = sqrt((m(1)+m(2)) / (a2*sqrt(a2)))
!
      do j = nbig + 1, nbod
        dx = x(1,j) - x(1,1)
        dy = x(2,j) - x(2,1)
        dz = x(3,j) - x(3,1)
        r = sqrt(dx*dx + dy*dy + dz*dz)
        dx = x(1,j) - x(1,2)
        dy = x(2,j) - x(2,2)
        dz = x(3,j) - x(3,2)
        d = sqrt(dx*dx + dy*dy + dz*dz)
!
        jac(j) = .5d0*(v(1,j)*v(1,j) + v(2,j)*v(2,j) + v(3,j)*v(3,j)) &
              - m(1)/r - m(2)/d - n*(x(1,j)*v(2,j) - x(2,j)*v(1,j))
      end do
!
!------------------------------------------------------------------------------
!
      return    
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MXX_SORT.FOR    (ErikSoft 24 May 1997)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Sorts an array X, of size N, using Shell's method. Also returns an array
! INDEX that gives the original index of every item in the sorted array X.
!
! N.B. The maximum array size is 29523.
! ===
!
!------------------------------------------------------------------------------
!
      subroutine mxx_sort (n,x,index)
!
      implicit none
!
! Input/Output
      integer n,index(n)
      real(8) x(n)
!
! Local
      integer i,j,k,l,m,inc,incarr(9),iy
      real(8) y
      data incarr/1,4,13,40,121,364,1093,3280,9841/
!
!------------------------------------------------------------------------------
!
      do i = 1, n
        index(i) = i
      end do
!
      m = 0
  10  m = m + 1
      if (incarr(m).lt.n) goto 10
      m = m - 1
!
      do i = m, 1, -1
        inc = incarr(i)
        do j = 1, inc
          do k = inc, n - j, inc
            y = x(j+k)
            iy = index(j+k)
            do l = j + k - inc, j, -inc
              if (x(l).le.y) goto 20
              x(l+inc) = x(l)
              index(l+inc) = index(l)
            end do
  20        x(l+inc) = y
            index(l+inc) = iy
          end do
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
!      MXX_SYNC.FOR    (ErikSoft   2 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Synchronizes the epochs of NBIG Big bodies (having a common epoch) and
! NBOD-NBIG Small bodies (possibly having differing epochs), for an 
! integration using MERCURY.
! The Small bodies are picked up in order starting with the one with epoch
! furthest from the time, TSTART, at which the main integration will begin
! producing output.
!
! N.B. The synchronization integrations use Everhart's RA15 routine.
! ---
!
!------------------------------------------------------------------------------
!
      subroutine mxx_sync (time,tstart,h0,tol,jcen,nbod,nbig,m,x,v,s, &
       rho,rceh,stat,id,epoch,ngf,opt,ngflag)
!
      implicit none
      include 'mercury95.inc'
!
! Input/Output
      integer nbod,nbig,ngflag,opt(8),stat(nbod)
      real(8) time,tstart,h0,tol,jcen(3),m(nbod),x(3,nbod),v(3,nbod)
      real(8) s(3,nbod),rceh(nbod),rho(nbod),epoch(nbod),ngf(4,nbod)
      character(8) id(nbod)
!
! Local
      integer j,k,l,nsml,nsofar,indx(NMAX),itemp,jtemp(NMAX)
      integer raflag,nce,ice(NMAX),jce(NMAX)
      real(8) temp,epsml(NMAX),rtemp(NMAX)
      real(8) h,hdid,tsmall,rphys(NMAX),rcrit(NMAX)
      character(8) ctemp(NMAX)
      external mfo_all
!
!------------------------------------------------------------------------------
!
! Reorder Small bodies by epoch so that ep(1) is furthest from TSTART
      nsml = nbod - nbig
      do j = nbig + 1, nbod
        epsml(j-nbig) = epoch(j)
      end do
      call mxx_sort (nsml,epsml,indx)
!
      if (abs(epsml(1)-tstart).lt.abs(epsml(nsml)-tstart)) then
        k = nsml + 1
        do j = 1, nsml / 2
          l = k - j
          temp = epsml(j)
          epsml(j) = epsml (l)
          epsml(l) = temp
          itemp = indx(j)
          indx(j) = indx (l)
          indx(l) = itemp
        end do
      end if
!
      do j = nbig + 1, nbod
        epoch(j) = epsml(j-nbig)
      end do
!
! Reorder the other arrays associated with each Small body
      do k = 1, 3
        do j = 1, nsml
          rtemp(j) = x(k,j+nbig)
        end do
        do j = 1, nsml
          x(k,j+nbig) = rtemp(indx(j))
        end do
        do j = 1, nsml
          rtemp(j) = v(k,j+nbig)
        end do
        do j = 1, nsml
          v(k,j+nbig) = rtemp(indx(j))
        end do
        do j = 1, nsml
          rtemp(j) = s(k,j+nbig)
        end do
        do j = 1, nsml
          s(k,j+nbig) = rtemp(indx(j))
        end do
      end do
!
      do j = 1, nsml
        rtemp(j) = m(j+nbig)
      end do
      do j = 1, nsml
        m(j+nbig) = rtemp(indx(j))
      end do
      do j = 1, nsml
        rtemp(j) = rceh(j+nbig)
      end do
      do j = 1, nsml
        rceh(j+nbig) = rtemp(indx(j))
      end do
      do j = 1, nsml
        rtemp(j) = rho(j+nbig)
      end do
      do j = 1, nsml
        rho(j+nbig) = rtemp(indx(j))
      end do
!
      do j = 1, nsml
        ctemp(j) = id(j+nbig)
        jtemp(j) = stat(j+nbig)
      end do
      do j = 1, nsml
        id(j+nbig) = ctemp(indx(j))
        stat(j+nbig) = jtemp(indx(j))
      end do
!
! Integrate Small bodies up to the same epoch
      h = h0
      tsmall = h0 * 1.d-12
      raflag = 0
!
      do j = nbig + 1, nbod
        nsofar = j - 1
        do while (abs(time-epoch(j)).gt.tsmall)
          temp = epoch(j) - time
          h = sign(max(min(abs(temp),abs(h)),tsmall),temp)
          call mdt_ra15 (time,h,hdid,tol,jcen,nsofar,nbig,m,x,v,s,rphys, &
           rcrit,ngf,stat,raflag,ngflag,opt,nce,ice,jce,mfo_all)
          time = time + hdid
        end do
        raflag = 1
      end do
!
!------------------------------------------------------------------------------
!
      return
      end
!

