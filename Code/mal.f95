!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MAL_HVAR.FOR    (ErikSoft   4 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Does an integration using a variable-timestep integration algorithm. The
! particular integrator routine is ONESTEP and the algorithm must use
! coordinates with respect to the central body.
!
! N.B. This routine is also called by the synchronisation routine mxx_sync,
! ===  in which case OPFLAG = -2. Beware when making changes involving OPFLAG.
!
!------------------------------------------------------------------------------
!
      subroutine mal_hvar (time,tstart,tstop,dtout,algor,h0,tol,jcen, &
       rcen,rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s,rho,rceh, &
       stat,id,ngf,opt,opflag,ngflag,outfile,dumpfile,mem,lmem,onestep)
!
      implicit none
      include 'mercury95.inc'
!
! Input/Output
      integer algor,nbod,nbig,stat(nbod),opt(8),opflag,ngflag,ndump,nfun
      integer lmem(NMESS)
      real(8) time,tstart,tstop,dtout,h0,tol,jcen(3),rcen,rmax
      real(8) en(3),am(3),cefac,m(nbod),xh(3,nbod),vh(3,nbod)
      real(8) s(3,nbod),rho(nbod),rceh(nbod),ngf(4,nbod)
      character(8) id(nbod)
      character(80) outfile(3),dumpfile(4),mem(NMESS)
!
! Local
      integer i,j,k,n,itmp,nhit,ihit(CMAX),jhit(CMAX),chit(CMAX)
      integer dtflag,ejflag,nowflag,stopflag,nstored,ce(NMAX)
      integer nclo,iclo(CMAX),jclo(CMAX),nce,ice(NMAX),jce(NMAX)
      real(8) tmp0,h,hdid,tout,tdump,tfun,tlog,tsmall,dtdump,dtfun
      real(8) thit(CMAX),dhit(CMAX),thit1,x0(3,NMAX),v0(3,NMAX)
      real(8) rce(NMAX),rphys(NMAX),rcrit(NMAX),a(NMAX)
      real(8) dclo(CMAX),tclo(CMAX),epoch(NMAX)
      real(8) ixvclo(6,CMAX),jxvclo(6,CMAX)
      external mfo_all,onestep
!
!------------------------------------------------------------------------------
!
! Initialize variables. DTFLAG = 0 implies first ever call to ONESTEP
      dtout  = abs(dtout)
      dtdump = abs(h0) * ndump
      dtfun  = abs(h0) * nfun
      dtflag = 0
      nstored = 0
      tsmall = h0 * 1.d-8
      h = h0
      do j = 2, nbod
        ce(j) = 0.d0
      end do
!
! Calculate close-encounter limits and physical radii for massive bodies
      call mce_init (tstart,algor,h0,jcen,rcen,rmax,cefac,nbod,nbig, &
       m,xh,vh,s,rho,rceh,rphys,rce,rcrit,id,opt,outfile(2),1)
!
! Set up time of next output, times of previous dump, log and periodic effect
      if (opflag.eq.-1) then
        tout = tstart
      else
        n = int (abs (time - tstart) / dtout) + 1
        tout = tstart  +  dtout * sign (dble(n), tstop - tstart)
        if ((tstop - tstart)*(tout - tstop).gt.0) tout = tstop
      endif
      tdump = time
      tfun  = time
      tlog  = time
!
!------------------------------------------------------------------------------
!
!  MAIN  LOOP  STARTS  HERE
!
 100  continue
!
! Is it time for output ?
      if (abs(tout-time).lt.abs(tsmall).and.opflag.ge.-1) then
!
! Beware: the integration may change direction at this point!!!!
        if (opflag.eq.-1) dtflag = 0
!
! Output data for all bodies
        call mio_out (time,jcen,rcen,rmax,nbod,nbig,m,xh,vh,s,rho, &
         stat,id,opt,opflag,algor,outfile(1))
        call mio_ce (time,tstart,rcen,rmax,nbod,nbig,m,stat,id, &
         0,iclo,jclo,opt,stopflag,tclo,dclo,ixvclo,jxvclo,mem,lmem, &
         outfile,nstored,0)
        tmp0 = tstop - tout
        tout = tout + sign( min( abs(tmp0), abs(dtout) ), tmp0 )
!
! Update the data dump files
        do j = 2, nbod
          epoch(j) = time
        end do
        call mio_dump (time,tstart,tstop,dtout,algor,h,tol,jcen,rcen, &
         rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s,rho,rceh,stat, &
         id,ngf,epoch,opt,opflag,dumpfile,mem,lmem)
        tdump = time
      endif
!
! If integration has finished return to the main part of programme
      if (abs(tstop-time).le.abs(tsmall).and.opflag.ne.-1) return
!
! Set the timestep
      if (opflag.eq.-1) tmp0 = tstart - time
      if (opflag.eq.-2) tmp0 = tstop  - time
      if (opflag.ge.0)  tmp0 = tout   - time
      h = sign ( max( min( abs(tmp0), abs(h) ), tsmall), tmp0 )
!
! Save the current coordinates and velocities
      call mco_iden (jcen,nbod,nbig,h,m,xh,vh,x0,v0)
!
! Advance one timestep
      call onestep (time,h,hdid,tol,jcen,nbod,nbig,m,xh,vh,s,rphys, &
       rcrit,ngf,stat,dtflag,ngflag,opt,nce,ice,jce,mfo_all)
      time = time + hdid
!
! Check if close encounters or collisions occurred
      nclo = 0
      call mce_stat (time,h,rcen,nbod,nbig,m,x0,v0,xh,vh,rce,rphys, &
       nclo,iclo,jclo,dclo,tclo,ixvclo,jxvclo,nhit,ihit,jhit, &
       chit,dhit,thit,thit1,nowflag,stat,outfile(3),mem,lmem)
!
!------------------------------------------------------------------------------
!
!  CLOSE  ENCOUNTERS
!
! If encounter minima occurred, output details and decide whether to stop
      if (nclo.gt.0.and.opflag.ge.-1) then
        itmp = 1
        if (nhit.ne.0) itmp = 0
        call mio_ce (time,tstart,rcen,rmax,nbod,nbig,m,stat,id,nclo, &
         iclo,jclo,opt,stopflag,tclo,dclo,ixvclo,jxvclo,mem,lmem, &
         outfile,nstored,itmp)
        if (stopflag.eq.1) return
      endif
!
!------------------------------------------------------------------------------
!
!  COLLISIONS
!
! If a collision occurred, output details and resolve the collision
      if (nhit.gt.0.and.opt(2).ne.0) then
        do k = 1, nhit
          if (chit(k).eq.1) then
            i = ihit(k)
            j = jhit(k)
            call mce_coll (thit(k),tstart,en(3),jcen,i,j,nbod,nbig,m,xh, &
             vh,s,rphys,stat,id,opt,mem,lmem,outfile(3))
          endif
        end do
!
! Remove lost objects, reset flags and recompute Hill and physical radii
        call mxx_elim (nbod,nbig,m,xh,vh,s,rho,rceh,rcrit,ngf,stat, &
         id,mem,lmem,outfile(3),itmp)
        dtflag = 1
        if (opflag.ge.0) opflag = 1
        call mce_init (tstart,algor,h0,jcen,rcen,rmax,cefac,nbod,nbig, &
         m,xh,vh,s,rho,rceh,rphys,rce,rcrit,id,opt,outfile(2),1)
      endif
!
!------------------------------------------------------------------------------
!
!  COLLISIONS  WITH  CENTRAL  BODY
!
! Check for collisions
      call mce_cent (time,hdid,rcen,jcen,2,nbod,nbig,m,x0,v0,xh,vh,nhit, &
       jhit,thit,dhit,algor,ngf,ngflag)
!
! Resolve the collisions
      if (nhit.gt.0) then
        do k = 1, nhit
          i = 1
          j = jhit(k)
          call mce_coll (thit(k),tstart,en(3),jcen,i,j,nbod,nbig,m,xh, &
           vh,s,rphys,stat,id,opt,mem,lmem,outfile(3))
        end do
!
! Remove lost objects, reset flags and recompute Hill and physical radii
        call mxx_elim (nbod,nbig,m,xh,vh,s,rho,rceh,rcrit,ngf,stat, &
         id,mem,lmem,outfile(3),itmp)
        dtflag = 1
        if (opflag.ge.0) opflag = 1
        call mce_init (tstart,algor,h0,jcen,rcen,rmax,cefac,nbod,nbig, &
         m,xh,vh,s,rho,rceh,rphys,rce,rcrit,id,opt,outfile(2),0)
      endif
!
!------------------------------------------------------------------------------
!
!  DATA  DUMP  AND  PROGRESS  REPORT
!
! Do the data dump
      if (abs(time-tdump).ge.abs(dtdump).and.opflag.ge.-1) then
        do j = 2, nbod
          epoch(j) = time
        end do
        call mio_ce (time,tstart,rcen,rmax,nbod,nbig,m,stat,id, &
         0,iclo,jclo,opt,stopflag,tclo,dclo,ixvclo,jxvclo,mem,lmem, &
         outfile,nstored,0)
        call mio_dump (time,tstart,tstop,dtout,algor,h,tol,jcen,rcen, &
         rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s,rho,rceh,stat, &
         id,ngf,epoch,opt,opflag,dumpfile,mem,lmem)
        tdump = time
      endif
!
! Write a progress report to the log file
      if (abs(time-tlog).ge.abs(dtdump).and.opflag.ge.0) then
        call mxx_en (jcen,nbod,nbig,m,xh,vh,s,en(2),am(2))
        call mio_log (time,tstart,en,am,opt,mem,lmem)
        tlog = time
      endif
!
!------------------------------------------------------------------------------
!
!  CHECK  FOR  EJECTIONS  AND  DO  OTHER  PERIODIC  EFFECTS
!
      if (abs(time-tfun).ge.abs(dtfun).and.opflag.ge.-1) then
!
! Recompute close encounter limits, to allow for changes in Hill radii
        call mce_hill (nbod,m,xh,vh,rce,a)
        do j = 2, nbod
          rce(j) = rce(j) * rceh(j)
        end do
!
! Check for ejections
        call mxx_ejec (time,tstart,rmax,en,am,jcen,2,nbod,nbig,m,xh,vh, &
         s,stat,id,opt,ejflag,outfile(3),mem,lmem)
!
! Remove lost objects, reset flags and recompute Hill and physical radii
        if (ejflag.ne.0) then
          call mxx_elim (nbod,nbig,m,xh,vh,s,rho,rceh,rcrit,ngf,stat, &
           id,mem,lmem,outfile(3),itmp)
          dtflag = 1
          if (opflag.ge.0) opflag = 1
          call mce_init (tstart,algor,h0,jcen,rcen,rmax,cefac,nbod,nbig, &
           m,xh,vh,s,rho,rceh,rphys,rce,rcrit,id,opt,outfile(2),0)
        endif
        tfun = time
      endif
!
! Go on to the next time step
      goto 100
!
!------------------------------------------------------------------------------
!
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MAL_HCON.FOR    (ErikSoft   28 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Does an integration using an integrator with a constant stepsize H.
! Input and output to this routine use coordinates XH, and velocities VH,
! with respect to the central body, but the integration algorithm uses
! its own internal coordinates X, and velocities V.
!
! The programme uses the transformation routines COORD and BCOORD to change
! to and from the internal coordinates, respectively.
!
!------------------------------------------------------------------------------
!
      subroutine mal_hcon (time,tstart,tstop,dtout,algor,h0,tol,jcen, &
       rcen,rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s,rho,rceh, &
       stat,id,ngf,opt,opflag,ngflag,outfile,dumpfile,mem,lmem,onestep, &
       coord,bcoord)
!
      implicit none
      include 'mercury95.inc'
!
! Input/Output
      integer algor,nbod,nbig,stat(nbod),opt(8),opflag,ngflag
      integer lmem(NMESS),ndump,nfun
      real(8) time,tstart,tstop,dtout,h0,tol,jcen(3),rcen,rmax
      real(8) en(3),am(3),cefac,m(nbod),xh(3,nbod),vh(3,nbod)
      real(8) s(3,nbod),rho(nbod),rceh(nbod),ngf(4,nbod)
      character(8) id(nbod)
      character(80) outfile(3),dumpfile(4),mem(NMESS)
!
! Local
      integer i,j,k,n,itmp,nclo,nhit,jhit(CMAX),iclo(CMAX),jclo(CMAX)
      integer dtflag,ejflag,stopflag,colflag,nstored
      real(8) x(3,NMAX),v(3,NMAX),xh0(3,NMAX),vh0(3,NMAX)
      real(8) rce(NMAX),rphys(NMAX),rcrit(NMAX),epoch(NMAX)
      real(8) hby2,tout,tmp0,tdump,tfun,tlog,dtdump,dtfun
      real(8) dclo(CMAX),tclo(CMAX),dhit(CMAX),thit(CMAX)
      real(8) ixvclo(6,CMAX),jxvclo(6,CMAX),a(NMAX)
      external onestep,coord,bcoord

! AVIMANDELL, test for 'confirmdumpfile'
      logical test
!
!------------------------------------------------------------------------------
!
! Initialize variables. DTFLAG = 0/2: first call ever/normal
      dtout  = abs(dtout)
      dtdump = abs(h0) * ndump
      dtfun  = abs(h0) * nfun
      dtflag = 0
      nstored = 0
      hby2 = 0.500001d0 * abs(h0)
!
! Calculate close-encounter limits and physical radii
      call mce_init (tstart,algor,h0,jcen,rcen,rmax,cefac,nbod,nbig, &
       m,xh,vh,s,rho,rceh,rphys,rce,rcrit,id,opt,outfile(2),1)
!
! Set up time of next output, times of previous dump, log and periodic effect
      if (opflag.eq.-1) then
        tout = tstart
      else
        n = int (abs (time-tstart) / dtout) + 1
        tout = tstart  +  dtout * sign (dble(n), tstop - tstart)
        if ((tstop-tstart)*(tout-tstop).gt.0) tout = tstop
      endif
      tdump = time
      tfun  = time
      tlog  = time
!
! Convert to internal coordinates and velocities
      call coord (time,jcen,nbod,nbig,h0,m,xh,vh,x,v,ngf,ngflag,opt)
!
!------------------------------------------------------------------------------
!
!  MAIN  LOOP  STARTS  HERE
!
 100  continue
!
! Is it time for output ?
      if (abs(tout-time).le.hby2.and.opflag.ge.-1) then
!
! Beware: the integration may change direction at this point!!!!
        if (opflag.eq.-1.and.dtflag.ne.0) dtflag = 1
!
! Convert to heliocentric coordinates and output data for all bodies
        call bcoord (time,jcen,nbod,nbig,h0,m,x,v,xh,vh,ngf,ngflag,opt)
        call mio_out (time,jcen,rcen,rmax,nbod,nbig,m,xh,vh,s,rho, &
         stat,id,opt,opflag,algor,outfile(1))
        call mio_ce (time,tstart,rcen,rmax,nbod,nbig,m,stat,id, &
         0,iclo,jclo,opt,stopflag,tclo,dclo,ixvclo,jxvclo,mem,lmem, &
         outfile,nstored,0)
        tmp0 = tstop - tout
        tout = tout + sign( min( abs(tmp0), abs(dtout) ), tmp0 )
!
! Update the data dump files
! AVIMANDELL
!  Data dump is commented; unnecessary at every output. Mistake?
!        do j = 2, nbod
!          epoch(j) = time
!        enddo
!
!        call mio_dump (time,tstart,tstop,dtout,algor,h0,tol,jcen,rcen,
!     %    rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s,rho,rceh,stat,
!     %    id,ngf,epoch,opt,opflag,dumpfile,mem,lmem)
!        tdump = time
      endif
!
! If integration has finished, convert to heliocentric coords and return
      if (abs(tstop-time).le.hby2.and.opflag.ge.0) then
        call bcoord (time,jcen,nbod,nbig,h0,m,x,v,xh,vh,ngf,ngflag,opt)
        return
      endif
!
! Make sure the integration is heading in the right direction
 150  continue
      tmp0 = tstop - time
      if (opflag.eq.-1) tmp0 = tstart - time
      h0 = sign (h0, tmp0)
!
! Save the current heliocentric coordinates and velocities
      if (algor.eq.1) then
        call mco_iden (jcen,nbod,nbig,h0,m,x,v,xh0,vh0)
      else
        call bcoord(time,jcen,nbod,nbig,h0,m,x,v,xh0,vh0,ngf,ngflag,opt)
      endif
      call onestep (time,tstart,h0,tol,rmax,en,am,jcen,rcen,nbod,nbig, &
       m,x,v,s,rphys,rcrit,rce,stat,id,ngf,algor,opt,dtflag,ngflag, &
       opflag,colflag,nclo,iclo,jclo,dclo,tclo,ixvclo,jxvclo,outfile, &
       mem,lmem)
      time = time + h0
!
!------------------------------------------------------------------------------
!
!  CLOSE  ENCOUNTERS
!
! If encounter minima occurred, output details and decide whether to stop
      if (nclo.gt.0.and.opflag.ge.-1) then
        itmp = 1
        if (colflag.ne.0) itmp = 0
        call mio_ce (time,tstart,rcen,rmax,nbod,nbig,m,stat,id,nclo, &
         iclo,jclo,opt,stopflag,tclo,dclo,ixvclo,jxvclo,mem,lmem, &
         outfile,nstored,itmp)
        if (stopflag.eq.1) return
      endif
!
!------------------------------------------------------------------------------
!
!  COLLISIONS
!
! If collisions occurred, output details and remove lost objects
      if (colflag.ne.0) then
!
! Reindex the surviving objects
        call bcoord (time,jcen,nbod,nbig,h0,m,x,v,xh,vh,ngf,ngflag,opt)
        call mxx_elim (nbod,nbig,m,xh,vh,s,rho,rceh,rcrit,ngf,stat, &
         id,mem,lmem,outfile(3),itmp)
!
! Reset flags, and calculate new Hill radii and physical radii
        dtflag = 1
        if (opflag.ge.0) opflag = 1
        call mce_init (tstart,algor,h0,jcen,rcen,rmax,cefac,nbod,nbig, &
         m,xh,vh,s,rho,rceh,rphys,rce,rcrit,id,opt,outfile(2),1)
        call coord (time,jcen,nbod,nbig,h0,m,xh,vh,x,v,ngf,ngflag,opt)
      endif
!
!------------------------------------------------------------------------------
!
!  COLLISIONS  WITH  CENTRAL  BODY
!
! Check for collisions with the central body
      if (algor.eq.1) then
        call mco_iden(jcen,nbod,nbig,h0,m,x,v,xh,vh)
      else
        call bcoord (time,jcen,nbod,nbig,h0,m,x,v,xh,vh,ngf,ngflag,opt)
      endif
      itmp = 2
      if (algor.eq.11.or.algor.eq.12) itmp = 3
      call mce_cent (time,h0,rcen,jcen,itmp,nbod,nbig,m,xh0,vh0,xh,vh, &
       nhit,jhit,thit,dhit,algor,ngf,ngflag)
!
! If something hit the central body, restore the coords prior to this step
      if (nhit.gt.0) then
        call mco_iden (jcen,nbod,nbig,h0,m,xh0,vh0,xh,vh)
        time = time - h0
!
! Merge the object(s) with the central body
        do k = 1, nhit
          i = 1
          j = jhit(k)
          call mce_coll (thit(k),tstart,en(3),jcen,i,j,nbod,nbig,m,xh, &
           vh,s,rphys,stat,id,opt,mem,lmem,outfile(3))
        end do
!
! Remove lost objects, reset flags and recompute Hill and physical radii
        call mxx_elim (nbod,nbig,m,xh,vh,s,rho,rceh,rcrit,ngf,stat, &
         id,mem,lmem,outfile(3),itmp)
        if (opflag.ge.0) opflag = 1
        dtflag = 1
        call mce_init (tstart,algor,h0,jcen,rcen,rmax,cefac,nbod,nbig, &
         m,xh,vh,s,rho,rceh,rphys,rce,rcrit,id,opt,outfile(2),0)
        if (algor.eq.1) then
          call mco_iden (jcen,nbod,nbig,h0,m,xh,vh,x,v)
        else
          call coord (time,jcen,nbod,nbig,h0,m,xh,vh,x,v,ngf,ngflag,opt)
        endif
!
! Redo that integration time step
        goto 150
      endif
!
!------------------------------------------------------------------------------
!
!  DATA  DUMP  AND  PROGRESS  REPORT
!
! AVIMANDELL
!  Added voluntary data dump when file 'confirmdumpfile' exists
!
! Convert to heliocentric coords and do the data dump
      inquire (file='confirmdumpfile', exist=test)

      if (test.and. &
       (abs(time-tdump).ge.abs(dtdump).and.opflag.ge.-1)) then
        call bcoord (time,jcen,nbod,nbig,h0,m,x,v,xh,vh,ngf,ngflag,opt)
        do j = 2, nbod
          epoch(j) = time
        end do
        call mio_ce (time,tstart,rcen,rmax,nbod,nbig,m,stat,id, &
         0,iclo,jclo,opt,stopflag,tclo,dclo,ixvclo,jxvclo,mem,lmem, &
         outfile,nstored,0)
        call mio_dump (time,tstart,tstop,dtout,algor,h0,tol,jcen,rcen, &
         rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s,rho,rceh,stat, &
         id,ngf,epoch,opt,opflag,dumpfile,mem,lmem)
        tdump = time
      endif
!
! Convert to heliocentric coords and write a progress report to the log file
      if (abs(time-tlog).ge.abs(dtdump).and.opflag.ge.0) then
        call bcoord (time,jcen,nbod,nbig,h0,m,x,v,xh,vh,ngf,ngflag,opt)
        call mxx_en (jcen,nbod,nbig,m,xh,vh,s,en(2),am(2))
        call mio_log (time,tstart,en,am,opt,mem,lmem)
        tlog = time
      endif
!
!------------------------------------------------------------------------------
!
!  CHECK  FOR  EJECTIONS  AND  DO  OTHER  PERIODIC  EFFECTS
!
      if (abs(time-tfun).ge.abs(dtfun).and.opflag.ge.-1) then
        if (algor.eq.1) then
          call mco_iden (jcen,nbod,nbig,h0,m,x,v,xh,vh)
        else
          call bcoord(time,jcen,nbod,nbig,h0,m,x,v,xh,vh,ngf,ngflag,opt)
        endif
!
! Recompute close encounter limits, to allow for changes in Hill radii
        call mce_hill (nbod,m,xh,vh,rce,a)
        do j = 2, nbod
          rce(j) = rce(j) * rceh(j)
        end do
!
! Check for ejections
        itmp = 2
        if (algor.eq.11.or.algor.eq.12) itmp = 3
        call mxx_ejec (time,tstart,rmax,en,am,jcen,itmp,nbod,nbig,m,xh, &
         vh,s,stat,id,opt,ejflag,outfile(3),mem,lmem)
!
! Remove ejected objects, reset flags, calculate new Hill and physical radii
        if (ejflag.ne.0) then
          call mxx_elim (nbod,nbig,m,xh,vh,s,rho,rceh,rcrit,ngf,stat, &
           id,mem,lmem,outfile(3),itmp)
          if (opflag.ge.0) opflag = 1
          dtflag = 1
          call mce_init (tstart,algor,h0,jcen,rcen,rmax,cefac,nbod,nbig, &
           m,xh,vh,s,rho,rceh,rphys,rce,rcrit,id,opt,outfile(2),0)
          if (algor.eq.1) then
            call mco_iden (jcen,nbod,nbig,h0,m,xh,vh,x,v)
          else
            call coord (time,jcen,nbod,nbig,h0,m,xh,vh,x,v,ngf,ngflag, &
             opt)
          endif
        endif
        tfun = time
      endif
!
! Go on to the next time step
      goto 100
!
!------------------------------------------------------------------------------
!
      end
!

