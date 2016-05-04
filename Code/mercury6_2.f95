!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MERCURY6_1.FOR    (ErikSoft   3 May 2002)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Mercury is a general-purpose N-body integration package for problems in
! celestial mechanics.
!
!------------------------------------------------------------------------------
! This package contains some subroutines taken from the Swift integration 
! package by H.F.Levison and M.J.Duncan (1994) Icarus, vol 108, pp18.
! Routines taken from Swift have names beginning with `drift' or `orbel'.
!
! The standard symplectic (MVS) algorithm is described in J.Widsom and
! M.Holman (1991) Astronomical Journal, vol 102, pp1528.
!
! The hybrid symplectic algorithm is described in J.E.Chambers (1999)
! Monthly Notices of the RAS, vol 304, pp793.
!
! RADAU is described in E.Everhart (1985) in ``The Dynamics of Comets:
! Their Origin and Evolution'' p185-202, eds. A.Carusi & G.B.Valsecchi,
! pub. Reidel.
!
! The Bulirsch-Stoer algorithms are described in W.H.Press et al. (1992)
! ``Numerical Recipes in Fortran'', pub. Cambridge.
!------------------------------------------------------------------------------
!
! #rjw# changes for f95 format:
! All tabs replaced with 4 spaces
! 'c' for comments replaced with '!'
! real*8 -> real(8), character*X -> character(X)
! mercury.inc -> mercury95.inc (and swift.inc. Switch comment chars in those)
! Updated line continuation symbols: '&' at end of unfinished line instead
! of or in addition to a character in 6th column of continuing line
!
!------------------------------------------------------------------------------
! User module by Dimitri Veras and R.J. Worth
!
! Adds Galactic tides and a dissipating gas cloud potential
! to the force calculations. Changes marked with either DGV or #rjw#,
! depending on author.

!------------------------------------------------------------------------------
! Performance improvement by Avi Mandell
!
! Some of the changes from Avi Mandell's version are included here, marked with 
! AVIMANDELL, including commenting out the output to the ce.out file because
! it's large and not useful to me, and commenting out writing output at every
! step, which causes it to run many times slower. 
! His other modifications are omitted. (He also wrote in a runtime printout, 
! explicit file closing, and forced migration.)

!------------------------------------------------------------------------------
! Bug fix courtesy of Alex Mustill
!
! The criteria for determining whether a body has collided with the central
! object worked poorly for moons/satellites. This adds the requirement that
! in order for a collision to be scored, the collision must occur during the
!  current timestep.
!
! The alteration is in the mce_cent subroutine and is marked with ##ajm##.

!------------------------------------------------------------------------------
! Bug fix found online
! Author: Karla de Souza Torres
!
! The variable array STAT is now explicitly initialized each time the 
! integrator (re-)starts. This prevents STAT from getting random values, 
! independently of compiler. The STAT variable is responsible for flagging 
! bodies to be deleted after events, and its non-initialization caused bodies
! to disappear in a non-physical manner in some cases.
!
! The alteration is made in the subroutine mio_in and is labelled with ##K##.
!
! Please, send your comments to karlchen79@gmail.com
!------------------------------------------------------------------------------
!
! Variables:
! ---------
!  M      = mass (in solar masses)
!  XH     = coordinates (x,y,z) with respect to the central body (AU)
!  VH     = velocities (vx,vy,vz) with respect to the central body (AU/day)
!  S      = spin angular momentum (solar masses AU^2/day)
!  RHO    = physical density (g/cm^3)
!  RCEH   = close-encounter limit (Hill radii)
!  STAT   = status (0 => alive, <>0 => to be removed)
!  ID     = name of the object (8 characters)
!  CE     = close encounter status
!  NGF    = (1-3) cometary non-gravitational (jet) force parameters
!   "     =  (4)  beta parameter for radiation pressure and P-R drag
!  EPOCH  = epoch of orbit (days)
!  NBOD  = current number of bodies (INCLUDING the central object)
!  NBIG  =    "       "    " big bodies (ones that perturb everything else)
!  TIME  = current epoch (days)
!  TOUT  = time of next output evaluation
!  TDUMP = time of next data dump
!  TFUN  = time of next periodic effect (e.g. next check for ejections)
!  H     = current integration timestep (days)
!  EN(1) = initial energy of the system
!  " (2) = current    "    "  "    "
!  " (3) = energy change due to collisions, ejections etc.
!  AM(1,2,3) = as above but for angular momentum
!
! Integration Parameters :
! ----------------------
!  ALGOR = 1  ->  Mixed-variable symplectic
!          2  ->  Bulirsch-Stoer integrator
!          3  ->         "           "      (conservative systems only)
!          4  ->  RA15 `radau' integrator
!          10 ->  Hybrid MVS/BS (democratic-heliocentric coords)
!          11 ->  Close-binary hybrid (close-binary coords)
!          12 ->  Wide-binary hybrid (wide-binary coords)
!
! TSTART = epoch of first required output (days)
! TSTOP  =   "      final required output ( "  )
! DTOUT  = data output interval           ( "  )
! DTDUMP = data-dump interval             ( "  )
! DTFUN  = interval for other periodic effects (e.g. check for ejections)
!  H0    = initial integration timestep (days)
!  TOL   = Integrator tolerance parameter (approx. error per timestep)
!  RMAX  = heliocentric distance at which objects are considered ejected (AU)
!  RCEN  = radius of central body (AU)
!  JCEN(1,2,3) = J2,J4,J6 for central body (units of RCEN^i for Ji)
!
! Options:
!  OPT(1) = close-encounter option (0=stop after an encounter, 1=continue)
!  OPT(2) = collision option (0=no collisions, 1=merge, 2=merge+fragment)
!  OPT(3) = time style (0=days 1=Greg.date 2/3=days/years w/respect to start)
!  OPT(4) = o/p precision (1,2,3 = 4,9,15 significant figures)
!  OPT(5) = < Not used at present >
!  OPT(6) = < Not used at present >
!  OPT(7) = apply post-Newtonian correction? (0=no, 1=yes)
!  OPT(8) = apply user-defined force routine mfo_user? (0=no, 1=yes)
!
! File variables :
! --------------
!  OUTFILE  (1) = osculating coordinates/velocities and masses
!     "     (2) = close encounter details
!     "     (3) = information file
!  DUMPFILE (1) = Big-body data
!     "     (2) = Small-body data
!     "     (3) = integration parameters
!     "     (4) = restart file
!
! Flags :
! -----
!  NGFLAG = do any bodies experience non-grav. forces?
!                            ( 0 = no non-grav forces)
!                              1 = cometary jets only
!                              2 = radiation pressure/P-R drag only
!                              3 = both
!  OPFLAG = integration mode (-2 = synchronising epochs)
!                             -1 = integrating towards start epoch
!                              0 = main integration, normal output
!                              1 = main integration, full output
!
!------------------------------------------------------------------------------
!
      implicit none
      include 'mercury95.inc'
!
      integer j,algor,nbod,nbig,opt(8),stat(NMAX),lmem(NMESS)
      integer opflag,ngflag,ndump,nfun
      real(8) m(NMAX),xh(3,NMAX),vh(3,NMAX),s(3,NMAX),rho(NMAX)
      real(8) rceh(NMAX),epoch(NMAX),ngf(4,NMAX),rmax,rcen,jcen(3)
      real(8) cefac,time,tstart,tstop,dtout,h0,tol,en(3),am(3)
      character(8) id(NMAX)
      character(80) outfile(3), dumpfile(4), mem(NMESS)
      external mdt_mvs, mdt_bs1, mdt_bs2, mdt_ra15, mdt_hy
      external mco_dh2h,mco_h2dh
      external mco_b2h,mco_h2b,mco_h2mvs,mco_mvs2h,mco_iden
!
      data opt/0,1,1,2,0,1,0,0/
!
!------------------------------------------------------------------------------
!
! Get initial conditions and integration parameters
      call mio_in (time,tstart,tstop,dtout,algor,h0,tol,rmax,rcen,jcen, &
       en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s,rho,rceh,stat,id, &
       epoch,ngf,opt,opflag,ngflag,outfile,dumpfile,lmem,mem)
!
! If this is a new integration, integrate all the objects to a common epoch.
      if (opflag.eq.-2) then
  20    open (23,file=outfile(3),status='old',access='append',err=20)
        write (23,'(/,a)') mem(55)(1:lmem(55))
        write (*,'(a)') mem(55)(1:lmem(55))
        call mxx_sync (time,tstart,h0,tol,jcen,nbod,nbig,m,xh,vh,s,rho, &
         rceh,stat,id,epoch,ngf,opt,ngflag)
        write (23,'(/,a,/)') mem(56)(1:lmem(56))
        write (*,'(a)') mem(56)(1:lmem(56))
        opflag = -1
        close (23)
      end if
!
! Main integration
      if (algor.eq.1) call mal_hcon (time,tstart,tstop,dtout,algor,h0, &
       tol,jcen,rcen,rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s, &
       rho,rceh,stat,id,ngf,opt,opflag,ngflag,outfile,dumpfile,mem, &
       lmem,mdt_mvs,mco_h2mvs,mco_mvs2h)
!
      if (algor.eq.9) call mal_hcon (time,tstart,tstop,dtout,algor,h0, &
       tol,jcen,rcen,rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s, &
       rho,rceh,stat,id,ngf,opt,opflag,ngflag,outfile,dumpfile,mem, &
       lmem,mdt_mvs,mco_iden,mco_iden)
!
      if (algor.eq.2) call mal_hvar (time,tstart,tstop,dtout,algor,h0, &
       tol,jcen,rcen,rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s, &
       rho,rceh,stat,id,ngf,opt,opflag,ngflag,outfile,dumpfile,mem, &
       lmem,mdt_bs1)
!
      if (algor.eq.3) call mal_hvar (time,tstart,tstop,dtout,algor,h0, &
       tol,jcen,rcen,rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s, &
       rho,rceh,stat,id,ngf,opt,opflag,ngflag,outfile,dumpfile,mem, &
       lmem,mdt_bs2)
!
      if (algor.eq.4) call mal_hvar (time,tstart,tstop,dtout,algor,h0, &
       tol,jcen,rcen,rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s, &
       rho,rceh,stat,id,ngf,opt,opflag,ngflag,outfile,dumpfile,mem, &
       lmem,mdt_ra15)
!
      if (algor.eq.10) call mal_hcon (time,tstart,tstop,dtout,algor,h0, &
       tol,jcen,rcen,rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s, &
       rho,rceh,stat,id,ngf,opt,opflag,ngflag,outfile,dumpfile,mem, &
       lmem,mdt_hy,mco_h2dh,mco_dh2h)
!
! Do a final data dump
      do j = 2, nbod
        epoch(j) = time
      end do
      call mio_dump (time,tstart,tstop,dtout,algor,h0,tol,jcen,rcen, &
       rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,xh,vh,s,rho,rceh,stat, &
       id,ngf,epoch,opt,opflag,dumpfile,mem,lmem)
!
! Calculate and record the overall change in energy and ang. momentum
  50  open  (23, file=outfile(3), status='old', access='append', &
       err=50)
      write (23,'(/,a)') mem(57)(1:lmem(57))
      call mxx_en (jcen,nbod,nbig,m,xh,vh,s,en(2),am(2))
!
      write (23,231) mem(58)(1:lmem(58)),  &
       abs((en(2) + en(3) - en(1)) / en(1))
      write (23,232) mem(59)(1:lmem(59)),  &
       abs((am(2) + am(3) - am(1)) / am(1))
      write (23,231) mem(60)(1:lmem(60)), abs(en(3) / en(1))
      write (23,232) mem(61)(1:lmem(61)), abs(am(3) / am(1))
      close (23)
      write (*,'(a)') mem(57)(1:lmem(57))
!
!------------------------------------------------------------------------------
!
 231  format (/,a,1p1e12.5)
 232  format (a,1p1e12.5)
      stop
      end
!

