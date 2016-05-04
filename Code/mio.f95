!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MIO_CE.FOR    (ErikSoft   1 March 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Writes details of close encounter minima to an output file, and decides how
! to continue the integration depending upon the close-encounter option
! chosen by the user. Close encounter details are stored until either 100
! have been accumulated, or a data dump is done, at which point the stored
! encounter details are also output.
!
! For each encounter, the routine outputs the time and distance of closest
! approach, the identities of the objects involved, and the output
! variables of the objects at this time. The output variables are:
! expressed as
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
      subroutine mio_ce (time,tstart,rcen,rmax,nbod,nbig,m,stat,id, &
       nclo,iclo,jclo,opt,stopflag,tclo,dclo,ixvclo,jxvclo,mem, &
       lmem,outfile,nstored,ceflush)
!
      implicit none
      include 'mercury95.inc'
!
! Input/Output
      integer nbod,nbig,opt(8),stat(nbod),lmem(NMESS),stopflag
      integer nclo,iclo(nclo),jclo(nclo),nstored,ceflush
      real(8) time,tstart,rcen,rmax,m(nbod),tclo(nclo),dclo(nclo)
      real(8) ixvclo(6,nclo),jxvclo(6,nclo)
      character(80) outfile(3),mem(NMESS)
      character(8) id(nbod)
!
! Local
      integer k,year,month
      real(8) tmp0,t1,rfac,fr,fv,theta,phi,vtheta,vphi
      character(80) c(200)
      character(38) fstop
      character(8) mio_fl2c, mio_re2c
      character(6) tstring
!
!------------------------------------------------------------------------------
!
      save c
!
! Scaling factor (maximum possible range) for distances
      rfac = log10 (rmax / rcen)
!
! Store details of each new close-encounter minimum
      do k = 1, nclo
        nstored = nstored + 1
        c(nstored)(1:8)   = mio_fl2c(tclo(k))
        c(nstored)(9:16)  = mio_re2c(dble(iclo(k)-1),0.d0,11239423.99d0)
        c(nstored)(12:19) = mio_re2c(dble(jclo(k)-1),0.d0,11239423.99d0)
        c(nstored)(15:22) = mio_fl2c(dclo(k))
!
        call mco_x2ov (rcen,rmax,m(1),0.d0,ixvclo(1,k),ixvclo(2,k), &
         ixvclo(3,k),ixvclo(4,k),ixvclo(5,k),ixvclo(6,k),fr,theta,phi, &
         fv,vtheta,vphi)
        c(nstored)(23:30) = mio_re2c (fr    , 0.d0, rfac)
        c(nstored)(27:34) = mio_re2c (theta , 0.d0, PI)
        c(nstored)(31:38) = mio_re2c (phi   , 0.d0, TWOPI)
        c(nstored)(35:42) = mio_re2c (fv    , 0.d0, 1.d0)
        c(nstored)(39:46) = mio_re2c (vtheta, 0.d0, PI)
        c(nstored)(43:50) = mio_re2c (vphi  , 0.d0, TWOPI)
!
        call mco_x2ov (rcen,rmax,m(1),0.d0,jxvclo(1,k),jxvclo(2,k), &
         jxvclo(3,k),jxvclo(4,k),jxvclo(5,k),jxvclo(6,k),fr,theta,phi, &
         fv,vtheta,vphi)
        c(nstored)(47:54) = mio_re2c (fr    , 0.d0, rfac)
        c(nstored)(51:58) = mio_re2c (theta , 0.d0, PI)
        c(nstored)(55:62) = mio_re2c (phi   , 0.d0, TWOPI)
        c(nstored)(59:66) = mio_re2c (fv    , 0.d0, 1.d0)
        c(nstored)(63:74) = mio_re2c (vtheta, 0.d0, PI)
        c(nstored)(67:78) = mio_re2c (vphi  , 0.d0, TWOPI)
      end do
!
! If required, output the stored close encounter details
! AVIMANDELL
!  *Commenting close encounter output for speed
!
      if (nstored.ge.100.or.ceflush.eq.0) then
!  10    open (22, file=outfile(2), status='old', access='append',err=10)
!        do k = 1, nstored
!          write (22,'(a1,a2,a70)') char(12),'6b',c(k)(1:70)
!        enddo
!        close (22)
        nstored = 0
      endif
!
! If new encounter minima have occurred, decide whether to stop integration
      stopflag = 0
      if (opt(1).eq.1.and.nclo.gt.0) then
  20    open (23, file=outfile(3), status='old', access='append',err=20)
! If time style is Gregorian date then...
        tmp0 = tclo(1)
        if (opt(3).eq.1) then
          fstop = '(5a,/,9x,a,i10,1x,i2,1x,f4.1)'
          call mio_jd2y (tmp0,year,month,t1)
          write (23,fstop) mem(121)(1:lmem(121)),mem(126) &
           (1:lmem(126)),id(iclo(1)),',',id(jclo(1)), &
           mem(71)(1:lmem(71)),year,month,t1
! Otherwise...
        else
          if (opt(3).eq.3) then
            tstring = mem(2)
            fstop = '(5a,/,9x,a,f14.3,a)'
            t1 = (tmp0 - tstart) / 365.25d0
          else
            tstring = mem(1)
            fstop = '(5a,/,9x,a,f14.1,a)'
            if (opt(3).eq.0) t1 = tmp0
            if (opt(3).eq.2) t1 = tmp0 - tstart
          end if
          write (23,fstop) mem(121)(1:lmem(121)),mem(126) &
           (1:lmem(126)),id(iclo(1)),',',id(jclo(1)), &
           mem(71)(1:lmem(71)),t1,tstring
        end if
        stopflag = 1
        close(23)
      end if
!
!------------------------------------------------------------------------------
!
      return
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MIO_DUMP.FOR    (ErikSoft   21 February 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Writes masses, coordinates, velocities etc. of all objects, and integration
! parameters, to dump files. Also updates a restart file containing other
! variables used internally by MERCURY.
!
!------------------------------------------------------------------------------
!
      subroutine mio_dump (time,tstart,tstop,dtout,algor,h0,tol,jcen, &
       rcen,rmax,en,am,cefac,ndump,nfun,nbod,nbig,m,x,v,s,rho,rceh, &
       stat,id,ngf,epoch,opt,opflag,dumpfile,mem,lmem)
!
      implicit none
      include 'mercury95.inc'
!
! Input/Output
      integer algor,nbod,nbig,stat(nbod),opt(8),opflag,ndump,nfun
      integer lmem(NMESS)
      real(8) time,tstart,tstop,dtout,h0,tol,rmax,en(3),am(3)
      real(8) jcen(3),rcen,cefac,m(nbod),x(3,nbod),v(3,nbod)
      real(8) s(3,nbod),rho(nbod),rceh(nbod),ngf(4,nbod),epoch(nbod)
      character(80) dumpfile(4),mem(NMESS)
      character(8) id(nbod)
!
! Local
      integer idp,i,j,k,len,j1,j2
      real(8) rhocgs,k_2,rcen_2,rcen_4,rcen_6,x0(3,NMAX),v0(3,NMAX)
      character(150) c
!
!------------------------------------------------------------------------------
!
      rhocgs = AU * AU * AU * K2 / MSUN
      k_2 = 1.d0 / K2
      rcen_2 = 1.d0 / (rcen * rcen)
      rcen_4 = rcen_2 * rcen_2
      rcen_6 = rcen_4 * rcen_2
!
! If using close-binary star, convert to user coordinates
!      if (algor.eq.11) call mco_h2ub (time,jcen,nbod,nbig,h0,m,x,v,
!     %   x0,v0)
!
! Dump to temporary files (idp=1) and real dump files (idp=2)
      do idp = 1, 2
!
! Dump data for the Big (i=1) and Small (i=2) bodies
        do i = 1, 2
          if (idp.eq.1) then
            if (i.eq.1) c(1:12) = 'big.tmp     '
            if (i.eq.2) c(1:12) = 'small.tmp   '
  20        open (31, file=c(1:12), status='unknown', err=20)
          else
  25        open (31, file=dumpfile(i), status='old', err=25)
          end if
!
! Write header lines, data style (and epoch for Big bodies)
          write (31,'(a)') mem(151+i)(1:lmem(151+i))
          if (i.eq.1) then
            j1 = 2
            j2 = nbig
          else
            j1 = nbig + 1
            j2 = nbod
          end if
          write (31,'(a)') mem(154)(1:lmem(154))
          write (31,'(a)') mem(155)(1:lmem(155))
          write (31,*) mem(156)(1:lmem(156)),'Cartesian'
          if (i.eq.1) write (31,*) mem(157)(1:lmem(157)),time
          write (31,'(a)') mem(155)(1:lmem(155))
!
! For each body...
          do j = j1, j2
            len = 37
            c(1:8) = id(j)
            write (c(9:37),'(1p,a3,e11.5,a3,e11.5)') ' r=',rceh(j), &
             ' d=',rho(j)/rhocgs
            if (m(j).gt.0) then
              write (c(len+1:len+25),'(a3,e22.15)') ' m=',m(j)*k_2
              len = len + 25
            end if
            do k = 1, 3
              if (ngf(k,j).ne.0) then
                write (c(len+1:len+16),'(a2,i1,a1,e12.5)') ' a',k,'=', &
                 ngf(k,j)
                len = len + 16
              end if
            end do
            if (ngf(4,j).ne.0) then
              write (c(len+1:len+15),'(a3,e12.5)') ' b=',ngf(4,j)
              len = len + 15
            end if
            write (31,'(a)') c(1:len)
            if (algor.eq.11) then
              write (31,312) x0(1,j), x0(2,j), x0(3,j)
              write (31,312) v0(1,j), v0(2,j), v0(3,j)
            else
              write (31,312) x(1,j), x(2,j), x(3,j)
              write (31,312) v(1,j), v(2,j), v(3,j)
            end if
            write (31,312) s(1,j)*k_2, s(2,j)*k_2, s(3,j)*k_2
          enddo
          close (31)
        end do
!
! Dump the integration parameters
  40    if (idp.eq.1) open (33,file='param.tmp',status='unknown',err=40)
  45    if (idp.eq.2) open (33, file=dumpfile(3), status='old', err=45)
!
! Important parameters
        write (33,'(a)') mem(151)(1:lmem(151))
        write (33,'(a)') mem(154)(1:lmem(154))
        write (33,'(a)') mem(155)(1:lmem(155))
        write (33,'(a)') mem(158)(1:lmem(158))
        write (33,'(a)') mem(155)(1:lmem(155))
        if (algor.eq.1) then
          write (33,*) mem(159)(1:lmem(159)),'MVS'
        else if (algor.eq.2) then
          write (33,*) mem(159)(1:lmem(159)),'BS'
        else if (algor.eq.3) then
          write (33,*) mem(159)(1:lmem(159)),'BS2'
        else if (algor.eq.4) then
          write (33,*) mem(159)(1:lmem(159)),'RADAU'
        else if (algor.eq.10) then
          write (33,*) mem(159)(1:lmem(159)),'HYBRID'
        else if (algor.eq.11) then
          write (33,*) mem(159)(1:lmem(159)),'CLOSE'
        else if (algor.eq.12) then
          write (33,*) mem(159)(1:lmem(159)),'WIDE'
        else
          write (33,*) mem(159)(1:lmem(159)),'0'
        end if
        write (33,*) mem(160)(1:lmem(160)),tstart
        write (33,*) mem(161)(1:lmem(161)),tstop
        write (33,*) mem(162)(1:lmem(162)),dtout
        write (33,*) mem(163)(1:lmem(163)),h0
        write (33,*) mem(164)(1:lmem(164)),tol
!
! Integration options
        write (33,'(a)') mem(155)(1:lmem(155))
        write (33,'(a)') mem(165)(1:lmem(165))
        write (33,'(a)') mem(155)(1:lmem(155))
        if (opt(1).eq.0) then
          write (33,'(2a)') mem(166)(1:lmem(166)),mem(5)(1:lmem(5))
        else
          write (33,'(2a)') mem(166)(1:lmem(166)),mem(6)(1:lmem(6))
        end if
        if (opt(2).eq.0) then
          write (33,'(2a)') mem(167)(1:lmem(167)),mem(5)(1:lmem(5))
          write (33,'(2a)') mem(168)(1:lmem(168)),mem(5)(1:lmem(5))
        else if (opt(2).eq.2) then
          write (33,'(2a)') mem(167)(1:lmem(167)),mem(6)(1:lmem(6))
          write (33,'(2a)') mem(168)(1:lmem(168)),mem(6)(1:lmem(6))
        else
          write (33,'(2a)') mem(167)(1:lmem(167)),mem(6)(1:lmem(6))
          write (33,'(2a)') mem(168)(1:lmem(168)),mem(5)(1:lmem(5))
        end if
        if (opt(3).eq.0.or.opt(3).eq.2) then
          write (33,'(2a)') mem(169)(1:lmem(169)),mem(1)(1:lmem(1))
        else
          write (33,'(2a)') mem(169)(1:lmem(169)),mem(2)(1:lmem(2))
        end if
        if (opt(3).eq.2.or.opt(3).eq.3) then
          write (33,'(2a)') mem(170)(1:lmem(170)),mem(6)(1:lmem(6))
        else
          write (33,'(2a)') mem(170)(1:lmem(170)),mem(5)(1:lmem(5))
        end if
        if (opt(4).eq.1) then
          write (33,'(2a)') mem(171)(1:lmem(171)),mem(7)(1:lmem(7))
        else if (opt(4).eq.3) then
          write (33,'(2a)') mem(171)(1:lmem(171)),mem(9)(1:lmem(9))
        else
          write (33,'(2a)') mem(171)(1:lmem(171)),mem(8)(1:lmem(8))
        end if
        write (33,'(a)') mem(172)(1:lmem(172))
        if (opt(7).eq.1) then
          write (33,'(2a)') mem(173)(1:lmem(173)),mem(6)(1:lmem(6))
        else
          write (33,'(2a)') mem(173)(1:lmem(173)),mem(5)(1:lmem(5))
        end if
        if (opt(8).eq.1) then
          write (33,'(2a)') mem(174)(1:lmem(174)),mem(6)(1:lmem(6))
        else
          write (33,'(2a)') mem(174)(1:lmem(174)),mem(5)(1:lmem(5))
        end if
!
! Infrequently-changed parameters
        write (33,'(a)') mem(155)(1:lmem(155))
        write (33,'(a)') mem(175)(1:lmem(175))
        write (33,'(a)') mem(155)(1:lmem(155))
        write (33,*) mem(176)(1:lmem(176)),rmax
        write (33,*) mem(177)(1:lmem(177)),rcen
        write (33,*) mem(178)(1:lmem(178)),m(1) * k_2
        write (33,*) mem(179)(1:lmem(179)),jcen(1) * rcen_2
        write (33,*) mem(180)(1:lmem(180)),jcen(2) * rcen_4
        write (33,*) mem(181)(1:lmem(181)),jcen(3) * rcen_6
        write (33,*) mem(182)(1:lmem(182))
        write (33,*) mem(183)(1:lmem(183))
        write (33,*) mem(184)(1:lmem(184)),cefac
        write (33,*) mem(185)(1:lmem(185)),ndump
        write (33,*) mem(186)(1:lmem(186)),nfun
        close (33)
!
! Create new version of the restart file
  60    if (idp.eq.1) open (35, file='restart.tmp', status='unknown', &
         err=60)
  65    if (idp.eq.2) open (35, file=dumpfile(4), status='old', err=65)
        write (35,'(1x,i2)') opflag
        write (35,*) en(1) * k_2
        write (35,*) am(1) * k_2
        write (35,*) en(3) * k_2
        write (35,*) am(3) * k_2
        write (35,*) s(1,1) * k_2
        write (35,*) s(2,1) * k_2
        write (35,*) s(3,1) * k_2
        close (35)
      end do
!
!------------------------------------------------------------------------------
!
 311  format (1x,a8,1x,a1,1p,e22.15,2(1x,e11.5))
 312  format (1p,3(1x,e22.15),1x,i8)
 313  format (1p,1x,e22.15,0p,2x,a)
 314  format (1x,a8,1x,a1,1p,e22.15,4(1x,e12.5),1x,e22.15,2(1x,e11.5))
      return
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MIO_FL2C.FOR    (ErikSoft  1 July 1998)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Converts a (floating point) real(8) variable X, into a character(8) ASCII 
! string, using the new format compression:
!
! X is first converted to base 224, and then each base 224 digit is converted 
! to an ASCII character, such that 0 -> character 32, 1 -> character 33...
! and 223 -> character 255.
! The first 7 characters in the string are used to store the mantissa, and the
! eighth character is used for the exponent.
!
! ASCII characters 0 - 31 (CTRL characters) are not used, because they
! cause problems when using some operating systems.
!
! N.B. X must lie in the range -1.e112 < X < 1.e112
! ===
!
!------------------------------------------------------------------------------
!
      function mio_fl2c (x)
!
      implicit none
!
! Input/Output
      real(8) x
      character(8) mio_fl2c
!
! Local
      integer ex
      real(8) ax,y
      character(8) mio_re2c
!
!------------------------------------------------------------------------------
!
      if (x.eq.0) then
        y = .5d0
      else
        ax = abs(x)
        ex = int(log10(ax))
        if (ax.ge.1) ex = ex + 1
        y = ax*(10.d0**(-ex))
        if (y.eq.1) then
          y = y * .1d0
          ex = ex + 1
        end if
        y = sign(y,x) *.5d0 + .5d0
      end if
!
      mio_fl2c(1:8) = mio_re2c (y, 0.d0, 1.d0)
      ex = ex + 112
      if (ex.gt.223) ex = 223
      if (ex.lt.0) ex = 0
      mio_fl2c(8:8) = char(ex+32)
!
!------------------------------------------------------------------------------
!
      return
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MIO_IN.FOR    (ErikSoft   4 May 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Reads names, masses, coordinates and velocities of all the bodies,
! and integration parameters for the MERCURY integrator package. 
! If DUMPFILE(4) exists, the routine assumes this is a continuation of
! an old integration, and reads all the data from the dump files instead
! of the input files.
!
! N.B. All coordinates are with respect to the central body!!
! ===
!
!------------------------------------------------------------------------------
!
      subroutine mio_in (time,tstart,tstop,dtout,algor,h0,tol,rmax,rcen, &
       jcen,en,am,cefac,ndump,nfun,nbod,nbig,m,x,v,s,rho,rceh,stat,id, &
       epoch,ngf,opt,opflag,ngflag,outfile,dumpfile,lmem,mem)
!
      implicit none
      include 'mercury95.inc'
!
! Input/Output
      integer algor,nbod,nbig,stat(NMAX),opt(8),opflag,ngflag
      integer lmem(NMESS),ndump,nfun
      real(8) time,tstart,tstop,dtout,h0,tol,rmax,rcen,jcen(3)
      real(8) en(3),am(3),m(NMAX),x(3,NMAX),v(3,NMAX),s(3,NMAX)
      real(8) rho(NMAX),rceh(NMAX),epoch(NMAX),ngf(4,NMAX),cefac
      character(80) outfile(3),dumpfile(4), mem(NMESS)
      character(8) id(NMAX)
!
! Local
      integer j,k,itmp,jtmp,informat,lim(2,10),nsub,year,month,lineno
      real(8) q,a,e,i,p,n,l,temp,tmp2,tmp3,rhocgs,t1,tmp4,tmp5,tmp6
!      real(8) v0(3,NMAX),x0(3,NMAX)
      logical test,oldflag,flag1,flag2
      character(1) c1
      character(3) c3,alg(60)
      character(80) infile(3),filename,c80
      character(150) string
!
!------------------------------------------------------------------------------
!
      data alg/'MVS','Mvs','mvs','mvs','mvs','BS ','Bs ','bs ','Bul', &
        'bul','BS2','Bs2','bs2','Bu2','bu2','RAD','Rad','rad','RA ', &
        'ra ','xxx','xxx','xxx','xxx','xxx','xxx','xxx','xxx','xxx', &
        'xxx','xxx','xxx','xxx','xxx','xxx','xxx','xxx','xxx','xxx', &
        'xxx','TES','Tes','tes','Tst','tst','HYB','Hyb','hyb','HY ', &
        'hy ','CLO','Clo','clo','CB ','cb ','WID','Wid','wid','WB ', &
        'wb '/
!
      rhocgs = AU * AU * AU * K2 / MSUN
      do j = 1, 80
        filename(j:j) = ' '
      end do
      do j = 1, 3
        infile(j)   = filename
        outfile(j)  = filename
        dumpfile(j) = filename
      end do
      dumpfile(4) = filename
!
! Read in output messages
      inquire (file='In/message.in', exist=test)
      if (.not.test) then
        write (*,'(/,2a)') ' ERROR: This file is needed to start', &
         ' the integration:  In/message.in'
        stop
      end if
      open (16, file='In/message.in', status='old')
  10  read (16,'(i3,1x,i2,1x,a80)',end=20) j,lmem(j),mem(j)
      goto 10
  20  close (16)
!
! Read in filenames and check for duplicate filenames
      inquire (file='In/files.in', exist=test)
      if (.not.test) call mio_err (6,mem(81),lmem(81),mem(88),lmem(88), &
       ' ',1,'In/files.in',8)
      open (15, file='In/files.in', status='old')
!
! Input files
      do j = 1, 3
        read (15,'(a150)') string
        call mio_spl (150,string,nsub,lim)
        infile(j)(1:(lim(2,1)-lim(1,1)+1)) = string(lim(1,1):lim(2,1))
        do k = 1, j - 1
          if (infile(j).eq.infile(k)) call mio_err (6,mem(81),lmem(81), &
           mem(89),lmem(89),infile(j),80,mem(86),lmem(86))
        end do
      end do
!
! Output files
      do j = 1, 3
        read (15,'(a150)') string
        call mio_spl (150,string,nsub,lim)
        outfile(j)(1:(lim(2,1)-lim(1,1)+1)) = string(lim(1,1):lim(2,1))
        do k = 1, j - 1
          if (outfile(j).eq.outfile(k)) call mio_err (6,mem(81), &
           lmem(81),mem(89),lmem(89),outfile(j),80,mem(86),lmem(86))
        end do
        do k = 1, 3
          if (outfile(j).eq.infile(k)) call mio_err (6,mem(81),lmem(81), &
           mem(89),lmem(89),outfile(j),80,mem(86),lmem(86))
        end do
      end do
!
! Dump files
      do j = 1, 4
        read (15,'(a150)') string
        call mio_spl (150,string,nsub,lim)
        dumpfile(j)(1:(lim(2,1)-lim(1,1)+1)) = string(lim(1,1):lim(2,1))
        do k = 1, j - 1
          if (dumpfile(j).eq.dumpfile(k)) call mio_err (6,mem(81), &
           lmem(81),mem(89),lmem(89),dumpfile(j),80,mem(86),lmem(86))
        end do
        do k = 1, 3
          if (dumpfile(j).eq.infile(k)) call mio_err (6,mem(81), &
           lmem(81),mem(89),lmem(89),dumpfile(j),80,mem(86),lmem(86))
        end do
        do k = 1, 3
          if (dumpfile(j).eq.outfile(k)) call mio_err (6,mem(81), &
           lmem(81),mem(89),lmem(89),dumpfile(j),80,mem(86),lmem(86))
        end do
      end do
      close (15)
!
! Find out if this is an old integration (i.e. does the restart file exist)
      inquire (file=dumpfile(4), exist=oldflag)
!
! Check if information file exists, and append a continuation message
      if (oldflag) then
        inquire (file=outfile(3), exist=test)
        if (.not.test) call mio_err (6,mem(81),lmem(81),mem(88), &
         lmem(88),' ',1,outfile(3),80)
 320    open(23,file=outfile(3),status='old',access='append',err=320)
      else
!
! If new integration, check information file doesn't exist, and then create it
        inquire (file=outfile(3), exist=test)
        if (test) call mio_err (6,mem(81),lmem(81),mem(87),lmem(87), &
         ' ',1,outfile(3),80)
 410    open(23, file = outfile(3), status = 'new', err=410)
      end if
!
!------------------------------------------------------------------------------
!
!  READ  IN  INTEGRATION  PARAMETERS
!
! Check if the file containing integration parameters exists, and open it
      filename = infile(3)
      if (oldflag) filename = dumpfile(3)
      inquire (file=filename, exist=test)
      if (.not.test) call mio_err (23,mem(81),lmem(81),mem(88),lmem(88), &
       ' ',1,filename,80)
  30  open  (13, file=filename, status='old', err=30)
!
! Read integration parameters
      lineno = 0
      do j = 1, 26
  40    lineno = lineno + 1
        read (13,'(a150)') string
        if (string(1:1).eq.')') goto 40
        call mio_spl (150,string,nsub,lim)
        c80(1:3) = '   '
        c80 = string(lim(1,nsub):lim(2,nsub))
        if (j.eq.1) then
          algor = 0
           do k = 1, 60
            if (c80(1:3).eq.alg(k)) algor = (k + 4) / 5
          end do
          if (algor.eq.0) call mio_err (23,mem(81),lmem(81),mem(98), &
           lmem(98),c80(lim(1,nsub):lim(2,nsub)),lim(2,nsub)- &
           lim(1,nsub)+1,mem(85),lmem(85))
        end if
        if (j.eq.2) read (c80,*,err=661) tstart
        if (j.eq.3) read (c80,*,err=661) tstop
        if (j.eq.4) read (c80,*,err=661) dtout
        if (j.eq.5) read (c80,*,err=661) h0
        if (j.eq.6) read (c80,*,err=661) tol
        c1 = c80(1:1)
        if (j.eq.7.and.(c1.eq.'y'.or.c1.eq.'Y')) opt(1) = 1
        if (j.eq.8.and.(c1.eq.'n'.or.c1.eq.'N')) opt(2) = 0
        if (j.eq.9.and.(c1.eq.'y'.or.c1.eq.'Y')) opt(2) = 2
        if (j.eq.10.and.(c1.eq.'d'.or.c1.eq.'D')) opt(3) = 0
        if (j.eq.11.and.(c1.eq.'y'.or.c1.eq.'Y')) opt(3) = opt(3) + 2
        if (j.eq.12) then
          if(c1.eq.'l'.or.c1.eq.'L') then
            opt(4) = 1
          else if (j.eq.12.and.(c1.eq.'m'.or.c1.eq.'M')) then
            opt(4) = 2
          else if (j.eq.12.and.(c1.eq.'h'.or.c1.eq.'H')) then
            opt(4) = 3
          else
            goto 661
          end if
        end if
        if (j.eq.15.and.(c1.eq.'y'.or.c1.eq.'Y')) opt(8) = 1
        if (j.eq.16) read (c80,*,err=661) rmax
        if (j.eq.17) read (c80,*,err=661) rcen
        if (j.eq.18) read (c80,*,err=661) m(1)
        if (j.eq.19) read (c80,*,err=661) jcen(1)
        if (j.eq.20) read (c80,*,err=661) jcen(2)
        if (j.eq.21) read (c80,*,err=661) jcen(3)
        if (j.eq.24) read (c80,*,err=661) cefac
        if (j.eq.25) read (c80,*,err=661) ndump
        if (j.eq.26) read (c80,*,err=661) nfun
      end do
      h0 = abs(h0)
      tol = abs(tol)
      rmax = abs(rmax)
      rcen = abs(rcen)
      cefac = abs(cefac)
      close (13)
!
! Change quantities for central object to suitable units
      m(1) = abs(m(1)) * K2
      jcen(1) = jcen(1) * rcen * rcen
      jcen(2) = jcen(2) * rcen * rcen * rcen * rcen
      jcen(3) = jcen(3) * rcen * rcen * rcen * rcen * rcen * rcen
      s(1,1) = 0.d0
      s(2,1) = 0.d0
      s(3,1) = 0.d0
!
! Make sure that RCEN isn't too small, since it is used to scale the output
! (Minimum value corresponds to a central body with density 100g/cm^3).
      temp = 1.1235d-3 * m(1) ** .333333333333333d0
      if (rcen.lt.temp) then
        rcen = temp
        write (13,'(/,2a)') mem(121)(1:lmem(121)),mem(131)(1:lmem(131))
      end if
!
!------------------------------------------------------------------------------
!
!  READ  IN  DATA  FOR  BIG  AND  SMALL  BODIES
!
      nbod = 1
      do j = 1, 2
        if (j.eq.2) nbig = nbod
!
! Check if the file containing data for Big bodies exists, and open it
        filename = infile(j)
        if (oldflag) filename = dumpfile(j)
        inquire (file=filename, exist=test)
        if (.not.test) call mio_err (23,mem(81),lmem(81),mem(88), &
         lmem(88),' ',1,filename,80)
 110    open (11, file=filename, status='old', err=110)
!
! Read data style
 120    read (11,'(a150)') string
        if (string(1:1).eq.')') goto 120
        call mio_spl (150,string,nsub,lim)
        c3 = string(lim(1,nsub):(lim(1,nsub)+2))
        if (c3.eq.'Car'.or.c3.eq.'car'.or.c3.eq.'CAR') then
          informat = 1
        else if (c3.eq.'Ast'.or.c3.eq.'ast'.or.c3.eq.'AST') then
          informat = 2
        else if (c3.eq.'Com'.or.c3.eq.'com'.or.c3.eq.'COM') then
          informat = 3
        else
          call mio_err (23,mem(81),lmem(81),mem(91),lmem(91),' ',1, &
           mem(82+j),lmem(82+j))
        end if
!
! Read epoch of Big bodies
        if (j.eq.1) then 
 125      read (11,'(a150)') string
          if (string(1:1).eq.')') goto 125
          call mio_spl (150,string,nsub,lim)
          read (string(lim(1,nsub):lim(2,nsub)),*,err=667) time
        end if
!
! Read information for each object
 130    read (11,'(a)',end=140) string
        if (string(1:1).eq.')') goto 130
        call mio_spl (150,string,nsub,lim)
        if (lim(1,1).eq.-1) goto 140
!
! Determine the name of the object
        nbod = nbod + 1
        if (nbod.gt.NMAX) call mio_err (23,mem(81),lmem(81),mem(90), &
         lmem(90),' ',1,mem(82),lmem(82))
!
        if ((lim(2,1)-lim(1,1)).gt.7) then
          write (23,'(/,3a)') mem(121)(1:lmem(121)), &
           mem(122)(1:lmem(122)),string( lim(1,1):lim(2,1) )
        end if
        id(nbod) = string( lim(1,1):min(7+lim(1,1),lim(2,1)) )
! Check if another object has the same name
        do k = 1, nbod - 1
          if (id(k).eq.id(nbod)) call mio_err (23,mem(81),lmem(81), &
           mem(103),lmem(103),id(nbod),8,' ',1)
        end do
!
! Default values of mass, close-encounter limit, density etc.
        m(nbod) = 0.d0
        rceh(nbod) = 1.d0
        rho(nbod) = rhocgs
        epoch(nbod) = time
        do k = 1, 4
          ngf(k,nbod) = 0.d0
        end do
!
! Read values of mass, close-encounter limit, density etc.
        do k = 3, nsub, 2
          c80 = string(lim(1,k-1):lim(2,k-1))
          read (string(lim(1,k):lim(2,k)),*,err=666) temp
          if (c80(1:1).eq.'m'.or.c80(1:1).eq.'M') then
            m(nbod) = temp * K2
          else if (c80(1:1).eq.'r'.or.c80(1:1).eq.'R') then
            rceh(nbod) = temp
          else if (c80(1:1).eq.'d'.or.c80(1:1).eq.'D') then
            rho(nbod) = temp * rhocgs
          else if (m(nbod).lt.0.or.rceh(nbod).lt.0.or.rho(nbod).lt.0) &
           then
            call mio_err (23,mem(81),lmem(81),mem(97),lmem(97),id(nbod), &
             8,mem(82+j),lmem(82+j))
          else if (c80(1:2).eq.'ep'.or.c80(1:2).eq.'EP'.or.c80(1:2) &
             .eq.'Ep') then
            epoch (nbod) = temp
          else if (c80(1:2).eq.'a1'.or.c80(1:2).eq.'A1') then
            ngf (1,nbod) = temp
          else if (c80(1:2).eq.'a2'.or.c80(1:2).eq.'A2') then
            ngf (2,nbod) = temp
          else if (c80(1:2).eq.'a3'.or.c80(1:2).eq.'A3') then
            ngf (3,nbod) = temp
          else if (c80(1:1).eq.'b'.or.c80(1:1).eq.'B') then
            ngf (4,nbod) = temp
          else
            goto 666
          end if
        end do
!
! If required, read Cartesian coordinates, velocities and spins of the bodies
        jtmp = 100
 135    read (11,'(a150)',end=666) string
        if (string(1:1).eq.')') goto 135
        backspace 11
        if (informat.eq.1) then
          read (11,*,err=666) x(1,nbod),x(2,nbod),x(3,nbod), &
           v(1,nbod),v(2,nbod),v(3,nbod),s(1,nbod),s(2,nbod),s(3,nbod)
        else
          read (11,*,err=666) a,e,i,p,n,l,s(1,nbod),s(2,nbod), &
           s(3,nbod)
          i = i * DR
          p = (p + n) * DR
          n = n * DR
          temp = m(nbod)  +  m(1)
!
! Alternatively, read Cometary or asteroidal elements
          if (informat.eq.3) then
            q = a
            a = q / (1.d0 - e)
            l = mod (sqrt(temp/(abs(a*a*a))) * (epoch(nbod) - l), TWOPI)
          else
            q = a * (1.d0 - e)
            l = l * DR
          end if
          if (algor.eq.11.and.nbod.ne.2) temp = temp + m(2)
          call mco_el2x (temp,q,e,i,p,n,l,x(1,nbod),x(2,nbod),x(3,nbod), &
           v(1,nbod),v(2,nbod),v(3,nbod))
        end if
!
        s(1,nbod) = s(1,nbod) * K2
        s(2,nbod) = s(2,nbod) * K2
        s(3,nbod) = s(3,nbod) * K2
!
        goto 130
 140    close (11)
      end do
!
! Set non-gravitational-forces flag, NGFLAG
      ngflag = 0
      do j = 2, nbod
        if (ngf(1,j).ne.0.or.ngf(2,j).ne.0.or.ngf(3,j).ne.0) then
          if (ngflag.eq.0) ngflag = 1
          if (ngflag.eq.2) ngflag = 3
        else if (ngf(4,j).ne.0) then
          if (ngflag.eq.0) ngflag = 2
          if (ngflag.eq.1) ngflag = 3
        end if
      end do
!
!------------------------------------------------------------------------------
!
!  IF  CONTINUING  AN  OLD  INTEGRATION
!
      if (oldflag) then
        if (opt(3).eq.1) then
          call mio_jd2y (time,year,month,t1)
          write (23,'(/,a,i10,i2,f8.5,/)') mem(62)(1:lmem(62)),year, &
           month,t1
        else if (opt(3).eq.3) then
          t1 = (time - tstart) / 365.25d0
          write (23,'(/,a,f18.7,a,/)') mem(62)(1:lmem(62)),t1, &
           mem(2)(1:lmem(2))
        else
          if (opt(3).eq.0) t1 = time
          if (opt(3).eq.2) t1 = time - tstart
          write (23,'(/,a,f18.5,a,/)') mem(62)(1:lmem(62)),t1, &
           mem(1)(1:lmem(1))
        end if
!
! Read in energy and angular momentum variables, and convert to internal units
 330    open (35, file=dumpfile(4), status='old', err=330)
          read (35,*) opflag
          read (35,*) en(1),am(1),en(3),am(3)
          en(1) = en(1) * K2
          en(3) = en(3) * K2
          am(1) = am(1) * K2
          am(3) = am(3) * K2
          read (35,*) s(1,1),s(2,1),s(3,1)
          s(1,1) = s(1,1) * K2
          s(2,1) = s(2,1) * K2
          s(3,1) = s(3,1) * K2
        close (35)
        if (opflag.eq.0) opflag = 1
!
!------------------------------------------------------------------------------
!
!  IF  STARTING  A  NEW  INTEGRATION
!
      else
        opflag = -2
!
! Write integration parameters to information file
        write (23,'(/,a)') mem(11)(1:lmem(11))
        write (23,'(a)') mem(12)(1:lmem(12))
        j = algor + 13
        write (23,'(/,2a)') mem(13)(1:lmem(13)),mem(j)(1:lmem(j))
        if (tstart.ge.1.d11.or.tstart.le.-1.d10) then
          write (23,'(/,a,1p,e19.13,a)') mem(26)(1:lmem(26)),tstart, &
         mem(1)(1:lmem(1))
        else
          write (23,'(/,a,f19.7,a)') mem(26)(1:lmem(26)),tstart, &
         mem(1)(1:lmem(1))
        end if
        if (tstop.ge.1.d11.or.tstop.le.-1.d10) then
          write (23,'(a,1p,e19.13)') mem(27)(1:lmem(27)),tstop
        else
          write (23,'(a,f19.7)') mem(27)(1:lmem(27)),tstop
        end if
        write (23,'(a,f15.3)') mem(28)(1:lmem(28)),dtout
        if (opt(4).eq.1) write (23,'(2a)') mem(40)(1:lmem(40)), &
         mem(7)(1:lmem(7))
        if (opt(4).eq.2) write (23,'(2a)') mem(40)(1:lmem(40)), &
         mem(8)(1:lmem(8))
        if (opt(4).eq.3) write (23,'(2a)') mem(40)(1:lmem(40)), &
         mem(9)(1:lmem(9))
!
        write (23,'(/,a,f10.3,a)') mem(30)(1:lmem(30)),h0, &
         mem(1)(1:lmem(1))
        write (23,'(a,1p1e10.4)') mem(31)(1:lmem(31)),tol
        write (23,'(a,1p1e10.4,a)') mem(32)(1:lmem(32)),m(1)/K2, &
         mem(3)(1:lmem(3))
        write (23,'(a,1p1e11.4)') mem(33)(1:lmem(33)),jcen(1)/rcen**2
        write (23,'(a,1p1e11.4)') mem(34)(1:lmem(34)),jcen(2)/rcen**4
        write (23,'(a,1p1e11.4)') mem(35)(1:lmem(35)),jcen(3)/rcen**6
        write (23,'(a,1p1e10.4,a)') mem(36)(1:lmem(36)),rmax, &
         mem (4)(1:lmem(4))
        write (23,'(a,1p1e10.4,a)') mem(37)(1:lmem(37)),rcen, &
         mem (4)(1:lmem(4))
!
        itmp = 5
        if (opt(2).eq.1.or.opt(2).eq.2) itmp = 6
        write (23,'(/,2a)') mem(41)(1:lmem(41)),mem(itmp)(1:lmem(itmp))
        itmp = 5
        if (opt(2).eq.2) itmp = 6
        write (23,'(2a)') mem(42)(1:lmem(42)),mem(itmp)(1:lmem(itmp))
        itmp = 5
        if (opt(7).eq.1) itmp = 6
        write (23,'(2a)') mem(45)(1:lmem(45)),mem(itmp)(1:lmem(itmp))
        itmp = 5
        if (opt(8).eq.1) itmp = 6
        write (23,'(2a)') mem(46)(1:lmem(46)),mem(itmp)(1:lmem(itmp))
!
! Check that element and close-encounter files don't exist, and create them
        do j = 1, 2
          inquire (file=outfile(j), exist=test)
          if (test) call mio_err (23,mem(81),lmem(81),mem(87),lmem(87), &
           ' ',1,outfile(j),80)
 430      open  (20+j, file=outfile(j), status='new', err=430)
          close (20+j)
        end do
!
! Check that dump files don't exist, and then create them
        do j = 1, 4
          inquire (file=dumpfile(j), exist=test)
          if (test) call mio_err (23,mem(81),lmem(81),mem(87),lmem(87), &
           ' ',1,dumpfile(j),80)
 450      open  (30+j, file=dumpfile(j), status='new', err=450)
          close (30+j)
        end do
!
! Write number of Big bodies and Small bodies to information file
        write (23,'(/,a,i4)') mem(38)(1:lmem(38)), nbig - 1
        write (23,'(a,i4)') mem(39)(1:lmem(39)), nbod - nbig
!
! Calculate initial energy and angular momentum and write to information file
        s(1,1) = 0.d0
        s(2,1) = 0.d0
        s(3,1) = 0.d0
        call mxx_en (jcen,nbod,nbig,m,x,v,s,en(1),am(1))
        write (23,'(//,a)') mem(51)(1:lmem(51))
        write (23,'(a)')    mem(52)(1:lmem(52))
        write (23,'(/,a,1p1e12.5,a)') mem(53)(1:lmem(53)),en(1)/K2, &
         mem(72)(1:lmem(72))
        write (23,'(a,1p1e12.5,a)')   mem(54)(1:lmem(54)),am(1)/K2, &
         mem(73)(1:lmem(73))
!
! Initialize lost energy and angular momentum
        en(3) = 0.d0
        am(3) = 0.d0
!
! Write warning messages if necessary
        if (tstop.lt.tstart) write (23,'(/,2a)') mem(121)(1:lmem(121)), &
         mem(123)(1:lmem(123))
        if (nbig.le.0) write (23,'(/,2a)') mem(121)(1:lmem(121)), &
         mem(124)(1:lmem(124))
        if (nbig.eq.nbod) write (23,'(/,2a)') mem(121)(1:lmem(121)), &
         mem(125)(1:lmem(125))
      end if
!
!------------------------------------------------------------------------------
!
!  CHECK  FOR  ATTEMPTS  TO  DO  INCOMPATIBLE  THINGS
!
! If using close-binary algorithm, set radius of central body to be no less
! than the periastron of binary star.
      if (algor.eq.11) then
        temp = m(1) + m(2)
        call mco_x2el (temp,x(1,2),x(2,2),x(3,2),v(1,2),v(2,2),v(3,2), &
         a,tmp2,tmp3,tmp4,tmp5,tmp6)
        rcen = max (rcen, a)
      end if
!
! Check if non-grav forces are being used with an incompatible algorithm
      if (ngflag.ne.0.and.(algor.eq.3.or.algor.eq.11.or.algor.eq.12)) &
       call mio_err (23,mem(81),lmem(81),mem(92),lmem(92),' ',1, &
       mem(85),lmem(85))
!
! Check if user-defined force routine is being used with wrong algorithm
      if (opt(8).eq.1.and.(algor.eq.11.or.algor.eq.12)) call mio_err &
       (23,mem(81),lmem(81),mem(93),lmem(93),' ',1,mem(85),lmem(85))
!
! Check whether MVS is being used to integrate massive Small bodies,
! or whether massive Small bodies have different epochs than Big bodies.
      flag1 = .false.
      flag2 = .false.
      do j = nbig + 1, nbod
        if (m(j).ne.0) then
          if (algor.eq.1) call mio_err (23,mem(81),lmem(81),mem(94), &
           lmem(94),' ',1,mem(85),lmem(85))
          flag1 = .true.
        end if
        if (epoch(j).ne.time) flag2 = .true.
      end do
      if (flag1.and.flag2) call mio_err (23,mem(81),lmem(81),mem(95), &
         lmem(95),' ',1,mem(84),lmem(84))
!
! Check if central oblateness is being used with close-binary algorithm
      if (algor.eq.11.and.(jcen(1).ne.0.or.jcen(2).ne.0.or.jcen(3) &
       .ne.0)) call mio_err (23,mem(81),lmem(81),mem(102),lmem(102), &
       ' ',1,mem(85),lmem(85))
!
! Check whether RCEN > RMAX or RMAX/RCEN is very large
      if (rcen.gt.rmax) call mio_err (23,mem(81),lmem(81),mem(105), &
       lmem(105),' ',1,mem(85),lmem(85))
      if (rmax/rcen.ge.1.d12) write (23,'(/,2a,/a)')  &
       mem(121)(1:lmem(121)),mem(106)(1:lmem(106)),mem(85)(1:lmem(85))
      close (23)
      return
!
! Error reading from the input file containing integration parameters
 661  write (c3,'(i3)') lineno
      call mio_err (23,mem(81),lmem(81),mem(99),lmem(99),c3,3, &
       mem(85),lmem(85))
!
! Error reading from the input file for Big or Small bodies
 666  call mio_err (23,mem(81),lmem(81),mem(100),lmem(100),id(nbod),8, &
       mem(82+j),lmem(82+j))
!
! Error reading epoch of Big bodies
 667  call mio_err (23,mem(81),lmem(81),mem(101),lmem(101),' ',1, &
       mem(83),lmem(83))
!
!------------------------------------------------------------------------------
! ##K##
! Always initialize the variable array STAT for all bodies with 0.

      do j = 2, nbod
        stat(j) = 0
      end do

! ##K##
!------------------------------------------------------------------------------
!
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MIO_LOG.FOR    (ErikSoft   25 February 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Writes a progress report to the log file (or the screen if you are running
! Mercury interactively).
!
!------------------------------------------------------------------------------
!
      subroutine mio_log (time,tstart,en,am,opt,mem,lmem)
!
      implicit none
      include 'mercury95.inc'
!
! Input/Output
      integer lmem(NMESS), opt(8)
      real(8) time, tstart, en(3), am(3)
      character(80) mem(NMESS)
!
! Local
      integer year, month
      real(8) tmp0, tmp1, t1
      character(38) flog
      character(6) tstring
!
!------------------------------------------------------------------------------
!
      if (opt(3).eq.0.or.opt(3).eq.2) then
        tstring = mem(1)
        flog = '(1x,a,f14.1,a,2(a,1p1e12.5))'
      else if (opt(3).eq.1) then
        flog = '(1x,a,i10,1x,i2,1x,f4.1,2(a,1p1e12.5))'
      else
        tstring = mem(2)
        flog = '(1x,a,f14.3,a,2(a,1p1e12.5))'
      end if
!
      tmp0 = 0.d0
      tmp1 = 0.d0
      if (en(1).ne.0) tmp0 = (en(2) + en(3) - en(1)) / abs(en(1))
      if (am(1).ne.0) tmp1 = (am(2) + am(3) - am(1)) / abs(am(1))
!
      if (opt(3).eq.1) then
        call mio_jd2y (time,year,month,t1)
        write (*,flog) mem(64)(1:lmem(64)), year, month, t1, &
         mem(65)(1:lmem(65)), tmp0,mem(66)(1:lmem(66)), tmp1
      else
        if (opt(3).eq.0) t1 = time
        if (opt(3).eq.2) t1 = time - tstart
        if (opt(3).eq.3) t1 = (time - tstart) / 365.25d0
        write (*,flog) mem(63)(1:lmem(63)), t1, tstring, &
         mem(65)(1:lmem(65)), tmp0, mem(66)(1:lmem(66)), tmp1
      end if
!
!------------------------------------------------------------------------------
!
      return
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MIO_OUT.FOR    (ErikSoft   13 February 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Writes output variables for each object to an output file. Each variable
! is scaled between the minimum and maximum possible values and then
! written in a compressed format using ASCII characters.
! The output variables are:
!  r = the radial distance
!  theta = polar angle
!  phi = azimuthal angle
!  fv = 1 / [1 + 2(ke/be)^2], where be and ke are the object's binding and
!                             kinetic energies. (Note that 0 < fv < 1).
!  vtheta = polar angle of velocity vector
!  vphi = azimuthal angle of the velocity vector
!
! If this is the first output (OPFLAG = -1), or the first output since the 
! number of the objects or their masses have changed (OPFLAG = 1), then 
! the names, masses and spin components of all the objects are also output.
!
! N.B. Each object's distance must lie between RCEN < R < RMAX
! ===  
!
!------------------------------------------------------------------------------
!
      subroutine mio_out (time,jcen,rcen,rmax,nbod,nbig,m,xh,vh,s,rho, &
       stat,id,opt,opflag,algor,outfile)
!
      implicit none
      include 'mercury95.inc'
!
! Input/Output
      integer nbod, nbig, stat(nbod), opt(8), opflag, algor
      real(8) time,jcen(3),rcen,rmax,m(nbod),xh(3,nbod),vh(3,nbod)
      real(8) s(3,nbod),rho(nbod)
      character(80) outfile
      character(8) id(nbod)
!
! Local
      integer k, len, nchar
      real(8) rhocgs,k_2,rfac,rcen_2,fr,fv,theta,phi,vtheta,vphi
      character(80) header,c(NMAX)
      character(8) mio_fl2c,mio_re2c
      character(5) fout
!
!------------------------------------------------------------------------------
!
      rhocgs = AU * AU * AU * K2 / MSUN
      k_2 = 1.d0 / K2
      rcen_2 = 1.d0 / (rcen * rcen)
!
! Scaling factor (maximum possible range) for distances
      rfac = log10 (rmax / rcen)
!
! Create the format list, FOUT, used when outputting the orbital elements
      if (opt(4).eq.1) nchar = 2
      if (opt(4).eq.2) nchar = 4
      if (opt(4).eq.3) nchar = 7
      len = 3  +  6 * nchar
      fout(1:5) = '(a  )'
      if (len.lt.10) write (fout(3:3),'(i1)') len
      if (len.ge.10) write (fout(3:4),'(i2)') len
!
! Open the orbital-elements output file
  10  open (21, file=outfile, status='old', access='append', err=10)
!
!------------------------------------------------------------------------------
!
!  SPECIAL  OUTPUT  PROCEDURE
!
! If this is a new integration or a complete output is required (e.g. because
! the number of objects has changed), then output object details & parameters.
      if (opflag.eq.-1.or.opflag.eq.1) then
!
! Compose a header line with time, number of objects and relevant parameters
        header(1:8)   = mio_fl2c (time)
        header(9:16)  = mio_re2c (dble(nbig - 1),   0.d0, 11239423.99d0)
        header(12:19) = mio_re2c (dble(nbod - nbig),0.d0, 11239423.99d0)
        header(15:22) = mio_fl2c (m(1) * k_2)
        header(23:30) = mio_fl2c (jcen(1) * rcen_2)
        header(31:38) = mio_fl2c (jcen(2) * rcen_2 * rcen_2)
        header(39:46) = mio_fl2c (jcen(3) * rcen_2 * rcen_2 * rcen_2)
        header(47:54) = mio_fl2c (rcen)
        header(55:62) = mio_fl2c (rmax)
!
! For each object, compress its index number, name, mass, spin components
! and density (some of these need to be converted to normal units).
        do k = 2, nbod
          c(k)(1:8) = mio_re2c (dble(k - 1), 0.d0, 11239423.99d0)
          c(k)(4:11) = id(k)
          c(k)(12:19) = mio_fl2c (m(k) * k_2)
          c(k)(20:27) = mio_fl2c (s(1,k) * k_2)
          c(k)(28:35) = mio_fl2c (s(2,k) * k_2)
          c(k)(36:43) = mio_fl2c (s(3,k) * k_2)
          c(k)(44:51) = mio_fl2c (rho(k) / rhocgs)
        end do
!
! Write compressed output to file
        write (21,'(a1,a2,i2,a62,i1)') char(12),'6a',algor,header(1:62), &
         opt(4)
        do k = 2, nbod
          write (21,'(a51)') c(k)(1:51)
        end do
      end if
!
!------------------------------------------------------------------------------
!
!  NORMAL  OUTPUT  PROCEDURE
!
! Compose a header line containing the time and number of objects
      header(1:8)   = mio_fl2c (time)
      header(9:16)  = mio_re2c (dble(nbig - 1),    0.d0, 11239423.99d0)
      header(12:19) = mio_re2c (dble(nbod - nbig), 0.d0, 11239423.99d0)
!
! Calculate output variables for each body and convert to compressed format
      do k = 2, nbod
        call mco_x2ov (rcen,rmax,m(1),m(k),xh(1,k),xh(2,k),xh(3,k), &
         vh(1,k),vh(2,k),vh(3,k),fr,theta,phi,fv,vtheta,vphi)
!
! Object's index number and output variables
        c(k)(1:8) = mio_re2c (dble(k - 1), 0.d0, 11239423.99d0)
        c(k)(4:11)                 = mio_re2c (fr,     0.d0, rfac)
        c(k)(4+  nchar:11+  nchar) = mio_re2c (theta,  0.d0, PI)
        c(k)(4+2*nchar:11+2*nchar) = mio_re2c (phi,    0.d0, TWOPI)
        c(k)(4+3*nchar:11+3*nchar) = mio_re2c (fv,     0.d0, 1.d0)
        c(k)(4+4*nchar:11+4*nchar) = mio_re2c (vtheta, 0.d0, PI)
        c(k)(4+5*nchar:11+5*nchar) = mio_re2c (vphi,   0.d0, TWOPI)
      end do
!
! Write compressed output to file
      write (21,'(a1,a2,a14)') char(12),'6b',header(1:14)
      do k = 2, nbod
        write (21,fout) c(k)(1:len)
      end do
!
      close (21)
      opflag = 0
!
!------------------------------------------------------------------------------
!
      return
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MIO_RE2C.FOR    (ErikSoft  27 June 1999)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Converts a real(8) variable X, where XMIN <= X < XMAX, into an ASCII string
! of 8 characters, using the new format compression: 
!
! X is first converted to base 224, and then each base 224 digit is converted 
! to an ASCII character, such that 0 -> character 32, 1 -> character 33...
! and 223 -> character 255.
!
! ASCII characters 0 - 31 (CTRL characters) are not used, because they
! cause problems when using some operating systems.
!
!------------------------------------------------------------------------------
!
      function mio_re2c (x,xmin,xmax)
!
      implicit none
!
! Input/output
      real(8) x,xmin,xmax
      character(8) mio_re2c
!
! Local
      integer j
      real(8) y,z
!
!------------------------------------------------------------------------------
!
      mio_re2c(1:8) = '        '
      y = (x - xmin) / (xmax - xmin)
!
      if (y.ge.1) then
        do j = 1, 8
          mio_re2c(j:j) = char(255)
        end do
      else if (y.gt.0) then
        z = y
        do j = 1, 8
          z = mod(z, 1.d0) * 224.d0
          mio_re2c(j:j) = char(int(z) + 32)
        end do
      end if
!
!------------------------------------------------------------------------------
!
      return
      end
!

