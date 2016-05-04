!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MCO_OV2X.FOR    (ErikSoft   28 February 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Converts output variables for an object to coordinates and velocities.
! The output variables are:
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
      subroutine mco_ov2x (rcen,rmax,mcen,m,fr,theta,phi,fv,vtheta,&
       vphi,x,y,z,u,v,w)
!
      implicit none
      include 'mercury95.inc'
!
! Input/Output
      real(8) rcen,rmax,mcen,m,x,y,z,u,v,w,fr,theta,phi,fv,vtheta,vphi
!
! Local
      real(8) r,v1,temp
!
!------------------------------------------------------------------------------
!
        r = rcen * 10.d0**fr
        temp = sqrt(.5d0*(1.d0/fv - 1.d0))
        v1 = sqrt(2.d0 * temp * (mcen + m) / r)
!
        x = r * sin(theta) * cos(phi)
        y = r * sin(theta) * sin(phi)
        z = r * cos(theta)
        u = v1 * sin(vtheta) * cos(vphi)
        v = v1 * sin(vtheta) * sin(vphi)
        w = v1 * cos(vtheta)
!
!------------------------------------------------------------------------------
!
      return
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MCE_SPIN.FOR    (ErikSoft  2 December 1999)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Calculates the spin rate (in rotations per day) for a fluid body given
! its mass, spin angular momentum and density. The routine assumes the
! body is a MacClaurin ellipsoid, whose axis ratio is defined by the
! quantity SS = SQRT(A^2/C^2 - 1), where A and C are the
! major and minor axes.
!
!------------------------------------------------------------------------------
!
      subroutine mce_spin (g,mass,spin,rho,rote)
!
      implicit none
      include 'mercury95.inc'
!
! Input/Output
      real(8) g,mass,spin,rho,rote
!
! Local
      integer k
      real(8) ss,s2,f,df,z,dz,tmp0,tmp1,t23
!
!------------------------------------------------------------------------------
!
      t23 = 2.d0 / 3.d0
      tmp1 = spin * spin / (2.d0 * PI * rho * g) &
          * ( 250.d0*PI*PI*rho*rho / (9.d0*mass**5) )**t23
!
! Calculate SS using Newton's method
      ss = 1.d0
      do k = 1, 20
        s2 = ss * ss
        tmp0 = (1.d0 + s2)**t23
        call m_sfunc (ss,z,dz)
        f = z * tmp0  -  tmp1
        df = tmp0 * ( dz  +  4.d0 * ss * z / (3.d0*(1.d0 + s2)) )
        ss = ss - f/df
      end do
!
      rote = sqrt(TWOPI * g * rho * z) / TWOPI
!
!------------------------------------------------------------------------------
!
      return
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MIO_AEI.FOR    (ErikSoft   31 January 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Creates a filename and opens a file to store aei information for an object.
! The filename is based on the name of the object.
!
!------------------------------------------------------------------------------
!
      subroutine mio_aei (id,extn,unitnum,header,lenhead,mem,lmem)
!
      implicit none
      include 'mercury95.inc'
!
! Input/Output
      integer unitnum,lenhead,lmem(NMESS)
      character(4) extn
      character(8) id
      character(250) header
      character(80) mem(NMESS)
!
! Local
      integer j,k,itmp,nsub,lim(2,4)
      logical test
      character(1) bad(5)
      character(250) filename
!
!------------------------------------------------------------------------------
!
      data bad/ '*', '/', '.', ':', '&'/
!
! Create a filename based on the object's name
      call mio_spl (8,id,nsub,lim)
      itmp = min(7,lim(2,1)-lim(1,1))
      filename(1:itmp+1) = id(1:itmp+1)
      filename(itmp+2:itmp+5) = extn
      do j = itmp + 6, 250
        filename(j:j) = ' '
      end do
!
! Check for inappropriate characters in the filename
      do j = 1, itmp + 1
        do k = 1, 5
          if (filename(j:j).eq.bad(k)) filename(j:j) = '_'
        end do
      end do
!
! If the file exists already, give a warning and don't overwrite it
      inquire (file=filename, exist=test)
      if (test) then
        write (*,'(/,3a)') mem(121)(1:lmem(121)),mem(87)(1:lmem(87)),&
         filename(1:80)
        unitnum = -1
      else
        open (unitnum, file=filename, status='new')
        write (unitnum, '(/,30x,a8,//,a)') id,header(1:lenhead)
      end if
!
!------------------------------------------------------------------------------
!
      return
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MCO_H2CB.FOR    (ErikSoft   2 November 2000)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Convert coordinates with respect to the central body to close-binary
! coordinates.
!
!------------------------------------------------------------------------------
!
      subroutine mco_h2cb (jcen,nbod,nbig,h,m,xh,vh,x,v)
!
      implicit none
!
! Input/Output
      integer nbod,nbig
      real(8) jcen(3),h,m(nbod),xh(3,nbod),vh(3,nbod),x(3,nbod),v(3,nbod)
!
! Local
      integer j
      real(8) msum,mvsum(3),temp,mbin,mbin_1,mtot_1
!
!------------------------------------------------------------------------------
!
      msum = 0.d0
      mvsum(1) = 0.d0
      mvsum(2) = 0.d0
      mvsum(3) = 0.d0
      mbin = m(1) + m(2)
      mbin_1 = 1.d0 / mbin
!
      x(1,2) = xh(1,2)
      x(2,2) = xh(2,2)
      x(3,2) = xh(3,2)
      temp = m(1) * mbin_1
      v(1,2) = temp * vh(1,2)
      v(2,2) = temp * vh(2,2)
      v(3,2) = temp * vh(3,2)
!
      do j = 3, nbod
        msum = msum + m(j)
        mvsum(1) = mvsum(1)  +  m(j) * vh(1,j)
        mvsum(2) = mvsum(2)  +  m(j) * vh(2,j)
        mvsum(3) = mvsum(3)  +  m(j) * vh(3,j)
      end do
      mtot_1 = 1.d0 / (msum + mbin)
      mvsum(1) = mtot_1 * (mvsum(1) + m(2)*vh(1,2))
      mvsum(2) = mtot_1 * (mvsum(2) + m(2)*vh(2,2))
      mvsum(3) = mtot_1 * (mvsum(3) + m(2)*vh(3,2))
!
      temp = m(2) * mbin_1
      do j = 3, nbod
        x(1,j) = xh(1,j)  -  temp * xh(1,2)
        x(2,j) = xh(2,j)  -  temp * xh(2,2)
        x(3,j) = xh(3,j)  -  temp * xh(3,2)
        v(1,j) = vh(1,j)  -  mvsum(1)
        v(2,j) = vh(2,j)  -  mvsum(2)
        v(3,j) = vh(3,j)  -  mvsum(3)
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
!      M_SFUNC.FOR     (ErikSoft  14 November 1998)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Calculates Z = [ (3 + S^2)arctan(S) - 3S ] / S^3 and its derivative DZ,
! for S > 0.
!
!------------------------------------------------------------------------------
!
      subroutine m_sfunc (s,z,dz)
!
      implicit none
!
! Input/Output
      real(8) s, z, dz
!
! Local
      real(8) s2,s4,s6,s8,a
!
!------------------------------------------------------------------------------
!
      s2 = s * s
!
      if (s.gt.1.d-2) then
        a  = atan(s)
        z  = ((3.d0 + s2)*a - 3.d0*s) / (s * s2)
        dz = (2.d0*s*a - 3.d0 + (3.d0+s2)/(1.d0+s2)) / (s * s2)&
          - 3.d0 * z / s
      else
        s4 = s2 * s2
        s6 = s2 * s4
        s8 = s4 * s4
        z  = - .1616161616161616d0*s8&
            + .1904761904761905d0*s6&
            - .2285714285714286d0*s4&
            + .2666666666666667d0*s2
        dz = s * (- 1.292929292929293d0*s6&
                 + 1.142857142857143d0*s4&
                 - 0.914285714285714d0*s2&
                 + 0.533333333333333d0)
      end if
!
!------------------------------------------------------------------------------
!
      return
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      M_FORMAT.FOR    (ErikSoft   31 January 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Makes an output format list and file header for the orbital-element files
! created by M_ELEM3.FOR
! Also identifies which orbital elements will be output for each object.
!
!------------------------------------------------------------------------------
!
      subroutine m_format (string,timestyle,nel,iel,fout,header,lenhead)
!
      implicit none
      include 'mercury95.inc'
!
! Input/Output
      integer timestyle,nel,iel(22),lenhead
      character(250) string,header,fout
!
! Local
      integer i,j,pos,nsub,lim(2,20),formflag,lenfout,f1,f2,itmp
      character(1) elcode(22)
      character(4) elhead(22)
!
!------------------------------------------------------------------------------
!
      data elcode/ 'a','e','i','g','n','l','p','q','b','x','y','z',&
       'u','v','w','r','f','m','o','s','d','c'/
      data elhead/ '  a ','  e ','  i ','peri','node','  M ','long',&
       '  q ','  Q ','  x ','  y ','  z ',' vx ',' vy ',' vz ','  r ',&
       '  f ','mass','oblq','spin','dens','comp'/
!
! Initialize header to a blank string
      do i = 1, 250
        header(i:i) = ' '
      end do
!
! Create part of the format list and header for the required time style
      if (timestyle.eq.0.or.timestyle.eq.2) then
        fout(1:9) = '(1x,f18.5'
        lenfout = 9
        header(1:19) = '    Time (days)    '
        lenhead = 19
      else if (timestyle.eq.1) then
        fout(1:21) = '(1x,i10,1x,i2,1x,f8.5'
        lenfout = 21
        header(1:23) = '    Year/Month/Day     '
        lenhead = 23
      else if (timestyle.eq.3) then
        fout(1:9) = '(1x,f18.7'
        lenfout = 9
        header(1:19) = '    Time (years)   '
        lenhead = 19
      end if
!
! Identify the required elements
      call mio_spl (250,string,nsub,lim)
      do i = 1, nsub
        do j = 1, 22
          if (string(lim(1,i):lim(1,i)).eq.elcode(j)) iel(i) = j
        end do
      end do
      nel = nsub
!
! For each element, see whether normal or exponential notation is required
      do i = 1, nsub
        formflag = 0
        do j = lim(1,i)+1, lim(2,i)
          if (formflag.eq.0) pos = j
          if (string(j:j).eq.'.') formflag = 1
          if (string(j:j).eq.'e') formflag = 2
        end do
!
! Create the rest of the format list and header
        if (formflag.eq.1) then
          read (string(lim(1,i)+1:pos-1),*) f1
          read (string(pos+1:lim(2,i)),*) f2
          write (fout(lenfout+1:lenfout+10),'(a10)') ',1x,f  .  '
          write (fout(lenfout+6:lenfout+7),'(i2)') f1
          write (fout(lenfout+9:lenfout+10),'(i2)') f2
          lenfout = lenfout + 10
        else if (formflag.eq.2) then
          read (string(lim(1,i)+1:pos-1),*) f1
          write (fout(lenfout+1:lenfout+16),'(a16)') ',1x,1p,e  .  ,0p'
          write (fout(lenfout+9:lenfout+10),'(i2)') f1
          write (fout(lenfout+12:lenfout+13),'(i2)') f1 - 7
          lenfout = lenfout + 16
        end if
        itmp = (f1 - 4) / 2
        header(lenhead+itmp+2:lenhead+itmp+5) = elhead(iel(i))
        lenhead = lenhead + f1 + 1
      end do
!
      lenfout = lenfout + 1
      fout(lenfout:lenfout) = ')'
!
!------------------------------------------------------------------------------
!
      return
      end
!

