!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MCO_EL2X.FOR    (ErikSoft  7 July 1999)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Calculates Cartesian coordinates and velocities given Keplerian orbital
! elements (for elliptical, parabolic or hyperbolic orbits).
!
! Based on a routine from Levison and Duncan's SWIFT integrator.
!
!  mu = grav const * (central + secondary mass)
!  q = perihelion distance
!  e = eccentricity
!  i = inclination                 )
!  p = longitude of perihelion !!! )   in
!  n = longitude of ascending node ) radians
!  l = mean anomaly                )
!
!  x,y,z = Cartesian positions  ( units the same as a )
!  u,v,w =     "     velocities ( units the same as sqrt(mu/a) )
!
!------------------------------------------------------------------------------
!
      subroutine mco_el2x (mu,q,e,i,p,n,l,x,y,z,u,v,w)
!
      implicit none
      include 'mercury95.inc'
!
! Input/Output
      real(8) mu,q,e,i,p,n,l,x,y,z,u,v,w
!
! Local
      real(8) g,a,ci,si,cn,sn,cg,sg,ce,se,romes,temp
      real(8) z1,z2,z3,z4,d11,d12,d13,d21,d22,d23
      real(8) mco_kep, orbel_fhybrid, orbel_zget
!
!------------------------------------------------------------------------------
!
! Change from longitude of perihelion to argument of perihelion
      g = p - n
!
! Rotation factors
      call mco_sine (i,si,ci)
      call mco_sine (g,sg,cg)
      call mco_sine (n,sn,cn)
      z1 = cg * cn
      z2 = cg * sn
      z3 = sg * cn
      z4 = sg * sn
      d11 =  z1 - z4*ci
      d12 =  z2 + z3*ci
      d13 = sg * si
      d21 = -z3 - z2*ci
      d22 = -z4 + z1*ci
      d23 = cg * si
!
! Semi-major axis
      a = q / (1.d0 - e)
!
! Ellipse
      if (e.lt.1.d0) then
        romes = sqrt(1.d0 - e*e)
        temp = mco_kep (e,l)
        call mco_sine (temp,se,ce)
        z1 = a * (ce - e)
        z2 = a * romes * se
        temp = sqrt(mu/a) / (1.d0 - e*ce)
        z3 = -se * temp
        z4 = romes * ce * temp
      else
! Parabola
        if (e.eq.1.d0) then
          ce = orbel_zget(l)
          z1 = q * (1.d0 - ce*ce)
          z2 = 2.d0 * q * ce
          z4 = sqrt(2.d0*mu/q) / (1.d0 + ce*ce)
          z3 = -ce * z4
        else
! Hyperbola
          romes = sqrt(e*e - 1.d0)
          temp = orbel_fhybrid(e,l)
          call mco_sinh (temp,se,ce)
          z1 = a * (ce - e)
          z2 = -a * romes * se
          temp = sqrt(mu/abs(a)) / (e*ce - 1.d0)
          z3 = -se * temp
          z4 = romes * ce * temp
        end if
      endif
!
      x = d11 * z1  +  d21 * z2
      y = d12 * z1  +  d22 * z2
      z = d13 * z1  +  d23 * z2
      u = d11 * z3  +  d21 * z4
      v = d12 * z3  +  d22 * z4
      w = d13 * z3  +  d23 * z4
!
!------------------------------------------------------------------------------
!
      return
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MCO_H2B.FOR    (ErikSoft   2 November 2000)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Converts coordinates with respect to the central body to barycentric
! coordinates.
!
!------------------------------------------------------------------------------
!
      subroutine mco_h2b (jcen,nbod,nbig,h,m,xh,vh,x,v)
!
      implicit none
!
! Input/Output
      integer nbod,nbig
      real(8) jcen(3),h,m(nbod),xh(3,nbod),vh(3,nbod),x(3,nbod),v(3,nbod)
!
! Local
      integer j
      real(8) mtot,temp
!
!------------------------------------------------------------------------------
!
      mtot = 0.d0
      x(1,1) = 0.d0
      x(2,1) = 0.d0
      x(3,1) = 0.d0
      v(1,1) = 0.d0
      v(2,1) = 0.d0
      v(3,1) = 0.d0
!
! Calculate coordinates and velocities of the central body
      do j = 2, nbod
        mtot = mtot  +  m(j)
        x(1,1) = x(1,1)  +  m(j) * xh(1,j)
        x(2,1) = x(2,1)  +  m(j) * xh(2,j)
        x(3,1) = x(3,1)  +  m(j) * xh(3,j)
        v(1,1) = v(1,1)  +  m(j) * vh(1,j)
        v(2,1) = v(2,1)  +  m(j) * vh(2,j)
        v(3,1) = v(3,1)  +  m(j) * vh(3,j)
      enddo
!
      temp = -1.d0 / (mtot + m(1))
      x(1,1) = temp * x(1,1)
      x(2,1) = temp * x(2,1)
      x(3,1) = temp * x(3,1)
      v(1,1) = temp * v(1,1)
      v(2,1) = temp * v(2,1)
      v(3,1) = temp * v(3,1)
!
! Calculate the barycentric coordinates and velocities
      do j = 2, nbod
        x(1,j) = xh(1,j) + x(1,1)
        x(2,j) = xh(2,j) + x(2,1)
        x(3,j) = xh(3,j) + x(3,1)
        v(1,j) = vh(1,j) + v(1,1)
        v(2,j) = vh(2,j) + v(2,1)
        v(3,j) = vh(3,j) + v(3,1)
      enddo
!
!------------------------------------------------------------------------------
!
      return
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MCO_H2J.FOR    (ErikSoft   2 November 2000)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Converts coordinates with respect to the central body to Jacobi coordinates.
! Note that the Jacobi coordinates of all small bodies are assumed to be the
! same as their coordinates with respect to the central body.
!
!------------------------------------------------------------------------------
!
      subroutine mco_h2j (jcen,nbod,nbig,h,m,xh,vh,x,v)
!
      implicit none
!
! Input/Output
      integer nbod,nbig
      real(8) jcen(3),h,m(nbig),xh(3,nbig),vh(3,nbig),x(3,nbig),v(3,nbig)
!
! Local
      integer j
      real(8) mtot, mx, my, mz, mu, mv, mw, temp
!
!------------------------------------------------------------------------------c
      mtot = m(2)
      x(1,2) = xh(1,2)
      x(2,2) = xh(2,2)
      x(3,2) = xh(3,2)
      v(1,2) = vh(1,2)
      v(2,2) = vh(2,2)
      v(3,2) = vh(3,2)
      mx = m(2) * xh(1,2)
      my = m(2) * xh(2,2)
      mz = m(2) * xh(3,2)
      mu = m(2) * vh(1,2)
      mv = m(2) * vh(2,2)
      mw = m(2) * vh(3,2)
!
      do j = 3, nbig - 1
        temp = 1.d0 / (mtot + m(1))
        mtot = mtot + m(j)
        x(1,j) = xh(1,j)  -  temp * mx
        x(2,j) = xh(2,j)  -  temp * my
        x(3,j) = xh(3,j)  -  temp * mz
        v(1,j) = vh(1,j)  -  temp * mu
        v(2,j) = vh(2,j)  -  temp * mv
        v(3,j) = vh(3,j)  -  temp * mw
        mx = mx  +  m(j) * xh(1,j)
        my = my  +  m(j) * xh(2,j)
        mz = mz  +  m(j) * xh(3,j)
        mu = mu  +  m(j) * vh(1,j)
        mv = mv  +  m(j) * vh(2,j)
        mw = mw  +  m(j) * vh(3,j)
      enddo
!
      if (nbig.gt.2) then
        temp = 1.d0 / (mtot + m(1))
        x(1,nbig) = xh(1,nbig)  -  temp * mx
        x(2,nbig) = xh(2,nbig)  -  temp * my
        x(3,nbig) = xh(3,nbig)  -  temp * mz
        v(1,nbig) = vh(1,nbig)  -  temp * mu
        v(2,nbig) = vh(2,nbig)  -  temp * mv
        v(3,nbig) = vh(3,nbig)  -  temp * mw
      end if
!
      do j = nbig + 1, nbod
        x(1,j) = xh(1,j)
        x(2,j) = xh(2,j)
        x(3,j) = xh(3,j)
        v(1,j) = vh(1,j)
        v(2,j) = vh(2,j)
        v(3,j) = vh(3,j)
      end do
!
!------------------------------------------------------------------------------
!
      return
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MCO_IDEN.FOR    (ErikSoft   2 November 2000)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Makes a new copy of a set of coordinates.
!
!------------------------------------------------------------------------------
!
      subroutine mco_iden (jcen,nbod,nbig,h,m,xh,vh,x,v)
!
      implicit none
!
! Input/Output
      integer nbod,nbig
      real(8) jcen(3),h,m(nbod),x(3,nbod),v(3,nbod),xh(3,nbod),vh(3,nbod)
!
! Local
      integer j
!
!------------------------------------------------------------------------------
!
      do j = 1, nbod
        x(1,j) = xh(1,j)
        x(2,j) = xh(2,j)
        x(3,j) = xh(3,j)
        v(1,j) = vh(1,j)
        v(2,j) = vh(2,j)
        v(3,j) = vh(3,j)
      enddo
!
!------------------------------------------------------------------------------
!
      return
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MCO_KEP.FOR    (ErikSoft  7 July 1999)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Solves Kepler's equation for eccentricities less than one.
! Algorithm from A. Nijenhuis (1991) Cel. Mech. Dyn. Astron. 51, 319-330.
!
!  e = eccentricity
!  l = mean anomaly      (radians)
!  u = eccentric anomaly (   "   )
!
!------------------------------------------------------------------------------
!
      function mco_kep (e,oldl)
      implicit none
!
! Input/Outout
      real(8) oldl,e,mco_kep
!
! Local
      real(8) l,pi,twopi,piby2,u1,u2,ome,sign
      real(8) x,x2,sn,dsn,z1,z2,z3,f0,f1,f2,f3
      real(8) p,q,p2,ss,cc
      logical flag,big,bigg
!
!------------------------------------------------------------------------------
!
      pi = 3.141592653589793d0
      twopi = 2.d0 * pi
      piby2 = .5d0 * pi
!
! Reduce mean anomaly to lie in the range 0 < l < pi
      if (oldl.ge.0) then
        l = mod(oldl, twopi)
      else
        l = mod(oldl, twopi) + twopi
      end if
      sign = 1.d0
      if (l.gt.pi) then
        l = twopi - l
        sign = -1.d0
      end if
!
      ome = 1.d0 - e
!
      if (l.ge..45d0.or.e.lt..55d0) then
!
! Regions A,B or C in Nijenhuis
! -----------------------------
!
! Rough starting value for eccentric anomaly
        if (l.lt.ome) then
          u1 = ome
        else
          if (l.gt.(pi-1.d0-e)) then
            u1 = (l+e*pi)/(1.d0+e)
          else
            u1 = l + e
          end if
        end if
!
! Improved value using Halley's method
        flag = u1.gt.piby2
        if (flag) then
          x = pi - u1
        else
          x = u1
        end if
        x2 = x*x
        sn = x*(1.d0 + x2*(-.16605 + x2*.00761) )
        dsn = 1.d0 + x2*(-.49815 + x2*.03805)
        if (flag) dsn = -dsn
        f2 = e*sn
        f0 = u1 - f2 - l
        f1 = 1.d0 - e*dsn
        u2 = u1 - f0/(f1 - .5d0*f0*f2/f1)
      else
!
! Region D in Nijenhuis
! ---------------------
!
! Rough starting value for eccentric anomaly
        z1 = 4.d0*e + .5d0
        p = ome / z1
        q = .5d0 * l / z1
        p2 = p*p
        z2 = exp( log( dsqrt( p2*p + q*q ) + q )/1.5 )
        u1 = 2.d0*q / ( z2 + p + p2/z2 )
!
! Improved value using Newton's method
        z2 = u1*u1
        z3 = z2*z2
        u2 = u1 - .075d0*u1*z3 / (ome + z1*z2 + .375d0*z3)
        u2 = l + e*u2*( 3.d0 - 4.d0*u2*u2 )
      end if
!
! Accurate value using 3rd-order version of Newton's method
! N.B. Keep cos(u2) rather than sqrt( 1-sin^2(u2) ) to maintain accuracy!
!
! First get accurate values for u2 - sin(u2) and 1 - cos(u2)
      bigg = (u2.gt.piby2)
      if (bigg) then
        z3 = pi - u2
      else
        z3 = u2
      end if
!
      big = (z3.gt.(.5d0*piby2))
      if (big) then
        x = piby2 - z3
      else
        x = z3
      end if
!
      x2 = x*x
      ss = 1.d0
      cc = 1.d0
!
      ss = x*x2/6.*(1. - x2/20.*(1. - x2/42.*(1. - x2/72.*(1. - &
        x2/110.*(1. - x2/156.*(1. - x2/210.*(1. - x2/272.)))))))
      cc =   x2/2.*(1. - x2/12.*(1. - x2/30.*(1. - x2/56.*(1. - &
        x2/ 90.*(1. - x2/132.*(1. - x2/182.*(1. - x2/240.*(1. - &
        x2/306.))))))))
!
      if (big) then
        z1 = cc + z3 - 1.d0
        z2 = ss + z3 + 1.d0 - piby2
      else
        z1 = ss
        z2 = cc
      end if
!
      if (bigg) then
        z1 = 2.d0*u2 + z1 - pi
        z2 = 2.d0 - z2
      end if
!
      f0 = l - u2*ome - e*z1
      f1 = ome + e*z2
      f2 = .5d0*e*(u2-z1)
      f3 = e/6.d0*(1.d0-z2)
      z1 = f0/f1
      z2 = f0/(f2*z1+f1)
      mco_kep = sign*( u2 + f0/((f3*z1+f2)*z2+f1) )
!
!------------------------------------------------------------------------------
!
      return
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MCO_SINE.FOR    (ErikSoft  17 April 1997)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Calculates sin and cos of an angle X (in radians).
!
!------------------------------------------------------------------------------
!
      subroutine mco_sine (x,sx,cx)
!
      implicit none
!
! Input/Output
      real(8) x,sx,cx
!
! Local
      real(8) pi,twopi
!
!------------------------------------------------------------------------------
!
      pi = 3.141592653589793d0
      twopi = 2.d0 * pi
!
      if (x.gt.0) then
        x = mod(x,twopi)
      else
        x = mod(x,twopi) + twopi
      end if
!
      cx = cos(x)
!
      if (x.gt.pi) then
        sx = -sqrt(1.d0 - cx*cx)
      else
        sx =  sqrt(1.d0 - cx*cx)
      end if
!
!------------------------------------------------------------------------------
!
      return
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MCO_SINH.FOR    (ErikSoft  12 June 1998)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Calculates sinh and cosh of an angle X (in radians)
!
!------------------------------------------------------------------------------
!
      subroutine mco_sinh (x,sx,cx)
!
      implicit none
!
! Input/Output
      real(8) x,sx,cx
!
!------------------------------------------------------------------------------
!
      sx = sinh(x)
      cx = sqrt (1.d0 + sx*sx)
!
!------------------------------------------------------------------------------
!
      return
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MCO_X2EL.FOR    (ErikSoft  20 February 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Calculates Keplerian orbital elements given relative coordinates and
! velocities, and mu = G times the sum of the masses.
!
! The elements are: q = perihelion distance
!                   e = eccentricity
!                   i = inclination
!                   p = longitude of perihelion (NOT argument of perihelion!!)
!                   n = longitude of ascending node
!                   l = mean anomaly (or mean longitude if e < 1.e-8)
!
!------------------------------------------------------------------------------
!
      subroutine mco_x2el (mu,x,y,z,u,v,w,q,e,i,p,n,l)
!
      implicit none
      include 'mercury95.inc'
!
! Input/Output
      real(8) mu,q,e,i,p,n,l,x,y,z,u,v,w
!
! Local
      real(8) hx,hy,hz,h2,h,v2,r,rv,s,true
      real(8) ci,to,temp,tmp2,bige,f,cf,ce
!
!------------------------------------------------------------------------------
!
      hx = y * w  -  z * v
      hy = z * u  -  x * w
      hz = x * v  -  y * u
      h2 = hx*hx + hy*hy + hz*hz
      v2 = u * u  +  v * v  +  w * w
      rv = x * u  +  y * v  +  z * w
      r = sqrt(x*x + y*y + z*z)
      h = sqrt(h2)
      s = h2 / mu
!
! Inclination and node
      ci = hz / h
      if (abs(ci).lt.1) then
        i = acos (ci)
        n = atan2 (hx,-hy)
        if (n.lt.0) n = n + TWOPI
      else
        if (ci.gt.0) i = 0.d0
        if (ci.lt.0) i = PI
        n = 0.d0
      end if
!
! Eccentricity and perihelion distance
      temp = 1.d0  +  s * (v2 / mu  -  2.d0 / r)
      if (temp.le.0) then
        e = 0.d0
      else
        e = sqrt (temp)
      end if
      q = s / (1.d0 + e)
!
! True longitude
      if (hy.ne.0) then
        to = -hx/hy
        temp = (1.d0 - ci) * to
        tmp2 = to * to
        true = atan2((y*(1.d0+tmp2*ci)-x*temp),(x*(tmp2+ci)-y*temp))
      else
        true = atan2(y * ci, x)
      end if
      if (ci.lt.0) true = true + PI
!
      if (e.lt.3.d-8) then
        p = 0.d0
        l = true
      else
        ce = (v2*r - mu) / (e*mu)
!
! Mean anomaly for ellipse
        if (e.lt.1) then
          if (abs(ce).gt.1) ce = sign(1.d0,ce)
          bige = acos(ce)
          if (rv.lt.0) bige = TWOPI - bige
          l = bige - e*sin(bige)
        else
!
! Mean anomaly for hyperbola
          if (ce.lt.1) ce = 1.d0
          bige = log( ce + sqrt(ce*ce-1.d0) )
          if (rv.lt.0) bige = - bige
          l = e*sinh(bige) - bige
        end if
!
! Longitude of perihelion
        cf = (s - r) / (e*r)
        if (abs(cf).gt.1) cf = sign(1.d0,cf)
        f = acos(cf)
        if (rv.lt.0) f = TWOPI - f
        p = true - f
        p = mod (p + TWOPI + TWOPI, TWOPI)
      end if
!
      if (l.lt.0.and.e.lt.1) l = l + TWOPI
      if (l.gt.TWOPI.and.e.lt.1) l = mod (l, TWOPI)
!
!------------------------------------------------------------------------------
!
      return
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MIO_C2FL.FOR    (ErikSoft   5 June 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Converts a character(8) ASCII string into a real(8) variable.
!
! N.B. X will lie in the range -1.e112 < X < 1.e112
! ===
!
!------------------------------------------------------------------------------
!
      function mio_c2fl (c)
!
      implicit none
!
! Input/Output
      real(8) mio_c2fl
      character(8) c
!
! Local
      real(8) x,mio_c2re
      integer ex
!
!------------------------------------------------------------------------------
!
      x = mio_c2re (c(1:8), 0.d0, 1.d0, 7)
      x = x * 2.d0 - 1.d0
      ex = mod(ichar(c(8:8)) + 256, 256) - 32 - 112
      mio_c2fl = x * (10.d0**dble(ex))
!
!------------------------------------------------------------------------------
!
      return
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MIO_C2RE.FOR    (ErikSoft   5 June 2001)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Converts an ASCII string into a real(8) variable X, where XMIN <= X < XMAX,
! using the new format compression:
!
! X is assumed to be made up of NCHAR base-224 digits, each one represented
! by a character in the ASCII string. Each digit is given by the ASCII
! number of the character minus 32.
! The first 32 ASCII characters (CTRL characters) are avoided, because they
! cause problems when using some operating systems.
!
!------------------------------------------------------------------------------
!
      function mio_c2re (c,xmin,xmax,nchar)
!
      implicit none
!
! Input/output
      integer nchar
      real(8) xmin,xmax,mio_c2re
      character(8) c
!
! Local
      integer j
      real(8) y
!
!------------------------------------------------------------------------------
!
      y = 0
      do j = nchar, 1, -1
        y = (y + dble(mod(ichar(c(j:j)) + 256, 256) - 32)) / 224.d0
      end do
!
      mio_c2re = xmin  +  y * (xmax - xmin)
!
!------------------------------------------------------------------------------
!
      return
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MIO_ERR.FOR    (ErikSoft  6 December 1999)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Writes out an error message and terminates Mercury.
!
!------------------------------------------------------------------------------
!
      subroutine mio_err (unit,s1,ls1,s2,ls2,s3,ls3,s4,ls4)
!
      implicit none
!
! Input/Output
      integer unit,ls1,ls2,ls3,ls4
      character(80) s1,s2,s3,s4
!
!------------------------------------------------------------------------------
!
      write (*,'(/,2a)') ' ERROR: Programme terminated. See information' &
       ,' file for details.'
!
      write (unit,'(/,3a,/,2a)') s1(1:ls1),s2(1:ls2),s3(1:ls3), &
       ' ',s4(1:ls4)
      stop
!
!------------------------------------------------------------------------------
!
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MIO_JD2Y.FOR    (ErikSoft  7 July 1999)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Converts from Julian day number to Julian/Gregorian Calendar dates, assuming
! the dates are those used by the English calendar.
!
! Algorithm taken from `Practical Astronomy with your calculator' (1988)
! by Peter Duffett-Smith, 3rd edition, C.U.P.
!
! Algorithm for negative Julian day numbers (Julian calendar assumed) by
! J. E. Chambers.
!
! N.B. The output date is with respect to the Julian Calendar on or before
! ===  4th October 1582, and with respect to the Gregorian Calendar on or 
!      after 15th October 1582.
!
!
!------------------------------------------------------------------------------
!
      subroutine mio_jd2y (jd0,year,month,day)
!
      implicit none
!
! Input/Output
      integer year,month
      real(8) jd0,day
!
! Local
      integer i,a,b,c,d,e,g
      real(8) jd,f,temp,x,y,z
!
!------------------------------------------------------------------------------
!
      if (jd0.le.0) goto 50
!
      jd = jd0 + 0.5d0
      i = sign( dint(dabs(jd)), jd )
      f = jd - 1.d0*i
!
! If on or after 15th October 1582
      if (i.gt.2299160) then
        temp = (1.d0*i - 1867216.25d0) / 36524.25d0
        a = sign( dint(dabs(temp)), temp )
        temp = .25d0 * a
        b = i + 1 + a - sign( dint(dabs(temp)), temp )
      else
        b = i
      end if
!
      c = b + 1524
      temp = (1.d0*c - 122.1d0) / 365.25d0
      d = sign( dint(dabs(temp)), temp )
      temp = 365.25d0 * d
      e = sign( dint(dabs(temp)), temp )
      temp = (c-e) / 30.6001d0
      g = sign( dint(dabs(temp)), temp )
!
      temp = 30.6001d0 * g
      day = 1.d0*(c-e) + f - 1.d0*sign( dint(dabs(temp)), temp )
!
      if (g.le.13) month = g - 1
      if (g.gt.13) month = g - 13
!
      if (month.gt.2) year = d - 4716
      if (month.le.2) year = d - 4715
!
      if (day.gt.32) then
        day = day - 32
        month = month + 1
      end if
!
      if (month.gt.12) then
        month = month - 12
        year = year + 1
      end if
      return
!
  50  continue
!
! Algorithm for negative Julian day numbers (Duffett-Smith doesn't work)
      x = jd0 - 2232101.5
      f = x - dint(x)
      if (f.lt.0) f = f + 1.d0
      y = dint(mod(x,1461.d0) + 1461.d0)
      z = dint(mod(y,365.25d0))
      month = int((z + 0.5d0) / 30.61d0)
      day = dint(z + 1.5d0 - 30.61d0*dble(month)) + f
      month = mod(month + 2, 12) + 1
!
      year = 1399 + int (x / 365.25d0)
      if (x.lt.0) year = year - 1
      if (month.lt.3) year = year + 1
!
!------------------------------------------------------------------------------
!
      return
      end
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!      MIO_SPL.FOR    (ErikSoft  14 November 1999)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Author: John E. Chambers
!
! Given a character string STRING, of length LEN bytes, the routine finds 
! the beginnings and ends of NSUB substrings present in the original, and 
! delimited by spaces. The positions of the extremes of each substring are 
! returned in the array DELIMIT.
! Substrings are those which are separated by spaces or the = symbol.
!
!------------------------------------------------------------------------------
!
      subroutine mio_spl (len,string,nsub,delimit)
!
      implicit none
!
! Input/Output
      integer len,nsub,delimit(2,100)
      character(1) string(len)
!
! Local
      integer j,k
      character(1) c
!
!------------------------------------------------------------------------------
!
      nsub = 0
      j = 0
      c = ' '
      delimit(1,1) = -1
!
! Find the start of string
  10  j = j + 1
      if (j.gt.len) goto 99
      c = string(j)
      if (c.eq.' '.or.c.eq.'=') goto 10
!
! Find the end of string
      k = j
  20  k = k + 1
      if (k.gt.len) goto 30
      c = string(k)
      if (c.ne.' '.and.c.ne.'=') goto 20
!
! Store details for this string
  30  nsub = nsub + 1
      delimit(1,nsub) = j
      delimit(2,nsub) = k - 1
!
      if (k.lt.len) then
        j = k
        goto 10
      end if
!
  99  continue
!
!------------------------------------------------------------------------------
!
      return
      end
!

