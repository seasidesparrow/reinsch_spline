      program reinsch
c version 1.1 - May 23, 2003
c version 1.0 - September 24, 2002

c---------------------------------------------------------------------------
c Reinsch.f -- a program to compute a smoothed spline through noisy data.
c Algorithm adopted from the ALGOL Language code in Numerische Mathematik
c volume 10, pp177-183 (1967) by Christian Reinsch.  This paper is often
c cited as a standard in fitting splines/smooth curves through noisy
c data.
c
c The program searches for a best-fitting third-order polynomial
c between the data points such that the global "least-squares" value falls
c below the user-set value of "s".  The coefficients of the polynomial are 
c calculated by the subroutine "reinschspline" and they are stored
c in the arrays a,b,c, and d.  There is one set of polynomial coefficients
c for each data point.  You model the data (see the makefit subroutine)
c by:
c
c	1) find the observed data points between which your desired
c	   interpolation point lies
c
c	2) Take the polynomial coefficients for the lower boundary of
c	   the interval
c
c	3) Compute the fitted y-value using the equations:
c
c	    h = xfit-xobs(i)
c
c	    yfit = f(xfit) = a(i) + (b(i)*h) + (c(i)*h*h) + (d(i)*h*h*h)
c
c The key parameter is the error on the data.  If the errors are not known
c or well estimated, the code will interpret every little wiggle as being
c real.
c
c As an example, see the work on R Sct (Buchler et al, ApJ 462, 489 (1996))
c
c The important arrays to read in are: time, mag, and (if known) error.
c In the reinschspline routine, time is "x", mag is "y", and error is "dy".
c Note that this version only reads in two columns of data in free, ascii
c format, and that the errors are set to the same value.  If the error per
c point is known and included within the input file, they can be read in
c by adding "error(i)" to the read statement in the readdat subroutine.
c
c Also, note that the maximum array sizes are set to 200000, and you are
c limited to reading in 200000 points.  This is trivially easy to change
c (just change 200000 to whatever you want and recompile), but make sure
c you change *every instance* of 200000 (and 200001) to the new value
c (and the new value+1).
c
c Writepoly is a tiny subroutine that prints the polynomial coefficients
c a,b,c, and d.  Uncomment the call statement to obtain this debugging
c information, otherwise, leave it be.
c
c Finally, DO NOT COMPILE WITH -O# OPTIMIZATION USING g77!  This produced
c a bug which made a bad fit.  Optimization really isn't necessary anyway.
c The majority of the time is taken by the makefit routine, not reinschspline.
c
c				-Matthew Templeton
c				 September 24, 2002
c
c Revised to automatically set "s" to (npts - sqrt(2*npts)) and tidy up
c the screen and file output.
c
c				-Matthew Templeton
c				 May 23, 2003
c---------------------------------------------------------------------------


c start of main program
      real*8 time(200000),mag(200000),error(200000)
      real*8 a(0:200001),b(0:200001),c(0:200001),d(0:200001)
      real*8 tout(200000),mout(200000)
      real*8 dt,dtlast,fnnew,t1new,tfnew,hh,dtin
      real*8 s
      integer npts,i

      write(6,100)
100   format(/,/,/,24h                 Reinsch,/,/,
     1       41hA program to compute a best-fitting curve,/,
     1       43hfor a set of data points with known errors.,/,/,
     1       41hTwo-column (time mag) data format assumed,/,/,/)

      call init(time,mag,error,npts)
      call queryf
      call readdat(time,mag,error,dtlast,npts)
      call setfit(s,error,npts)
      call queryt(time,t1new,tfnew,dtlast,dtin,nnew,npts)
      call maket(t1new,tout,mout,dtin,nnew)
      call reinschspline(time,mag,error,1,npts,s,a,b,c,d)
c     call writepoly(a,b,c,d,npts)
      call makefit(time,tout,mout,a,b,c,d,nnew,npts,error)

c end of main program
      end

c------------------------------------------------------------------------

      subroutine reinschspline(x,y,dy,n1,n2,s,a,b,c,d)
      
      real*8 x(200000),y(200000),dy(200000)
      integer n1,n2
      real*8 s
      integer i,m1,m2
      real*8 e,f,f2,g,h,p
      real*8 v(0:200001),r(0:200001),r1(0:200001),r2(0:200001)
      real*8 t(0:200001),t1(0:200001),u(0:200001)
      real*8 a(0:200001),b(0:200001),c(0:200001),d(0:200001)
      
      e=0.d0

      do i=0,200001
       v(i)=0.d0
       r(i)=0.d0
       r1(i)=0.d0
       r2(i)=0.d0
       t(i)=0.d0
       t1(i)=0.d0
       u(i)=0.d0
       a(i)=0.d0
       b(i)=0.d0
       c(i)=0.d0
       d(i)=0.d0
      enddo

      m1=n1-1
      m2=n2+1
      r(m1)=0.d0
      r(n1)=0.d0
      r1(n2)=0.d0
      r2(n2)=0.d0
      r2(m2)=0.d0
      u(m1)=0.d0
      u(n1)=0.d0
      u(n2)=0.d0
      u(m2)=0.d0
      p=0.d0
      m1=n1+1
      m2=n2-1
      h=x(m1)-x(n1)
      f=(y(m1)-y(n1))/h

      do i=m1,m2
       g=h
       h=x(i+1)-x(i)
       e=f
       f=(y(i+1)-y(i))/h
       a(i)=f-e
       t(i)=2.d0*(g+h)/3.d0
       t1(i)=h/3.d0
       r2(i)=dy(i-1)/g
       r(i)=dy(i+1)/h
       r1(i)=(-1.d0*dy(i)/g) - (dy(i)/h)
      enddo

      do i=m1,m2
       b(i)=(r(i)*r(i))+(r1(i)*r1(i))+(r2(i)*r2(i))
       c(i)=(r(i)*r1(i+1))+(r1(i)*r2(i+1))
       d(i)=(r(i)*r2(i+2))
      enddo

      f2=-1.d0*s

      do j=1,100
       do i=m1,m2
        r1(i-1)=f*r(i-1)
        r2(i-2)=g*r(i-2)
        r(i)=1.d0/((p*b(i))+t(i)-(f*r1(i-1))-(g*r2(i-2)))
        u(i)=a(i)-(r1(i-1)*u(i-1))-(r2(i-2)*u(i-2))
        f=(p*c(i))+t1(i)-(h*r1(i-1))
        g=h
        h=d(i)*p
       enddo
      
       do i=m2,m1,-1
        u(i)=(r(i)*u(i))-(r1(i)*u(i+1))-(r2(i)*u(i+2))
        e=0.d0
        h=0.d0
       enddo
 
       do i=n1,m2
        g=h
        h=(u(i+1)-u(i))/(x(i+1)-x(i))
        v(i)=(h-g)*dy(i)*dy(i)
        e=e+(v(i)*(h-g))
       enddo
 
       g=-1.d0*h*dy(n2)*dy(n2)
       v(n2)=g
       e=e-(g*h)
       g=f2
       f2=e*p*p
       if(f2.ge.s.and.f2.le.g) then
        do i=n1,n2
         a(i)=y(i)-(p*v(i))
         c(i)=u(i)
        enddo
        do i=n1,n2
         h=x(i+1)-x(i)
         d(i)=(c(i+1)-c(i))/(3.d0*h)
         b(i)=((a(i+1)-a(i))/h) - ((h*d(i))+c(i))*h
        enddo
        return
       endif
       f=0.d0
       h=(v(m1)-v(n1))/(x(m1)-x(n1))
       do i=m1,m2
        g=h
        h=(v(i+1)-v(i))/(x(i+1)-x(i))
        g=h-g-(r1(i-1)*r(i-1))-(r2(i-2)*r(i-2))
        f=f+(g*r(i)*g)
        r(i)=g
       enddo
       h=e-(p*f)
       if(h.lt.0.d0) then
        do i=n1,n2
         a(i)=y(i)-(p*v(i))
         c(i)=u(i)
        enddo
        do i=n1,n2
         h=x(i+1)-x(i)
         d(i)=(c(i+1)-c(i))/(3.d0*h)
         b(i)=((a(i+1)-a(i))/h) - ((h*d(i))+c(i))*h
        enddo
        return
       endif
       p=p+((s-f2)/((dsqrt(dabs(s/e))+p)*h))
c      print*,p,h,f2,g,s
      enddo
      write(6,909)
909   format(61hThis fit attempt reached 100 iterations, and is probably
     1 bad.,/,13hStopping now.,/,/)
      stop
      return
      end

c--------------------------------------------------------------------------

      subroutine queryt(time,t1new,tfnew,dtlast,dtin,nnew,npts)
      real*8 time(200000)
      real*8 t1new,tfnew,dtin,dtlast
      integer nnew,npts
      real*8 fnnew

      write(6,*),'The start and end times of your data are:'
      write(6,*),'Start:', time(1)
      write(6,*),'End:', time(npts)
      write(6,*),'Minimum data spacing is:', dtlast
      write(6,*),'Choose a new data spacing:'
      read*,dtin
      t1new=time(1)+(dtin/2.d0)
      tfnew=time(npts)-(dtin/2.d0)

      fnnew=(tfnew-t1new)/dtin
      nnew=idint(fnnew)+1
      return
      end

c--------------------------------------------------------------------------

      subroutine maket(t1new,tout,mout,dtin,nnew)
      real*8 t1new,tout(200000),mout(200000)
      real*8 dtin
      integer nnew,i

      tout(1)=t1new
      mout(1)=0.d0
      do i=2,nnew
       tout(i)=tout(i-1)+dtin
       mout(i)=0.d0 
      enddo
      return
      end

c--------------------------------------------------------------------------

      subroutine makefit(time,tout,mout,a,b,c,d,nnew,npts,error)
      real*8 time(200000),tout(200000),mout(200000)
      real*8 error(200000)
      integer nnew,npts,i,j
      real*8 a(0:200001),b(0:200001),c(0:200001),d(0:200001)

      write(6,606) nnew
606   format(48hFit search complete. Writing best-fit curve now.,/,
     1       16hWriting fit for ,i5,8h points.)

      do i=1,nnew
       do j=2,npts
        if(tout(i).ge.time(j-1).and.tout(i).lt.time(j)) then
         hh=tout(i)-time(j-1)
         mout(i)=a(j-1)+(b(j-1)*hh)+(c(j-1)*hh*hh)+(d(j-1)*hh*hh*hh)
         write(2,404) tout(i),mout(i),error(i)
         go to 808
        endif
       enddo
808    continue
      enddo
      close(2)

c404   format(f13.4,2(2x,f8.3))
404   format(f12.4,2(1x,f5.2))
      return
      end

c------------------------------------------------------------------------
      
      subroutine queryf
      character*80 infile,outfile
      write(6,*) 'Input filename?'
      read*,infile
      write(6,*) 'Output filename?'
      read*,outfile
      open(unit=1,file=infile,status="old")
      open(unit=2,file=outfile,status="unknown")
      return
      end

c------------------------------------------------------------------------

      subroutine init(time,mag,error,npts)
      real*8 time(200000),mag(200000),error(200000)
      integer npts,i
      do i=1,200000
       time(i)=0.d0
       mag(i)=0.d0
       error(i)=0.d0
      enddo
      npts=0
      return
      end

c------------------------------------------------------------------------

      subroutine readdat(time,mag,error,dtlast,npts)
      real*8 time(200000),mag(200000),error(200000)
      real*8 dtlast,dt
      integer npts,i
      dtlast=1.d20
      do i=1,200000
       read(1,*,end=999) time(i),mag(i)
       npts=npts+1
       if(i.gt.1) then
        dt=time(i)-time(i-1)
        if(dt.lt.dtlast) dtlast=dt
       endif
      enddo
999   continue
      close(1)
      return
      end

c------------------------------------------------------------------------
      
      subroutine setfit(s,error,npts)
      real*8 s
      real*8 error(200000)
      s=dfloat(npts)-dsqrt(2.d0*dfloat(npts))
      write(6,*)'What is your estimate of the error per point?'
      read*,ein
      do i=1,200000
       error(i)=ein
      enddo
      return
      end

c------------------------------------------------------------------------

      subroutine writepoly(a,b,c,d,npts)
      real*8 a(0:200001),b(0:200001),c(0:200001),d(0:200001)
      integer npts,i
      do i=1,npts
       write(64,909) a(i),b(i),c(i),d(i)
      enddo
909   format(3(f12.4,1x),f12.4)
      return
      end
