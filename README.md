# reinsch_spline
Fortran code that fits a third-order spline through noisy data; based on Reinsch (1967)

## Copy of code header
```
c version 1.1 - May 23, 2003
c version 1.0 - September 24, 2002

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
```
