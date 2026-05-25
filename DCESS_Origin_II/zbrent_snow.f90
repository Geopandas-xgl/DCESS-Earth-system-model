FUNCTION zbrent_snow(func,x1,x2,tol,par1,par2,par3)
use parameters, only: wp
IMPLICIT NONE

REAL(wp), INTENT(IN) :: x1,x2,tol,par1,par2,par3
REAL(wp) :: zbrent_snow
INTERFACE
	FUNCTION func(x,c,d,e)
		use parameters, only: wp
		IMPLICIT NONE
		REAL(wp), INTENT(IN) :: x,c,d,e
		REAL(wp) :: func
	END FUNCTION func
END INTERFACE
INTEGER, PARAMETER :: ITMAX=100
REAL(wp), PARAMETER :: EPS=epsilon(x1)
!Using Brent’s method, find the root of a function func known to lie between x1 and x2 .
!The root, returned as zbrent , will be refined until its accuracy is tol .
!Parameters: Maximum allowed number of iterations, and machine floating-point precision.
INTEGER :: iter
REAL(wp) :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
a=x1
b=x2
fa=func(a,par1,par2,par3)
fb=func(b,par1,par2,par3)
if ((fa > 0.0 .and. fb > 0.0) .or. (fa < 0.0 .and. fb < 0.0)) then
	WRITE(*,*) "root must be bracketed for zbrent"
	zbrent_snow=0.0
	RETURN 
end if
c=b
fc=fb
do iter=1,ITMAX
	if ((fb > 0.0 .and. fc > 0.0) .or. (fb < 0.0 .and. fc < 0.0)) then
		c=a
		!Rename a, b, c and adjust bounding interval d.
		fc=fa
		d=b-a
		e=d
	end if
	if (abs(fc) < abs(fb)) then
		a=b
		b=c
		c=a
		fa=fb
		fb=fc
		fc=fa
	end if
	tol1=2.0_wp*EPS*abs(b)+0.5_wp*tol
	!Convergence check.
	xm=0.5_wp*(c-b)
	if (abs(xm) <= tol1 .or. fb == 0.0) then
		zbrent_snow=b
		RETURN
	end if
	if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
	s=fb/fa
	!Attempt inverse quadratic interpolation.
	if (a == c) then
		p=2.0_wp*xm*s
		q=1.0_wp-s
	else
		q=fa/fc
		r=fb/fc
		p=s*(2.0_wp*xm*q*(q-r)-(b-a)*(r-1.0_wp))
		q=(q-1.0_wp)*(r-1.0_wp)*(s-1.0_wp)
	end if
	if (p > 0.0) q=-q
		!Check whether in bounds.
		p=abs(p)
		if (2.0_wp*p < min(3.0_wp*xm*q-abs(tol1*q),abs(e*q))) then
			e=d
			!Accept interpolation.
			d=p/q
		else
			d=xm
			!Interpolation failed; use bisection.
			e=d
		end if
	else
		!Bounds decreasing too slowly; use bisection.
		d=xm
		e=d
	end if
	a=b
	!Move last best guess to a.
	fa=fb
	b=b+merge(d,sign(tol1,xm), abs(d) > tol1 )
	!Evaluate new trial root.
	fb=func(b,par1,par2,par3)
end do
WRITE(*,*) "zbrent: exceeded maximum iterations"
zbrent_snow=b
END FUNCTION zbrent_snow
