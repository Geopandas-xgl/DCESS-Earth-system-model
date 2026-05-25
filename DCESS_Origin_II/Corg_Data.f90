!SHARED DATA IN CORR[1-7]W7A_EQ2_KIND FUNTIONS
MODULE CORG_DATA
	use parameters, only: wp

	!Number Parameters
	REAL(wp), PARAMETER:: zro=0.0_wp,half=0.5_wp,one=1.0_wp,two=2.0_wp,four=4.0_wp
  !Input Shared Parameter of Corg[1-7]w7a_eq2_kind functions
	REAL(wp) :: w,ddb,jo2,beta,ffcc,fi(7),ff(7),del(7),iinioo2,zb(8)
  !Output Shared Parameters
	REAL(wp) :: s(28),rrm(8)
	REAL(wp) :: cc1,cc2,cc3,cc4,cc5,cc6,cc7
	REAL(wp) :: dd1,dd2,dd3,dd4,dd5,dd6,dd7
	REAL(wp) :: x1(16),x2(16),x3(16),x4(16),x5(16),x6(16),x7(16)
 	!TEMPORAL DATA FOR LACPACK LIBRARY 
	!MATRIX in BANDED FORMAT
	REAL(wp) :: a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15
	REAL(wp) :: a16,a17,a18,a19,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29,a30
	REAL(wp) :: a31,a32,a33,a34,a35,a36,a37,a38,a39,a40,a41,a42,a43,a44,a45
	REAL(wp) :: a46,a47,a48,a49,a50,a51,a52,a53,a54,a55,a56,a57,a58,a59,a60
	REAL(wp) :: bandedaa(7,16)
	INTEGER  :: info,ipiv(16)
	REAL(wp) :: rrm1,rrm2,rrm3,rrm4,rrm5,rrm6,rrm7,rrm8

END MODULE CORG_DATA
