FUNCTION corg4a7w_eq(zo) result(res)
! USE PARVAL_N2O, ONLY: kr,kii
use parameters, only: wp
USE CORG_DATA

IMPLICIT NONE
REAL(wp), intent(in) :: zo
REAL(wp)             :: res
!REAL(wp)::aa(16,16),a(69)
!INTEGER::i,j
!------------------------------------------------------------------------------------------------

  a1 =(w-ddb*s(1))*(1._wp-fi(1))
  a2 =(w-ddb*s(2))*(1._wp-fi(1))
  a3 =exp(s(1)*zb(2))
  a4 =exp(s(2)*zb(2))
  a5 =-a3
  a6 =-a4
  a7 =(w-ddb*s(1))*a3*(1.0_wp-fi(1))
  a8 =(w-ddb*s(2))*a4*(1.0_wp-fi(1))
  a9 =(w-ddb*s(3))*a5*(1.0_wp-fi(2))
  a10=(w-ddb*s(4))*a6*(1.0_wp-fi(2))
  a11=exp(s(3)*zb(3))
  a12=exp(s(4)*zb(3))
  a13=-a11
  a14=-a12
  a15=(w-ddb*s(3))*a11*(1.0_wp-fi(2))
  a16=(w-ddb*s(4))*a12*(1.0_wp-fi(2))
  a17=(w-ddb*s(5))*a13*(1.0_wp-fi(3))
  a18=(w-ddb*s(6))*a14*(1.0_wp-fi(3))
  a19=exp(s(5)*zb(4))
  a20=exp(s(6)*zb(4))
  a21=-a19
  a22=-a20
  a23=(w-ddb*s(5))*a19*(1.0_wp-fi(3))
  a24=(w-ddb*s(6))*a20*(1.0_wp-fi(3))
  a25=(w-ddb*s(7))*a21*(1.0_wp-fi(4))
  a26=(w-ddb*s(8))*a22*(1.0_wp-fi(4))
  a27=exp(s(7)*zo)
  a28=exp(s(8)*zo)
  a29=-exp(s(17)*zo)
  a30=-exp(s(18)*zo)
  a31=(w-ddb*s(7))*a27
  a32=(w-ddb*s(8))*a28
  a33=(w-ddb*s(17))*a29
  a34=(w-ddb*s(18))*a30
  a35=exp(s(17)*zb(5))
  a36=exp(s(18)*zb(5))
  a37=-a35
  a38=-a36
  a39=(w-ddb*s(17))*a35*(1.0_wp-fi(4))
  a40=(w-ddb*s(18))*a36*(1.0_wp-fi(4))
  a41=(w-ddb*s(19))*a37*(1.0_wp-fi(5))
  a42=(w-ddb*s(20))*a38*(1.0_wp-fi(5))
  a43=exp(s(19)*zb(6))
  a44=exp(s(20)*zb(6))
  a45=-a43
  a46=-a44
  a47=(w-ddb*s(19))*a45*(1.0_wp-fi(5))
  a48=(w-ddb*s(20))*a46*(1.0_wp-fi(5))
  a49=(w-ddb*s(21))*a43*(1.0_wp-fi(6))
  a50=(w-ddb*s(20))*a44*(1.0_wp-fi(6))
  a51=exp(s(21)*zb(7))
  a52=exp(s(22)*zb(7))
  a53=-a51
  a54=-a52
  a55=(w-ddb*s(21))*a53*(1.0_wp-fi(6))
  a56=(w-ddb*s(22))*a54*(1.0_wp-fi(6))
  a57=(w-ddb*s(23))*a51*(1.0_wp-fi(7))
  a58=(w-ddb*s(24))*a52*(1.0_wp-fi(7))
  a59=(ddb*s(23))*exp(s(23)*zb(8))
  a60=(ddb*s(24))*exp(s(24)*zb(8))

  bandedaa=transpose(reshape(&
     (/zro,zro,zro,zro,zro,zro,zro,zro,zro,zro,zro,zro,zro,zro,zro,zro,&
       zro,zro,zro,zro,zro,zro,zro,zro,zro,zro,zro,zro,zro,zro,zro,zro,&
       zro,zro,zro,a6 ,zro,a14,zro,a22,zro,a30,zro,a38,zro,a46,zro,a54,&
       zro,a2 ,a5 ,a10,a13,a18,a21,a26,a29,a34,a37,a42,a45,a50,a53,a58,&
       a1 ,a4 ,a9 ,a12,a17,a20,a25,a28,a33,a36,a41,a44,a49,a52,a57,a60,&
       a3 ,a8 ,a11,a16,a19,a24,a27,a32,a35,a40,a43,a48,a51,a56,a59,zro,&
       a7 ,zro,a15,zro,a23,zro,a31,zro,a39,zro,a47,zro,a55,zro,zro,zro/),(/16,7/)))

  x4 = (/ffcc,zro,zro,zro,zro,zro,zro,zro,zro,zro,zro,zro,zro,zro,zro,zro/)
	!aa matrix
	!DO i=1,16
	!	DO j=1,16
	!		aa(i,j)=0.0_kr
	!	END DO
	!END DO
	!aa(1,1:2)=a(1:2)
	!DO i=2,14,2
	!		aa(i  ,i-1:i+2)=a(4*i-5:4*i-2)
	!		aa(i+1,i-1:i+2)=a(4*i-1:4*i+2)
	!END DO
	!aa(16,15:16)=a(59:60)
	!call dgesv(16,1,aa,16, ipiv,x4,16,info )
  call DGBSV(16,2,2,1,bandedaa,7,ipiv,x4,16,info)

  rrm1=(1.0_wp-fi(1))*jo2*(x4(1)*(exp(s(1)*zb(2))-1.0_wp)/s(1)+&
                           x4(2)*(exp(s(2)*zb(2))-1.0_wp)/s(2))

  rrm2=(1.0_wp-fi(2))*jo2*(x4(3)*(exp(s(3)*zb(3))-exp(s(3)*zb(2)))/s(3)+&
                           x4(4)*(exp(s(4)*zb(3))-exp(s(4)*zb(2)))/s(4))

  rrm3=(1.0_wp-fi(3))*jo2*(x4(5)*(exp(s(5)*zb(4))-exp(s(5)*zb(3)))/s(5)+&
                           x4(6)*(exp(s(6)*zb(4))-exp(s(6)*zb(3)))/s(6))

	rrm4=(1.0_wp-fi(4))*jo2*(x4(7)*(exp(s(7)*zo)-exp(s(7)*zb(4)))/s(7)+&
                           x4(8)*(exp(s(8)*zo)-exp(s(8)*zb(4)))/s(8)+&
                     beta*(x4(9)*(exp(s(17)*zb(5))-exp(s(17)*zo))/s(17)+&
                          x4(10)*(exp(s(18)*zb(5))-exp(s(18)*zo))/s(18)))

	rrm5=(1.0_wp-fi(5))*beta*jo2*(x4(11)*(exp(s(19)*zb(6))-exp(s(19)*zb(5)))/s(19)+&
                               x4(12)*(exp(s(20)*zb(6))-exp(s(20)*zb(5)))/s(20))
	rrm6=(1.0_wp-fi(6))*beta*jo2*(x4(13)*(exp(s(21)*zb(7))-exp(s(21)*zb(6)))/s(21)+&
                               x4(14)*(exp(s(22)*zb(7))-exp(s(22)*zb(6)))/s(22))
  rrm7=(1.0_wp-fi(7))*beta*jo2*(x4(15)*(exp(s(23)*zb(8))-exp(s(23)*zb(7)))/s(23)+&
                               x4(16)*(exp(s(24)*zb(8))-exp(s(24)*zb(7)))/s(24))

	rrm8=(1.0_wp-fi(4))*jo2*(beta*(x4(9)*(exp(s(17)*zb(5))-exp(s(17)*zo))/s(17)+&
                            x4(10)*(exp(s(18)*zb(5))-exp(s(18)*zo))/s(18)))+&
                            rrm5+rrm6+rrm7

  rrm=(/ rrm1, rrm2, rrm3, rrm4, rrm5, rrm6, rrm7, rrm8 /)

  cc1 = iinioo2 - del(1)*(x4(1)/(s(1)**2)+x4(2)/(s(2)**2));

  dd4 = - del(4)*((x4(7)/s(7)*exp(s(7)*zo)+x4(8)/s(8)*exp(s(8)*zo)&
 +beta*(x4(9)*(exp(s(17)*zb(5))-exp(s(17)*zo))/s(17)+x4(10)*(exp(s(18)*zb(5))-exp(s(18)*zo))/s(18)&
 +del(5)/del(4)*(x4(11)*(exp(s(19)*zb(6))-exp(s(19)*zb(5)))/s(19)+x4(12)*(exp(s(20)*zb(6))-exp(s(20)*zb(5)))/s(20))&
 +del(6)/del(4)*(x4(13)*(exp(s(21)*zb(7))-exp(s(21)*zb(6)))/s(21)+x4(14)*(exp(s(22)*zb(7))-exp(s(22)*zb(6)))/s(22))&
 +del(7)/del(4)*(x4(15)*(exp(s(23)*zb(8))-exp(s(23)*zb(7)))/s(23)+x4(16)*(exp(s(24)*zb(8))-exp(s(24)*zb(7)))/s(24)))))

  dd3 = fi(4)*ff(3)/(fi(3)*ff(4))*(dd4 + del(4)*(x4(7)*exp(s(7)*zb(4))/s(7)+x4(8)*exp(s(8)*zb(4))/s(8)))-&
      del(3)*(x4(5)*exp(s(5)*zb(4))/s(5)+x4(6)*exp(s(6)*zb(4))/s(6))

  dd2 = fi(3)*ff(2)/(fi(2)*ff(3))*(dd3 + del(3)*(x4(5)*exp(s(5)*zb(3))/s(5)+x4(6)*exp(s(6)*zb(3))/s(6)))-&
      del(2)*(x4(3)*exp(s(3)*zb(3))/s(3)+x4(4)*exp(s(4)*zb(3))/s(4))

  dd1 = fi(2)*ff(1)/(fi(1)*ff(2))*(dd2 + del(2)*(x4(3)*exp(s(3)*zb(2))/s(3)+x4(4)*exp(s(4)*zb(2))/s(4)))-&
      del(1)*(x4(1)*exp(s(1)*zb(2))/s(1)+x4(2)*exp(s(2)*zb(2))/s(2))


  cc2 = cc1 - (dd2-dd1)*zb(2)-&
      del(2)*(x4(3)*exp(s(3)*zb(2))/(s(3)**2)+x4(4)*exp(s(4)*zb(2))/(s(4)**2))+&
      del(1)*(x4(1)*exp(s(1)*zb(2))/(s(1)**2)+x4(2)*exp(s(2)*zb(2))/(s(2)**2))

  cc3 = cc2 - (dd3-dd2)*zb(3)-&
      del(3)*(x4(5)*exp(s(5)*zb(3))/(s(5)**2)+x4(6)*exp(s(6)*zb(3))/(s(6)**2))+&
      del(2)*(x4(3)*exp(s(3)*zb(3))/(s(3)**2)+x4(4)*exp(s(4)*zb(3))/(s(4)**2))

  cc4 = cc3 - (dd4-dd3)*zb(4)-&
      del(4)*(x4(7)*exp(s(7)*zb(4))/(s(7)**2)+x4(8)*exp(s(8)*zb(4))/(s(8)**2))+&
      del(3)*(x4(5)*exp(s(5)*zb(4))/(s(5)**2)+x4(6)*exp(s(6)*zb(4))/(s(6)**2))

!  rres4 = -(cc4 + dd4*zo + del(4)*(x4(7)*exp(s(7)*zo)/(s(7)**2)+x4(8)*exp(s(8)*zo)/(s(8)**2)))
  res = -(cc4 + dd4*zo + del(4)*(x4(7)*exp(s(7)*zo)/(s(7)**2)+x4(8)*exp(s(8)*zo)/(s(8)**2)))

END FUNCTION corg4a7w_eq
