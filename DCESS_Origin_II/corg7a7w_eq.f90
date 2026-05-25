FUNCTION corg7a7w_eq(zo) result(res)
! USE PARVAL_N2O, ONLY: kr,kii
use parameters, only: wp
USE CORG_DATA

IMPLICIT NONE
REAL(wp), intent(in) :: zo
REAL(wp)             :: res
!---------------------------------------------------------------------------------------------

  a1=(w-ddb*s(1))*(1.0_wp-fi(1))
  a2=(w-ddb*s(2))*(1.0_wp-fi(1))
  a3=exp(s(1)*zb(2))
  a4=exp(s(2)*zb(2))
  a5=-exp(s(3)*zb(2))
  a6=-exp(s(4)*zb(2))
  a7=(w-ddb*s(1))*exp(s(1)*zb(2))*(1.0_wp-fi(1))
  a8=(w-ddb*s(2))*exp(s(2)*zb(2))*(1.0_wp-fi(1))
  a9=-(w-ddb*s(3))*exp(s(3)*zb(2))*(1.0_wp-fi(2))
  a10=-(w-ddb*s(4))*exp(s(4)*zb(2))*(1.0_wp-fi(2))
  a11=exp(s(3)*zb(3))
  a12=exp(s(4)*zb(3))
  a13=-exp(s(5)*zb(3))
  a14=-exp(s(6)*zb(3))
  a15=(w-ddb*s(3))*exp(s(3)*zb(3))*(1.0_wp-fi(2))
  a16=(w-ddb*s(4))*exp(s(4)*zb(3))*(1.0_wp-fi(2))
  a17=-(w-ddb*s(5))*exp(s(5)*zb(3))*(1.0_wp-fi(3))
  a18=-(w-ddb*s(6))*exp(s(6)*zb(3))*(1.0_wp-fi(3))
  a19=exp(s(5)*zb(4))
  a20=exp(s(6)*zb(4))
  a21=-exp(s(7)*zb(4))
  a22=-exp(s(8)*zb(4))
  a23=(w-ddb*s(5))*exp(s(5)*zb(4))*(1.0_wp-fi(3))
  a24=(w-ddb*s(6))*exp(s(6)*zb(4))*(1.0_wp-fi(3))
  a25=-(w-ddb*s(7))*exp(s(7)*zb(4))*(1.0_wp-fi(4))
  a26=-(w-ddb*s(8))*exp(s(8)*zb(4))*(1.0_wp-fi(4))
  a27=exp(s(7)*zb(5))
  a28=exp(s(8)*zb(5))
  a29=-exp(s(9)*zb(5))
  a30=-exp(s(10)*zb(5))
  a31=(w-ddb*s(7))*exp(s(7)*zb(5))*(1.0_wp-fi(4))
  a32=(w-ddb*s(8))*exp(s(8)*zb(5))*(1.0_wp-fi(4))
  a33=-(w-ddb*s(9))*exp(s(9)*zb(5))*(1.0_wp-fi(5))
  a34=-(w-ddb*s(10))*exp(s(10)*zb(5))*(1.0_wp-fi(5))
  a35=exp(s(9)*zb(6))
  a36=exp(s(10)*zb(6))
  a37=-exp(s(25)*zb(6))
  a38=-exp(s(26)*zb(6))
  a39=(w-ddb*s(9))*exp(s(9)*zb(6))*(1.0_wp-fi(5))
  a40=(w-ddb*s(10))*exp(s(10)*zb(6))*(1.0_wp-fi(5))
  a41=-(w-ddb*s(25))*exp(s(25)*zb(6))*(1.0_wp-fi(6))
  a42=-(w-ddb*s(26))*exp(s(26)*zb(6))*(1.0_wp-fi(6))
  a43=exp(s(25)*zb(7))
  a44=exp(s(26)*zb(7))
  a45=-exp(s(27)*zb(7))
  a46=-exp(s(28)*zb(7))
  a47=-(w-ddb*s(25))*exp(s(25)*zb(7))*(1.0_wp-fi(6))
  a48=-(w-ddb*s(26))*exp(s(26)*zb(7))*(1.0_wp-fi(6))
  a49=(w-ddb*s(27))*exp(s(27)*zb(7))*(1.0_wp-fi(7))
  a50=(w-ddb*s(28))*exp(s(28)*zb(7))*(1.0_wp-fi(7))
  a51=exp(s(27)*zo)
  a52=exp(s(28)*zo)
  a53=-exp(s(23)*zo)
  a54=-exp(s(24)*zo)
  a55=-(w-ddb*s(27))*exp(s(27)*zo)
  a56=-(w-ddb*s(28))*exp(s(28)*zo)
  a57=(w-ddb*s(23))*exp(s(23)*zo)
  a58=(w-ddb*s(24))*exp(s(24)*zo)
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

  x7 = (/ffcc,zro,zro,zro,zro,zro,zro,zro,zro,zro,zro,zro,zro,zro,zro,zro/)

  call DGBSV(16,2,2,1,bandedaa,7,ipiv,x7,16,info)

 rrm1=(1._wp-fi(1))*jo2*(x7(1)*(exp(s(1)*zb(2))-1)/s(1)+&
                         x7(2)*(exp(s(2)*zb(2))-1)/s(2))
 rrm2=(1._wp-fi(2))*jo2*(x7(3)*(exp(s(3)*zb(3))-exp(s(3)*zb(2)))/s(3)+&
                     x7(4)*(exp(s(4)*zb(3))-exp(s(4)*zb(2)))/s(4))
 rrm3=(1._wp-fi(3))*jo2*(x7(5)*(exp(s(5)*zb(4))-exp(s(5)*zb(3)))/s(5)+&
                    x7(6)*(exp(s(6)*zb(4))-exp(s(6)*zb(3)))/s(6))
 rrm4=(1._wp-fi(4))*jo2*(x7(7)*(exp(s(7)*zb(5))-exp(s(7)*zb(4)))/s(7)+&
                     x7(8)*(exp(s(8)*zb(5))-exp(s(8)*zb(4)))/s(8))
 rrm5=(1._wp-fi(5))*jo2*(x7(9)*(exp(s(9)*zb(6))-exp(s(9)*zb(5)))/s(9)+&
                     x7(10)*(exp(s(10)*zb(6))-exp(s(10)*zb(5)))/s(10))
 rrm6=(1._wp-fi(6))*jo2*(x7(11)*(exp(s(25)*zb(7))-exp(s(25)*zb(6)))/s(25)+&
                     x7(12)*(exp(s(26)*zb(7))-exp(s(26)*zb(6)))/s(26))
 rrm7=(1._wp-fi(7))*jo2*(x7(13)*(exp(s(27)*zo)-exp(s(27)*zb(7)))/s(27)+&
                     x7(14)*(exp(s(28)*zo)-exp(s(28)*zb(7)))/s(28)+&
                     beta*(x7(15)*(exp(s(23)*zb(8))-exp(s(23)*zo))/s(23)+&
                     x7(16)*(exp(s(24)*zb(8))-exp(s(24)*zo))/s(24)))
 rrm8=(1._wp-fi(7))*jo2*(beta*(x7(15)*(exp(s(23)*zb(8))-exp(s(23)*zo))/s(23)+&
                     x7(16)*(exp(s(24)*zb(8))-exp(s(24)*zo))/s(24)))
 rrm = (/rrm1, rrm2, rrm3, rrm4, rrm5, rrm6, rrm7, rrm8/)

  cc1 = iinioo2 - del(1)*(x7(1)/(s(1)**2)+x7(2)/(s(2)**2))

  dd7 = -del(7)*((x7(13)/s(27)*exp(s(27)*zo)+x7(14)/s(28)*exp(s(28)*zo)&
   +beta*(x7(15)*(exp(s(23)*zb(8))-exp(s(23)*zo))/s(23)+x7(16)*(exp(s(24)*zb(8))-exp(s(24)*zo))/s(24))))

  dd6 = fi(7)*ff(6)/(fi(6)*ff(7))*(dd7 + del(7)*(x7(13)*exp(s(27)*zb(7))/s(27)+x7(14)*exp(s(28)*zb(7))/s(28)))-&
      del(6)*(x7(11)*exp(s(25)*zb(7))/s(9)+x7(12)*exp(s(26)*zb(7))/s(26))

  dd5 = fi(6)*ff(5)/(fi(5)*ff(6))*(dd6 + del(6)*(x7(11)*exp(s(25)*zb(6))/s(25)+x7(12)*exp(s(26)*zb(6))/s(26)))-&
      del(5)*(x7(9)*exp(s(9)*zb(6))/s(9)+x7(10)*exp(s(10)*zb(6))/s(10))

  dd4 = fi(5)*ff(4)/(fi(4)*ff(5))*(dd5 + del(5)*(x7(9)*exp(s(9)*zb(5))/s(9)+x7(10)*exp(s(10)*zb(5))/s(10)))-&
      del(4)*(x7(7)*exp(s(7)*zb(5))/s(7)+x7(8)*exp(s(8)*zb(5))/s(8))

  dd3 = fi(4)*ff(3)/(fi(3)*ff(4))*(dd4 + del(4)*(x7(7)*exp(s(7)*zb(4))/s(7)+x7(8)*exp(s(8)*zb(4))/s(8)))-&
      del(3)*(x7(5)*exp(s(5)*zb(4))/s(5)+x7(6)*exp(s(6)*zb(4))/s(6))

  dd2 = fi(3)*ff(2)/(fi(2)*ff(3))*(dd3 + del(3)*(x7(5)*exp(s(5)*zb(3))/s(5)+x7(6)*exp(s(6)*zb(3))/s(6)))-&
      del(2)*(x7(3)*exp(s(3)*zb(3))/s(3)+x7(4)*exp(s(4)*zb(3))/s(4))

  dd1 = fi(2)*ff(1)/(fi(1)*ff(2))*(dd2 + del(2)*(x7(3)*exp(s(3)*zb(2))/s(3)+x7(4)*exp(s(4)*zb(2))/s(4)))-&
      del(1)*(x7(1)*exp(s(1)*zb(2))/s(1)+x7(2)*exp(s(2)*zb(2))/s(2))

  cc2 = cc1 - (dd2-dd1)*zb(2)-&
      del(2)*(x7(3)*exp(s(3)*zb(2))/(s(3)**2)+x7(4)*exp(s(4)*zb(2))/(s(4)**2))+&
      del(1)*(x7(1)*exp(s(1)*zb(2))/(s(1)**2)+x7(2)*exp(s(2)*zb(2))/(s(2)**2))

  cc3 = cc2 - (dd3-dd2)*zb(3)-&
      del(3)*(x7(5)*exp(s(5)*zb(3))/(s(5)**2)+x7(6)*exp(s(6)*zb(3))/(s(6)**2))+&
      del(2)*(x7(3)*exp(s(3)*zb(3))/(s(3)**2)+x7(4)*exp(s(4)*zb(3))/(s(4)**2))

  cc4 =  cc3 - (dd4-dd3)*zb(4)-&
      del(4)*(x7(7)*exp(s(7)*zb(4))/(s(7)**2)+x7(8)*exp(s(8)*zb(4))/(s(8)**2))+&
      del(3)*(x7(5)*exp(s(5)*zb(4))/(s(5)**2)+x7(6)*exp(s(6)*zb(4))/(s(6)**2))

  cc5 =  cc4 - (dd5-dd4)*zb(5)-&
      del(5)*(x7(9)*exp(s(9)*zb(5))/(s(9)**2)+x7(10)*exp(s(10)*zb(5))/(s(10)**2))+&
      del(4)*(x7(7)*exp(s(7)*zb(5))/(s(7)**2)+x7(8)*exp(s(8)*zb(5))/(s(8)**2))

  cc6 =  cc5 - (dd6-dd5)*zb(6)-&
      del(6)*(x7(11)*exp(s(25)*zb(6))/(s(25)**2)+x7(12)*exp(s(26)*zb(6))/(s(26)**2))+&
      del(5)*(x7(9)*exp(s(9)*zb(6))/(s(9)**2)+x7(10)*exp(s(10)*zb(6))/(s(10)**2))

  cc7 =  cc6 - (dd7-dd6)*zb(7)-&
      del(7)*(x7(13)*exp(s(27)*zb(7))/(s(27)**2)+x7(14)*exp(s(28)*zb(7))/(s(28)**2))+&
      del(6)*(x7(11)*exp(s(25)*zb(7))/(s(25)**2)+x7(12)*exp(s(26)*zb(7))/(s(26)**2))


	res = -(cc7 + dd7*zo + del(7)*(x7(13)*exp(s(27)*zo)/(s(27)**2)+x7(14)*exp(s(28)*zo)/(s(28)**2)))

END FUNCTION corg7a7w_eq
