FUNCTION corg1a7w_eq(zo) result(res)
! USE PARVAL_N2O, ONLY: kr,kii
use parameters, only: wp
USE CORG_DATA

IMPLICIT NONE
REAL(wp), intent(in) :: zo
REAL(wp)             :: res
!------------------------------------------------------------------------------------------------------------
  a1=(w-ddb*s(1))*(1._wp-fi(1))
  a2=(w-ddb*s(2))*(1._wp-fi(1))
  a3=exp(s(1)*zo)
  a4=exp(s(2)*zo)
  a5=-exp(s(11)*zo)
  a6=-exp(s(12)*zo)
  a7=-(w-ddb*s(1))*exp(s(1)*zo)
  a8=-(w-ddb*s(2))*exp(s(2)*zo)
  a9=(w-ddb*s(11))*exp(s(11)*zo)
  a10=(w-ddb*s(12))*exp(s(12)*zo)
  a11=exp(s(11)*zb(2))
  a12=exp(s(12)*zb(2))
  a13=-exp(s(13)*zb(2))
  a14=-exp(s(14)*zb(2))
  !a15=-(Db*s(11))*exp(s(11)*zb(2))*(1-fi(1))
  !a16=-(Db*s(12))*exp(s(12)*zb(2))*(1-fi(1))
  !a17=(Db*s(13))*exp(s(13)*zb(2))*(1-fi(2))
  !a18=(Db*s(14))*exp(s(14)*zb(2))*(1-fi(2))
  a15=-(w-ddb*s(11))*exp(s(11)*zb(2))*(1._wp-fi(1))
  a16=-(w-ddb*s(12))*exp(s(12)*zb(2))*(1._wp-fi(1))
  a17=(w-ddb*s(13))*exp(s(13)*zb(2))*(1._wp-fi(2))
  a18=(w-ddb*s(14))*exp(s(14)*zb(2))*(1._wp-fi(2))
  a19=exp(s(13)*zb(3))
  a20=exp(s(14)*zb(3))
  a21=-exp(s(15)*zb(3))
  a22=-exp(s(16)*zb(3))
  a23=-(w-ddb*s(13))*exp(s(13)*zb(3))*(1._wp-fi(2))
  a24=-(w-ddb*s(14))*exp(s(14)*zb(3))*(1._wp-fi(2))
  a25=(w-ddb*s(15))*exp(s(15)*zb(3))*(1._wp-fi(3))
  a26=(w-ddb*s(16))*exp(s(16)*zb(3))*(1._wp-fi(3))
  a27=exp(s(15)*zb(4))
  a28=exp(s(16)*zb(4))
  a29=-exp(s(17)*zb(4))
  a30=-exp(s(18)*zb(4))
  a31=-(w-ddb*s(15))*exp(s(15)*zb(4))*(1._wp-fi(3))
  a32=-(w-ddb*s(16))*exp(s(16)*zb(4))*(1._wp-fi(3))
  a33=(w-ddb*s(17))*exp(s(17)*zb(4))*(1._wp-fi(4))
  a34=(w-ddb*s(18))*exp(s(18)*zb(4))*(1._wp-fi(4))
  a35=exp(s(17)*zb(5))
  a36=exp(s(18)*zb(5))
  a37=-exp(s(19)*zb(5))
  a38=-exp(s(20)*zb(5))
  a39=-(w-ddb*s(17))*exp(s(17)*zb(5))*(1._wp-fi(4))
  a40=-(w-ddb*s(18))*exp(s(18)*zb(5))*(1._wp-fi(4))
  a41=(w-ddb*s(19))*exp(s(19)*zb(5))*(1._wp-fi(5))
  a42=(w-ddb*s(20))*exp(s(20)*zb(5))*(1._wp-fi(5))
  a43=exp(s(19)*zb(6))
  a44=exp(s(20)*zb(6))
  a45=-exp(s(21)*zb(6))
  a46=-exp(s(22)*zb(6))
  a47=-(w-ddb*s(19))*exp(s(19)*zb(6))*(1._wp-fi(5))
  a48=-(w-ddb*s(20))*exp(s(20)*zb(6))*(1._wp-fi(5))
  a49=(w-ddb*s(21))*exp(s(21)*zb(6))*(1._wp-fi(6))
  a50=(w-ddb*s(20))*exp(s(20)*zb(6))*(1._wp-fi(6))
  a51=exp(s(21)*zb(7))
  a52=exp(s(22)*zb(7))
  a53=-exp(s(23)*zb(7))
  a54=-exp(s(24)*zb(7))
  a55=-(w-ddb*s(21))*exp(s(21)*zb(7))*(1._wp-fi(6))
  a56=-(w-ddb*s(22))*exp(s(22)*zb(7))*(1._wp-fi(6))
  a57=(w-ddb*s(23))*exp(s(23)*zb(7))*(1._wp-fi(7))
  a58=(w-ddb*s(24))*exp(s(24)*zb(7))*(1._wp-fi(7))
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

  x1 = (/ffcc,zro,zro,zro,zro,zro,zro,zro,zro,zro,zro,zro,zro,zro,zro,zro/)

  call DGBSV(16,2,2,1,bandedaa,7,ipiv,x1,16,info)

                       
rrm1=(1._wp-fi(1))*jo2*(x1(1)*(exp(s(1)*zo)-1._wp)/s(1)+&
          x1(2)*(exp(s(2)*zo)-1._wp)/s(2)+beta*(x1(3)*(exp(s(11)*zb(2))-&
          exp(s(11)*zo))/s(11)+x1(4)*(exp(s(12)*zb(2))-&
          exp(s(12)*zo))/s(12)))
          
rrm2=(1._wp-fi(2))*beta*jo2*(x1(5)*(exp(s(13)*zb(3))-exp(s(13)*zb(2)))/s(13)+&
                        x1(6)*(exp(s(14)*zb(3))-exp(s(14)*zb(2)))/s(14))    
rrm3=(1._wp-fi(3))*beta*jo2*(x1(7)*(exp(s(15)*zb(4))-exp(s(15)*zb(3)))/s(15)+&
                        x1(8)*(exp(s(16)*zb(4))-exp(s(16)*zb(3)))/s(16))
rrm4=(1._wp-fi(4))*beta*jo2*(x1(9)*(exp(s(17)*zb(5))-exp(s(17)*zb(4)))/s(17)+&
                        x1(10)*(exp(s(18)*zb(5))-exp(s(18)*zb(4)))/s(18))
rrm5=(1._wp-fi(5))*beta*jo2*(x1(11)*(exp(s(19)*zb(6))-exp(s(19)*zb(5)))/s(19)+&
                         x1(12)*(exp(s(20)*zb(6))-exp(s(20)*zb(5)))/s(20))
rrm6=(1._wp-fi(6))*beta*jo2*(x1(13)*(exp(s(21)*zb(7))-exp(s(21)*zb(6)))/s(21)+&
                         x1(14)*(exp(s(22)*zb(7))-exp(s(22)*zb(6)))/s(22))
rrm7=(1._wp-fi(7))*beta*jo2*(x1(15)*(exp(s(23)*zb(8))-exp(s(23)*zb(7)))/s(23)+&
                         x1(16)*(exp(s(24)*zb(8))-exp(s(24)*zb(7)))/s(24))
rrm8=(1._wp-fi(1))*jo2*(beta*(x1(3)*(exp(s(11)*zb(2))-exp(s(11)*zo))/s(11)+&
                         x1(4)*(exp(s(12)*zb(2))-exp(s(12)*zo))/s(12)))+&
                         rrm2+rrm3+rrm4+rrm5+rrm6+rrm7
                         
rrm = (/rrm1,rrm2,rrm3,rrm4,rrm5,rrm6,rrm7,rrm8/)

 cc1 = iinioo2 - del(1)*(x1(1)/(s(1)**2)+x1(2)/(s(2)**2))
 
 dd1 = -del(1)*((x1(1)/s(1)*exp(s(1)*zo)+x1(2)/s(2)*exp(s(2)*zo)+&
           beta*(x1(3)*(exp(s(11)*zb(2))-exp(s(11)*zo))/s(11)+&
           x1(4)*(exp(s(12)*zb(2))-exp(s(12)*zo))/s(12)+&
           del(2)/del(1)*(x1(5)*(exp(s(13)*zb(3))-exp(s(13)*zb(2)))/s(13)+&
           x1(6)*(exp(s(14)*zb(3))-exp(s(14)*zb(2)))/s(14))+&
           del(3)/del(1)*(x1(7)*(exp(s(15)*zb(4))-exp(s(15)*zb(3)))/s(15)+&
           x1(8)*(exp(s(16)*zb(4))-exp(s(16)*zb(3)))/s(16))+&
           del(4)/del(1)*(x1(9)*(exp(s(17)*zb(5))-exp(s(17)*zb(4)))/s(17)+&
           x1(10)*(exp(s(18)*zb(5))-exp(s(18)*zb(4)))/s(18))+&
           del(5)/del(1)*(x1(11)*(exp(s(19)*zb(6))-exp(s(19)*zb(5)))/s(19)+&
           x1(12)*(exp(s(20)*zb(6))-exp(s(20)*zb(5)))/s(20))+&
           del(6)/del(1)*(x1(13)*(exp(s(21)*zb(7))-exp(s(21)*zb(6)))/s(21)+&
           x1(14)*(exp(s(22)*zb(7))-exp(s(22)*zb(6)))/s(22))+&
           del(7)/del(1)*(x1(15)*(exp(s(23)*zb(8))-exp(s(23)*zb(7)))/s(23)+&
           x1(16)*(exp(s(24)*zb(8))-exp(s(24)*zb(7)))/s(24)))))
           
           
! rres1 = -(cc1 + dd1*zo + del(1)*(x1(1)*exp(s(1)*zo)/(s(1)**2)+x1(2)*exp(s(2)*zo)/(s(2)**2)))
 res = -(cc1 + dd1*zo + del(1)*(x1(1)*exp(s(1)*zo)/(s(1)**2)+x1(2)*exp(s(2)*zo)/(s(2)**2)))
END FUNCTION corg1a7w_eq                                       
