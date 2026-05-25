SUBROUTINE EVALUATE_S_CORG()
	use parameters, only: wp
	USE CORG_DATA, ONLY: s,w,ddb,beta,jo2
	IMPLICIT NONE
	s(1) =   0.5_wp*(w/ddb + (w**2/(ddb**2)+4._wp*jo2/ddb)**0.5)
	s(2) =   0.5_wp*(w/ddb - (w**2/(ddb**2)+4._wp*jo2/ddb)**0.5)
	s(27)=s(1); s(25)=s(1); s(9) =s(1); s(7)=s(1); s(5)=s(1); s(3)=s(1)
	s(28)=s(2); s(26)=s(2);	s(10)=s(2); s(8)=s(2); s(6)=s(2); s(4)=s(2)

	s(11) = 0.5_wp*(w/ddb + (w**2/(ddb**2)+4._wp*beta*jo2/ddb)**0.5)
	s(12) = 0.5_wp*(w/ddb - (w**2/(ddb**2)+4._wp*beta*jo2/ddb)**0.5)
	s(24)=s(12); s(22)=s(12); s(20)=s(12); s(18)=s(12); s(16)=s(12); s(14)=s(12)
	s(23)=s(11); s(21)=s(11); s(19)=s(11); s(17)=s(11); s(15)=s(11); s(13)=s(11)

END SUBROUTINE EVALUATE_S_CORG
