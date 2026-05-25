subroutine step_rk4(at,lbn,lbs,hna,mna,ena,esa,msa,hnp,mnp,enp,esp,msp,arc,soc,slso,&
                kat,klbn,klbs,&
                khna,kmna,kena,kesa,kmsa,khnp,kmnp,kenp,kesp,kmsp,karc,ksoc,kslso,f,&
                nat,nlbn,nlbs,nhna,nmna,nena,nesa,nmsa,nhnp,nmnp,nenp,nesp,nmsp,narc,nsoc,nslso)
use parameters, only: wp,nta,nto,nlb,n,nab
implicit none
real(wp), intent(in) , dimension(nta,nab) :: at,kat
real(wp), intent(in) , dimension(nto,n) :: hna,mna,ena,esa,msa,hnp,mnp,enp,esp,msp,arc,soc,&
                                           khna,kmna,kena,kesa,kmsa,khnp,kmnp,kenp,kesp,kmsp,karc,ksoc
real(wp), intent(in) , dimension(nto) :: slso,kslso
real(wp), intent(in) , dimension(nlb) :: lbn,lbs,klbn,klbs
real(wp), intent(in) :: f

real(wp), intent(out) , dimension(nta,nab) :: nat
real(wp), intent(out) , dimension(nto,n) :: nhna,nmna,nena,nesa,nmsa,nhnp,nmnp,nenp,nesp,nmsp,narc,nsoc
real(wp), intent(out) , dimension(nto) :: nslso
real(wp), intent(out) , dimension(nlb) :: nlbn,nlbs

nat    = at  + kat *f
nlbn   = lbn + klbn*f
nlbs   = lbs + klbs*f

narc   = arc + karc  *f

nhna  = hna  + khna *f
nmna  = mna  + kmna *f
nena  = ena  + kena *f
nesa  = esa  + kesa *f
nmsa  = msa  + kmsa *f

nhnp  = hnp  + khnp *f
nmnp  = mnp  + kmnp *f
nenp  = enp  + kenp *f
nesp  = esp  + kesp *f
nmsp  = msp  + kmsp *f

nsoc   = soc   + ksoc  *f
nslso  = slso  + kslso *f


end subroutine step_rk4