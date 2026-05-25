module dimensions
use parameters, only: wp,n,nto,nta,nlb,nz_all,nob,nab

implicit none

real(wp) :: inso(19,14),oc_pres(n,7),shl_dat(nz_all,3),olf_bios(1801,2),&
            gaa(nob,n),gap(nob,n),gaar(n),gaso(n)
real(wp), dimension(nto,n):: hlna,mlna,llna,llsa,mlsa,hlnp,mlnp,llnp,llsp,mlsp,soc,arc
real(wp), dimension(nto) :: slso
real(wp), dimension(nob,n) :: pcala,pcalp,porga,porgp
real(wp), dimension(n) :: pcalar,pcalso,porgar,porgso
real(wp) :: at(nta,nab),atam(nta,nab),lbn(nlb),lbs(nlb),taam(nab)
real(wp), dimension(nob,n):: dwpcala,dwpcalssa,dwporga,dwporgssa,fima,wseda,&
                           dwpcalp,dwpcalssp,dwporgp,dwporgssp,fimp,wsedp
real(wp), dimension(n) :: dwpcalar,dwpcalssar,dwporgar,dwporgssar,fimar,wsedar,&
                          dwpcalso,dwpcalssso,dwporgso,dwporgssso,fimso,wsedso

real(wp), dimension(nto,n) :: arctmp,hlnatmp,mlnatmp,llnatmp,llsatmp,mlsatmp,hlnptmp,mlnptmp,llnptmp,llsptmp,mlsptmp,soctmp,&
                              narc,nhlna,nmlna,nllna,nllsa,nmlsa,nhlnp,nmlnp,nllnp,nllsp,nmlsp,nsoc,&
                              arcam,hlnaam,mlnaam,llnaam,llsaam,mlsaam,hlnpam,mlnpam,llnpam,llspam,mlspam,socam
real(wp) :: attmp(nta,nab),lbntmp(nlb),lbstmp(nlb),slsotmp(nto)
integer :: i,cc,c,tend,dtout,iz,fout,k
real(wp) :: t,h,ko,dt,fstp,ntp,mon

real(wp), dimension(nto,n) :: k1arc,k1hlna,k1mlna,k1llna,k1llsa,k1mlsa,k1hlnp,k1mlnp,k1llnp,k1llsp,k1mlsp,k1soc,&
                              k2arc,k2hlna,k2mlna,k2llna,k2llsa,k2mlsa,k2hlnp,k2mlnp,k2llnp,k2llsp,k2mlsp,k2soc,&
                              k3arc,k3hlna,k3mlna,k3llna,k3llsa,k3mlsa,k3hlnp,k3mlnp,k3llnp,k3llsp,k3mlsp,k3soc,&
                              k4arc,k4hlna,k4mlna,k4llna,k4llsa,k4mlsa,k4hlnp,k4mlnp,k4llnp,k4llsp,k4mlsp,k4soc
real(wp) :: k1at(nta,nab),k1lbn(nlb),k1lbs(nlb),k1slso(nto),nat(nta,nab),nlbn(nlb),nlbs(nlb),nslso(nto)
real(wp) :: k2at(nta,nab),k2lbn(nlb),k2lbs(nlb),k2slso(nto),&
            k3at(nta,nab),k3lbn(nlb),k3lbs(nlb),k3slso(nto),&
            k4at(nta,nab),k4lbn(nlb),k4lbs(nlb),k4slso(nto)
real(wp) :: aoif_at(nob),aoif_pa(nob),aoif_ar,aoif_so
real(wp) :: qsom,qarm,qm(nob)
real(wp),dimension(nob) :: co2a,co3a,omea,co2p,co3p,omep
real(wp) :: co2ar,co3ar,omear,co2so,co3so,omeso
real(wp), dimension(nto,n) :: srar,srso,srhlna,srmlna,srllna,srllsa,srmlsa,srhlnp,srmlnp,srllnp,srllsp,srmlsp
real(wp) :: rcarat(nob),rcarpa(nob),rcarar,rcarso,rorgat(nob),rorgpa(nob),rorgar,rorgso
real(wp) :: npar,npso,npat(nob),nppa(nob),npartmp,npsotmp,npatmp(nob),npptmp(nob),nparam,npsoam,npaam(nob),nppam(nob)
real(wp) :: lbnam(nlb),lbsam(nlb),slsoam(nto),psi(nto,n)
! Save data
real(wp), allocatable, dimension(:,:,:) :: shlna,smlna,sllna,sllsa,smlsa,shlnp,smlnp,sllnp,sllsp,smlsp,ssoc,sarc,&
sdwpcala,sdwpcalp,sdwporga,sdwporgp,sfima,sfimp,swseda,swsedp,sat
real(wp), allocatable, dimension(:,:) :: slbn,slbs,sslso,sdwpcalar,sdwpcalso,sdwporgar,sdwporgso,sfimar,sfimso,swsedar,swsedso
real(wp), allocatable, dimension(:) :: st

! netcdf
integer :: ncid,dimt,dimoz,dimot,dimat,dimob,dimab,dimlt,dimoda,dimodp,dimek

end module dimensions