! Danish Center for Earth System Science (DCESS) model (DCESS II v1.0)
! Authors: Esteban Fernández and Gary Shaffer
! Contact: est.paleo@gmail.com


program rk4

use parameters, only: wp,n,nto,nta,nab,nob,nlb,nz_all,sy,b_north,fbr,bsh,fdh,fdl,fde
use dimensions
use netcdf

implicit none
CHARACTER(LEN = 100), parameter :: name_out = 'Out.nc', name_res = 'Res.nc'
CHARACTER(LEN = 100) :: name_var
real(wp) :: ini,fin,vms,vmn,vmg
real(wp) :: q_atl(n,6),q_pac(n,5),ekman(n,2)
CALL cpu_time(ini)

!-----------------------------
! Control parameters for the time integration 
tend  = 100                 ! End-time of integration (yr)
dtout = 10                      ! Output interval (yr)
fout  = tend/dtout
h     = 1._wp/24._wp        ! Timestep (times per year)
ko    = 1._wp/h             ! Multiple of timestep at which CO3 and time dependent sediment model are calculated
mon   = 1._wp/h
!-----------------------------

! Load data files
call input_data(inso,oc_pres,shl_dat,olf_bios,gaa,gap,gaar,gaso)

! Load initial conditions
call init_cond(at,atam,lbn,lbs,hlna,mlna,llna,llsa,mlsa,hlnp,mlnp,llnp,llsp,mlsp,soc,arc,slso,&
               pcala,pcalp,pcalar,pcalso,porga,porgp,porgar,porgso,&
               dwpcala,dwpcalssa,dwporga,dwporgssa,fima,wseda,dwpcalp,dwpcalssp,dwporgp,dwporgssp,fimp,wsedp,&
               dwpcalar,dwpcalssar,dwporgar,dwporgssar,fimar,wsedar,&
               dwpcalso,dwpcalssso,dwporgso,dwporgssso,fimso,wsedso)

taam  = atam(1,:)

! We prepare variables for saving data
allocate(st(fout))
allocate(sat(nta,nab,fout))
allocate(slbn(fout,nlb))
allocate(slbs(fout,nlb))

allocate(shlna(nto,n,fout));allocate(smlna(nto,n,fout));allocate(sllna(nto,n,fout));
allocate(sllsa(nto,n,fout));allocate(smlsa(nto,n,fout));allocate(shlnp(nto,n,fout));
allocate(smlnp(nto,n,fout));allocate(sllnp(nto,n,fout));allocate(sllsp(nto,n,fout));
allocate(smlsp(nto,n,fout));allocate(ssoc(nto,n,fout));allocate(sarc(nto,n,fout));allocate(sslso(fout,nto))

allocate(sdwpcala(nob,n,fout));allocate(sdwpcalp(nob,n,fout));allocate(sdwpcalar(fout,n));allocate(sdwpcalso(fout,n))
allocate(sdwporgp(nob,n,fout));allocate(sdwporga(nob,n,fout));allocate(sdwporgar(fout,n));allocate(sdwporgso(fout,n))

allocate(sfima(nob,n,fout ));allocate(sfimp(nob,n,fout ));allocate(sfimar(fout,n));allocate(sfimso(fout,n))
allocate(swseda(nob,n,fout));allocate(swsedp(nob,n,fout));allocate(swsedar(fout,n));allocate(swsedso(fout,n))

! at(5,:) = 560.e-6_wp


! Integration
t     = 0._wp
dt    = h*sy            ! Timestep (sec)
k = 1
dotend: do c=1,tend

attmp    = 0._wp
lbntmp   = 0._wp
lbstmp   = 0._wp

arctmp   = 0._wp
hlnatmp  = 0._wp
mlnatmp  = 0._wp
llnatmp  = 0._wp
llsatmp  = 0._wp
mlsatmp  = 0._wp

hlnptmp  = 0._wp
mlnptmp  = 0._wp
llnptmp  = 0._wp
llsptmp  = 0._wp
mlsptmp  = 0._wp

slsotmp = 0._wp
soctmp  = 0._wp

npartmp = 0._wp
npsotmp = 0._wp
npatmp  = 0._wp
npptmp  = 0._wp

doyr: do cc=1,int(sy/dt)
    t=t+dt

    call odeexp(t,at,taam,lbn,lbs,hlna,mlna,llna,llsa,mlsa,hlnp,mlnp,llnp,llsp,mlsp,arc,soc,slso,&
                porga,porgp,porgar,porgso,pcala,pcalp,pcalar,pcalso,&
                inso,oc_pres,shl_dat,olf_bios,gaa,gap,gaar,gaso,&
                k1at,k1lbn,k1lbs,&
                k1arc,k1hlna,k1mlna,k1llna,k1llsa,k1mlsa,k1hlnp,k1mlnp,k1llnp,k1llsp,k1mlsp,k1soc,k1slso)

    fstp = dt/2._wp
    call step_rk4(at,lbn,lbs,hlna,mlna,llna,llsa,mlsa,hlnp,mlnp,llnp,llsp,mlsp,arc,soc,slso,&
                  k1at,k1lbn,k1lbs,&
                  k1hlna,k1mlna,k1llna,k1llsa,k1mlsa,k1hlnp,k1mlnp,k1llnp,k1llsp,k1mlsp,k1arc,k1soc,k1slso,fstp,&
                  nat,nlbn,nlbs,nhlna,nmlna,nllna,nllsa,nmlsa,nhlnp,nmlnp,nllnp,nllsp,nmlsp,narc,nsoc,nslso)

    call odeexp(t+fstp,nat,taam,nlbn,nlbs,nhlna,nmlna,nllna,nllsa,nmlsa,nhlnp,nmlnp,nllnp,nllsp,nmlsp,narc,nsoc,nslso,&
                porga,porgp,porgar,porgso,pcala,pcalp,pcalar,pcalso,&
                inso,oc_pres,shl_dat,olf_bios,gaa,gap,gaar,gaso,&
                k2at,k2lbn,k2lbs,&
                k2arc,k2hlna,k2mlna,k2llna,k2llsa,k2mlsa,k2hlnp,k2mlnp,k2llnp,k2llsp,k2mlsp,k2soc,k2slso)

    fstp = dt/2._wp
    call step_rk4(at,lbn,lbs,hlna,mlna,llna,llsa,mlsa,hlnp,mlnp,llnp,llsp,mlsp,arc,soc,slso,&
                  k2at,k2lbn,k2lbs,&
                  k2hlna,k2mlna,k2llna,k2llsa,k2mlsa,k2hlnp,k2mlnp,k2llnp,k2llsp,k2mlsp,k2arc,k2soc,k2slso,fstp,&
                  nat,nlbn,nlbs,nhlna,nmlna,nllna,nllsa,nmlsa,nhlnp,nmlnp,nllnp,nllsp,nmlsp,narc,nsoc,nslso)

    call odeexp(t+fstp,nat,taam,nlbn,nlbs,nhlna,nmlna,nllna,nllsa,nmlsa,nhlnp,nmlnp,nllnp,nllsp,nmlsp,narc,nsoc,nslso,&
                porga,porgp,porgar,porgso,pcala,pcalp,pcalar,pcalso,&
                inso,oc_pres,shl_dat,olf_bios,gaa,gap,gaar,gaso,&
                k3at,k3lbn,k3lbs,&
                k3arc,k3hlna,k3mlna,k3llna,k3llsa,k3mlsa,k3hlnp,k3mlnp,k3llnp,k3llsp,k3mlsp,k3soc,k3slso)

    fstp = dt
    call step_rk4(at,lbn,lbs,hlna,mlna,llna,llsa,mlsa,hlnp,mlnp,llnp,llsp,mlsp,arc,soc,slso,&
                  k3at,k3lbn,k3lbs,&
                  k3hlna,k3mlna,k3llna,k3llsa,k3mlsa,k3hlnp,k3mlnp,k3llnp,k3llsp,k3mlsp,k3arc,k3soc,k3slso,fstp,&
                  nat,nlbn,nlbs,nhlna,nmlna,nllna,nllsa,nmlsa,nhlnp,nmlnp,nllnp,nllsp,nmlsp,narc,nsoc,nslso)

    call odeexp(t+fstp,nat,taam,nlbn,nlbs,nhlna,nmlna,nllna,nllsa,nmlsa,nhlnp,nmlnp,nllnp,nllsp,nmlsp,narc,nsoc,nslso,&
                porga,porgp,porgar,porgso,pcala,pcalp,pcalar,pcalso,&
                inso,oc_pres,shl_dat,olf_bios,gaa,gap,gaar,gaso,&
                k4at,k4lbn,k4lbs,&
                k4arc,k4hlna,k4mlna,k4llna,k4llsa,k4mlsa,k4hlnp,k4mlnp,k4llnp,k4llsp,k4mlsp,k4soc,k4slso)

    AT    = AT    + 1._wp/6._wp*dt*( k1AT    + 2._wp*k2AT    + 2._wp*k3AT    + k4AT    )
    LBN   = LBN   + 1._wp/6._wp*dt*( k1LBN   + 2._wp*k2LBN   + 2._wp*k3LBN   + k4LBN   )
    LBS   = LBS   + 1._wp/6._wp*dt*( k1LBS   + 2._wp*k2LBS   + 2._wp*k3LBS   + k4LBS   )

    Arc   = Arc   + 1._wp/6._wp*dt*( k1Arc   + 2._wp*k2Arc   + 2._wp*k3Arc   + k4Arc   )

    HLNA  = HLNA  + 1._wp/6._wp*dt*( k1HLNA  + 2._wp*k2HLNA  + 2._wp*k3HLNA  + k4HLNA  )
    MLNA  = MLNA  + 1._wp/6._wp*dt*( k1MLNA  + 2._wp*k2MLNA  + 2._wp*k3MLNA  + k4MLNA  )
    LLNA  = LLNA  + 1._wp/6._wp*dt*( k1LLNA  + 2._wp*k2LLNA  + 2._wp*k3LLNA  + k4LLNA  )
    LLSA  = LLSA  + 1._wp/6._wp*dt*( k1LLSA  + 2._wp*k2LLSA  + 2._wp*k3LLSA  + k4LLSA  )
    MLSA  = MLSA  + 1._wp/6._wp*dt*( k1MLSA  + 2._wp*k2MLSA  + 2._wp*k3MLSA  + k4MLSA  )

    HLNP  = HLNP  + 1._wp/6._wp*dt*( k1HLNP  + 2._wp*k2HLNP  + 2._wp*k3HLNP  + k4HLNP  )
    MLNP  = MLNP  + 1._wp/6._wp*dt*( k1MLNP  + 2._wp*k2MLNP  + 2._wp*k3MLNP  + k4MLNP  )
    LLNP  = LLNP  + 1._wp/6._wp*dt*( k1LLNP  + 2._wp*k2LLNP  + 2._wp*k3LLNP  + k4LLNP  )
    LLSP  = LLSP  + 1._wp/6._wp*dt*( k1LLSP  + 2._wp*k2LLSP  + 2._wp*k3LLSP  + k4LLSP  )
    MLSP  = MLSP  + 1._wp/6._wp*dt*( k1MLSP  + 2._wp*k2MLSP  + 2._wp*k3MLSP  + k4MLSP  )

    SOc   = SOc   + 1._wp/6._wp*dt*( k1SOc   + 2._wp*k2SOc   + 2._wp*k3SOc   + k4SOc   )
    SlSO  = SlSO  + 1._wp/6._wp*dt*( k1SlSO  + 2._wp*k2SlSO  + 2._wp*k3SlSO  + k4SlSO  )

    ! Ice-line constraints
    at(10,1) = min(at(10,1),b_north )   ! North Atlantic/Arctic
    at(10,2) = min(at(10,2),fbr)        ! North Pacific
    at(10,6) = min(at(10,6),bsh)        ! Southern Ocean

    ! at(5,:) = 560.e-6_wp


    ATtmp    = ATtmp    + AT
    LBNtmp   = LBNtmp   + LBN
    LBStmp   = LBStmp   + LBS

    Arctmp   = Arctmp   + Arc

    HLNAtmp  = HLNAtmp  + HLNA
    MLNAtmp  = MLNAtmp  + MLNA
    LLNAtmp  = LLNAtmp  + LLNA
    LLSAtmp  = LLSAtmp  + LLSA
    MLSAtmp  = MLSAtmp  + MLSA

    HLNPtmp  = HLNPtmp  + HLNP
    MLNPtmp  = MLNPtmp  + MLNP
    LLNPtmp  = LLNPtmp  + LLNP
    LLSPtmp  = LLSPtmp  + LLSP
    MLSPtmp  = MLSPtmp  + MLSP

    SOctmp   = SOctmp   + SOc
    SlSOtmp  = SlSOtmp  + SlSO

    ATam   = ATtmp/mon
    LBNam  = LBNtmp/mon
    LBSam  = LBStmp/mon
    
    Arcam  = Arctmp/mon
    
    HLNAam = HLNAtmp/mon
    MLNAam = MLNAtmp/mon
    LLNAam = LLNAtmp/mon
    LLSAam = LLSAtmp/mon
    MLSAam = MLSAtmp/mon
    
    HLNPam = HLNPtmp/mon
    MLNPam = MLNPtmp/mon
    LLNPam = LLNPtmp/mon
    LLSPam = LLSPtmp/mon
    MLSPam = MLSPtmp/mon
    
    SOcam  = SOctmp/mon
    SlSOam = SlSOtmp/mon

    call pre_sed(at,taam,hlna,mlna,llna,llsa,mlsa,hlnp,mlnp,llnp,llsp,mlsp,arc,soc,slso,&
                 aoif_at,aoif_pa,aoif_ar,aoif_so,rorgat,rorgpa,rorgar,rorgso,&
                 co2a,co3a,omea,co2p,co3p,omep,co2ar,co3ar,omear,co2so,co3so,omeso)

    call ocean_rad(t,inso,qsom,qarm,qm)

    call oce_srcsnk(hlna,mlna,llna,llsa,mlsa,hlnp,mlnp,llnp,llsp,mlsp,arc,soc,&
                aoif_at,aoif_pa,aoif_ar,aoif_so,gaa,gap,gaar,gaso,&
                omea,omep,omear,omeso,co2a,co2p,co2ar,co2so,&
                co3a,co3p,co3ar,co3so,porga,porgp,porgar,porgso,&
                pcala,pcalp,pcalar,pcalso,qsom,qarm,qm,&
                srar,srso,srhlna,srmlna,srllna,srllsa,srmlsa,srhlnp,srmlnp,srllnp,srllsp,srmlsp)

    call river_input(at,srhlna,srmlna,srllna,srllsa,srmlsa,srhlnp,srmlnp,srllnp,srllsp,srmlsp,srar,srso,&
                     rorgat,rorgpa,rorgar,rorgso,rcarat,rcarpa,rcarar,rcarso)

    call new_prod_sed(srhlna,srmlna,srllna,srllsa,srmlsa,srhlnp,srmlnp,srllnp,srllsp,srmlsp,srar,srso,&
                      rorgat,rorgpa,rorgar,rorgso,npat,nppa,npar,npso)

    npartmp = npartmp + npar
    npsotmp = npsotmp + npso
    npatmp  = npatmp  + npat
    npptmp  = npptmp  + nppa
    
    nparam = npartmp/mon
    npsoam = npsotmp/mon
    npaam  = npatmp/mon
    nppam  = npptmp/mon
end do doyr
taam=atam(1,:)

! ============================================================
! Sediment model using annual mean values
! ============================================================

! Arctic
psi = arcam
call smtdcorg(psi,nparam,dwpcalssar,dwporgssar,&
              dwpcalar,dwporgar,fimar,wsedar,pcalar,porgar,ko,h)

! Southern Ocean
psi = socam
call smtdcorg(psi,npsoam,dwpcalssso,dwporgssso,&
              dwpcalso,dwporgso,fimso,wsedso,pcalso,porgso,ko,h)

! hlna
iz  = 1
psi = hlnaam
call smtdcorg(psi,npaam(iz),dwpcalssa(iz,:),dwporgssa(iz,:),&
              dwpcala(iz,:),dwporga(iz,:),fima(iz,:),wseda(iz,:),pcala(iz,:),porga(iz,:),ko,h)

! mlna
iz  = 2
psi = mlnaam
call smtdcorg(psi,npaam(iz),dwpcalssa(iz,:),dwporgssa(iz,:),&
              dwpcala(iz,:),dwporga(iz,:),fima(iz,:),wseda(iz,:),pcala(iz,:),porga(iz,:),ko,h)

! llna
iz  = 3
psi = llnaam
call smtdcorg(psi,npaam(iz),dwpcalssa(iz,:),dwporgssa(iz,:),&
              dwpcala(iz,:),dwporga(iz,:),fima(iz,:),wseda(iz,:),pcala(iz,:),porga(iz,:),ko,h)
              
! llsa
iz  = 4
psi = llsaam
call smtdcorg(psi,npaam(iz),dwpcalssa(iz,:),dwporgssa(iz,:),&
              dwpcala(iz,:),dwporga(iz,:),fima(iz,:),wseda(iz,:),pcala(iz,:),porga(iz,:),ko,h)

! mlsa
iz  = 5
psi = mlsaam
call smtdcorg(psi,npaam(iz),dwpcalssa(iz,:),dwporgssa(iz,:),&
              dwpcala(iz,:),dwporga(iz,:),fima(iz,:),wseda(iz,:),pcala(iz,:),porga(iz,:),ko,h)

! hlnp
iz  = 1
psi = hlnpam
call smtdcorg(psi,nppam(iz),dwpcalssp(iz,:),dwporgssp(iz,:),&
              dwpcalp(iz,:),dwporgp(iz,:),fimp(iz,:),wsedp(iz,:),pcalp(iz,:),porgp(iz,:),ko,h)

! mlnp
iz  = 2
psi = mlnpam
call smtdcorg(psi,nppam(iz),dwpcalssp(iz,:),dwporgssp(iz,:),&
              dwpcalp(iz,:),dwporgp(iz,:),fimp(iz,:),wsedp(iz,:),pcalp(iz,:),porgp(iz,:),ko,h)

! llnp
iz  = 3
psi = llnpam
call smtdcorg(psi,nppam(iz),dwpcalssp(iz,:),dwporgssp(iz,:),&
              dwpcalp(iz,:),dwporgp(iz,:),fimp(iz,:),wsedp(iz,:),pcalp(iz,:),porgp(iz,:),ko,h)
              
! llsp
iz  = 4
psi = llspam
call smtdcorg(psi,nppam(iz),dwpcalssp(iz,:),dwporgssp(iz,:),&
              dwpcalp(iz,:),dwporgp(iz,:),fimp(iz,:),wsedp(iz,:),pcalp(iz,:),porgp(iz,:),ko,h)

! mlsp
iz  = 5
psi = mlspam
call smtdcorg(psi,nppam(iz),dwpcalssp(iz,:),dwporgssp(iz,:),&
              dwpcalp(iz,:),dwporgp(iz,:),fimp(iz,:),wsedp(iz,:),pcalp(iz,:),porgp(iz,:),ko,h)

! End sediment

! Save data
if (mod(c, dtout) == 0) then

write (*,'(A)') 'Saving'

st(k)      = t/sy
sat(:,:,k) = atam
slbn(k,:)  = lbnam
slbs(k,:)  = lbsam

sarc(:,:,k) = arcam
ssoc(:,:,k) = socam
sslso(k,:)  = slsoam

shlna(:,:,k) = hlnaam
smlna(:,:,k) = mlnaam
sllna(:,:,k) = llnaam
sllsa(:,:,k) = llsaam
smlsa(:,:,k) = mlsaam

shlnp(:,:,k) = hlnpam
smlnp(:,:,k) = mlnpam
sllnp(:,:,k) = llnpam
sllsp(:,:,k) = llspam
smlsp(:,:,k) = mlspam

sdwpcalar(k,:)  = dwpcalar
sdwpcalso(k,:)  = dwpcalso
sdwpcala(:,:,k) = dwpcala
sdwpcalp(:,:,k) = dwpcalp

sdwporgar(k,:)  = dwporgar
sdwporgso(k,:)  = dwporgso
sdwporga(:,:,k) = dwporga
sdwporgp(:,:,k) = dwporgp

sfimar(k,:)  = fimar
sfimso(k,:)  = fimso
sfima(:,:,k) = fima
sfimp(:,:,k) = fimp

swsedar(k,:)  = wsedar
swsedso(k,:)  = wsedso
swseda(:,:,k) = wseda
swsedp(:,:,k) = wsedp

k=k+1

end if

write (*,'(A)') '------------------------------------------------------------'
write (*,'(A,I10,A,I10,A)') 'Integration at time: ', int(t/sy), '   (tend = ', int(tend), ')'

vmn = atam(1,1)*(1._wp-sin(fdh)) + atam(1,2)*(sin(fdh)-sin(fdl)) + atam(1,3)*(sin(fdl)-sin(fde))
vms = atam(1,6)*(1._wp-sin(fdh)) + atam(1,5)*(sin(fdh)-sin(fdl)) + atam(1,4)*(sin(fdl)-sin(fde))
vmg  = .5_wp*(vmN+vmS)
print*, vms,vmn,vmg


end do dotend ! end main loop

write (*,'(A)') '------------------------------------------------------------'
write (*,'(A)') 'End of simulation'
write (*,'(A)') '------------------------------------------------------------'

CALL cpu_time(fin)

PRINT '("Time execution =",f20.5," hours.")',(fin-ini)/3600._wp

! =========================================================================
! Save netcdf files
! =========================================================================
! Save restart file
CALL CHECK( NF90_CREATE(name_res, NF90_NETCDF4, ncid) )
CALL CHECK( NF90_DEF_DIM(ncid, 'oc_layers', n    ,dimoz) )
CALL CHECK( NF90_DEF_DIM(ncid, 'oc_trc'   , nto  ,dimot) )
CALL CHECK( NF90_DEF_DIM(ncid, 'oc_boxn'  , nob  ,dimob) )
CALL CHECK( NF90_DEF_DIM(ncid, 'at_trc'   , nta  ,dimat) )
CALL CHECK( NF90_DEF_DIM(ncid, 'at_boxn'  , nab  ,dimab) )
CALL CHECK( NF90_DEF_DIM(ncid, 'lb_trc'   , nlb  ,dimlt) )
CALL CHECK(NF90_CLOSE(ncid))

! save 1D data
name_var = 'lbn'       ; call save_1d(name_res,name_var,lbn       ,nlb,dimlt)
name_var = 'lbs'       ; call save_1d(name_res,name_var,lbs       ,nlb,dimlt)
name_var = 'pcalar'    ; call save_1d(name_res,name_var,pcalar    ,n  ,dimoz)
name_var = 'pcalso'    ; call save_1d(name_res,name_var,pcalso    ,n  ,dimoz)
name_var = 'porgar'    ; call save_1d(name_res,name_var,porgar    ,n  ,dimoz)
name_var = 'porgso'    ; call save_1d(name_res,name_var,porgso    ,n  ,dimoz)
name_var = 'fimar'     ; call save_1d(name_res,name_var,fimar     ,n  ,dimoz)
name_var = 'fimso'     ; call save_1d(name_res,name_var,fimso     ,n  ,dimoz)
name_var = 'wsedar'    ; call save_1d(name_res,name_var,wsedar    ,n  ,dimoz)
name_var = 'wsedso'    ; call save_1d(name_res,name_var,wsedso    ,n  ,dimoz)
name_var = 'slso'      ; call save_1d(name_res,name_var,slso      ,nto,dimot)
name_var = 'dwpcalar'  ; call save_1d(name_res,name_var,dwpcalar  ,n  ,dimoz)
name_var = 'dwporgar'  ; call save_1d(name_res,name_var,dwporgar  ,n  ,dimoz)
name_var = 'dwpcalso'  ; call save_1d(name_res,name_var,dwpcalso  ,n  ,dimoz)
name_var = 'dwporgso'  ; call save_1d(name_res,name_var,dwporgso  ,n  ,dimoz)
name_var = 'dwpcalssar'; call save_1d(name_res,name_var,dwpcalssar,n  ,dimoz)
name_var = 'dwporgssar'; call save_1d(name_res,name_var,dwporgssar,n  ,dimoz)
name_var = 'dwpcalssso'; call save_1d(name_res,name_var,dwpcalssso,n  ,dimoz)
name_var = 'dwporgssso'; call save_1d(name_res,name_var,dwporgssso,n  ,dimoz)

! save 2D data
name_var = 'at'       ; call save_2d(name_res,name_var,at       ,nta,nab,dimat,dimab)
name_var = 'atam'     ; call save_2d(name_res,name_var,atam     ,nta,nab,dimat,dimab)
name_var = 'arc'      ; call save_2d(name_res,name_var,arc      ,nto,n  ,dimot,dimoz)
name_var = 'soc'      ; call save_2d(name_res,name_var,soc      ,nto,n  ,dimot,dimoz)
name_var = 'hlna'     ; call save_2d(name_res,name_var,hlna     ,nto,n  ,dimot,dimoz)
name_var = 'mlna'     ; call save_2d(name_res,name_var,mlna     ,nto,n  ,dimot,dimoz)
name_var = 'llna'     ; call save_2d(name_res,name_var,llna     ,nto,n  ,dimot,dimoz)
name_var = 'llsa'     ; call save_2d(name_res,name_var,llsa     ,nto,n  ,dimot,dimoz)
name_var = 'mlsa'     ; call save_2d(name_res,name_var,mlsa     ,nto,n  ,dimot,dimoz)
name_var = 'hlnp'     ; call save_2d(name_res,name_var,hlnp     ,nto,n  ,dimot,dimoz)
name_var = 'mlnp'     ; call save_2d(name_res,name_var,mlnp     ,nto,n  ,dimot,dimoz)
name_var = 'llnp'     ; call save_2d(name_res,name_var,llnp     ,nto,n  ,dimot,dimoz)
name_var = 'llsp'     ; call save_2d(name_res,name_var,llsp     ,nto,n  ,dimot,dimoz)
name_var = 'mlsp'     ; call save_2d(name_res,name_var,mlsp     ,nto,n  ,dimot,dimoz)
name_var = 'pcala'    ; call save_2d(name_res,name_var,pcala    ,nob,n  ,dimob,dimoz)
name_var = 'pcalp'    ; call save_2d(name_res,name_var,pcalp    ,nob,n  ,dimob,dimoz)
name_var = 'porga'    ; call save_2d(name_res,name_var,porga    ,nob,n  ,dimob,dimoz)
name_var = 'porgp'    ; call save_2d(name_res,name_var,porgp    ,nob,n  ,dimob,dimoz)
name_var = 'fima'     ; call save_2d(name_res,name_var,fima     ,nob,n  ,dimob,dimoz)
name_var = 'fimp'     ; call save_2d(name_res,name_var,fimp     ,nob,n  ,dimob,dimoz)
name_var = 'wseda'    ; call save_2d(name_res,name_var,wseda    ,nob,n  ,dimob,dimoz)
name_var = 'wsedp'    ; call save_2d(name_res,name_var,wsedp    ,nob,n  ,dimob,dimoz)
name_var = 'dwpcala'  ; call save_2d(name_res,name_var,dwpcala  ,nob,n  ,dimob,dimoz)
name_var = 'dwpcalp'  ; call save_2d(name_res,name_var,dwpcalp  ,nob,n  ,dimob,dimoz)
name_var = 'dwporga'  ; call save_2d(name_res,name_var,dwporga  ,nob,n  ,dimob,dimoz)
name_var = 'dwporgp'  ; call save_2d(name_res,name_var,dwporgp  ,nob,n  ,dimob,dimoz)
name_var = 'dwporgssa'; call save_2d(name_res,name_var,dwporgssa,nob,n  ,dimob,dimoz)
name_var = 'dwporgssp'; call save_2d(name_res,name_var,dwporgssp,nob,n  ,dimob,dimoz)
name_var = 'dwpcalssa'; call save_2d(name_res,name_var,dwpcalssa,nob,n  ,dimob,dimoz)
name_var = 'dwpcalssp'; call save_2d(name_res,name_var,dwpcalssp,nob,n  ,dimob,dimoz)

PRINT *, "*** SUCCESS WRITING NETCDF FILE: ", name_res

! Save output model
! Call ocean circulation fluxes
call ocean_circulation(atam,hlnaam,mlnaam,llnaam,llsaam,mlsaam,hlnpam,mlnpam,llnpam,llspam,mlspam,socam,arcam,slsoam,&
                       inso,oc_pres,shl_dat,&
                       q_atl,q_pac,ekman)

CALL CHECK( NF90_CREATE(name_out, NF90_NETCDF4, ncid) )
CALL CHECK( NF90_DEF_DIM(ncid, 'time'      , fout ,dimt ) )
CALL CHECK( NF90_DEF_DIM(ncid, 'oc_layers' , n    ,dimoz) )
CALL CHECK( NF90_DEF_DIM(ncid, 'oc_trc'    , nto  ,dimot) )
CALL CHECK( NF90_DEF_DIM(ncid, 'oc_boxn'   , nob  ,dimob) )
CALL CHECK( NF90_DEF_DIM(ncid, 'at_trc'    , nta  ,dimat) )
CALL CHECK( NF90_DEF_DIM(ncid, 'at_boxn'   , nab  ,dimab) )
CALL CHECK( NF90_DEF_DIM(ncid, 'lb_trc'    , nlb  ,dimlt) )
CALL CHECK( NF90_DEF_DIM(ncid, 'oc_div_atl', 6    ,dimoda) )
CALL CHECK( NF90_DEF_DIM(ncid, 'oc_div_pac', 5    ,dimodp) )
CALL CHECK( NF90_DEF_DIM(ncid, 'ekman_flxs', 2    ,dimek) )
CALL CHECK(NF90_CLOSE(ncid))

! save 1D data
name_var = 'st'; call save_1d(name_out,name_var,st,fout,dimt)

! save 2D data
name_var = 'slbn'   ; call save_2d(name_out,name_var,slbn   ,fout,nlb,dimt,dimlt)
name_var = 'slbs'   ; call save_2d(name_out,name_var,slbs   ,fout,nlb,dimt,dimlt)
name_var = 'sdwpcalar'; call save_2d(name_out,name_var,sdwpcalar,fout,n  ,dimt,dimoz)
name_var = 'sdwpcalso'; call save_2d(name_out,name_var,sdwpcalso,fout,n  ,dimt,dimoz)
name_var = 'sdwporgar'; call save_2d(name_out,name_var,sdwporgar,fout,n  ,dimt,dimoz)
name_var = 'sdwporgso'; call save_2d(name_out,name_var,sdwporgso,fout,n  ,dimt,dimoz)
name_var = 'sfimar' ; call save_2d(name_out,name_var,sfimar ,fout,n  ,dimt,dimoz)
name_var = 'sfimso' ; call save_2d(name_out,name_var,sfimso ,fout,n  ,dimt,dimoz)
name_var = 'swsedar'; call save_2d(name_out,name_var,swsedar,fout,n  ,dimt,dimoz)
name_var = 'swsedso'; call save_2d(name_out,name_var,swsedso,fout,n  ,dimt,dimoz)
name_var = 'sslso'  ; call save_2d(name_out,name_var,sslso  ,fout,nto,dimt,dimot)
name_var = 'q_atl'  ; call save_2d(name_out,name_var,q_atl  ,n   ,6  ,dimoz,dimoda)
name_var = 'q_pac'  ; call save_2d(name_out,name_var,q_pac  ,n   ,5  ,dimoz,dimodp)
name_var = 'ekman'  ; call save_2d(name_out,name_var,ekman  ,n   ,2  ,dimoz,dimek )

! save 3D data
name_var = 'sat'; call save_3d(name_out,name_var,sat,nta,nab,fout,dimat,dimab,dimt)

name_var = 'sarc'  ; call save_3d(name_out,name_var,sarc  ,nto,n,fout,dimot,dimoz,dimt)
name_var = 'ssoc'  ; call save_3d(name_out,name_var,ssoc  ,nto,n,fout,dimot,dimoz,dimt)
name_var = 'shlna' ; call save_3d(name_out,name_var,shlna ,nto,n,fout,dimot,dimoz,dimt)
name_var = 'smlna' ; call save_3d(name_out,name_var,smlna ,nto,n,fout,dimot,dimoz,dimt)
name_var = 'sllna' ; call save_3d(name_out,name_var,sllna ,nto,n,fout,dimot,dimoz,dimt)
name_var = 'sllsa' ; call save_3d(name_out,name_var,sllsa ,nto,n,fout,dimot,dimoz,dimt)
name_var = 'smlsa' ; call save_3d(name_out,name_var,smlsa ,nto,n,fout,dimot,dimoz,dimt)
name_var = 'shlnp' ; call save_3d(name_out,name_var,shlnp ,nto,n,fout,dimot,dimoz,dimt)
name_var = 'smlnp' ; call save_3d(name_out,name_var,smlnp ,nto,n,fout,dimot,dimoz,dimt)
name_var = 'sllnp' ; call save_3d(name_out,name_var,sllnp ,nto,n,fout,dimot,dimoz,dimt)
name_var = 'sllsp' ; call save_3d(name_out,name_var,sllsp ,nto,n,fout,dimot,dimoz,dimt)
name_var = 'smlsp' ; call save_3d(name_out,name_var,smlsp ,nto,n,fout,dimot,dimoz,dimt)
name_var = 'sdwpcala'; call save_3d(name_out,name_var,sdwpcala,nob,n,fout,dimob,dimoz,dimt)
name_var = 'sdwpcalp'; call save_3d(name_out,name_var,sdwpcalp,nob,n,fout,dimob,dimoz,dimt)
name_var = 'sdwporga'; call save_3d(name_out,name_var,sdwporga,nob,n,fout,dimob,dimoz,dimt)
name_var = 'sdwporgp'; call save_3d(name_out,name_var,sdwporgp,nob,n,fout,dimob,dimoz,dimt)
name_var = 'sfima' ; call save_3d(name_out,name_var,sfima ,nob,n,fout,dimob,dimoz,dimt)
name_var = 'sfimp' ; call save_3d(name_out,name_var,sfimp ,nob,n,fout,dimob,dimoz,dimt)
name_var = 'swseda'; call save_3d(name_out,name_var,swseda,nob,n,fout,dimob,dimoz,dimt)
name_var = 'swsedp'; call save_3d(name_out,name_var,swsedp,nob,n,fout,dimob,dimoz,dimt)



PRINT *, "*** SUCCESS WRITING NETCDF FILE: ", name_out


deallocate(st);deallocate(sat);deallocate(slbn);deallocate(slbs)

deallocate(shlna);deallocate(smlna);deallocate(sllna);deallocate(sllsa);deallocate(smlsa)
deallocate(shlnp);deallocate(smlnp);deallocate(sllnp);deallocate(sllsp);deallocate(smlsp)
deallocate(ssoc);deallocate(sarc);deallocate(sslso)

deallocate(sdwpcala);deallocate(sdwpcalp);deallocate(sdwpcalar);deallocate(sdwpcalso)
deallocate(sdwporgp);deallocate(sdwporga);deallocate(sdwporgar);deallocate(sdwporgso)

deallocate(sfima);deallocate(sfimp);deallocate(sfimar);deallocate(sfimso)
deallocate(swseda);deallocate(swsedp);deallocate(swsedar);deallocate(swsedso)

!=====================================================================
! End model
!=====================================================================
end program rk4

include 'Input_Data.f90'
include 'Init_Cond.f90'
include 'OdeExp.f90'
include 'Pre_Sed.f90'
include 'New_Prod_Sed.f90'
include 'River_input.f90'
include 'SMtdCorg.f90'
include 'Step_RK4.f90'
include 'Ocean_Circulation.f90'
include 'netcdf_check.f90'
include 'Save_1D.f90'
include 'Save_2D.f90'
include 'Save_3D.f90'
!------------------------------
! 
!  Copyright © 2024 Danish Center for Earth System Science
!  Permission is hereby granted, free of charge, to any person obtaining a copy of this software 
!  and associated documentation files (the "Software"), to deal in the Software without restriction, 
!  including without limitation the rights to use, copy, and modify copies of the Software, and 
!  to permit persons to whom the Software is furnished to do so, subject to the following conditions:
 
!  The above copyright notice and this permission notice shall be included in all copies or 
!  substantial portions of the Software.
 
!  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING 
!  BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
!  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
!  DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
!  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

