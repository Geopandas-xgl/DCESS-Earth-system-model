module parameters
! Set global parameters used in the model.
implicit none
    
integer, parameter :: wp = selected_real_kind(15, 307)

! Useful constans
real(wp), parameter :: pi = acos(-1.0_wp)
real(wp), parameter :: sy = 31536000._wp                        ! Seconds per year (s/yr)
real(wp), parameter :: rc = 4.e6_wp                             ! Scaled heat capacity (J/m3/K) 
real(wp), parameter :: rva = 1._wp/(1.7676e20_wp/2._wp),&       ! Inverse atmospheric volume (one hemisphere) meassure [atm/mol]
                       mgt = 8.326e13_wp,&                      ! Moles C per GtC
                       avg = 6.024e23_wp                        ! Avogadro number [atoms/mol]
real(wp), parameter :: e_r = 6.371e6_wp                         ! Earth radius (m)
real(wp), parameter :: g   = 9.82_wp                            ! Gravity (m/s/s)
real(wp), parameter :: dr  = pi/180._wp                         ! degree to radians



integer, parameter :: nta = 10,&                                ! number of atmospheric tracers
                      nto = 10,&                                ! number of oceanic tracers
                      nlb = 36,&                                ! number of land biosphere tracers
                      nab = 6 ,&                                ! number of atmospheric boxes
                      nob = 5                                   ! number of ocean boxes (per basin, Atl/Pac)

!========================================================================
! Geometry
!========================================================================
! Ocean and Atmosphere boundaries (degrees)
real(wp), parameter :: b_north = 90._wp*dr          ! North
real(wp), parameter :: fbr     = 65._wp*dr          ! Bering Strait
real(wp), parameter :: fdh     = 55._wp*dr          ! high latitudes
real(wp), parameter :: fdl     = 35._wp*dr          ! low latitudes
real(wp), parameter :: fde     = 0._wp *dr          ! equator
real(wp), parameter :: b_south = 69._wp*dr          ! Southern Ocean limit
real(wp), parameter :: bsh     = 70._wp*dr          ! Southern Ocean shelf limit (Antarctica)


! Ocean vertical geometry parameters
integer, parameter ::  n    = 55,&                  ! number of ocean layers
                       n_ds = 10,&                  ! "Denmark Strait" depth
                       n_dp = 20,&                  ! Drake Passage sill depth
                       shl  = 5                     ! shelf level

real(wp), parameter :: d  = 100._wp,&               ! vertical resolution (m)
                       dm = 100._wp                 ! mixing layer depth (m)
integer :: i
real(wp), parameter, dimension(n)  :: zmid = (/(real(i,wp)*d+d/2._wp, i=0,n-1)/) ! centered layer depths (m)

! Center of ocean boxes in radians
real(wp), parameter :: ly_ar = ( b_north - fbr     )/2._wp + fbr, &
                       ly_n  = ( fbr     - fdh     )/2._wp + fdh, &
                       ly_mn = ( fdh     - fdl     )/2._wp + fdl, &
                       ly_en = ( fdl     - fde     )/2._wp + fde, &
                       ly_es = ( fde     + fdl     )/2._wp - fdl, &
                       ly_ms = (-fdl     + fdh     )/2._wp - fdh, &
                       ly_s  = (-fdh     + b_south )/2._wp - b_south
! Ocean widths
real(wp), parameter :: w_55np = 80._wp /360._wp, &
                       w_35np = 100._wp/360._wp, &
                       w_eqp  = 200._wp/360._wp, &
                       w_35sp = 220._wp/360._wp, &
                       w_55sp = 260._wp/360._wp, &
                       w_65na = 60._wp /360._wp, &
                       w_55na = 70._wp /360._wp, &
                       w_35na = 60._wp /360._wp, &
                       w_eqa  = 60._wp /360._wp, &
                       w_35sa = 70._wp /360._wp, &
                       w_55sa = 90._wp /360._wp

! Catchment area factors for each ocean box (drainage areas)
real(wp), parameter :: c_ar  = 1._wp      , &
                       c_na  = 0.375_wp   , &
                       c_mna = 0.38_wp    , &
                       c_msa = 0.25_wp    , &
                       c_sa  = 0.30_wp    , &
                       c_np  = 0.275_wp   , &
                       c_mnp = 1._wp-c_mna, &
                       c_msp = 1._wp-c_msa, &
                       c_sp  = 1._wp-c_sa , &
                       c_arx = 1._wp-(c_na + c_np)

! Areas (m2)
real(wp), parameter :: h_area = 2._wp*pi*e_r**2                 ! Hemispheric Area
! Atmospheric areas
real(wp), parameter :: aa_h = h_area*(1._wp    - sin(fdh))      ! High latitude (55 N/S - 90 N/S)
real(wp), parameter :: aa_m = h_area*(sin(fdh) - sin(fdl))      ! Mid latitude  (35 N/S - 55 N/S)
real(wp), parameter :: aa_e = h_area*sin(fdl)                   ! Low latitude  (Eq     - 35 N/S)


! Ocean to land area fractions
real(wp), parameter :: olf_ar  = 0.5903_wp,&
                       olf_na  = 0.2414_wp,&
                       olf_mna = 0.2223_wp,&
                       olf_ena = 0.1936_wp,&
                       olf_esa = 0.1559_wp,&
                       olf_msa = 0.2336_wp,&
                       olf_np  = 0.1395_wp,&
                       olf_mnp = 0.2771_wp,&
                       olf_enp = 0.5022_wp,&
                       olf_esp = 0.6231_wp,&
                       olf_msp = 0.7283_wp
! Ocean areas (m2)
real(wp), parameter :: ao_ar  =  olf_ar*h_area*(sin(b_north) - sin(fbr)    ),& ! Arctic
                       ao_so  =         h_area*(sin(b_south) - sin(fdh)    ),& ! Southern Ocean
                       ao_sh  =         h_area*(sin(bsh)     - sin(b_south)),& ! SO shelf
                       ao_na  = olf_na *h_area*(sin(fbr)     - sin(fdh)    ),& ! hna
                       ao_mna = olf_mna*h_area*(sin(fdh)     - sin(fdl)    ),& ! mna
                       ao_ena = olf_ena*h_area*(sin(fdl)     - sin(fde)    ),& ! ena
                       ao_esa = olf_esa*h_area*(sin(fdl)     - sin(fde)    ),& ! esa
                       ao_msa = olf_msa*h_area*(sin(fdh)     - sin(fdl)    ),& ! msa
                       ao_np  = olf_np *h_area*(sin(fbr)     - sin(fdh)    ),& ! hnp
                       ao_mnp = olf_mnp*h_area*(sin(fdh)     - sin(fdl)    ),& ! mnp
                       ao_enp = olf_enp*h_area*(sin(fdl)     - sin(fde)    ),& ! enp
                       ao_esp = olf_esp*h_area*(sin(fdl)     - sin(fde)    ),& ! esp
                       ao_msp = olf_msp*h_area*(sin(fdh)     - sin(fdl)    )   ! msp
!========================================================================
! End geometry
!========================================================================

! ==========================================================================
! Atmosphere
! ==========================================================================
real(wp), parameter :: ko = 40._wp,&                ! Ocean Atmosphere energy exchange
                       tsnow = 0._wp,&              ! Temperature for snowline 
                       tn2o = 150._wp               ! Atmospheric lifetime N2O
! Pre-industrial GHG
real(wp), parameter :: pco2int = 280.e-6_wp         ! Carbon dioxide
real(wp), parameter :: pch4int = 0.72e-6_wp         ! Methane
real(wp), parameter :: pn2oint = 0.27e-6_wp         ! Nitrous oxide


! PI atmospheric temperature for weathering rates
real(wp), parameter :: tamnpi = 11.35_wp,&          ! mn
                       taenpi = 24.76_wp,&          ! en
                       taespi = 22.77_wp,&          ! es
                       tamspi = 9.09_wp             ! ms

! Atmospheric volume meassure [mol/atm]
real(wp), parameter :: rvah   = 1.7676e20_wp/2._wp,&
                       rvahp  = rvah*(1._wp    - sin(fdh)),&
                       rvahm  = rvah*(sin(fdh) - sin(fdl)),&
                       rvahe  = rvah*sin(fdl)
! ==========================================================================
! End Atmosphere
! ==========================================================================


! ==========================================================================
! Ocean
! ==========================================================================
real(wp), parameter :: rho0 = 1028._wp,&          ! Density scale (kg/m3)
                       rq   = 1.e-4_wp,&          ! Flow resistance (1/s)
                       to_f = -1.9_wp             ! Freezing point temperature (°C)

! Horizontal diffusion
real(wp), parameter :: kh_s = 2.5e4_wp,&        ! surface diffusion coefficient (m2/s)
                       kh_d = 1.e3_wp ,&        ! deep diffusion coefficient (m2/s)
                       z_kh = 200._wp           ! e-folding length scale (m)

! Ocean volume fluxes (Sverdrup, Sv)
real(wp), parameter :: brvf  = 1.e6_wp    ,&            ! Volume flux across Bering Strait
                       evfn  = 0.e6_wp    ,&
                       evfs  = 30.e6_wp   ,&            ! Northward Southern Ocean
                       evfsa = 0.4_wp*evfs,&
                       evfsp = 0.6_wp*evfs,&
                       fw_z  = 0.15e6_wp                ! North-equatorial Atl-Pac zonal freshwater flux

! Downsloping flow
real(wp), parameter :: shl_dz = 2._wp,&                 ! Vertical resolution for mixing (m)
                       hgt0   = 100._wp                 ! Initial plume height (m)

integer, parameter :: nz_bot = int((n       - shl    )*d/shl_dz + 1),&
                      nz_all = int((n*d     - 0      )  /shl_dz + 1),&
                      nz_mid = int((zmid(n) - zmid(1))  /shl_dz + 1),&
                      nz_upp = int((d*shl   - 0      )  /shl_dz + 1)

! Shelf volume fluxes 
real(wp), parameter :: fice     = 0.13e6_wp,&           ! sea-ice export (Sv)
                       fgla_tot = 0.07e6_wp,&           ! Total glacial (freshwater) flux (Sv)
                       fgla_sh  = 0.7_wp*fgla_tot,&     ! Direct to shelf
                       fgla_so  = 0.3_wp*fgla_tot       ! Export as iceberg

! Sea ice/glacier properties
real(wp), parameter :: sice = 5._wp,&                   ! Salintiy of sea ice           ! (g kg-1)
                       sgla = 0._wp,&                   ! Salintiy of glacier water     ! (g kg-1)
                       o18fg = -40.e-3_wp               ! O18 glacier water             ! (permil)

! Biogeochemistry
real(wp), parameter :: r13pdb  = 0.0112372_wp,&                     ! 'Pee-Dee Belemnite' standard for the 13C/12C ratio
                       r14oas  = 1.176e-12_wp,&                     ! 'Oxalic Acid Standard' for the 14C/12C ratio
                       dc14    = 3.84e-12_wp                        ! C14 decay rate [s-1]

! Carbon System
real(wp), parameter :: s0      =  34.8_wp,&             ! Mean Salinity
                       ca0     =  0.01028_wp            ! Ocean mean calcium concentration, modern value 0.01028

real(wp), parameter :: knut    = 0.13_wp,&                      ! Remineralization rate (nutrient)
                       kcar    = 0.11_wp,&                      ! Remineralization rate (carbon)
                       lmdcar  = 1._wp/2500._wp,&               ! Water column dissolution scale for CaCO3 (m-1)
                       rpm     = 0.35_wp                        ! maximum rain ratio  [-] Maier-Reimer 93
! Redfield ratios
real(wp), parameter :: rcp     = 106._wp ,&                     ! Carbon to phosphate, organic matter [-]
                       rcd     = 1._wp   ,&                     ! Carbonate to DIC [-]
                       rda     = -16._wp ,&                     ! DIC to ALK [-]
                       rca     = 2._wp   ,&                     ! Carbonate to ALK [-]
                       rno     = -32._wp ,&                     ! Oxygen to nitrogen, organic matter [-],    
                       rcop    = -118._wp,&                     ! C+H to phosphate, organic matter [-], 
                       mino2   = 3.e-3_wp                       ! Min. O2 conc. for oxidation of POM, implies onset of DN [mol/m3]

! Limitation factor for new production (-)
real(wp), parameter :: lfar    = 0.3_wp,&                                       ! Arctic
                       lfso    = 0.15_wp                                        ! Southern Ocean
real(wp), parameter, dimension(5) :: lf = (/1._wp,1._wp,1._wp,1._wp,0.15_wp/)   ! [hn mn en es ms]
! ==========================================================================
! End Ocean
! ==========================================================================


! ==========================================================================
! Land Biosphere
! ==========================================================================
real(wp), parameter :: co2fer  = 0.37_wp,&              ! Fertilization factor based on C4MIP results, Friedlingstein 2006
                       q10     = 2._wp,&                ! Q10 temperature dependence for respiration rate (ocean and land) and weathering rates, st. val 2.0
                       q10_met = 2._wp,&                ! Q10 for terrestrial methane release, st. val 2.0
                       mdts    = 9.5_wp,&               ! Initial methane decay time scale (yr)
                       fch4    = 2.285971176296663_wp,&
                       fn2o    = 0.663607278017540_wp,& 
                       ta_tn0  = 24.8749_wp,&           ! Tropical initial atmosphere temperature (°C)
                       ta_gn0  = 20.0187_wp,&           ! Grassland initial atmosphere temperature (°C)
                       ta_bn0  = 7.5835_wp,&            ! Extratropical initial atmosphere temperature (°C)
                       ta_ts0  = 23.9617_wp,&
                       ta_gs0  = 18.3717_wp,&
                       ta_bs0  = 9.6277_wp

real(wp), parameter :: olftn0   = 0.7704_wp,&       ! Initial tropical forest ocean-land fraction
                       olfgn0   = 0.6497_wp,&       ! Initial grassland ocean-land fraction
                       olfbn0   = 0.4835_wp,&       ! Initial extratropical forest ocean-land fraction
                       olfts0   = 0.7706_wp,&
                       olfgs0   = 0.7957_wp,&
                       olfbs0   = 0.9641_wp,&       
                       ftropn0  = 0.1969_wp,&       ! Initial tropical forest line in rad
                       fgrassn0 = 0.6592_wp,&       ! Initial grassland line in rad
                       fsnown0  = 0.9749_wp,&       ! Initial snow line in rad
                       ftrops0  = 0.3431_wp,&
                       fgrasss0 = 0.6219_wp

! Land Areas (m2)
real(wp), parameter :: labfn0 = (1._wp-olfbn0)*h_area*(sin(fsnown0)  - sin(fgrassn0)),&     ! Area BF PI NH
                       lagsn0 = (1._wp-olfgn0)*h_area*(sin(fgrassn0) - sin(ftropn0) ),&     ! Area GSD PI Nh
                       latfn0 = (1._wp-olftn0)*h_area*(sin(ftropn0)  - sin(fde)     ),&     ! Area TF PI NH
                       latfs0 = (1._wp-olfts0)*h_area*(sin(ftrops0)  - sin(fde)     ),&     ! Area TF PI SH
                       lagss0 = (1._wp-olfgs0)*h_area*(sin(fgrasss0) - sin(ftrops0) ),&     ! Area GSD PI SH
                       labfs0 = (1._wp-olfbs0)*h_area*(sin(fdh)      - sin(fgrasss0))       ! Area BF PI SH

! Pre-industrial area factors
real(wp), parameter :: fac_tfn0 = latfn0/(latfn0 + latfs0),&
                       fac_gsn0 = lagsn0/(lagsn0 + lagss0),&
                       fac_bfn0 = labfn0/(labfn0 + labfs0),&
                       fac_tfs0 = latfs0/(latfn0 + latfs0),&
                       fac_gss0 = lagss0/(lagsn0 + lagss0),&
                       fac_bfs0 = labfs0/(labfn0 + labfs0)

! Global reservoirs
real(wp), parameter :: gr0_t   = 30._wp,&           ! Pre-Industrial(PI) leafy biomass on tropical forest, GtC
                       wo0_t   = 270._wp,&          ! PI woody biomass on tropical forest, GtC
                       li0_t   = 16._wp,&           ! PI litter biomass on tropical forest, GtC
                       sl0_t   = 200._wp,&          ! PI soil biomass on tropical forest, GtC
                       np0_t   = 25._wp,&           ! PI primary production on tropical forest,  GtC/yr
                       gr0_g   = 50._wp,&           ! Pre-Industrial(PI) leafy biomas s on Grassland, GtC
                       wo0_g   = 50._wp,&           ! PI woody biomass on Grassland, GtC
                       li0_g   = 40._wp,&           ! PI litter biomass on Grassland, GtC
                       sl0_g   = 500._wp,&          ! PI soil biomass on Grassland, GtC
                       np0_g   = 20._wp,&           ! PI primary production on Grassland,  GtC/yr
                       gr0_b   = 20._wp,&           ! Pre-Industrial(PI) leafy biomass on Grassland, GtC
                       wo0_b   = 180._wp,&          ! PI woody biomass on Grassland, GtC
                       li0_b   = 64._wp,&           ! PI litter biomass on Grassland, GtC
                       sl0_b   = 800._wp,&          ! PI soil biomass on Grassland, GtC
                       np0_b   = 15._wp             ! PI primary production on boreal forest,  GtC/yr

real(wp), parameter :: gr_tn0   = fac_tfn0*gr0_t,&
                       wo_tn0   = fac_tfn0*wo0_t,&
                       li_tn0   = fac_tfn0*li0_t,&
                       sl_tn0   = fac_tfn0*sl0_t,&
                       gr_ts0   = fac_tfs0*gr0_t,&
                       wo_ts0   = fac_tfs0*wo0_t,&
                       li_ts0   = fac_tfs0*li0_t,&
                       sl_ts0   = fac_tfs0*sl0_t,&
                       
                       gr_gn0   = fac_gsn0*gr0_g,&
                       wo_gn0   = fac_gsn0*wo0_g,&
                       li_gn0   = fac_gsn0*li0_g,&
                       sl_gn0   = fac_gsn0*sl0_g,&
                       gr_gs0   = fac_gss0*gr0_g,&
                       wo_gs0   = fac_gss0*wo0_g,&
                       li_gs0   = fac_gss0*li0_g,&
                       sl_gs0   = fac_gss0*sl0_g,&

                       gr_bn0   = fac_bfn0*gr0_b,&
                       wo_bn0   = fac_bfn0*wo0_b,&
                       li_bn0   = fac_bfn0*li0_b,&
                       sl_bn0   = fac_bfn0*sl0_b,&
                       gr_bs0   = fac_bfs0*gr0_b,&
                       wo_bs0   = fac_bfs0*wo0_b,&
                       li_bs0   = fac_bfs0*li0_b,&
                       sl_bs0   = fac_bfs0*sl0_b

! New production
real(wp), parameter :: np_tn0 = fac_tfn0*np0_t,&
                       np_gn0 = fac_gsn0*np0_g,&
                       np_bn0 = fac_bfn0*np0_b,&
                       np_ts0 = fac_tfs0*np0_t,&
                       np_gs0 = fac_gss0*np0_g,&
                       np_bs0 = fac_bfs0*np0_b

! PI CH4 and N2O production
real(wp), parameter :: lbmp0 = fch4*0.86_wp*pch4int/(rva*mgt*mdts),&    ! PI land biosphere methane production, GtC/yr
                       lbnp0 = fn2o        *pn2oint/(rva*tn2o*sy)       ! PA N2O production, mol/s

real(wp), parameter :: lbmp_tn0 = 0.5_wp*fac_tfn0*lbmp0,&
                       lbmp_gn0 = 0._wp,&
                       lbmp_bn0 = 0.5_wp*fac_bfn0*lbmp0,&
                       lbmp_ts0 = 0.5_wp*fac_tfs0*lbmp0,&
                       lbmp_gs0 = 0._wp,&
                       lbmp_bs0 = 0.5_wp*fac_bfs0*lbmp0,&
                       lbnp_tn0 = fac_tfn0*lbnp0,&
                       lbnp_gn0 = fac_gsn0*lbnp0,&
                       lbnp_bn0 = fac_bfn0*lbnp0,&
                       lbnp_ts0 = fac_tfs0*lbnp0,&
                       lbnp_gs0 = fac_gss0*lbnp0,&
                       lbnp_bs0 = fac_bfs0*lbnp0
! ==========================================================================
! End Land Biosphere
! ==========================================================================

! ==========================================================================
! External forcing
! ==========================================================================
! River input factors
real(wp), dimension(5), parameter :: friva  = (/0.0811_wp,0.1686_wp,0.1103_wp,0.1632_wp,0.0054_wp/),&   ! Pacific
                                     frivp  = (/0.0184_wp,0.0562_wp,0.1589_wp,0.1059_wp,0.0562_wp/)     ! Atlantic
real(wp), parameter :: frivar = 0.0758_wp,&                                                             ! Arctic
                       frivso = 0._wp                                                                   ! Southern Ocean

real(wp), parameter :: bcarpi     = 5.6410e5_wp,&               ! global carbonate burial (from calibration steps)
                       borgpi     = 2.5126e3_wp,&               ! global PO4 burial
                       bcorgpi    = 2.7934e5_wp,&               ! global organic carbon burial
                       d13ccarpi  = 2.2206e-3_wp,&              ! d13C carbonate
                       d13corgcpi = -22.0611e-3_wp,&            ! d13C organic carbon
                       d13volpi   = -0.005_wp,&                 ! d13C volcanic emissions
                       gamma_sil  = 0.85_wp,&                   ! carbonate to silicate weathering ratio
                       volpi      = gamma_sil/(1._wp+gamma_sil)*bcarpi*(d13ccarpi - d13corgcpi)/(d13volpi - d13corgcpi) ! volcanic input
real(wp), dimension(6), parameter :: fvol = (/0.12_wp,0.24_wp,0.33_wp,0.24_wp,0.07_wp,0._wp/)        ! meridional volcanoes distribution
! ==========================================================================
! End External Forcing
! ==========================================================================

! ==========================================================================
! Ocean Sediment
! ==========================================================================
real(wp), parameter :: ncf      = 0.3_wp,&          ! Open ocean, non-Calcite flux to sediment,Low Lat., st. val 0.3 [g/(cm2*kyr)]
                       caf      = 20._wp,&          ! Amplification factor for NCF at the "coast"
                       kcalcite = 0.0015_wp,&       ! CaCO3 dissolution rate coefficient, st. val 0.001 [1/day];
                       mmca     = 100._wp,&         ! Molar weight of CaCO3 [g/mol]
                       rhom     = 2.7_wp            ! Density of mineral fraction in sediments [g/cm3]
! ==========================================================================
! End Ocean Sediment
! ==========================================================================

end module parameters