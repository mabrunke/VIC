  subroutine vic2clmtype(nlevgrnd, rec, nrec, adspinup, init_state, yr, mo, day, &
	secs, jday, yr1, jday1, dt, lat, lon, begg, endg, begc, endc, begp, &
	endp, num_soilc, num_soilp, vic_cn) BIND(C, &
	name='vic2clmtype')

! 03/06/2014  Added sending of CN variables for inclusion in VIC 
!             state files.                                          MAB
!
! !USES:
   use, intrinsic :: ISO_C_BINDING
   use shr_kind_mod, only: r8 => shr_kind_r8
   use shr_const_mod, only: SHR_CONST_PI, SHR_CONST_TKFRZ
   use clmtype
   use clm_varpar, only : max_bands, max_layers, max_nodes, max_veg, mxpft, numrad
   use shr_orb_mod
   use SurfaceAlbedoMod, only : SurfaceAlbedo
   use SurfaceRadiationMod, only : SurfaceRadiation
   use CanopyFluxesMod, only : CanopyFluxes
   use CNEcosystemDynMod, only : CNEcosystemDyn
   use CNAnnualUpdateMod, only : CNAnnualUpdate
   use CNBalanceCheckMod, only : CBalanceCheck, NBalanceCheck
!
! !ARGUMENTS
   implicit none

   integer, intent(in) :: nlevgrnd ! Number of below ground layers
   integer, intent(in) :: yr       ! Current year
   integer, intent(in) :: mo       ! Current month
   integer, intent(in) :: day      ! Current day
   integer, intent(in) :: secs     ! Current time of the day in seconds
   real(r8), intent(in) :: jday     ! Julian day
   integer, intent(in) :: yr1      ! Next time step year
   real(r8), intent(in) :: jday1   ! Next time step Julian day
   real(r8), intent(in) :: dt      ! Time step in seconds
   integer, intent(in) :: rec   ! Record number
   integer, intent(in) :: nrec  ! Total number of records
!   integer, intent(in) :: nspinup ! Year in spin-up
   integer, intent(in) :: adspinup ! Flag for AD spin-up
   integer , intent(in)    :: init_state  ! flag to determine where initial state comes from
   real(r8), intent(in) :: lat  ! Latitude (degrees N)
   real(r8), intent(in) :: lon  ! Longitude (degrees E)
   integer, intent(in) :: begg	! Beginning grid cell index = 1
   integer, intent(in) :: endg  ! Ending grid cell index = 1
   integer, intent(inout) :: begc  ! Beginning precip. dist. index = 1
   integer, intent(inout) :: endc  ! Ending precip. dist. index = Ndist
   integer, intent(inout) :: begp  ! Beginning veg. type index = 1
   integer, intent(inout) :: endp  ! Ending veg. type index = Nveg
   integer, intent(in) :: num_soilc ! number of soil columns in filter
   integer, intent(in) :: num_soilp ! number of soil pfts in filter

    TYPE, BIND(C) :: vic_cn_data_type
   	! atmospheric quantities
   	REAL(C_DOUBLE) :: Tair                               ! air temperature
   	REAL(C_DOUBLE) :: vp                                 ! vapor pressure
   	REAL(C_DOUBLE) :: vpd                                ! vapor pressure depression
   	REAL(C_DOUBLE) :: psfc                               ! surface pressure
   	REAL(C_DOUBLE) :: lwrad                              ! longwave radiation
   	REAL(C_DOUBLE) :: swrad                              ! shortwave radiation
   	REAL(C_DOUBLE) :: precip                             ! precipitation
   	REAL(C_DOUBLE) :: swrd(numrad)                            ! direct shortwave radiation
   	REAL(C_DOUBLE) :: swri(numrad)                            ! diffuse shortwave radiation
   	REAL(C_DOUBLE) :: alb                                ! albedo
   	REAL(C_DOUBLE) :: t2m(0:mxpft)                            ! 2-m air temperature
   	REAL(C_DOUBLE) :: t_soisno(max_nodes+2)                        ! soil/snow temperature
   	! Soil properties
   	REAL(C_DOUBLE) :: z(max_nodes+2)                               ! Node depth
   	REAL(C_DOUBLE) :: dz(max_nodes+2)                              ! Node thickness
   	REAL(C_DOUBLE) :: z0                                 ! Surface roughness
   	REAL(C_DOUBLE) :: z0s                                ! Snow roughness
   	! Vegetation characteristics
   	REAL(C_DOUBLE) :: Tveg(0:mxpft)                           ! vegetation temperature
   	REAL(C_DOUBLE) :: rveg(0:mxpft)                           ! vegetation resistance
   	REAL(C_DOUBLE) :: zov(0:mxpft)                            ! vegetation roughness
   	REAL(C_DOUBLE) :: displ(0:mxpft)                          ! displacement height
   	REAL(C_DOUBLE) :: fwet(0:mxpft)                           ! wet veg fraction
   	! Soil hydrology
   	REAL(C_DOUBLE) :: baseflow                           ! baseflow
   	REAL(C_DOUBLE) :: moist(max_nodes)                           ! soil moisture
   	REAL(C_DOUBLE) :: ice(max_nodes+2)                               ! soil ice
   	REAL(C_DOUBLE) :: rootfr(max_nodes,0:mxpft)                    ! root fraction
   	REAL(C_DOUBLE) :: bsw(2,max_nodes)                           ! Clapp-Hornberger coefficient
   	REAL(C_DOUBLE) :: sucsat(2,max_nodes)                        ! saturated suction
   	REAL(C_DOUBLE) :: soisuc(2,max_nodes)                        ! soil matrix potential
   	REAL(C_DOUBLE) :: snowdep                            ! snow depth
        ! PFT-level ecophysiological variables
        REAL(C_DOUBLE) :: LAI(0:mxpft)                   ! leaf area index
        REAL(C_DOUBLE) :: dormant_flag(0:mxpft)          ! dormancy flag
        REAL(C_DOUBLE) :: days_active(0:mxpft)           ! # days since last dormancy
        REAL(C_DOUBLE) :: onset_flag(0:mxpft)            ! onset flag
        REAL(C_DOUBLE) :: onset_counter(0:mxpft)         ! onset days counter
        REAL(C_DOUBLE) :: onset_gddflag(0:mxpft)         ! onset flag for growing deg day sum
        REAL(C_DOUBLE) :: onset_fdd(0:mxpft)             ! onset freezing deg day counter
        REAL(C_DOUBLE) :: onset_gdd(0:mxpft)             ! onset growing degree days
        REAL(C_DOUBLE) :: onset_swi(0:mxpft)             ! onset soil water index
        REAL(C_DOUBLE) :: offset_flag(0:mxpft)           ! offset flag
        REAL(C_DOUBLE) :: offset_counter(0:mxpft)        ! offset days counter
        REAL(C_DOUBLE) :: offset_fdd(0:mxpft)            ! offset freezing deg day counter
        REAL(C_DOUBLE) :: offset_swi(0:mxpft)            ! offset soil water index
        REAL(C_DOUBLE) :: lgsf(0:mxpft)                  ! long growing season factor
        REAL(C_DOUBLE) :: bglfr(0:mxpft)                 ! background litterfall rate (1/s)
        REAL(C_DOUBLE) :: bgtr(0:mxpft)                  ! background transfer growth rate (1/s)
        REAL(C_DOUBLE) :: dayl(0:mxpft)                  ! daylength (s)
        REAL(C_DOUBLE) :: prev_dayl(0:mxpft)             ! daylength at previous timestep (s)
        REAL(C_DOUBLE) :: annavg_t2m(0:mxpft)            ! annual average 2-m air temperature (K)
        REAL(C_DOUBLE) :: tempavg_t2m(0:mxpft)           ! temporary average 2-m air temperature (K)
        REAL(C_DOUBLE) :: gpp2(0:mxpft)                  ! GPP flux before downregulation (g C/m^2/s)
        REAL(C_DOUBLE) :: availc(0:mxpft)                ! C flux available for allocation (g C/m^2/s)
        REAL(C_DOUBLE) :: xsmrpool_recover(0:mxpft)      ! C flux assigned to recovery (g C/m^2/s)
        REAL(C_DOUBLE) :: alloc_pnow(0:mxpft)            ! fraction of current allocation as new growth
        REAL(C_DOUBLE) :: c_allometry(0:mxpft)           ! C allocation index
        REAL(C_DOUBLE) :: n_allometry(0:mxpft)           ! N allocation index
        REAL(C_DOUBLE) :: plant_ndemand(0:mxpft)         ! N flux required to support GPP (g N/m^2/s)
        REAL(C_DOUBLE) :: tempsum_potential_gpp(0:mxpft) ! temporary annual sum of potential GPP
        REAL(C_DOUBLE) :: annsum_potential_gpp(0:mxpft)  ! annuals sum of potential GPP
        REAL(C_DOUBLE) :: tempmax_retransn(0:mxpft)      ! temporary annual max of retrans N pool (g N/m^2)
        REAL(C_DOUBLE) :: annmax_retransn(0:mxpft)       ! annual max of retransloc N pool (g N/m^2)
        REAL(C_DOUBLE) :: avail_retransn(0:mxpft)        ! N flux avail for retransloc (g N/m^2/s)
        REAL(C_DOUBLE) :: plant_nalloc(0:mxpft)          ! total allocated N flux (g N/m^2/s)
        REAL(C_DOUBLE) :: plant_calloc(0:mxpft)          ! total allocated C flux (g C/m^2/s)
        REAL(C_DOUBLE) :: excess_cflux(0:mxpft)          ! C flux not allocated (g C/m^2/s)
        REAL(C_DOUBLE) :: downreg(0:mxpft)               ! fract reduction in GPP due to N limit
        REAL(C_DOUBLE) :: prev_leafc_to_litter(0:mxpft)  ! previous leaf C litterfall (g C/m^2/s)
        REAL(C_DOUBLE) :: prev_frootc_to_litter(0:mxpft) ! previous froot C litterfall (g C/m^2/s)
        REAL(C_DOUBLE) :: tempsum_npp(0:mxpft)           ! temporary annual sum of NPP (g C/m^2/yr)
        REAL(C_DOUBLE) :: annsum_npp(0:mxpft)            ! annual sum of NPP (g C/m^2/yr)
        REAL(C_DOUBLE) :: gpp(0:mxpft)                   ! gross primary production (g C/m^2/s)
        REAL(C_DOUBLE) :: npp(0:mxpft)                   ! net primary production (g C/m^2/s)
        REAL(C_DOUBLE) :: ar(0:mxpft)                    ! autotrophic respiration (g C/m^2/s)
   	REAL(C_DOUBLE) :: gr(0:mxpft)                    ! growth respiration (g C/m^2/s)
   	REAL(C_DOUBLE) :: mr(0:mxpft)                    ! maintenance respiration (g C/m^2/s)
   	REAL(C_DOUBLE) :: leaf_mr(0:mxpft)               ! leaf maintenance respiration (g C/m^2/s)
   	REAL(C_DOUBLE) :: litfall(0:mxpft)               ! litterfall (g C/m^2/s)
   	REAL(C_DOUBLE) :: rs(0:mxpft)                    ! stomatal resistance
   	REAL(C_DOUBLE) :: par(0:mxpft)                   ! absorbed PAR
        ! PFT-level carbon state
        REAL(C_DOUBLE) :: leafc(0:mxpft)                 ! leaf C (g C/m^2)
        REAL(C_DOUBLE) :: leafc_storage(0:mxpft)         ! leaf C storage (g C/m^2)
        REAL(C_DOUBLE) :: leafc_xfer(0:mxpft)            ! leaf C transfer (g C/m^2)
        REAL(C_DOUBLE) :: frootc(0:mxpft)                ! fine root C (g C/m^2)
        REAL(C_DOUBLE) :: frootc_storage(0:mxpft)        ! fine root C storage (g C/m^2)
        REAL(C_DOUBLE) :: frootc_xfer(0:mxpft)           ! fine root C transfer (g C/m^2)
        REAL(C_DOUBLE) :: livestemc(0:mxpft)             ! live stem C (g C/m^2)
        REAL(C_DOUBLE) :: livestemc_storage(0:mxpft)     ! live stem C storage (g C/m^2)
        REAL(C_DOUBLE) :: livestemc_xfer(0:mxpft)        ! live stem C transfer (g C/m^2)
        REAL(C_DOUBLE) :: deadstemc(0:mxpft)             ! dead stem C (g C/m^2)
        REAL(C_DOUBLE) :: deadstemc_storage(0:mxpft)     ! dead stem C storage (g C/m^2)
        REAL(C_DOUBLE) :: deadstemc_xfer(0:mxpft)        ! dead stem C transfer (g C/m^2)
        REAL(C_DOUBLE) :: livecrootc(0:mxpft)            ! live coarse root C (g C/m^2)
        REAL(C_DOUBLE) :: livecrootc_storage(0:mxpft)    ! live coarse root C storage (g C/m^2)
        REAL(C_DOUBLE) :: livecrootc_xfer(0:mxpft)       ! live coarse root C transfer (g C/m^2)
        REAL(C_DOUBLE) :: deadcrootc(0:mxpft)            ! dead coarse root C (g C/m^2)
        REAL(C_DOUBLE) :: deadcrootc_storage(0:mxpft)    ! dead coarse root C storage (g C/m^2)
        REAL(C_DOUBLE) :: deadcrootc_xfer(0:mxpft)       ! dead coarse root C transfer (g C/m^2)
        REAL(C_DOUBLE) :: gresp_storage(0:mxpft)         ! growth respiration storage (g C/m^2)
        REAL(C_DOUBLE) :: gresp_xfer(0:mxpft)            ! growth respiration transfer (g C/m^2)
        REAL(C_DOUBLE) :: cpool(0:mxpft)                 ! temporary photosynthate C pool (g C/m^2)
        REAL(C_DOUBLE) :: xsmrpool(0:mxpft)              ! abstract C pool to meet excess MR demand (g C/m^2)
        REAL(C_DOUBLE) :: pft_ctrunc(0:mxpft)            ! PFT-level sink for C truncation (g C/m^2)
        REAL(C_DOUBLE) :: totvegc(0:mxpft)               ! total vegetation C (g C/m^2)
        REAL(C_DOUBLE) :: woodc(0:mxpft)                 ! wood C (g C/m^2)
    	REAL(C_DOUBLE) :: fpsn(0:mxpft)                  ! photosynthesis
   	REAL(C_DOUBLE) :: ci(0:mxpft)                    ! intracellular CO2
       ! PFT-level nitrogen state
        REAL(C_DOUBLE) :: leafn(0:mxpft)                 ! leaf N (g N/m^2)
        REAL(C_DOUBLE) :: leafn_storage(0:mxpft)         ! leaf N storage (g N/m^2)
        REAL(C_DOUBLE) :: leafn_xfer(0:mxpft)            ! leaf N transfer (g N/m^2)
        REAL(C_DOUBLE) :: frootn(0:mxpft)                ! fine root N (g N/m^2)
        REAL(C_DOUBLE) :: frootn_storage(0:mxpft)        ! fine root N storage (g N/m^2)
        REAL(C_DOUBLE) :: frootn_xfer(0:mxpft)           ! fine root N transfer (g N/m^2)
        REAL(C_DOUBLE) :: livestemn(0:mxpft)             ! live stem N (g N/m^2)
        REAL(C_DOUBLE) :: livestemn_storage(0:mxpft)     ! live stem N storage (g N/m^2)
        REAL(C_DOUBLE) :: livestemn_xfer(0:mxpft)        ! live stem N transfer (g N/m^2)
        REAL(C_DOUBLE) :: deadstemn(0:mxpft)             ! dead stem N (g N/m^2)
        REAL(C_DOUBLE) :: deadstemn_storage(0:mxpft)     ! dead stem N storage (g N/m^2)
        REAL(C_DOUBLE) :: deadstemn_xfer(0:mxpft)        ! dead stem N transfer (g N/m^2)
        REAL(C_DOUBLE) :: livecrootn(0:mxpft)            ! live coarse root N (g N/m^2)
        REAL(C_DOUBLE) :: livecrootn_storage(0:mxpft)    ! live coarse root N storage (g N/m^2)
        REAL(C_DOUBLE) :: livecrootn_xfer(0:mxpft)       ! live coarse root N transfer (g N/m^2)
        REAL(C_DOUBLE) :: deadcrootn(0:mxpft)            ! dead coarse root N (g N/m^2)
        REAL(C_DOUBLE) :: deadcrootn_storage(0:mxpft)    ! dead coarse root N storage (g N/m^2)
        REAL(C_DOUBLE) :: deadcrootn_xfer(0:mxpft)       ! dead coarse root N transfer (g N/m^2)
        REAL(C_DOUBLE) :: retransn(0:mxpft)              ! retranslocated N (g N/m^2)

        REAL(C_DOUBLE) :: npool(0:mxpft)                 ! temporary photosynthate N pool (g N/m^2)
        REAL(C_DOUBLE) :: pft_ntrunc(0:mxpft)            ! PFT-level sink for N truncation (g N/m^2)
        ! column (band) physical state
        REAL(C_DOUBLE) :: decl                      ! solar declination angle (radians)
        REAL(C_DOUBLE) :: fpi                       ! fraction of potential immobilization
        REAL(C_DOUBLE) :: fpg                       ! fraction of potential GPP
        REAL(C_DOUBLE) :: annsum_counter            ! seconds since last ann accumulation turnover
        REAL(C_DOUBLE) :: cannsum_npp               ! annual sum of NPP, averaged from PFT-level (g C/m^2/yr)
        REAL(C_DOUBLE) :: cannavg_t2m               ! annual avg. of 2-m air temperature, averaged from PFT-level (K)
        REAL(C_DOUBLE) :: watfc(max_nodes)                  ! volumetric soil water at field capacity
        REAL(C_DOUBLE) :: me                        ! moisture of extinction
        REAL(C_DOUBLE) :: fire_prob                 ! daily fire probability
        REAL(C_DOUBLE) :: mean_fire_prob            ! e-folding mean of daily fire prob.
        REAL(C_DOUBLE) :: fireseasonl               ! annual fire season length (days)
        REAL(C_DOUBLE) :: farea_burned              ! timestep fractional area burned
        REAL(C_DOUBLE) :: ann_farea_burned          ! annual total fract. area burned
        REAL(C_DOUBLE) :: hr                        ! heterotrophic respiration (g C/m^2/s)
   	REAL(C_DOUBLE) :: lithr                     ! litter heterotrophic resp. (g C/m^2/s)
        REAL(C_DOUBLE) :: nee                       ! net ecosystem exchange (g C/m^2/s)
        REAL(C_DOUBLE) :: nep                       ! net ecosystem production (g C/m^2/s)
        ! column (band) carbon state
        REAL(C_DOUBLE) :: cwdc                      ! coarse woody debris C (g C/m^2)
        REAL(C_DOUBLE) :: litr1c                    ! litter labile C (g C/m^2)
        REAL(C_DOUBLE) :: litr2c                    ! litter cellulose C (g C/m^2)
        REAL(C_DOUBLE) :: litr3c                    ! litter lignin C (g C/m^2)
        REAL(C_DOUBLE) :: soil1c                    ! fastest soil organic matter C
        REAL(C_DOUBLE) :: soil2c                    ! medium soil organic matter C
        REAL(C_DOUBLE) :: soil3c                    ! slow soil organic matter C
        REAL(C_DOUBLE) :: soil4c                    ! slowest soil organic matter C
        REAL(C_DOUBLE) :: seedc                     ! column-lev pool for seeding new PFTs
        REAL(C_DOUBLE) :: col_ctrunc                ! column-lev sink for C truncation
        REAL(C_DOUBLE) :: totlitc                   ! total litter C (g C/m^2)
        REAL(C_DOUBLE) :: totsomc                   ! total soil organic C (g C/m^2)
        REAL(C_DOUBLE) :: totcolc                   ! total column C (g C/m^2)
        REAL(C_DOUBLE) :: prod10c                   ! wood product C pool, 10-yr lifespan (g C/m^2)
        REAL(C_DOUBLE) :: prod100c                  ! wood product C pool, 100-yr lifespan (g C/m^2)
        ! column (band) nitrogen state
        REAL(C_DOUBLE) :: cwdn                      ! coarse woody debris N (g N/m^2)
        REAL(C_DOUBLE) :: litr1n                    ! litter labile N (g N/m^2)
        REAL(C_DOUBLE) :: litr2n                    ! litter cellulose N (g N/m^2)
        REAL(C_DOUBLE) :: litr3n                    ! litter lignin N (g N/m^2)
        REAL(C_DOUBLE) :: soil1n                    ! fastest soil organic matter N
        REAL(C_DOUBLE) :: soil2n                    ! medium soil organic matter N
        REAL(C_DOUBLE) :: soil3n                    ! slow soil organic matter N
        REAL(C_DOUBLE) :: soil4n                    ! slowest soil organic matter N
        REAL(C_DOUBLE) :: sminn                     ! soil mineral N (g N/m^2)
        REAL(C_DOUBLE) :: seedn                     ! column-lev pool for seeding new PFTs
        REAL(C_DOUBLE) :: col_ntrunc                ! column-lev sink for N truncation
        REAL(C_DOUBLE) :: totcoln                   ! total column N (g N/m^2)
        REAL(C_DOUBLE) :: prod10n                   ! wood product N pool, 10-yr lifespan (g N/m^2)
        REAL(C_DOUBLE) :: prod100n                  ! wood product N pool, 100-yr lifespan (g N/m^2)
    END TYPE vic_cn_data_type

   type(c_ptr), value :: vic_cn
   type(vic_cn_data_type), dimension(:), pointer :: cn

   real(r8), pointer :: z(:,:)          ! layer depth (m)
   real(r8), pointer :: dz(:,:)         ! layer thickness (m)
   real(r8), pointer :: qflx_drain(:)   ! sub-surface runoff (mm H20/s)
   real(r8), pointer :: h2osoi_liq(:,:) ! liquid water (kg / m2)
   real(r8), pointer :: h2osoi_ice(:,:) ! ice lens (kg / m2)
   real(r8), pointer :: t_soisno(:,:)     ! soil temperature (K)
   real(r8), pointer :: t_ref2m(:)      ! 2 m height temperature (K) 
   real(r8), pointer :: t_veg(:)        ! vegetation temperature (K)
   real(r8), pointer :: fwet(:)       ! wet vegetation fraction
   real(r8), pointer :: snowdp(:)       ! snow depth (m)
   real(r8), pointer :: rootfr(:,:)     ! fraction of roots in ea. soil layer
   real(r8), pointer :: psisat(:,:)     ! soil water at saturation (MPa) for CN
   real(r8), pointer :: soilpsi(:,:)    ! soil water potential (MPa) for CN
   real(r8), pointer :: bsw2(:,:)        ! Clapp-Hornberger coefficient for CN
   real(r8), pointer :: bsw(:,:)        ! Clapp-Hornberger coefficient
   real(r8), pointer :: sucsat(:,:)     ! minimum soil suction (mm)
   real(r8), pointer :: forc_pbot(:)    ! atmospheric pressure (Pa)
   real(r8), pointer :: forc_t(:)       ! atmospheric temperatuer (K)
   real(r8), pointer :: forc_q(:)       ! atmos. spec. hum. (kg/kg)
   real(r8), pointer :: forc_vp(:)      ! atmos. vapor pressure (Pa)
   real(r8), pointer :: qg(:)           ! surface spec. hum. (kg/kg)
   real(r8), pointer :: thm(:)          ! forc_t + 0.0098 * forc_hgt_t_pft
   real(r8), pointer :: forc_hgt_u(:)   ! wind forcing height (m)
   real(r8), pointer :: forc_hgt_t(:)   ! temperature forcing height (m)
   real(r8), pointer :: forc_hgt_q(:)   ! humidity forcing height (m)
   real(r8), pointer :: forc_hgt_u_pft(:) ! forc_hgt_u + z0m + d
   real(r8), pointer :: forc_hgt_t_pft(:) ! forc_hgt_t + z0m + d
   real(r8), pointer :: forc_hgt_q_pft(:) ! forc_hgt_q + z0m + d
   real(r8), pointer :: forc_pco2(:)    ! CO2 partial pressure (Pa)
   real(r8), pointer :: forc_po2(:)     ! O2 partial pressure (Pa)
   real(r8), pointer :: forc_lwrad(:)   ! longwave radiation (W/m^2)
   real(r8), pointer :: forc_solar(:)   ! shortwave radiation (W/m^2)
   real(r8), pointer :: forc_solad(:,:) ! direct SW radiation (W/m^2)
   real(r8), pointer :: forc_solai(:,:) ! diffuse SW radiation (W/m^2)
   real(r8), pointer :: soil1c(:)       ! soil organic matter C (fast pool)
   real(r8), pointer :: soil2c(:)       ! soil organic matter C (medium pool)
   real(r8), pointer :: soil3c(:)       ! soil organic matter C (slow pool)
   real(r8), pointer :: soil4c(:)       ! soil organic matter C (slowest pool)
   real(r8), pointer :: litr1c(:)       ! litter labile C
   real(r8), pointer :: litr2c(:)       ! litter cellulose C
   real(r8), pointer :: litr3c(:)       ! litter lignin C
   real(r8), pointer :: cwdc(:)         ! coarse woody debris C
   real(r8), pointer :: leafc(:)        ! leaf C
   real(r8), pointer :: frootc(:)       ! fine root C
   real(r8), pointer :: livestemc(:)    ! live stem C
   real(r8), pointer :: deadstemc(:)    ! dead stem C
   real(r8), pointer :: livecrootc(:)   ! live coarse root C
   real(r8), pointer :: deadcrootc(:)   ! dead coarse root C
   real(r8), pointer :: woodc(:)        ! wood C
   real(r8), pointer :: totvegc(:)      ! total vegetation C (g C/m^2)
   real(r8), pointer :: totlitc(:)      ! total litter C (g C/m^2)
   real(r8), pointer :: totsomc(:)      ! total SOM C (g C/m^2)
   real(r8), pointer :: soil1n(:)       ! soil organic matter N (fast pool)
   real(r8), pointer :: soil2n(:)       ! soil organic matter N (medium pool)
   real(r8), pointer :: soil3n(:)       ! soil organic matter N (slow pool)
   real(r8), pointer :: soil4n(:)       ! soil organic matter N (slowest pool)
   real(r8), pointer :: sminn(:)        ! soil mineral N
   real(r8), pointer :: litr1n(:)       ! litter labile N
   real(r8), pointer :: litr2n(:)       ! litter cellulose N
   real(r8), pointer :: litr3n(:)       ! litter lignin N
   real(r8), pointer :: cwdn(:)         ! coarse woody debris N
   real(r8), pointer :: leafn(:)        ! leaf N
   real(r8), pointer :: frootn(:)       ! fine root N
   real(r8), pointer :: livestemn(:)    ! live stem N
   real(r8), pointer :: deadstemn(:)    ! dead stem N
   real(r8), pointer :: livecrootn(:)   ! live coarse root N
   real(r8), pointer :: deadcrootn(:)   ! dead coarse root N
   real(r8), pointer :: cpool_to_leafc(:)
   real(r8), pointer :: gpp2(:)         ! GPP before downregulation
   real(r8), pointer :: gpp(:)          ! GPP
   real(r8), pointer :: npp(:)          ! NPP
   real(r8), pointer :: fpsn(:)         ! photosynthesis
   real(r8), pointer :: psnsun(:)       ! sunlit photosynthesis
   real(r8), pointer :: psnsha(:)       ! shaded photosynthesis
   real(r8), pointer :: leaf_mr(:)      ! leaf maintenance respiration
   real(r8), pointer :: mr(:)           ! maintenance respiration
   real(r8), pointer :: gr(:)           ! growth respiration
   real(r8), pointer :: ar(:)           ! autotrophic respiration
   real(r8), pointer :: hr(:)           ! heterotrophic respiration
   real(r8), pointer :: lithr(:)        ! litter hetero respiration
   real(r8), pointer :: nee(:)          ! NEE
   real(r8), pointer :: nep(:)          ! NEP
   real(r8), pointer :: dormant_flag(:) ! dormant flag
   real(r8), pointer :: days_active(:)  ! # days since dormancy
   real(r8), pointer :: onset_flag(:)   ! onset flag
   real(r8), pointer :: onset_counter(:) ! onset counter
   real(r8), pointer :: onset_gddflag(:) ! onset flag for growing deg days
   real(r8), pointer :: onset_fdd(:)    ! onset freezing deg days counter
   real(r8), pointer :: onset_gdd(:)    ! onset growing deg days counter
   real(r8), pointer :: onset_swi(:)    ! onset soil water index
   real(r8), pointer :: offset_flag(:)  ! offset flag
   real(r8), pointer :: offset_counter(:) ! onset counter
   real(r8), pointer :: offset_fdd(:)   ! onset freezing deg days counter
   real(r8), pointer :: offset_swi(:)   ! onset soil water index
   real(r8), pointer :: lgsf(:)         ! long growing season factor
   real(r8), pointer :: bglfr(:)        ! background litterfall rate
   real(r8), pointer :: bgtr(:)         ! background transfer growth rate
   real(r8), pointer :: dayl(:)         ! daylength
   real(r8), pointer :: prev_dayl(:)    ! previous timestep daylength
   real(r8), pointer :: annavg_t2m(:)   ! annual avg 2-m air temperature
   real(r8), pointer :: tempavg_t2m(:)  ! temp avg 2-ma air temperature
   real(r8), pointer :: availc(:)       ! C flux for allocation
   real(r8), pointer :: xsmrpool_recover(:) ! C flux for recovery
   real(r8), pointer :: alloc_pnow(:)   ! fract of alloc for new growth
   real(r8), pointer :: c_allometry(:)  ! C allocation index
   real(r8), pointer :: n_allometry(:)  ! N allocation index
   real(r8), pointer :: plant_ndemand(:) ! N flux for initial GPP
   real(r8), pointer :: tempsum_potential_gpp(:) ! temp ann sum of pot GPP
   real(r8), pointer :: annsum_potential_gpp(:) ! ann sum of potential GPP
   real(r8), pointer :: tempmax_retransn(:) ! temp ann max retransloc N
   real(r8), pointer :: annmax_retransn(:) ! ann max of retranslocated N
   real(r8), pointer :: avail_retransn(:) ! N flux avail for retranslocation
   real(r8), pointer :: plant_nalloc(:) ! total allocated N flux
   real(r8), pointer :: plant_calloc(:) ! total allocated C flux
   real(r8), pointer :: excess_cflux(:) ! C flux not allocated
   real(r8), pointer :: downreg(:)      ! fract reduct in GPP from N limit
   real(r8), pointer :: prev_leafc_to_litter(:) ! prev leaf C to litterfall
   real(r8), pointer :: prev_frootc_to_litter(:) ! prev fine root C to litter
   real(r8), pointer :: tempsum_npp(:)  ! temp ann sum of NPP
   real(r8), pointer :: annsum_npp(:)   ! annual sum of NPP
   real(r8), pointer :: leafc_storage(:) ! leaf C storage
   real(r8), pointer :: leafc_xfer(:)   ! leaf C transfer
   real(r8), pointer :: frootc_storage(:) ! fine root C storage
   real(r8), pointer :: frootc_xfer(:)  ! fine root C transfer
   real(r8), pointer :: livestemc_storage(:) ! live stem C storage
   real(r8), pointer :: livestemc_xfer(:) ! live stem C transfer
   real(r8), pointer :: deadstemc_storage(:) ! dead stem C storage
   real(r8), pointer :: deadstemc_xfer(:) ! dead stem C transfer
   real(r8), pointer :: livecrootc_storage(:) ! live coarse root C storage
   real(r8), pointer :: livecrootc_xfer(:) ! live coarse root C transfer
   real(r8), pointer :: deadcrootc_storage(:) ! dead coarse root C storage
   real(r8), pointer :: deadcrootc_xfer(:) ! dead coarse root C transfer
   real(r8), pointer :: gresp_storage(:) ! growth respiration storage
   real(r8), pointer :: gresp_xfer(:)   ! growth respiration transfer
   real(r8), pointer :: cpool(:)        ! temp photosynthate C pool
   real(r8), pointer :: xsmrpool(:)     ! abstract C pool for excess MR demand
   real(r8), pointer :: pft_ctrunc(:)   ! PFT-lev sink for C truncation
   real(r8), pointer :: leafn_storage(:) ! leaf N storage
   real(r8), pointer :: leafn_xfer(:)   ! leaf N transfer
   real(r8), pointer :: frootn_storage(:) ! fine root N storage
   real(r8), pointer :: frootn_xfer(:)  ! fine root N transfer
   real(r8), pointer :: livestemn_storage(:) ! live stem N storage
   real(r8), pointer :: livestemn_xfer(:) ! live stem N transfer
   real(r8), pointer :: deadstemn_storage(:) ! dead stem N storage
   real(r8), pointer :: deadstemn_xfer(:) ! dead stem N transfer
   real(r8), pointer :: livecrootn_storage(:) ! live coarse root N storage
   real(r8), pointer :: livecrootn_xfer(:) ! live coarse root N transfer
   real(r8), pointer :: deadcrootn_storage(:) ! dead coarse root N storage
   real(r8), pointer :: deadcrootn_xfer(:) ! dead coarse root N transfer
   real(r8), pointer :: retransn(:)     ! retranslocated N
   real(r8), pointer :: npool(:)        ! temp photosynthate N pool
   real(r8), pointer :: pft_ntrunc(:)   ! PFT-lev sink for N truncation
   real(r8), pointer :: fpi(:)          ! fract of potential immobilization
   real(r8), pointer :: fpg(:)          ! fract of potential GPP
   real(r8), pointer :: annsum_counter(:) ! ann sum turnover counter
   real(r8), pointer :: cannsum_npp(:)  ! col-avg ann sum of NPP
   real(r8), pointer :: cannsum_t2m(:)  ! col-avg ann sum of 2-m air temp
   real(r8), pointer :: watfc(:,:)      ! vol soil water at field capacity
   real(r8), pointer :: me(:)           ! moisture of extinction
   real(r8), pointer :: fire_prob(:)    ! daily fire probability
   real(r8), pointer :: mean_fire_prob(:) ! mean of daily fire probability
   real(r8), pointer :: fireseasonl(:)  ! fire season length
   real(r8), pointer :: farea_burned(:) ! fract area burned
   real(r8), pointer :: ann_farea_burned(:) ! ann sum of fract area burned
   real(r8), pointer :: seedc(:)        ! C pool for seeding new PFTs
   real(r8), pointer :: col_ctrunc(:)   ! col-lev sink for C truncation
   real(r8), pointer :: totcolc(:)      ! total column C
   real(r8), pointer :: prod10c(:)      ! 10-yr wood product C
   real(r8), pointer :: prod100c(:)     ! 100-yr wood product C
   real(r8), pointer :: seedn(:)        ! N pool for seeding new PFTs
   real(r8), pointer :: col_ntrunc(:)   ! col-lev sink for N truncation
   real(r8), pointer :: totcoln(:)      ! total column N
   real(r8), pointer :: prod10n(:)      ! 10-yr wood product N
   real(r8), pointer :: prod100n(:)     ! 100-yr wood product N
   real(r8), pointer :: tlai(:)         ! total leaf area index
   real(r8), pointer :: elai(:)         ! effective leaf area index
   real(r8), pointer :: litfall(:)      ! litterfall
   integer, pointer :: ivt(:)          ! PFT index
   integer, pointer :: ncol(:)          ! column index
   real(r8), pointer :: latrad(:)       ! latitude (radians)
   real(r8), pointer :: lonrad(:)       ! longitude (radians)
   real(r8), pointer :: latdeg(:)       ! latitude (degrees)
   real(r8), pointer :: londeg(:)       ! longitude (degrees)
   real(r8), pointer :: lat_a(:)        ! "atm" latitude (radians)
   real(r8), pointer :: lon_a(:)        ! "atm" longitude (radians)
   real(r8), pointer :: latdeg_a(:)     ! "atm" latitude (degrees)
   real(r8), pointer :: londeg_a(:)     ! "atm" longitude (degrees)
   real(r8), pointer :: coszen(:)       ! cosine of zenith angle
   real(r8), pointer :: decl(:)       ! declination angle (radians)
   real(r8), pointer :: albgrd(:,:)     ! direct surface reflectance
   real(r8), pointer :: albgri(:,:)     ! diffuse surface reflectance
   real(r8), pointer :: rb(:)           ! canopy resistance
   real(r8), pointer :: rssun(:)        ! sunlit stomatal resistance
   real(r8), pointer :: rssha(:)        ! shaded stomatal resistance
   real(r8), pointer :: cisun(:)        ! sunlit intracellular CO2
   real(r8), pointer :: cisha(:)        ! shaded intracellular CO2
   real(r8), pointer :: parsun(:)       ! sunlit absorbed PAR
   real(r8), pointer :: parsha(:)       ! shaded absorbed PAR
   real(r8), pointer :: laisun(:)       ! sunlit LAI
   real(r8), pointer :: laisha(:)       ! shaded LAI
   real(r8), pointer :: displa(:)       ! displacement height
   real(r8), pointer :: z0mv(:)         ! roughness length over vegetation
   real(r8), pointer :: z0hv(:)
   real(r8), pointer :: z0qv(:)
   real(r8), pointer :: z0mg(:)         ! roughness length over ground
   real(r8), pointer :: z0hg(:)
   real(r8), pointer :: z0qg(:)
   integer, pointer :: frac_veg_nosno(:)

   integer g, c, p, k, r, doy
   integer, parameter :: nlevsno = 2
   real(r8) :: vg                       ! surface vapor pressure (Pa)
   real(r8) :: eccen                    ! orbital eccentricity
   real(r8) :: obliq                    ! obliquity (degrees)
   real(r8) :: mvelp                    ! moving vernal equinox longitude
   real(r8) :: obliqr                   ! obliquity (radians)
   real(r8) :: lambm0                   ! mean lon of perihelion at vernal equinox (radians)
   real(r8) :: mvelpp                   ! moving vernal equinox lon + pi (rads)
   real(r8) :: delta                    ! declination angle (radians)
   real(r8) :: eccf                     ! Earth-sun distance factor
   real(r8) :: cosz                     ! cosine of zenith angle
   real(r8) :: eccen1                    ! orbital eccentricity
   real(r8) :: obliq1                    ! obliquity (degrees)
   real(r8) :: mvelp1                    ! moving vernal equinox longitude
   real(r8) :: obliqr1                   ! obliquity (radians)
   real(r8) :: lambm01                   ! mean lon of perihelion at vernal equinox (radians)
   real(r8) :: mvelpp1                   ! moving vernal equinox lon + pi (rads)
   real(r8) :: delta1                    ! declination angle (radians)
   real(r8) :: eccf1                     ! Earth-sun distance factor

   call c_f_pointer(vic_cn, cn, [max_bands])

   ! Assign local pointers to derived type arrays
   ncol       => clm3%g%c%p%column
   ivt        => clm3%g%c%p%itype
   z          => clm3%g%c%cps%z
   dz         => clm3%g%c%cps%dz
   frac_veg_nosno => clm3%g%c%p%pps%frac_veg_nosno
   qflx_drain => clm3%g%c%cwf%qflx_drain
   h2osoi_liq => clm3%g%c%cws%h2osoi_liq
   h2osoi_ice => clm3%g%c%cws%h2osoi_ice
   t_soisno   => clm3%g%c%ces%t_soisno
   t_ref2m    => clm3%g%c%p%pes%t_ref2m
   t_veg      => clm3%g%c%p%pes%t_veg
   fwet       => clm3%g%c%p%pps%fwet
   snowdp     => clm3%g%c%cps%snowdp
   rootfr     => clm3%g%c%p%pps%rootfr
   psisat     => clm3%g%c%cps%psisat
   soilpsi    => clm3%g%c%cps%soilpsi
   bsw2        => clm3%g%c%cps%bsw2
   sucsat     => clm3%g%c%cps%sucsat
   bsw        => clm3%g%c%cps%bsw
   watfc      => clm3%g%c%cps%watfc
   forc_t     => clm_a2l%forc_t
   forc_q     => clm_a2l%forc_q
   forc_pbot  => clm_a2l%forc_pbot
   forc_vp    => clm_a2l%forc_vp
   qg         => clm3%g%c%cws%qg
   thm        => clm3%g%c%p%pes%thm
   forc_hgt_u => clm_a2l%forc_hgt_u
   forc_hgt_t => clm_a2l%forc_hgt_t
   forc_hgt_q => clm_a2l%forc_hgt_q
   forc_hgt_u_pft => clm3%g%c%p%pps%forc_hgt_u_pft
   forc_hgt_t_pft => clm3%g%c%p%pps%forc_hgt_t_pft
   forc_hgt_q_pft => clm3%g%c%p%pps%forc_hgt_q_pft
   rb         => clm3%g%c%p%pps%rb1
   tlai       => clm3%g%c%p%pps%tlai
   elai       => clm3%g%c%p%pps%elai
   displa     => clm3%g%c%p%pps%displa
   z0mv       => clm3%g%c%p%pps%z0mv
   z0hv       => clm3%g%c%p%pps%z0hv
   z0qv       => clm3%g%c%p%pps%z0qv
   z0mg       => clm3%g%c%cps%z0mg
   z0hg       => clm3%g%c%cps%z0hg
   z0qg       => clm3%g%c%cps%z0qg
   forc_pco2  => clm_a2l%forc_pco2
   forc_po2   => clm_a2l%forc_po2
   forc_lwrad => clm_a2l%forc_lwrad
   forc_solar => clm_a2l%forc_solar
   forc_solad => clm_a2l%forc_solad
   forc_solai => clm_a2l%forc_solai
   latrad     => clm3%g%lat
   lonrad     => clm3%g%lon
   latdeg     => clm3%g%latdeg
   londeg     => clm3%g%londeg
   lat_a      => clm3%g%lat_a
   lon_a      => clm3%g%lon_a
   latdeg_a   => clm3%g%latdeg_a
   londeg_a   => clm3%g%londeg_a
   coszen     => clm3%g%c%cps%coszen
   albgrd     => clm3%g%c%cps%albgrd
   albgri     => clm3%g%c%cps%albgri
   decl       => clm3%g%c%cps%decl
   fpsn       => clm3%g%c%p%pcf%fpsn
   psnsun     => clm3%g%c%p%pcf%psnsun
   psnsha     => clm3%g%c%p%pcf%psnsha

   ! Assign VIC data to local pointers

   print *, rec, yr, mo, day, secs
   do g = begg, endg

!	print *, g

   ! Added by MAB, 8/27/13

     latdeg(g) = lat
     latdeg_a(g) = lat
     if(lon .lt. 0.) then
       londeg(g) = 360._r8 + lon
     else
       londeg(g) = lon
     endif
     londeg_a(g) = londeg(g)

     latrad(g) = latdeg(g) * SHR_CONST_PI / 180.
     lonrad(g) = londeg(g) * SHR_CONST_PI / 180.
     lat_a(g) = latrad(g)
     lon_a(g) = lonrad(g)

     call shr_orb_params(yr, eccen, obliq, mvelp, obliqr, lambm0, mvelpp)
     call shr_orb_decl(jday, eccen, mvelpp, lambm0, obliqr, delta, eccf)
     cosz = shr_orb_cosz(jday, latrad(g), lonrad(g), delta)
     do c = begc, endc
       coszen((g - 1) * endc + c) = cosz
       decl((g - 1) * endc + c) = delta
     end do

   ! Added by MAB, 10/8/13

     if(yr1 .gt. 0 .and. jday1 .gt. 0.) then
       call shr_orb_params(yr1, eccen1, obliq1, mvelp1, obliqr1, lambm01, &
	mvelpp1)
       call shr_orb_decl(jday1, eccen1, mvelpp1, lambm01, obliqr1, delta1, &
	eccf1)
     end if

   ! Added by MAB, 8/29/13

     forc_t(g) = cn(1)%tair + 273.16
     forc_q(g) = 0.622_r8 * cn(1)%vp / cn(1)%psfc
     forc_vp(g) = cn(1)%vp
     vg = cn(1)%vp + cn(1)%vpd
     do c = begc, endc
!	print *, c
       qg((g - 1) * endc + c) = 0.622_r8 * vg / cn(1)%psfc
     end do
	print *, cn(1)%tair, cn(1)%vp, cn(1)%vpd, cn(1)%psfc

   ! Added by MAB, 10/11/13

     forc_hgt_u(g) = 10._r8
     forc_hgt_t(g) = 2._r8
     forc_hgt_q(g) = 2._r8

   ! Added by MAB, 8/14/13

     forc_pbot(g) = cn(1)%psfc
     forc_pco2(g) = cn(1)%psfc * 3.8e-4
     forc_po2(g) = cn(1)%psfc * 0.2095
     forc_lwrad(g) = cn(1)%lwrad
     forc_solar(g) = cn(1)%swrad
     do r = 1, numrad
       forc_solad(g,r) = cn(1)%swrd(r)
       forc_solai(g,r) = cn(1)%swri(r)
     end do

   do c = begc, endc

!	print *, c
     qflx_drain((g - 1) * endc + c) = cn(c)%baseflow
     snowdp((g - 1) * endc + c) = cn(c)%snowdep

     ! Added by MAB, 10/11/13
     if(snowdp((g - 1) * endc + c) .gt. 0.) then
       z0mg((g - 1) * endc + c) = cn(c)%z0s
     else
       z0mg((g - 1) * endc + c) = cn(c)%z0
     end if
     z0hg((g - 1) * endc + c) = z0mv(c)
     z0qg((g - 1) * endc + c) = z0mv(c)

     do r = 1, numrad
       albgrd((g - 1) * endc + c,r) = cn(1)%alb
       albgri((g - 1) * endc + c,r) = cn(1)%alb
     end do
!	print *, cn(1)%baseflow, cn(1)%snowdep, cn(1)%z0s, cn(1)%z0, cn(1)%alb

     do k = -nlevsno + 1, nlevgrnd  ! -nlevsno+1:nlevgrnd
         z((g - 1) * endc + c,k) = cn(c)%z(k + nlevsno)
         dz((g - 1) * endc + c,k) = cn(c)%dz(k + nlevsno)
         h2osoi_liq((g - 1) * endc + c,k) = cn(c)%moist(k + nlevsno)
         h2osoi_ice((g - 1) * endc + c,k) = cn(c)%ice(k + nlevsno)
         t_soisno((g - 1) * endc + c,k) = cn(c)%t_soisno(k + nlevsno) + &
		SHR_CONST_TKFRZ
     end do

     do k = 1, nlevgrnd
         bsw((g - 1) * endc + c,k) = cn(c)%bsw(1,k)
         sucsat((g - 1) * endc + c,k) = cn(c)%sucsat(1,k)
         bsw2((g - 1) * endc + c,k) = cn(c)%bsw(2,k)
         psisat((g - 1) * endc + c,k) = cn(c)%sucsat(2,k)
         soilpsi((g - 1) * endc + c,k) = cn(c)%soisuc(2,k)
         watfc((g - 1) * endc + c,k) = cn(c)%watfc(k)
     end do

     ! How to find the right p's for g != 1?  (g - 1) * endc * endp?

     do p = begp, endp
         t_ref2m((g - 1) * endc * endp + (c - 1) * endp + p) = cn(1)%t2m(ivt(p)) + 273.16
         do k = 1, nlevgrnd
           rootfr((g - 1) * endc * endp + (c - 1) * endp + p,k) = cn(c)%rootfr(k,ivt(p))
!	   if(ivt(p) .lt. 10) then
!		print *, p, ivt(p), k, rootfr((g - 1) * endc * endp + (c - 1) * endp + p, k), cn(c)%rootfr(k,ivt(p))
!	   end if
         end do
         t_veg((g - 1) * endc * endp + (c - 1) * endp + p) = cn(c)%Tveg(ivt(p)) + &
		273.16
         fwet((g - 1) * endc * endp + (c - 1) * endp + p) = cn(c)%fwet(ivt(p))
         rb((g - 1) * endc * endp + (c - 1) * endp + p) = cn(c)%rveg(ivt(p))

         ! Divide up photosynthesis, all photosynthesis goes to sunlit leaves
         ! for now, added by MAB, 12/10/14
!         fpsn((g - 1) * endc * endp + (c - 1) * endp + p) = &
!		photosynth(ivt(p),c)
!         psnsun((g - 1) * endc * endp + (c - 1) * endp + p) = &
!		photosynth(ivt(p),c)
!         psnsha((g - 1) * endc * endp + (c - 1) * endp + p) = 0._r8
!         print *, c, p, fpsn((g - 1) * endc * endp + (c - 1) * endp + p), &
!		psnsun((g - 1) * endc * endp + (c - 1) * endp + p), &
!		psnsha((g - 1) * endc * endp + (c - 1) * endp + p)

         ! Added by MAB, 10/11/13
!         if(rec == 0 .and. nspinup == 0) then
         if(rec == 0) then
           tlai((g - 1) * endc * endp + (c - 1) * endp + p) = cn(c)%LAI(ivt(p))
         end if
         displa((g - 1) * endc * endp + (c - 1) * endp + p) = cn(c)%displ(ivt(p))
         z0mv((g - 1) * endc * endp + (c - 1) * endp + p) = cn(c)%zov(ivt(p))
         z0hv((g - 1) * endc * endp + (c - 1) * endp + p) = cn(c)%zov(ivt(p))
         z0qv((g - 1) * endc * endp + (c - 1) * endp + p) = cn(c)%zov(ivt(p))
!	print *, p, ivt(p), cn(c)%LAI(ivt(p)), cn(c)%displ(ivt(p)), cn(c)%zov(ivt(p))
         if(frac_veg_nosno((g - 1) * endc * endp + (c - 1) * endp + p) == 0) &
		then
           forc_hgt_u_pft((g - 1) * endc * endp + (c - 1) * endp + p) = &
		forc_hgt_u(1) + z0mg(c) + &
		displa((g - 1) * endc * endp + (c - 1) * endp + p)
           forc_hgt_t_pft((g - 1) * endc * endp + (c - 1) * endp + p) = &
		forc_hgt_t(1) + z0mg(c) + &
		displa((g - 1) * endc * endp + (c - 1) * endp + p)
           forc_hgt_q_pft((g - 1) * endc * endp + (c - 1) * endp + p) = &
		forc_hgt_q(1) + z0mg(c) + &
		displa((g - 1) * endc * endp + (c - 1) * endp + p)
         else
           forc_hgt_u_pft((g - 1) * endc * endp + (c - 1) * endp + p) = &
		forc_hgt_u(1) + &
		z0mv((g - 1) * endc * endp + (c - 1) * endp + p) + &
		displa((g - 1) * endc * endp + (c - 1) * endp + p)
           forc_hgt_t_pft((g - 1) * endc * endp + (c - 1) * endp + p) = &
		forc_hgt_t(1) + &
		z0mv((g - 1) * endc * endp + (c - 1) * endp + p) + &
		displa((g - 1) * endc * endp + (c - 1) * endp + p)
           forc_hgt_q_pft((g - 1) * endc * endp + (c - 1) * endp + p) = &
		forc_hgt_q(1) + &
		z0mv((g - 1) * endc * endp + (c - 1) * endp + p) + &
		displa((g - 1) * endc * endp + (c - 1) * endp + p)
         end if
         thm((g - 1) * endc * endp + (c - 1) * endp + p) = forc_t(1) + &
		0.0098_r8 * &
		forc_hgt_t_pft((g - 1) * endc * endp + (c - 1) * endp + p)
       end do
     end do
   end do

   begc = 1
   endc = num_soilc
   begp = 1
   endp = num_soilc * num_soilp

   soil1c => clm3%g%c%ccs%soil1c
   soil2c => clm3%g%c%ccs%soil2c
   soil3c => clm3%g%c%ccs%soil3c
   soil4c => clm3%g%c%ccs%soil4c
   litr1c => clm3%g%c%ccs%litr1c
   litr2c => clm3%g%c%ccs%litr2c
   litr3c => clm3%g%c%ccs%litr3c
   cwdc => clm3%g%c%ccs%cwdc
   leafc => clm3%g%c%p%pcs%leafc
   frootc => clm3%g%c%p%pcs%frootc
   livestemc => clm3%g%c%p%pcs%livestemc
   deadstemc => clm3%g%c%p%pcs%deadstemc
   livecrootc => clm3%g%c%p%pcs%livecrootc
   deadcrootc => clm3%g%c%p%pcs%deadcrootc
   woodc => clm3%g%c%p%pcs%woodc
   totvegc => clm3%g%c%p%pcs%totvegc
   totlitc => clm3%g%c%ccs%totlitc
   totsomc => clm3%g%c%ccs%totsomc
   soil1n => clm3%g%c%cns%soil1n
   soil2n => clm3%g%c%cns%soil2n
   soil3n => clm3%g%c%cns%soil3n
   soil4n => clm3%g%c%cns%soil4n
   sminn => clm3%g%c%cns%sminn
   litr1n => clm3%g%c%cns%litr1n
   litr2n => clm3%g%c%cns%litr2n
   litr3n => clm3%g%c%cns%litr3n
   cwdn => clm3%g%c%cns%cwdn
   leafn => clm3%g%c%p%pns%leafn
   frootn => clm3%g%c%p%pns%frootn
   livestemn => clm3%g%c%p%pns%livestemn
   deadstemn => clm3%g%c%p%pns%deadstemn
   livecrootn => clm3%g%c%p%pns%livecrootn
   deadcrootn => clm3%g%c%p%pns%deadcrootn
   cpool_to_leafc => clm3%g%c%p%pcf%cpool_to_leafc
   gpp2 => clm3%g%c%p%pepv%gpp
   gpp => clm3%g%c%p%pcf%gpp
   npp => clm3%g%c%p%pcf%npp
   ar => clm3%g%c%p%pcf%ar
   gr => clm3%g%c%p%pcf%gr
   mr => clm3%g%c%p%pcf%mr
   leaf_mr => clm3%g%c%p%pcf%leaf_mr
   hr => clm3%g%c%ccf%hr
   nee => clm3%g%c%ccf%nee
   nep => clm3%g%c%ccf%nep
   tlai       => clm3%g%c%p%pps%tlai
   elai       => clm3%g%c%p%pps%elai
   dormant_flag => clm3%g%c%p%pepv%dormant_flag
   days_active => clm3%g%c%p%pepv%days_active
   onset_flag => clm3%g%c%p%pepv%onset_flag
   onset_counter => clm3%g%c%p%pepv%onset_counter
   onset_gddflag => clm3%g%c%p%pepv%onset_gddflag
   onset_fdd => clm3%g%c%p%pepv%onset_fdd
   onset_gdd => clm3%g%c%p%pepv%onset_gdd
   onset_swi => clm3%g%c%p%pepv%onset_swi
   offset_flag => clm3%g%c%p%pepv%offset_flag
   offset_counter => clm3%g%c%p%pepv%offset_counter
   offset_fdd => clm3%g%c%p%pepv%offset_fdd
   offset_swi => clm3%g%c%p%pepv%offset_swi
   lgsf => clm3%g%c%p%pepv%lgsf
   bglfr => clm3%g%c%p%pepv%bglfr
   bgtr => clm3%g%c%p%pepv%bgtr
   dayl => clm3%g%c%p%pepv%dayl
   prev_dayl => clm3%g%c%p%pepv%prev_dayl
   annavg_t2m => clm3%g%c%p%pepv%annavg_t2m
   tempavg_t2m => clm3%g%c%p%pepv%tempavg_t2m
   availc => clm3%g%c%p%pepv%availc
   xsmrpool_recover => clm3%g%c%p%pepv%xsmrpool_recover
   alloc_pnow => clm3%g%c%p%pepv%alloc_pnow
   c_allometry => clm3%g%c%p%pepv%c_allometry
   n_allometry => clm3%g%c%p%pepv%n_allometry
   plant_ndemand => clm3%g%c%p%pepv%plant_ndemand
   tempsum_potential_gpp => clm3%g%c%p%pepv%tempsum_potential_gpp
   annsum_potential_gpp => clm3%g%c%p%pepv%annsum_potential_gpp
   tempmax_retransn => clm3%g%c%p%pepv%tempmax_retransn
   annmax_retransn => clm3%g%c%p%pepv%annmax_retransn
   avail_retransn => clm3%g%c%p%pepv%avail_retransn
   plant_nalloc => clm3%g%c%p%pepv%plant_nalloc
   plant_calloc => clm3%g%c%p%pepv%plant_calloc
   excess_cflux => clm3%g%c%p%pepv%excess_cflux
   downreg => clm3%g%c%p%pepv%downreg
   prev_leafc_to_litter => clm3%g%c%p%pepv%prev_leafc_to_litter
   prev_frootc_to_litter => clm3%g%c%p%pepv%prev_frootc_to_litter
   tempsum_npp => clm3%g%c%p%pepv%tempsum_npp
   annsum_npp => clm3%g%c%p%pepv%annsum_npp
   leafc_storage => clm3%g%c%p%pcs%leafc_storage
   leafc_xfer => clm3%g%c%p%pcs%leafc_xfer
   frootc_storage => clm3%g%c%p%pcs%frootc_storage
   frootc_xfer => clm3%g%c%p%pcs%frootc_xfer
   livestemc_storage => clm3%g%c%p%pcs%livestemc_storage
   livestemc_xfer => clm3%g%c%p%pcs%livestemc_xfer
   deadstemc_storage => clm3%g%c%p%pcs%deadstemc_storage
   deadstemc_xfer => clm3%g%c%p%pcs%deadstemc_xfer
   livecrootc_storage => clm3%g%c%p%pcs%livecrootc_storage
   livecrootc_xfer => clm3%g%c%p%pcs%livecrootc_xfer
   deadcrootc_storage => clm3%g%c%p%pcs%deadcrootc_storage
   deadcrootc_xfer => clm3%g%c%p%pcs%deadcrootc_xfer
   gresp_storage => clm3%g%c%p%pcs%gresp_storage
   gresp_xfer => clm3%g%c%p%pcs%gresp_xfer
   cpool => clm3%g%c%p%pcs%cpool
   xsmrpool => clm3%g%c%p%pcs%xsmrpool
   pft_ctrunc => clm3%g%c%p%pcs%pft_ctrunc
   leafn_storage => clm3%g%c%p%pns%leafn_storage
   leafn_xfer => clm3%g%c%p%pns%leafn_xfer
   frootn_storage => clm3%g%c%p%pns%frootn_storage
   frootn_xfer => clm3%g%c%p%pns%frootn_xfer
   livestemn_storage => clm3%g%c%p%pns%livestemn_storage
   livestemn_xfer => clm3%g%c%p%pns%livestemn_xfer
   deadstemn_storage => clm3%g%c%p%pns%deadstemn_storage
   deadstemn_xfer => clm3%g%c%p%pns%deadstemn_xfer
   livecrootn_storage => clm3%g%c%p%pns%livecrootn_storage
   livecrootn_xfer => clm3%g%c%p%pns%livecrootn_xfer
   deadcrootn_storage => clm3%g%c%p%pns%deadcrootn_storage
   deadcrootn_xfer => clm3%g%c%p%pns%deadcrootn_xfer
   retransn => clm3%g%c%p%pns%retransn
   npool => clm3%g%c%p%pns%npool
   pft_ntrunc => clm3%g%c%p%pns%pft_ntrunc
!   decl => clm3%g%c%cps%decl
   fpi => clm3%g%c%cps%fpi
   fpg => clm3%g%c%cps%fpg
   annsum_counter => clm3%g%c%cps%annsum_counter
   cannsum_npp => clm3%g%c%cps%cannsum_npp
   cannsum_t2m => clm3%g%c%cps%cannavg_t2m
   me => clm3%g%c%cps%me
   fire_prob => clm3%g%c%cps%fire_prob
   mean_fire_prob => clm3%g%c%cps%mean_fire_prob
   fireseasonl => clm3%g%c%cps%fireseasonl
   farea_burned => clm3%g%c%cps%farea_burned
   ann_farea_burned => clm3%g%c%cps%ann_farea_burned
   seedc => clm3%g%c%ccs%seedc
   col_ctrunc => clm3%g%c%ccs%col_ctrunc
   totcolc => clm3%g%c%ccs%totcolc
   prod10c => clm3%g%c%ccs%prod10c
   prod100c => clm3%g%c%ccs%prod100c
   seedn => clm3%g%c%cns%seedn
   col_ntrunc => clm3%g%c%cns%col_ntrunc
   totcoln => clm3%g%c%cns%totcoln
   prod10n => clm3%g%c%cns%prod10n
   prod100n => clm3%g%c%cns%prod100n

   if(rec .eq. 0 .and. init_state .eq. 1) then

    do c = begc, endc
     soil1c(c) = cn(c)%soil1c
     soil2c(c) = cn(c)%soil2c
     soil3c(c) = cn(c)%soil3c
     soil4c(c) = cn(c)%soil4c
     litr1c(c) = cn(c)%litr1c
     litr2c(c) = cn(c)%litr2c
     litr3c(c) = cn(c)%litr3c
     cwdc(c) = cn(c)%cwdc
     totlitc(c) = cn(c)%totlitc
     totsomc(c) = cn(c)%totsomc
     soil1n(c) = cn(c)%soil1n
     soil2n(c) = cn(c)%soil2n
     soil3n(c) = cn(c)%soil3n
     soil4n(c) = cn(c)%soil4n
     sminn(c) = cn(c)%sminn
     litr1n(c) = cn(c)%litr1n
     litr2n(c) = cn(c)%litr2n
     litr3n(c) = cn(c)%litr3n
     cwdn(c) = cn(c)%cwdn
     hr(c) = cn(c)%hr
     nee(c) = cn(c)%nee
     nep(c) = cn(c)%nep
!     decl(c) = cn(c)%decl
     fpi(c) = cn(c)%fpi
     fpg(c) = cn(c)%fpg
     annsum_counter(c) = cn(c)%annsum_counter
     cannsum_npp(c) = cn(c)%cannsum_npp
     cannsum_t2m(c) = cn(c)%cannavg_t2m
     do k = 1, nlevgrnd
       watfc(c,k) = cn(c)%watfc(k)
     end do
     me(c) = cn(c)%me
     fire_prob(c) = cn(c)%fire_prob
     mean_fire_prob(c) = cn(c)%mean_fire_prob
     fireseasonl(c) = cn(c)%fireseasonl
     farea_burned(c) = cn(c)%farea_burned
     ann_farea_burned(c) = cn(c)%ann_farea_burned
     seedc(c) = cn(c)%seedc
     col_ctrunc(c) = cn(c)%col_ctrunc
     totcolc(c) = cn(c)%totcolc
     prod10c(c) = cn(c)%prod10c
     prod100c(c) = cn(c)%prod100c
     seedn(c) = cn(c)%seedn
     col_ntrunc(c) = cn(c)%col_ntrunc
     totcoln(c) = cn(c)%totcoln
     prod10n(c) = cn(c)%prod10n
     prod100n(c) = cn(c)%prod100n
   end do

   do p = begp, endp
     c = ncol(p)
       elai(p) = cn(c)%LAI(ivt(p))
       leafc(p) = cn(c)%leafc(ivt(p))
       frootc(p) = cn(c)%frootc(ivt(p))
       livestemc(p) = cn(c)%livestemc(ivt(p))
       deadstemc(p) = cn(c)%deadstemc(ivt(p))
       livecrootc(p) = cn(c)%livecrootc(ivt(p))
       deadcrootc(p) = cn(c)%deadcrootc(ivt(p))
       woodc(p) = cn(c)%woodc(ivt(p))
       totvegc(p) = cn(c)%totvegc(ivt(p))
       leafn(p) = cn(c)%leafn(ivt(p))
       frootn(p) = cn(c)%frootn(ivt(p))
       livestemn(p) = cn(c)%livestemn(ivt(p))
       deadstemn(p) = cn(c)%deadstemn(ivt(p))
       livecrootn(p) = cn(c)%livecrootn(ivt(p))
       deadcrootn(p) = cn(c)%deadcrootn(ivt(p))
       gpp2(p) = cn(c)%gpp2(ivt(p))
       gpp(p) = cn(c)%gpp(ivt(p))
       npp(p) = cn(c)%npp(ivt(p))
       ar(p) = cn(c)%ar(ivt(p))
       gr(p) = cn(c)%gr(ivt(p))
       mr(p) = cn(c)%mr(ivt(p))
       leaf_mr(p) = cn(c)%leaf_mr(ivt(p))
       dormant_flag(p) = cn(c)%dormant_flag(ivt(p))
       days_active(p) = cn(c)%days_active(ivt(p))
       onset_flag(p) = cn(c)%onset_flag(ivt(p))
       onset_counter(p) = cn(c)%onset_counter(ivt(p))
       onset_gddflag(p) = cn(c)%onset_gddflag(ivt(p))
       onset_fdd(p) = cn(c)%onset_fdd(ivt(p))
       onset_gdd(p) = cn(c)%onset_gdd(ivt(p))
       onset_swi(p) = cn(c)%onset_swi(ivt(p))
       offset_flag(p) = cn(c)%offset_flag(ivt(p))
       offset_counter(p) = cn(c)%offset_counter(ivt(p))
       offset_fdd(p) = cn(c)%offset_fdd(ivt(p))
       offset_swi(p) = cn(c)%offset_swi(ivt(p))
       lgsf(p) = cn(c)%lgsf(ivt(p))
       bglfr(p) = cn(c)%bglfr(ivt(p))
       bgtr(p) = cn(c)%bgtr(ivt(p))
       dayl(p) = cn(c)%dayl(ivt(p))
       prev_dayl(p) = cn(c)%prev_dayl(ivt(p))
       annavg_t2m(p) = cn(c)%annavg_t2m(ivt(p))
       tempavg_t2m(p) = cn(c)%tempavg_t2m(ivt(p))
       availc(p) = cn(c)%availc(ivt(p))
       xsmrpool_recover(p) = cn(c)%xsmrpool_recover(ivt(p))
       alloc_pnow(p) = cn(c)%alloc_pnow(ivt(p))
       c_allometry(p) = cn(c)%c_allometry(ivt(p))
       n_allometry(p) = cn(c)%n_allometry(ivt(p))
       plant_ndemand(p) = cn(c)%plant_ndemand(ivt(p))
       tempsum_potential_gpp(p) = cn(c)%tempsum_potential_gpp(ivt(p))
       annsum_potential_gpp(p) = cn(c)%annsum_potential_gpp(ivt(p))
       tempmax_retransn(p) = cn(c)%tempmax_retransn(ivt(p))
       annmax_retransn(p) = cn(c)%annmax_retransn(ivt(p))
       avail_retransn(p) = cn(c)%avail_retransn(ivt(p))
       plant_nalloc(p) = cn(c)%plant_nalloc(ivt(p))
       plant_calloc(p) = cn(c)%plant_calloc(ivt(p))
       excess_cflux(p) = cn(c)%excess_cflux(ivt(p))
       downreg(p) = cn(c)%downreg(ivt(p))
       prev_leafc_to_litter(p) = cn(c)%prev_leafc_to_litter(ivt(p))
       prev_frootc_to_litter(p) = cn(c)%prev_frootc_to_litter(ivt(p))
       tempsum_npp(p) = cn(c)%tempsum_npp(ivt(p))
       annsum_npp(p) = cn(c)%annsum_npp(ivt(p))
       leafc_storage(p) = cn(c)%leafc_storage(ivt(p))
       leafc_xfer(p) = cn(c)%leafc_xfer(ivt(p))
       frootc_storage(p) = cn(c)%frootc_storage(ivt(p))
       frootc_xfer(p) = cn(c)%frootc_xfer(ivt(p))
       livestemc_storage(p) = cn(c)%livestemc_storage(ivt(p))
       livestemc_xfer(p) = cn(c)%livestemc_xfer(ivt(p))
       deadstemc_storage(p) = cn(c)%deadstemc_storage(ivt(p))
       deadstemc_xfer(p) = cn(c)%deadstemc_xfer(ivt(p))
       livecrootc_storage(p) = cn(c)%livecrootc_storage(ivt(p))
       livecrootc_xfer(p) = cn(c)%livecrootc_xfer(ivt(p))
       deadcrootc_storage(p) = cn(c)%deadcrootc_storage(ivt(p))
       deadcrootc_xfer(p) = cn(c)%deadcrootc_xfer(ivt(p))
       gresp_storage(p) = cn(c)%gresp_storage(ivt(p))
       gresp_xfer(p) = cn(c)%gresp_xfer(ivt(p))
       cpool(p) = cn(c)%cpool(ivt(p))
       xsmrpool(p) = cn(c)%xsmrpool(ivt(p))
       pft_ctrunc(p) = cn(c)%pft_ctrunc(ivt(p))
       leafn_storage(p) = cn(c)%leafn_storage(ivt(p))
       leafn_xfer(p) = cn(c)%leafn_xfer(ivt(p))
       frootn_storage(p) = cn(c)%frootn_storage(ivt(p))
       frootn_xfer(p) = cn(c)%frootn_xfer(ivt(p))
       livestemn_storage(p) = cn(c)%livestemn_storage(ivt(p))
       livestemn_xfer(p) = cn(c)%livestemn_xfer(ivt(p))
       deadstemn_storage(p) = cn(c)%deadstemn_storage(ivt(p))
       deadstemn_xfer(p) = cn(c)%deadstemn_xfer(ivt(p))
       livecrootn_storage(p) = cn(c)%livecrootn_storage(ivt(p))
       livecrootn_xfer(p) = cn(c)%livecrootn_xfer(ivt(p))
       deadcrootn_storage(p) = cn(c)%deadcrootn_storage(ivt(p))
       deadcrootn_xfer(p) = cn(c)%deadcrootn_xfer(ivt(p))
       retransn(p) = cn(c)%retransn(ivt(p))
       npool(p) = cn(c)%npool(ivt(p))
       pft_ntrunc(p) = cn(c)%pft_ntrunc(ivt(p))
!	print *, init_state, cn(c)%leafc(ivt(p))
    end do
   end if

   if(jday1 .gt. 0.) then
     call SurfaceAlbedo(begg, endg, begc, endc, begp, endp, jday1, delta1)
   end if
   call SurfaceRadiation(yr, mo, day, secs, dt, begp, endp)
   call CanopyFluxes(dt, nlevgrnd, begg, endg, begc, endc, begp, endp)

   doy = int(jday)
   call CNEcosystemDyn(adspinup, nlevgrnd, rec, yr, doy, dt, begg, endg, &
	begc, endc, begp, endp, num_soilc, num_soilp)
   call CNAnnualUpdate(dt, yr, begc, endc, begp, endp, num_soilc, num_soilp)

   ! Check the carbon and nitrogen balance
   call CBalanceCheck(dt, begc, endc, num_soilc)
   call NBalanceCheck(dt, begc, endc, num_soilc)

   soil1c => clm3%g%c%ccs%soil1c
   soil2c => clm3%g%c%ccs%soil2c
   soil3c => clm3%g%c%ccs%soil3c
   soil4c => clm3%g%c%ccs%soil4c
   litr1c => clm3%g%c%ccs%litr1c
   litr2c => clm3%g%c%ccs%litr2c
   litr3c => clm3%g%c%ccs%litr3c
   cwdc => clm3%g%c%ccs%cwdc
   leafc => clm3%g%c%p%pcs%leafc
   frootc => clm3%g%c%p%pcs%frootc
   livestemc => clm3%g%c%p%pcs%livestemc
   deadstemc => clm3%g%c%p%pcs%deadstemc
   livecrootc => clm3%g%c%p%pcs%livecrootc
   deadcrootc => clm3%g%c%p%pcs%deadcrootc
   woodc => clm3%g%c%p%pcs%woodc
   totvegc => clm3%g%c%p%pcs%totvegc
   totlitc => clm3%g%c%ccs%totlitc
   totsomc => clm3%g%c%ccs%totsomc
   soil1n => clm3%g%c%cns%soil1n
   soil2n => clm3%g%c%cns%soil2n
   soil3n => clm3%g%c%cns%soil3n
   soil4n => clm3%g%c%cns%soil4n
   sminn => clm3%g%c%cns%sminn
   litr1n => clm3%g%c%cns%litr1n
   litr2n => clm3%g%c%cns%litr2n
   litr3n => clm3%g%c%cns%litr3n
   cwdn => clm3%g%c%cns%cwdn
   leafn => clm3%g%c%p%pns%leafn
   frootn => clm3%g%c%p%pns%frootn
   livestemn => clm3%g%c%p%pns%livestemn
   deadstemn => clm3%g%c%p%pns%deadstemn
   livecrootn => clm3%g%c%p%pns%livecrootn
   deadcrootn => clm3%g%c%p%pns%deadcrootn
   cpool_to_leafc => clm3%g%c%p%pcf%cpool_to_leafc
   gpp2 => clm3%g%c%p%pepv%gpp
   gpp => clm3%g%c%p%pcf%gpp
   npp => clm3%g%c%p%pcf%npp
   mr => clm3%g%c%p%pcf%mr
   leaf_mr => clm3%g%c%p%pcf%leaf_mr
   gr => clm3%g%c%p%pcf%gr
   ar => clm3%g%c%p%pcf%ar
   hr => clm3%g%c%ccf%hr
   lithr => clm3%g%c%ccf%lithr
   nee => clm3%g%c%ccf%nee
   nep => clm3%g%c%ccf%nep
   tlai       => clm3%g%c%p%pps%tlai
   elai       => clm3%g%c%p%pps%elai
   dormant_flag => clm3%g%c%p%pepv%dormant_flag
   days_active => clm3%g%c%p%pepv%days_active
   onset_flag => clm3%g%c%p%pepv%onset_flag
   onset_counter => clm3%g%c%p%pepv%onset_counter
   onset_gddflag => clm3%g%c%p%pepv%onset_gddflag
   onset_fdd => clm3%g%c%p%pepv%onset_fdd
   onset_gdd => clm3%g%c%p%pepv%onset_gdd
   onset_swi => clm3%g%c%p%pepv%onset_swi
   offset_flag => clm3%g%c%p%pepv%offset_flag
   offset_counter => clm3%g%c%p%pepv%offset_counter
   offset_fdd => clm3%g%c%p%pepv%offset_fdd
   offset_swi => clm3%g%c%p%pepv%offset_swi
   lgsf => clm3%g%c%p%pepv%lgsf
   bglfr => clm3%g%c%p%pepv%bglfr
   bgtr => clm3%g%c%p%pepv%bgtr
   dayl => clm3%g%c%p%pepv%dayl
   prev_dayl => clm3%g%c%p%pepv%prev_dayl
   annavg_t2m => clm3%g%c%p%pepv%annavg_t2m
   tempavg_t2m => clm3%g%c%p%pepv%tempavg_t2m
   availc => clm3%g%c%p%pepv%availc
   xsmrpool_recover => clm3%g%c%p%pepv%xsmrpool_recover
   alloc_pnow => clm3%g%c%p%pepv%alloc_pnow
   c_allometry => clm3%g%c%p%pepv%c_allometry
   n_allometry => clm3%g%c%p%pepv%n_allometry
   plant_ndemand => clm3%g%c%p%pepv%plant_ndemand
   tempsum_potential_gpp => clm3%g%c%p%pepv%tempsum_potential_gpp
   annsum_potential_gpp => clm3%g%c%p%pepv%annsum_potential_gpp
   tempmax_retransn => clm3%g%c%p%pepv%tempmax_retransn
   annmax_retransn => clm3%g%c%p%pepv%annmax_retransn
   avail_retransn => clm3%g%c%p%pepv%avail_retransn
   plant_nalloc => clm3%g%c%p%pepv%plant_nalloc
   plant_calloc => clm3%g%c%p%pepv%plant_calloc
   excess_cflux => clm3%g%c%p%pepv%excess_cflux
   downreg => clm3%g%c%p%pepv%downreg
   prev_leafc_to_litter => clm3%g%c%p%pepv%prev_leafc_to_litter
   prev_frootc_to_litter => clm3%g%c%p%pepv%prev_frootc_to_litter
   tempsum_npp => clm3%g%c%p%pepv%tempsum_npp
   annsum_npp => clm3%g%c%p%pepv%annsum_npp
   leafc_storage => clm3%g%c%p%pcs%leafc_storage
   leafc_xfer => clm3%g%c%p%pcs%leafc_xfer
   frootc_storage => clm3%g%c%p%pcs%frootc_storage
   frootc_xfer => clm3%g%c%p%pcs%frootc_xfer
   livestemc_storage => clm3%g%c%p%pcs%livestemc_storage
   livestemc_xfer => clm3%g%c%p%pcs%livestemc_xfer
   deadstemc_storage => clm3%g%c%p%pcs%deadstemc_storage
   deadstemc_xfer => clm3%g%c%p%pcs%deadstemc_xfer
   livecrootc_storage => clm3%g%c%p%pcs%livecrootc_storage
   livecrootc_xfer => clm3%g%c%p%pcs%livecrootc_xfer
   deadcrootc_storage => clm3%g%c%p%pcs%deadcrootc_storage
   deadcrootc_xfer => clm3%g%c%p%pcs%deadcrootc_xfer
   gresp_storage => clm3%g%c%p%pcs%gresp_storage
   gresp_xfer => clm3%g%c%p%pcs%gresp_xfer
   cpool => clm3%g%c%p%pcs%cpool
   xsmrpool => clm3%g%c%p%pcs%xsmrpool
   pft_ctrunc => clm3%g%c%p%pcs%pft_ctrunc
   leafn_storage => clm3%g%c%p%pns%leafn_storage
   leafn_xfer => clm3%g%c%p%pns%leafn_xfer
   frootn_storage => clm3%g%c%p%pns%frootn_storage
   frootn_xfer => clm3%g%c%p%pns%frootn_xfer
   livestemn_storage => clm3%g%c%p%pns%livestemn_storage
   livestemn_xfer => clm3%g%c%p%pns%livestemn_xfer
   deadstemn_storage => clm3%g%c%p%pns%deadstemn_storage
   deadstemn_xfer => clm3%g%c%p%pns%deadstemn_xfer
   livecrootn_storage => clm3%g%c%p%pns%livecrootn_storage
   livecrootn_xfer => clm3%g%c%p%pns%livecrootn_xfer
   deadcrootn_storage => clm3%g%c%p%pns%deadcrootn_storage
   deadcrootn_xfer => clm3%g%c%p%pns%deadcrootn_xfer
   retransn => clm3%g%c%p%pns%retransn
   npool => clm3%g%c%p%pns%npool
   pft_ntrunc => clm3%g%c%p%pns%pft_ntrunc
   decl => clm3%g%c%cps%decl
   fpi => clm3%g%c%cps%fpi
   fpg => clm3%g%c%cps%fpg
   annsum_counter => clm3%g%c%cps%annsum_counter
   cannsum_npp => clm3%g%c%cps%cannsum_npp
   cannsum_t2m => clm3%g%c%cps%cannavg_t2m
   me => clm3%g%c%cps%me
   fire_prob => clm3%g%c%cps%fire_prob
   mean_fire_prob => clm3%g%c%cps%mean_fire_prob
   fireseasonl => clm3%g%c%cps%fireseasonl
   farea_burned => clm3%g%c%cps%farea_burned
   ann_farea_burned => clm3%g%c%cps%ann_farea_burned
   seedc => clm3%g%c%ccs%seedc
   col_ctrunc => clm3%g%c%ccs%col_ctrunc
   totcolc => clm3%g%c%ccs%totcolc
   prod10c => clm3%g%c%ccs%prod10c
   prod100c => clm3%g%c%ccs%prod100c
   seedn => clm3%g%c%cns%seedn
   col_ntrunc => clm3%g%c%cns%col_ntrunc
   totcoln => clm3%g%c%cns%totcoln
   prod10n => clm3%g%c%cns%prod10n
   prod100n => clm3%g%c%cns%prod100n
   litfall => clm3%g%c%p%pcf%litfall
   cisun => clm3%g%c%p%pps%cisun
   cisha => clm3%g%c%p%pps%cisha
   rssun => clm3%g%c%p%pps%rssun
   rssha => clm3%g%c%p%pps%rssha
   parsun => clm3%g%c%p%pef%parsun
   parsha => clm3%g%c%p%pef%parsha
   laisun => clm3%g%c%p%pps%laisun
   laisha => clm3%g%c%p%pps%laisha

   do c = begc, endc
     cn(c)%soil1c = soil1c(c)
     cn(c)%soil2c = soil2c(c)
     cn(c)%soil3c = soil3c(c)
     cn(c)%soil4c = soil4c(c)
     cn(c)%litr1c = litr1c(c)
     cn(c)%litr2c = litr2c(c)
     cn(c)%litr3c = litr3c(c)
     cn(c)%cwdc = cwdc(c)
     cn(c)%totlitc = totlitc(c)
     cn(c)%totsomc = totsomc(c)
     cn(c)%soil1n = soil1n(c)
     cn(c)%soil2n = soil2n(c)
     cn(c)%soil3n = soil3n(c)
     cn(c)%soil4n = soil4n(c)
     cn(c)%sminn = sminn(c)
     cn(c)%litr1n = litr1n(c)
     cn(c)%litr2n = litr2n(c)
     cn(c)%litr3n = litr3n(c)
     cn(c)%cwdn = cwdn(c)
     cn(c)%hr = hr(c)
     cn(c)%lithr = lithr(c)
     cn(c)%nee = nee(c)
     cn(c)%nep = nep(c)
     cn(c)%decl = decl(c)
     cn(c)%fpi = fpi(c)
     cn(c)%fpg = fpg(c)
     cn(c)%annsum_counter = annsum_counter(c)
     cn(c)%cannsum_npp = cannsum_npp(c)
     cn(c)%cannavg_t2m = cannsum_t2m(c)
     do k = 1, nlevgrnd
       cn(c)%watfc(k) = watfc(c,k)
     end do
     cn(c)%me = me(c)
     cn(c)%fire_prob = fire_prob(c)
     cn(c)%mean_fire_prob = mean_fire_prob(c)
     cn(c)%fireseasonl = fireseasonl(c)
     cn(c)%farea_burned = farea_burned(c)
     cn(c)%ann_farea_burned = ann_farea_burned(c)
     cn(c)%seedc = seedc(c)
     cn(c)%col_ctrunc = col_ctrunc(c)
     cn(c)%totcolc = totcolc(c)
     cn(c)%prod10c = prod10c(c)
     cn(c)%prod100c = prod100c(c)
     cn(c)%seedn = seedn(c)
     cn(c)%col_ntrunc = col_ntrunc(c)
     cn(c)%totcoln = totcoln(c)
     cn(c)%prod10n = prod10n(c)
     cn(c)%prod100n = prod100n(c)
!	print *, c, cn(c)%totlitc, cn(c)%soil2c, cn(c)%soil3c, totlitc(c), soil2c(c), soil3c(c)
   end do

   do p = begp, endp
     c = ncol(p)
       cn(c)%LAI(ivt(p)) = elai(p)
       cn(c)%leafc(ivt(p)) = leafc(p)
       cn(c)%frootc(ivt(p)) = frootc(p)
       cn(c)%livestemc(ivt(p)) = livestemc(p)
       cn(c)%deadstemc(ivt(p)) = deadstemc(p)
       cn(c)%livecrootc(ivt(p)) = livecrootc(p)
       cn(c)%deadcrootc(ivt(p)) = deadcrootc(p)
       cn(c)%woodc(ivt(p)) = woodc(p)
       cn(c)%totvegc(ivt(p)) = totvegc(p)
       cn(c)%leafn(ivt(p)) = leafn(p)
       cn(c)%frootn(ivt(p)) = frootn(p)
       cn(c)%livestemn(ivt(p)) = livestemn(p)
       cn(c)%deadstemn(ivt(p)) = deadstemn(p)
       cn(c)%livecrootn(ivt(p)) = livecrootn(p)
       cn(c)%deadcrootn(ivt(p)) = deadcrootn(p)
       cn(c)%gpp2(ivt(p)) = gpp2(p)
       cn(c)%gpp(ivt(p)) = gpp(p)
       cn(c)%npp(ivt(p)) = npp(p)
       cn(c)%leaf_mr(ivt(p)) = leaf_mr(p)
       cn(c)%mr(ivt(p)) = mr(p)
       cn(c)%gr(ivt(p)) = gr(p)
       cn(c)%ar(ivt(p)) = ar(p)
       cn(c)%dormant_flag(ivt(p)) = dormant_flag(p)
       cn(c)%days_active(ivt(p)) = days_active(p)
       cn(c)%onset_flag(ivt(p)) = onset_flag(p)
       cn(c)%onset_counter(ivt(p)) = onset_counter(p)
       cn(c)%onset_gddflag(ivt(p)) = onset_gddflag(p)
       cn(c)%onset_fdd(ivt(p)) = onset_fdd(p)
       cn(c)%onset_gdd(ivt(p)) = onset_gdd(p)
       cn(c)%onset_swi(ivt(p)) = onset_swi(p)
       cn(c)%offset_flag(ivt(p)) = offset_flag(p)
       cn(c)%offset_counter(ivt(p)) = offset_counter(p)
       cn(c)%offset_fdd(ivt(p)) = offset_fdd(p)
       cn(c)%offset_swi(ivt(p)) = offset_swi(p)
       cn(c)%lgsf(ivt(p)) = lgsf(p)
       cn(c)%bglfr(ivt(p)) = bglfr(p)
       cn(c)%bgtr(ivt(p)) = bgtr(p)
       cn(c)%dayl(ivt(p)) = dayl(p)
       cn(c)%prev_dayl(ivt(p)) = prev_dayl(p)
       cn(c)%annavg_t2m(ivt(p)) = annavg_t2m(p)
       cn(c)%tempavg_t2m(ivt(p)) = tempavg_t2m(p)
       cn(c)%availc(ivt(p)) = availc(p)
       cn(c)%xsmrpool_recover(ivt(p)) = xsmrpool_recover(p)
       cn(c)%alloc_pnow(ivt(p)) = alloc_pnow(p)
       cn(c)%c_allometry(ivt(p)) = c_allometry(p)
       cn(c)%n_allometry(ivt(p)) = n_allometry(p)
       cn(c)%plant_ndemand(ivt(p)) = plant_ndemand(p)
       cn(c)%tempsum_potential_gpp(ivt(p)) = tempsum_potential_gpp(p)
       cn(c)%annsum_potential_gpp(ivt(p)) = annsum_potential_gpp(p)
       cn(c)%tempmax_retransn(ivt(p)) = tempmax_retransn(p)
       cn(c)%annmax_retransn(ivt(p)) = annmax_retransn(p)
       cn(c)%avail_retransn(ivt(p)) = avail_retransn(p)
       cn(c)%plant_nalloc(ivt(p)) = plant_nalloc(p)
       cn(c)%plant_calloc(ivt(p)) = plant_calloc(p)
       cn(c)%excess_cflux(ivt(p)) = excess_cflux(p)
       cn(c)%downreg(ivt(p)) = downreg(p)
       cn(c)%prev_leafc_to_litter(ivt(p)) = prev_leafc_to_litter(p)
       cn(c)%prev_frootc_to_litter(ivt(p)) = prev_frootc_to_litter(p)
       cn(c)%tempsum_npp(ivt(p)) = tempsum_npp(p)
       cn(c)%annsum_npp(ivt(p)) = annsum_npp(p)
       cn(c)%leafc_storage(ivt(p)) = leafc_storage(p)
       cn(c)%leafc_xfer(ivt(p)) = leafc_xfer(p)
       cn(c)%frootc_storage(ivt(p)) = frootc_storage(p)
       cn(c)%frootc_xfer(ivt(p)) = frootc_xfer(p)
       cn(c)%livestemc_storage(ivt(p)) = livestemc_storage(p)
       cn(c)%livestemc_xfer(ivt(p)) = livestemc_xfer(p)
       cn(c)%deadstemc_storage(ivt(p)) = deadstemc_storage(p)
       cn(c)%deadstemc_xfer(ivt(p)) = deadstemc_xfer(p)
       cn(c)%livecrootc_storage(ivt(p)) = livecrootc_storage(p)
       cn(c)%livecrootc_xfer(ivt(p)) = livecrootc_xfer(p)
       cn(c)%deadcrootc_storage(ivt(p)) = deadcrootc_storage(p)
       cn(c)%deadcrootc_xfer(ivt(p)) = deadcrootc_xfer(p)
       cn(c)%gresp_storage(ivt(p)) = gresp_storage(p)
       cn(c)%gresp_xfer(ivt(p)) = gresp_xfer(p)
       cn(c)%cpool(ivt(p)) = cpool(p)
       cn(c)%xsmrpool(ivt(p)) = xsmrpool(p)
       cn(c)%pft_ctrunc(ivt(p)) = pft_ctrunc(p)
       cn(c)%leafn_storage(ivt(p)) = leafn_storage(p)
       cn(c)%leafn_xfer(ivt(p)) = leafn_xfer(p)
       cn(c)%frootn_storage(ivt(p)) = frootn_storage(p)
       cn(c)%frootn_xfer(ivt(p)) = frootn_xfer(p)
       cn(c)%livestemn_storage(ivt(p)) = livestemn_storage(p)
       cn(c)%livestemn_xfer(ivt(p)) = livestemn_xfer(p)
       cn(c)%deadstemn_storage(ivt(p)) = deadstemn_storage(p)
       cn(c)%deadstemn_xfer(ivt(p)) = deadstemn_xfer(p)
       cn(c)%livecrootn_storage(ivt(p)) = livecrootn_storage(p)
       cn(c)%livecrootn_xfer(ivt(p)) = livecrootn_xfer(p)
       cn(c)%deadcrootn_storage(ivt(p)) = deadcrootn_storage(p)
       cn(c)%deadcrootn_xfer(ivt(p)) = deadcrootn_xfer(p)
       cn(c)%retransn(ivt(p)) = retransn(p)
       cn(c)%npool(ivt(p)) = npool(p)
       cn(c)%pft_ntrunc(ivt(p)) = pft_ntrunc(p)
       cn(c)%litfall(ivt(p)) = litfall(p) * dt
       cn(c)%fpsn(ivt(p)) = fpsn(p)
       cn(c)%ci(ivt(p)) = cisun(p) * laisun(p) + cisha(p) * laisha(p)
       cn(c)%rs(ivt(p)) = rssun(p) * laisun(p) + rssha(p) * laisha(p)
       cn(c)%par(ivt(p)) = parsun(p) * laisun(p) + parsha(p) * laisha(p)
!	print *, p, ivt(p), elai(p), cn(c)%LAI(ivt(p))
	print *, p, ivt(p), leafc(p), cn(c)%leafc(ivt(p))
   end do

end subroutine vic2clmtype
