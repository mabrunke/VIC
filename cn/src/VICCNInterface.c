#include <stdio.h>
#include <stdlib.h>
#include <vic_def.h>
#include <vic_run.h>
#include <string.h>
#define MAX_DIST 2
#define MAX_PFT 20
#define NUMRAD 2

extern void vic2clmtype(int *nlevgrnd, int *rec, int *nrec, int *adspinup, \
			int *init_state, int *yr, int *mo, int *day,	\
			 int *secs, double *jday, int *yrnxt, double *jdynxt, \
			 double *dt, double *lat, double *lon, \
			 int *begg, int *endg, int *begc,      \
			 int *endc, int *begp, int *endp, \
			 int *num_soilc, int *num_soilp, \
			 cn_data_struct *cn);

int VICCNInterface(int                 rec,
                   int                 Nveg,
		   int                 Npft,
		   dmy_struct          *dmy, 
		   global_param_struct *global_param,
                   atmos_data_struct   *atmos,
		   all_vars_struct    *all_vars,
		   soil_con_struct     *soil_con,
		   veg_con_struct      *veg_con,
                   veg_lib_struct      *veg_lib,
                   cn_data_struct      *cn)

/**********************************************************************
	VICCNInterface	Michael Brunke		July 23, 2014

Prepares VIC data for use in CN.

Modifications:
  24-Jul-2014 Replaced dist by bands.                                   MAB
  15-Oct-2014 Added veg_var and veg_hist structures as input.           MAB

**********************************************************************/

{

  extern option_struct options;

  double z[MAX_BANDS][MAX_NODES+2];             /* Layer depth (m) */
  double dz[MAX_BANDS][MAX_NODES+2];            /* Layer thickness (m) */
  double baseflow[MAX_BANDS];                   /* Baseflow */
  double moist[MAX_BANDS][MAX_NODES+2];          /* Soil moisture */
  double ice[MAX_BANDS][MAX_NODES+2];            /* Soil ice lens */
  double t_soisno[MAX_BANDS][MAX_NODES+2];     /* Snow/soil temperature */
  double t2m[MAX_CN_PFTS];                         /* 2-m air temperature */
  double snowdep[MAX_BANDS];                    /* Snow depth */
  double fwet[MAX_BANDS][MAX_CN_PFTS+1];            /* Wet veg. fraction */
  double rootfr[MAX_BANDS][MAX_CN_PFTS][MAX_NODES];/* Root fraction in ea. layer */
  double psisat[2][MAX_BANDS][MAX_NODES];      /* Saturated matric potential */
  double soipsi[2][MAX_BANDS][MAX_NODES];         /* Matric potential */
  double coeff[2][MAX_BANDS][MAX_NODES];          /* Clapp-Hornberger coeff. */
  double rveg[MAX_BANDS][MAX_CN_PFTS+1];           /* Leaf resistance */
  double *rcan;                               /* Canopy resistance */
  double LAI[MAX_BANDS][MAX_CN_PFTS+1];             /* Leaf area index */
  double soilcfast[MAX_BANDS];                  /* Fast soil C pool */
  double soilcmid[MAX_BANDS];                   /* Medium soil C pool */
  double soilcslo1[MAX_BANDS];                  /* Slow soil C pool */
  double soilcslo2[MAX_BANDS];                  /* Slowest soil C pool */
  double litrlabc[MAX_BANDS];                   /* Litter labile C pool */
  double litrcellc[MAX_BANDS];                  /* Litter cellulose C pool */
  double litrligc[MAX_BANDS];                   /* Litter lignin C pool */
  double cwdc[MAX_BANDS];                       /* Coarse woody debris C pool */
  double leafc[MAX_BANDS][MAX_CN_PFTS+1];           /* Leaf C pool */
  double frootc[MAX_BANDS][MAX_CN_PFTS+1];          /* Fine root C pool */
  double livestemc[MAX_BANDS][MAX_CN_PFTS+1];       /* Live stem C pool */
  double deadstemc[MAX_BANDS][MAX_CN_PFTS+1];       /* Dead stem C pool */
  double livecrootc[MAX_BANDS][MAX_CN_PFTS+1];      /* Live coarse root C pool */
  double deadcrootc[MAX_BANDS][MAX_CN_PFTS+1];      /* Dead coarse root C pool */
  double woodc[MAX_BANDS][MAX_CN_PFTS+1];           /* Wood C pool */
  double totvegc[MAX_BANDS][MAX_CN_PFTS+1];         /* Total vegetation C */
  double totlitc[MAX_BANDS];                    /* Total litter C */
  double totsomc[MAX_BANDS];                    /* Total SOM C */
  double soilnfast[MAX_BANDS];                  /* Fast soil N pool */
  double soilnmid[MAX_BANDS];                   /* Medium soil N pool */
  double soilnslo1[MAX_BANDS];                  /* Slow soil N pool */
  double soilnslo2[MAX_BANDS];                  /* Slowest soil N pool */
  double soilminn[MAX_BANDS];                   /* Mineral soil N pool */
  double litrlabn[MAX_BANDS];                   /* Litter labile N pool */
  double litrcelln[MAX_BANDS];                  /* Litter cellulose N pool */
  double litrlign[MAX_BANDS];                   /* Litter lignin N pool */
  double cwdn[MAX_BANDS];                       /* Coarse woody debris N pool */
  double leafn[MAX_BANDS][MAX_CN_PFTS+1];           /* Leaf N pool */
  double frootn[MAX_BANDS][MAX_CN_PFTS+1];          /* Fine root N pool */
  double livestemn[MAX_BANDS][MAX_CN_PFTS+1];       /* Live stem N pool */
  double deadstemn[MAX_BANDS][MAX_CN_PFTS+1];       /* Dead stem N pool */
  double livecrootn[MAX_BANDS][MAX_CN_PFTS+1];      /* Live coarse root N pool */
  double deadcrootn[MAX_BANDS][MAX_CN_PFTS+1];      /* Dead coarse root N pool */
  double totlitn[MAX_BANDS];                    /* Total litter N */
  double gpp2[MAX_BANDS][MAX_CN_PFTS+1];            /* GPP before downregulation */
  double gpp[MAX_BANDS][MAX_CN_PFTS+1];             /* Gross primary production */
  double npp[MAX_BANDS][MAX_CN_PFTS+1];             /* Net primary production */
  double darkr[MAX_BANDS][MAX_CN_PFTS+1];           /* Leaf maint respiration */
  double mr[MAX_BANDS][MAX_CN_PFTS+1];              /* Maintenance respiration */
  double gr[MAX_BANDS][MAX_CN_PFTS+1];              /* Growth respiration */
  double ar[MAX_BANDS][MAX_CN_PFTS+1];              /* Autotrophic respiration */
  double hr[MAX_BANDS];                         /* Heterotrophic respiration */
  double lithr[MAX_BANDS];                      /* Litter hetero respiration */
  double nee[MAX_BANDS];                        /* Net ecosystem exchange */
  double nep[MAX_BANDS];                        /* Net ecosystem production */

  double dormant_flag[MAX_BANDS][MAX_CN_PFTS+1];    /* Dormancy flag */
  double days_active[MAX_BANDS][MAX_CN_PFTS+1];     /* # days since last dormancy */
  double onset_flag[MAX_BANDS][MAX_CN_PFTS+1];      /* Onset flag */
  double onset_counter[MAX_BANDS][MAX_CN_PFTS+1];   /* Onset days counter */
  double onset_gddflag[MAX_BANDS][MAX_CN_PFTS+1];  /* Onset flag for grow deg sum */
  double onset_fdd[MAX_BANDS][MAX_CN_PFTS+1];       /* Onset freeze days counter */
  double onset_gdd[MAX_BANDS][MAX_CN_PFTS+1];       /* Onset growing deg days */
  double onset_swi[MAX_BANDS][MAX_CN_PFTS+1];       /* Onset soil water index */
  double offset_flag[MAX_BANDS][MAX_CN_PFTS+1];      /* Offset flag */
  double offset_counter[MAX_BANDS][MAX_CN_PFTS+1];   /* Offset days counter */
  double offset_fdd[MAX_BANDS][MAX_CN_PFTS+1];       /* Offset freeze days counter */
  double offset_swi[MAX_BANDS][MAX_CN_PFTS+1];       /* Offset soil water index */
  double lgsf[MAX_BANDS][MAX_CN_PFTS+1];            /* Long growing season factor */
  double bglfr[MAX_BANDS][MAX_CN_PFTS+1];           /* Background litterfall rate */
  double bgtr[MAX_BANDS][MAX_CN_PFTS+1];            /* Background transfer growth */
  double dayl[MAX_BANDS][MAX_CN_PFTS+1];            /* Daylength (s) */
  double prev_dayl[MAX_BANDS][MAX_CN_PFTS+1];       /* Previous daylength (s) */
  double annavg_t2m[MAX_BANDS][MAX_CN_PFTS+1];      /* Annual avg 2m air temp. */
  double tempavg_t2m[MAX_BANDS][MAX_CN_PFTS+1];     /* Temp. avg. 2m air temp. */
  double availc[MAX_BANDS][MAX_CN_PFTS+1];          /* C flux avail for alloc */
  double xsmrpool_recover[MAX_BANDS][MAX_CN_PFTS+1];/* C flx assigned to recovery */
  double alloc_pnow[MAX_BANDS][MAX_CN_PFTS+1];    /* Fract of alloc as new growth */
  double c_allometry[MAX_BANDS][MAX_CN_PFTS+1];     /* C allocation index */
  double n_allometry[MAX_BANDS][MAX_CN_PFTS+1];     /* N allocation index */
  double plant_ndemand[MAX_BANDS][MAX_CN_PFTS+1];   /* N flux to support GPP */
  double tempsum_potential_gpp[MAX_BANDS][MAX_CN_PFTS+1];/* Temp ann sum of pot GPP */
  double annsum_potential_gpp[MAX_BANDS][MAX_CN_PFTS+1];/* Ann sum of pot GPP */
  double tempmax_retransn[MAX_BANDS][MAX_CN_PFTS+1];/* Temp ann max retrans N */
  double annmax_retransn[MAX_BANDS][MAX_CN_PFTS+1]; /* Ann max retrans N pool */
  double avail_retransn[MAX_BANDS][MAX_CN_PFTS+1];  /* N flux from retrans pool */
  double plant_nalloc[MAX_BANDS][MAX_CN_PFTS+1];    /* Total allocated N flux */
  double plant_calloc[MAX_BANDS][MAX_CN_PFTS+1];    /* Total allocated C flux */
  double excess_cflux[MAX_BANDS][MAX_CN_PFTS+1];    /* C flux not allocated */
  double downreg[MAX_BANDS][MAX_CN_PFTS+1];        /* Fract reduct GPP from N lim */
  double prev_leafc_to_litter[MAX_BANDS][MAX_CN_PFTS+1];/* Prev leaf C litterfall */
  double prev_frootc_to_litter[MAX_BANDS][MAX_CN_PFTS+1];/* Prev froot C litterfall */
  double tempsum_npp[MAX_BANDS][MAX_CN_PFTS+1];     /* Temp ann sum of NPP */
  double annsum_npp[MAX_BANDS][MAX_CN_PFTS+1];      /* Annual sum of NPP */
  double leafc_storage[MAX_BANDS][MAX_CN_PFTS+1];   /* Leaf C storage */
  double leafc_xfer[MAX_BANDS][MAX_CN_PFTS+1];      /* Leaf C transfer */
  double frootc_storage[MAX_BANDS][MAX_CN_PFTS+1];  /* Fine root C storage */
  double frootc_xfer[MAX_BANDS][MAX_CN_PFTS+1];     /* Fine root C transfer */
  double livestemc_storage[MAX_BANDS][MAX_CN_PFTS+1];/* Live stem C storage */
  double livestemc_xfer[MAX_BANDS][MAX_CN_PFTS+1];  /* Live stem C transfer */
  double deadstemc_storage[MAX_BANDS][MAX_CN_PFTS+1]; /* Dead stem C storage */
  double deadstemc_xfer[MAX_BANDS][MAX_CN_PFTS+1];  /* Dead stem C transfer */
  double livecrootc_storage[MAX_BANDS][MAX_CN_PFTS+1];/* Live coarse root C storage */
  double livecrootc_xfer[MAX_BANDS][MAX_CN_PFTS+1]; /* Live coarse root C transfer */
  double deadcrootc_storage[MAX_BANDS][MAX_CN_PFTS+1];/* Dead coarse root C storage */
  double deadcrootc_xfer[MAX_BANDS][MAX_CN_PFTS+1]; /* Dead coarse root C transfer */
  double gresp_storage[MAX_BANDS][MAX_CN_PFTS+1];   /* Growth respiration storage */
  double gresp_xfer[MAX_BANDS][MAX_CN_PFTS+1];      /* Growth respiration transfer */
  double cpool[MAX_BANDS][MAX_CN_PFTS+1];           /* Temp photosynthate C pool */
  double xsmrpool[MAX_BANDS][MAX_CN_PFTS+1];        /* Abstract C pool */
  double pft_ctrunc[MAX_BANDS][MAX_CN_PFTS+1];      /* PFT sink for C truncation */
  double leafn_storage[MAX_BANDS][MAX_PFT+1];   /* Leaf N storage */
  double leafn_xfer[MAX_BANDS][MAX_PFT+1];      /* Leaf N storage */
  double frootn_storage[MAX_BANDS][MAX_PFT+1];  /* Fine root N storage */
  double frootn_xfer[MAX_BANDS][MAX_PFT+1];     /* Fine root N transfer */
  double livestemn_storage[MAX_BANDS][MAX_PFT+1]; /* Live stem N storage */
  double livestemn_xfer[MAX_BANDS][MAX_PFT+1];  /* Live stem N transfer */
  double deadstemn_storage[MAX_BANDS][MAX_PFT+1]; /* Dead stem N storage */
  double deadstemn_xfer[MAX_BANDS][MAX_PFT+1];  /* Dead stem N transfer */
  double livecrootn_storage[MAX_BANDS][MAX_PFT+1]; /* Live coarse root N storage */
  double livecrootn_xfer[MAX_BANDS][MAX_PFT+1]; /* Live coarse root N transfer */
  double deadcrootn_storage[MAX_BANDS][MAX_PFT+1]; /* Dead coarse root N storage */
  double deadcrootn_xfer[MAX_BANDS][MAX_PFT+1]; /* Dead coarse root N transfer */
  double retransn[MAX_BANDS][MAX_PFT+1];        /* Retranslocated N */
  double npool[MAX_BANDS][MAX_PFT+1];           /* Temp. photosynthate N pool */
  double pft_ntrunc[MAX_BANDS][MAX_PFT+1];      /* PFT sink for N truncation */
  double decl[MAX_BANDS];                       /* Solar declination angle */
  double fpi[MAX_BANDS];                        /* Fract. pot. immobilization */
  double fpg[MAX_BANDS];                        /* Fract. potential GPP */
  double annsum_counter[MAX_BANDS];             /* Secs since last ann accum. turnover */
  double cannsum_npp[MAX_BANDS];                /* Annual sum of NPP */
  double cannavg_t2m[MAX_BANDS];                /* Annual avg. of 2-m air temp. */
  double watfc[MAX_BANDS][MAX_NODES];           /* Volumetric soil moisture at field capcity */
  double me[MAX_BANDS];                         /* Moisture of extinction */
  double fire_prob[MAX_BANDS];                  /* Daily fire probability */
  double mean_fire_prob[MAX_BANDS];             /* E-fold mean daily fire prob. */
  double fireseasonl[MAX_BANDS];                /* Ann. fire season length */
  double farea_burned[MAX_BANDS];               /* Timestep fract. area burned */
  double ann_farea_burned[MAX_BANDS];           /* Ann. tot. fract. area burned */
  double seedc[MAX_BANDS];                      /* Col-lev pool for seeding new PFTs */
  double col_ctrunc[MAX_BANDS];                 /* Col-lev sink for C trunc */
  double totcolc[MAX_BANDS];                    /* Total column C */
  double prod10c[MAX_BANDS];                    /* Wood product C pool, 10-yr lifespan */
  double prod100c[MAX_BANDS];                   /* Wood product C pool, 100-yr lifespan */
  double seedn[MAX_BANDS];                      /* Col-lev pool for seeding new PFTs */
  double col_ntrunc[MAX_BANDS];                 /* Col-lev sink for N trunc */
  double totcoln[MAX_BANDS];                    /* Total column N */
  double prod10n[MAX_BANDS];                    /* Wood product N pool, 10-yr lifespan */
  double prod100n[MAX_BANDS];                   /* Wood product N pool, 100-yr lifespan */
  double litfall[MAX_BANDS][MAX_PFT+1];         /* Litterfall */
  double fpsn[MAX_BANDS][MAX_PFT+1];             /* Photosynthesis */
  double ci[MAX_BANDS][MAX_PFT+1];               /* Intracellular C02 */
  double rc[MAX_BANDS][MAX_PFT+1];               /* Canopy stomatal resistance */
  double apar[MAX_BANDS][MAX_PFT+1];            /* absorbed PAR */

  double Tair;                                 /* Air temperature */
  double *Tcan;                                /* Canopy temperature */
  double Tveg[MAX_BANDS][MAX_PFT+1];            /* Leaf temperature */
  double psfc;                                 /* Surface pressure */
  double vp;                                   /* Atmos. vapor pressure */
  double vpd;                                  /* Vapor pressure deficit */
  double lwrad;                                /* LW radiation */
  double swrad;                                /* SW radiation */
  double swrd[NUMRAD];                         /* direct SW radiation */
  double swri[NUMRAD];                         /* diffuse SW radiation */
  double alb[MAX_BANDS];                        /* surface albedo */
  double zo, zos, zov[MAX_BANDS][MAX_PFT+1];      /* roughness lengths */
  double displ[MAX_BANDS][MAX_PFT+1];             /* displacement height */
  double lat;
  double lon;

  int i;
  int iveg;
  int band;
  int Nbands;
  int lidx, nidx;
  int Nlayers;
  int Nnodes;
  int Nrecs;
  int adspinup;
  int lbg;
  int ubg;
  int lbc;
  int ubc;
  int lbp;
  int ubp;
  int num_soilc;
  int num_soilp;
  int init_state;
  int year, month, day, hour, secs, yrnxt;
  double jday, jdynxt;
  int err = 0;
  int veg_class;
  double factor;
  double precip;
  double sucsat[2][MAX_LAYERS];      /* Saturated matric potential in layers */
  double soisuc[2][MAX_LAYERS];      /* Soil matric potential */
  double theta[MAX_LAYERS];          /* Volumetric water content */
  double watsat[MAX_LAYERS];         /* Saturated volumetric water content */
  double psand;                      /* % sand */
  double pclay;                      /* % clay */
  double bsw[2][MAX_LAYERS];         /* Clapp-Hornberger coefficient */
  double ksat;                       /* Saturated water conductivity */
  double watfrac[MAX_LAYERS];        /* Water content at field capacity */
  double fsat;                       /* Saturation fraction */
  double Lsum;                        /* cumulative depth of moisture layer */
  double fsum[MAX_PFT + 1];
  char PAST_BOTTOM;

  cell_data_struct **cell;
  energy_bal_struct **energy;
  snow_data_struct  **snow;
  veg_var_struct **veg;

  double dt;
  double b;

  /* Set local pointers */
  cell = all_vars->cell;
  energy = all_vars->energy;
  snow = all_vars->snow;
  veg = all_vars->veg_var;

  if(options.CARBON == CN_ADECOMP)
    adspinup = 1;
  else
    adspinup = 0;
  if(options.INIT_STATE)
    init_state = 1;
  else
    init_state = 0;

  /* Set number of snow bands */
  Nbands = options.SNOW_BAND;

  /* Set number of soil layers */
  Nlayers = options.Nlayer;

  /* Set number of nodes */
  Nnodes = options.Nnode;

  /* Set number of records */
  Nrecs = global_param->nrecs;

  /* Convert VIC time step to seconds */
  dt = (double) global_param->dt * 3600.0;

  /* Get current date from dmy */
  year = dmy[rec].year;
  month = dmy[rec].month;
  day = dmy[rec].day;
  hour = dmy[rec].hour;
  secs = hour * 3600.0;
  jday = dmy[rec].day_in_year;

  /* Get next date from dmy, MAB 10/8/13 */
  if(rec + 1 < Nrecs)
    {
    yrnxt = dmy[rec+1].year;
    jdynxt = dmy[rec+1].day_in_year;
    }
  else
    {
    yrnxt = -1;
    jdynxt = -1.0;
    }

  /* Get latitude/longitude, MAB 8/27/13 */
  lat = (double) soil_con->lat;
  lon = (double) soil_con->lng;

  /* Place atmospheric quantities into CN structure, MAB 1/23/15 */
  get_CNatmos(lat, atmos, cn);

  get_CNsoil(Nbands, Nnodes, Nveg, soil_con, all_vars, cn);

  get_CNsoiltherm(Nbands, Nnodes, Nveg, all_vars, cn);

  get_CNsoilhydro(Nbands, Nnodes, Nlayers, Nveg, Npft, soil_con, veg_con, \
		  all_vars, cn);

  get_CNveg(Nbands, Nnodes, atmos->prec[NR], veg_con, veg_lib, all_vars, cn);

  /* This will need to be adapted when implemented into RASM or when VIC is */
  /* run regionally. */

  lbg = 1;
  ubg = 1;
  lbc = 1;
  ubc = Nbands;
  lbp = 1;
  ubp = Npft;
  num_soilc = Nbands;
  num_soilp = Npft;

  /* use_cn_data_struct((size_t) Nbands, (size_t) Nnodes, (size_t) Npft, cn); */
  vic2clmtype(&Nnodes, &rec, &Nrecs, &adspinup, &init_state, &year, &month, \
	       &day, &secs, &jday, &yrnxt, &jdynxt, &dt, &lat, &lon, &lbg, \
	       &ubg, &lbc, &ubc, &lbp, &ubp, &num_soilc, &num_soilp, cn);

  for(band = 0; band < Nbands; band++)
    {
      for(iveg = 0; iveg < Nveg; iveg++)

	{

	  /* printf("%d %f %f %f\n", cell[iveg][band].CLitter, cell[iveg][band].CInter, cell[iveg][band].CSlow); */
	  cell[iveg][band].CLitter = cn[band].totlitc;
	  cell[iveg][band].CInter = cn[band].soil2c;
	  cell[iveg][band].CSlow = cn[band].soil3c;
	  cell[iveg][band].RhLitter = cn[band].lithr;
	  cell[iveg][band].RhTot = cn[band].hr;
	  /* printf("%d %f %f %f %f %f %f\n", cn[band].totlitc, cn[band].soil2c, cn[band].soil3c, cell[iveg][band].CLitter, cell[iveg][band].CInter, cell[iveg][band].CSlow); */

	  veg[iveg][band].LAI = pft2cov(cn[band].LAI, veg_con[iveg].veg_class);
	  veg[iveg][band].aPAR = pft2cov(cn[band].par, \
					 veg_con[iveg].veg_class);
	  veg[iveg][band].Ci = pft2cov(cn[band].ci, veg_con[iveg].veg_class);
	  veg[iveg][band].rc = pft2cov(cn[band].rs, veg_con[iveg].veg_class);
	  veg[iveg][band].GPP = pft2cov(cn[band].gpp, veg_con[iveg].veg_class);
	  veg[iveg][band].Rphoto = pft2cov(cn[band].fpsn,		\
	     veg_con[iveg].veg_class);
	  veg[iveg][band].Rdark = pft2cov(cn[band].leaf_mr, \
					   veg_con[iveg].veg_class);
	  veg[iveg][band].Rmaint = pft2cov(cn[band].mr, \
					   veg_con[iveg].veg_class);
	  veg[iveg][band].Rgrowth = pft2cov(cn[band].gr, \
					    veg_con[iveg].veg_class);
	  veg[iveg][band].Raut = pft2cov(cn[band].ar, veg_con[iveg].veg_class);
	  veg[iveg][band].NPP = pft2cov(cn[band].npp, veg_con[iveg].veg_class);
	  veg[iveg][band].Litterfall = pft2cov(cn[band].litfall, \
					       veg_con[iveg].veg_class);
	  veg[iveg][band].AnnualNPP = pft2cov(cn[band].annsum_npp, \
					      veg_con[iveg].veg_class);

	}

    }	

  return(err);

}
