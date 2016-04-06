/******************************************************************************
 * @section DESCRIPTION
 *
 * Header file for vic_run routines
 *
 * @section LICENSE
 *
 * The Variable Infiltration Capacity (VIC) macroscale hydrological model
 * Copyright (C) 2014 The Land Surface Hydrology Group, Department of Civil
 * and Environmental Engineering, University of Washington.
 *
 * The VIC model is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 *****************************************************************************/

#ifndef SYN_VEG_H
#define SYN_VEG_H

/******************************************************************
 * @brief   This structure stores output from CN, 
 *          added by MAB, 7/23/14
 *****************************************************************/
 typedef struct {
   /* atmospheric quantities */
   double Tair;                               /* air temperature */
   double vp;                                 /* vapor pressure */
   double vpd;                                /* vapor pressure depression */
   double psfc;                               /* surface pressure */
   double lwrad;                              /* longwave radiation */
   double swrad;                              /* shortwave radiation */
   double precip;                             /* precipitation */
   double swrd[NUMRAD];                    /* direct shortwave radiation */
   double swri[NUMRAD];                    /* diffuse shortwave radiation */
   double alb;                                /* albedo */
   double t2m[MAX_PFT+1];                   /* 2-m air temperature */
   double t_soisno[MAX_NODES+2];              /* soil/snow temperature */
   /* Soil properties */
   double z[MAX_NODES+2];                     /* Node depth */
   double dz[MAX_NODES+2];                    /* Node thickness */
   double z0;                                 /* Surface roughness */
   double z0s;                                /* Snow roughness */
   /* Vegetation characteristics */
   double Tveg[MAX_PFT+1];                  /* vegetation temperature */
   double rveg[MAX_PFT+1];                  /* vegetation resistance */
   double zov[MAX_PFT+1];                   /* vegetation roughness */
   double displ[MAX_PFT+1];                 /* displacement height */
   double fwet[MAX_PFT+1];                  /* wet veg fraction */
   /* Soil hydrology */
   double baseflow;                           /* baseflow */
   double moist[MAX_NODES];                   /* soil moisture */
   double ice[MAX_NODES+2];                     /* soil ice */
   double rootfr[MAX_PFT+1][MAX_NODES];     /* root fraction */
   double bsw[MAX_NODES][2];                  /* Clapp-Hornberger coefficient */
   double sucsat[MAX_NODES][2];               /* saturated suction */
   double soisuc[MAX_NODES][2];               /* soil matrix potential */
   double snowdep;                            /* snow depth */
   /* PFT-level ecophysiological variables */
   double LAI[MAX_PFT+1];                   /* leaf area index */
   double dormant_flag[MAX_PFT+1];          /* dormancy flag */
   double days_active[MAX_PFT+1];           /* # days since last dormancy */
   double onset_flag[MAX_PFT+1];            /* onset flag */
   double onset_counter[MAX_PFT+1];         /* onset days counter */
   double onset_gddflag[MAX_PFT+1];         /* onset flag for growing deg day sum */
   double onset_fdd[MAX_PFT+1];             /* onset freezing deg day counter */
   double onset_gdd[MAX_PFT+1];             /* onset growing degree days */
   double onset_swi[MAX_PFT+1];             /* onset soil water index */
   double offset_flag[MAX_PFT+1];            /* offset flag */
   double offset_counter[MAX_PFT+1];         /* offset days counter */
   double offset_fdd[MAX_PFT+1];             /* offset freezing deg day counter */
   double offset_swi[MAX_PFT+1];             /* offset soil water index */
   double lgsf[MAX_PFT+1];                  /* long growing season factor */
   double bglfr[MAX_PFT+1];                 /* background litterfall rate (1/s) */
   double bgtr[MAX_PFT+1];                  /* background transfer growth rate (1/s) */
   double dayl[MAX_PFT+1];                  /* daylength (s) */
   double prev_dayl[MAX_PFT+1];             /* daylength at previous timestep (s) */
   double annavg_t2m[MAX_PFT+1];            /* annual average 2-m air temperature (K) */
   double tempavg_t2m[MAX_PFT+1];           /* temporary average 2-m air temperature (K) */
   double gpp2[MAX_PFT+1];                   /* GPP flux before downregulation (g C/m^2/s) */
   double availc[MAX_PFT+1];                /* C flux available for allocation (g C/m^2/s) */
   double xsmrpool_recover[MAX_PFT+1];      /* C flux assigned to recovery (g C/m^2/s) */
   double alloc_pnow[MAX_PFT+1];            /* fraction of current allocation as new growth */
   double c_allometry[MAX_PFT+1];           /* C allocation index */
   double n_allometry[MAX_PFT+1];           /* N allocation index */
   double plant_ndemand[MAX_PFT+1];         /* N flux required to support GPP (g N/m^2/s) */
   double tempsum_potential_gpp[MAX_PFT+1]; /* temporary annual sum of potential GPP */
   double annsum_potential_gpp[MAX_PFT+1];  /* annuals sum of potential GPP */
   double tempmax_retransn[MAX_PFT+1];      /* temporary annual max of retrans N pool (g N/m^2) */
   double annmax_retransn[MAX_PFT+1];       /* annual max of retransloc N pool (g N/m^2) */
   double avail_retransn[MAX_PFT+1];        /* N flux avail for retransloc (g N/m^2/s) */
   double plant_nalloc[MAX_PFT+1];          /* total allocated N flux (g N/m^2/s) */
   double plant_calloc[MAX_PFT+1];          /* total allocated C flux (g C/m^2/s) */
   double excess_cflux[MAX_PFT+1];          /* C flux not allocated (g C/m^2/s) */
   double downreg[MAX_PFT+1];               /* fract reduction in GPP due to N limit */
   double prev_leafc_to_litter[MAX_PFT+1];  /* previous leaf C litterfall (g C/m^2/s) */
   double prev_frootc_to_litter[MAX_PFT+1]; /* previous froot C litterfall (g C/m^2/s) */
   double tempsum_npp[MAX_PFT+1];           /* temporary annual sum of NPP (g C/m^2/yr) */
   double annsum_npp[MAX_PFT+1];            /* annual sum of NPP (g C/m^2/yr) */
   double gpp[MAX_PFT+1];                   /* gross primary production (g C/m^2/s) */
   double npp[MAX_PFT+1];                   /* net primary production (g C/m^2/s) */
   double ar[MAX_PFT+1];                    /* autotrophic respiration (g C/m^2/s) */
   double gr[MAX_PFT+1];                    /* growth respiration (g C/m^2/s) */
   double mr[MAX_PFT+1];                    /* maintenance respiration (g C/m^2/s) */
   double leaf_mr[MAX_PFT+1];               /* leaf maintenance respiration (g C/m^2/s) */
   double litfall[MAX_PFT+1];               /* litterfall (g C/m^2/s) */
   double rs[MAX_PFT+1];                    /* stomatal resistance */
   double par[MAX_PFT+1];                   /* absorbed PAR */
   /* PFT-level carbon state */
   double leafc[MAX_PFT+1];                 /* leaf C (g C/m^2) */
   double leafc_storage[MAX_PFT+1];         /* leaf C storage (g C/m^2) */
   double leafc_xfer[MAX_PFT+1];            /* leaf C transfer (g C/m^2) */
   double frootc[MAX_PFT+1];                /* fine root C (g C/m^2) */
   double frootc_storage[MAX_PFT+1];        /* fine root C storage (g C/m^2) */
   double frootc_xfer[MAX_PFT+1];           /* fine root C transfer (g C/m^2) */
   double livestemc[MAX_PFT+1];             /* live stem C (g C/m^2) */
   double livestemc_storage[MAX_PFT+1];     /* live stem C storage (g C/m^2) */
   double livestemc_xfer[MAX_PFT+1];        /* live stem C transfer (g C/m^2) */
   double deadstemc[MAX_PFT+1];             /* dead stem C (g C/m^2) */
   double deadstemc_storage[MAX_PFT+1];     /* dead stem C storage (g C/m^2) */
   double deadstemc_xfer[MAX_PFT+1];        /* dead stem C transfer (g C/m^2) */
   double livecrootc[MAX_PFT+1];            /* live coarse root C (g C/m^2) */
   double livecrootc_storage[MAX_PFT+1];    /* live coarse root C storage (g C/m^2) */
   double livecrootc_xfer[MAX_PFT+1];       /* live coarse root C transfer (g C/m^2) */
   double deadcrootc[MAX_PFT+1];            /* dead coarse root C (g C/m^2) */
   double deadcrootc_storage[MAX_PFT+1];    /* dead coarse root C storage (g C/m^2) */
   double deadcrootc_xfer[MAX_PFT+1];       /* dead coarse root C transfer (g C/m^2) */
   double gresp_storage[MAX_PFT+1];         /* growth respiration storage (g C/m^2) */
   double gresp_xfer[MAX_PFT+1];            /* growth respiration transfer (g C/m^2) */
   double cpool[MAX_PFT+1];                 /* temporary photosynthate C pool (g C/m^2) */
   double xsmrpool[MAX_PFT+1];              /* abstract C pool to meet excess MR demand (g C/m^2) */
   double pft_ctrunc[MAX_PFT+1];            /* PFT-level sink for C truncation (g C/m^2) */
   double totvegc[MAX_PFT+1];               /* total vegetation C (g C/m^2) */
   double woodc[MAX_PFT+1];                 /* wood C (g C/m^2) */
   double fpsn[MAX_PFT+1];                  /* photosynthesis */
   double ci[MAX_PFT+1];                    /* intracellular CO2 */
   /* PFT-level nitrogen state */
   double leafn[MAX_PFT+1];                 /* leaf N (g N/m^2) */
   double leafn_storage[MAX_PFT+1];         /* leaf N storage (g N/m^2) */
   double leafn_xfer[MAX_PFT+1];            /* leaf N transfer (g N/m^2) */
   double frootn[MAX_PFT+1];                /* fine root N (g N/m^2) */
   double frootn_storage[MAX_PFT+1];        /* fine root N storage (g N/m^2) */
   double frootn_xfer[MAX_PFT+1];           /* fine root N transfer (g N/m^2) */
   double livestemn[MAX_PFT+1];             /* live stem N (g N/m^2) */
   double livestemn_storage[MAX_PFT+1];     /* live stem N storage (g N/m^2) */
   double livestemn_xfer[MAX_PFT+1];        /* live stem N transfer (g N/m^2) */
   double deadstemn[MAX_PFT+1];             /* dead stem N (g N/m^2) */
   double deadstemn_storage[MAX_PFT+1];     /* dead stem N storage (g N/m^2) */
   double deadstemn_xfer[MAX_PFT+1];        /* dead stem N transfer (g N/m^2) */
   double livecrootn[MAX_PFT+1];            /* live coarse root N (g N/m^2) */
   double livecrootn_storage[MAX_PFT+1];    /* live coarse root N storage (g N/m^2) */
   double livecrootn_xfer[MAX_PFT+1];       /* live coarse root N transfer (g N/m^2) */
   double deadcrootn[MAX_PFT+1];            /* dead coarse root N (g N/m^2) */
   double deadcrootn_storage[MAX_PFT+1];    /* dead coarse root N storage (g N/m^2) */
   double deadcrootn_xfer[MAX_PFT+1];       /* dead coarse root N transfer (g N/m^2) */
   double retransn[MAX_PFT+1];              /* retranslocated N (g N/m^2) */

   double npool[MAX_PFT+1];                 /* temporary photosynthate N pool (g N/m^2) */
   double pft_ntrunc[MAX_PFT+1];            /* PFT-level sink for N truncation (g N/m^2) */
   /* column (band) physical state */
   double decl;                      /* solar declination angle (radians) */
   double fpi;                       /* fraction of potential immobilization */
   double fpg;                       /* fraction of potential GPP */
   double annsum_counter;            /* seconds since last ann accumulation turnover */
   double cannsum_npp;               /* annual sum of NPP, averaged from PFT-level (g C/m^2/yr) */
   double cannavg_t2m;               /* annual avg. of 2-m air temperature, averaged from PFT-level (K) */
   double watfc[MAX_NODES];          /* volumetric soil water at field capacity */
   double me;                        /* moisture of extinction */
   double fire_prob;                 /* daily fire probability */
   double mean_fire_prob;            /* e-folding mean of daily fire prob. */
   double fireseasonl;               /* annual fire season length (days) */
   double farea_burned;              /* timestep fractional area burned */
   double ann_farea_burned;          /* annual total fract. area burned */
   double hr;                    /* heterotrophic respiration (g C/m^2/s) */
   double lithr;                 /* litter heterotrophic resp. (g C/m^2/s) */
   double nee;                   /* net ecosystem exchange (g C/m^2/s) */
   double nep;                   /* net ecosystem production (g C/m^2/s) */
   /* column (band) carbon state */
   double cwdc;                      /* coarse woody debris C (g C/m^2) */
   double litr1c;                    /* litter labile C (g C/m^2) */
   double litr2c;                    /* litter cellulose C (g C/m^2) */
   double litr3c;                    /* litter lignin C (g C/m^2) */
   double soil1c;                    /* fastest soil organic matter C */
   double soil2c;                    /* medium soil organic matter C */
   double soil3c;                    /* slow soil organic matter C */
   double soil4c;                    /* slowest soil organic matter C */
   double seedc;                     /* column-lev pool for seeding new PFTs */
   double col_ctrunc;                /* column-lev sink for C truncation */
   double totlitc;                   /* total litter C (g C/m^2) */
   double totsomc;                   /* total soil organic C (g C/m^2) */
   double totcolc;                   /* total column C (g C/m^2) */
   double prod10c;                   /* wood product C pool, 10-yr lifespan (g C/m^2) */
   double prod100c;                  /* wood product C pool, 100-yr lifespan (g C/m^2) */
   /* column (band) nitrogen state */
   double cwdn;                      /* coarse woody debris N (g N/m^2) */
   double litr1n;                    /* litter labile N (g N/m^2) */
   double litr2n;                    /* litter cellulose N (g N/m^2) */
   double litr3n;                    /* litter lignin N (g N/m^2) */
   double soil1n;                    /* fastest soil organic matter N */
   double soil2n;                    /* medium soil organic matter N */
   double soil3n;                    /* slow soil organic matter N */
   double soil4n;                    /* slowest soil organic matter N */
   double sminn;                     /* soil mineral N (g N/m^2) */
   double seedn;                     /* column-lev pool for seeding new PFTs */
   double col_ntrunc;                /* column-lev sink for N truncation */
   double totcoln;                   /* total column N (g N/m^2) */
   double prod10n;                   /* wood product N pool, 10-yr lifespan (g N/m^2) */
   double prod100n;                   /* wood product N pool, 100-yr lifespan (g N/m^2) */
} cn_data_struct;

cn_data_struct *alloc_cn(int, int);
void cov2pft(double, int, double, double []);
void get_CNatmos(double, atmos_data_struct *, cn_data_struct *);
void get_CNsoil(int, int, int, soil_con_struct *, all_vars_struct *, cn_data_struct *);
void get_CNsoilhydro(int, int, int, int, int, soil_con_struct *, veg_con_struct *, all_vars_struct *, cn_data_struct *);
void get_CNsoiltherm(int, int, int, all_vars_struct *, cn_data_struct *);
void get_CNveg(int, int, double, veg_con_struct *, veg_lib_struct *, all_vars_struct *, cn_data_struct *);
double pft2cov(double [], int);
void synveg_alloc(void);
void synveg_finalize(void);
void synveg_init(void);
void synveg_init_output(void);
void synveg_restore(void);
void synveg_run(void);
void synveg_start(void);
void synveg_write(void);
void VICCNInterface(int, int, int, int, int, int, int, int, dmy_struct *, global_param_struct *, atmos_data_struct *, all_vars_struct *, soil_con_struct *, veg_con_struct *, veg_lib_struct *, cn_data_struct *);

extern void clm_initialize2_(double *, int *, int *, int *, int *, int *, \
			     int *, int *, int [], int [], double [], \
			     int []);
extern void init_clmpftwt_(int *, int *, int *, int *, int *, double *);

extern void vic2clmtype(int *nlevgrnd, int *rec, int *nrec, int *adspinup, \
			int *init_state, int *yr, int *mo, int *day,	\
			 int *secs, double *jday, int *yrnxt, double *jdynxt, \
			 double *dt, double *lat, double *lon, \
			 int *begg, int *endg, int *begc,      \
			 int *endc, int *begp, int *endp, \
			int *num_soilc, int *num_soilp, cn_data_struct *cn);

#endif
