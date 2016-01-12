/******************************************************************************
 * @section DESCRIPTION
 *
 * Initilization library.
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

#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_cesm.h>

/******************************************************************************
 * @brief    Initialize soil con sructure.
 *****************************************************************************/
void
initialize_soil_con(soil_con_struct *soil_con)
{
    extern option_struct options;
    size_t               i;
    size_t               j;


    soil_con->FS_ACTIVE = 0;
    soil_con->gridcel = -1;

    soil_con->AlbedoPar = 0.;
    soil_con->elevation = 0.;
    soil_con->lat = 0.;
    soil_con->lng = 0.;
    soil_con->time_zone_lng = 0.;

    soil_con->annual_prec = 0.;
    soil_con->aspect = 0.;
    soil_con->avg_temp = 0.;
    soil_con->avgJulyAirTemp = 0.;
    soil_con->b_infilt = 0.;
    soil_con->c = 0.;
    soil_con->cell_area = 0.;
    soil_con->dp = 0.;
    soil_con->Ds = 0.;
    soil_con->Dsmax = 0.;
    soil_con->ehoriz = 0.;
    soil_con->frost_slope = 0.;
    soil_con->max_infil = 0.;
    soil_con->max_snow_distrib_slope = 0.;
    soil_con->rough = 0.;
    soil_con->slope = 0.;
    soil_con->snow_rough = 0.;
    soil_con->whoriz = 0.;
    soil_con->Ws = 0.;

    for (i = 0; i < MAX_LAYERS; i++) {
        soil_con->bubble[i] = 0.;
        soil_con->bulk_density[i] = 0.;
        soil_con->bulk_dens_min[i] = 0.;
        soil_con->bulk_dens_org[i] = 0.;
        soil_con->depth[i] = 0.;
        soil_con->expt[i] = 0.;
        soil_con->init_moist[i] = 0.;
        soil_con->Ksat[i] = 0.;
        soil_con->max_moist[i] = 0.;
        soil_con->phi_s[i] = 0.;
        soil_con->porosity[i] = 0.;
        soil_con->quartz[i] = 0.;
        soil_con->organic[i] = 0.;
        soil_con->resid_moist[i] = 0.;
        soil_con->soil_density[i] = 0.;
        soil_con->soil_dens_min[i] = 0.;
        soil_con->soil_dens_org[i] = 0.;
        soil_con->Wcr[i] = 0.;
        soil_con->Wpwp[i] = 0.;
    }

    for (i = 0; i < MAX_NODES; i++) {
        soil_con->alpha[i] = 0.;
        soil_con->beta[i] = 0.;
        soil_con->bubble_node[i] = 0.;
        soil_con->dz_node[i] = 0.;
        soil_con->Zsum_node[i] = 0.;
        soil_con->expt_node[i] = 0.;
        soil_con->gamma[i] = 0.;
        soil_con->max_moist_node[i] = 0.;
    }

    for (i = 0; i < MAX_FROST_AREAS; i++) {
        soil_con->frost_fract[i] = 0.;
    }

    for (i = 0; i < options.SNOW_BAND; i++) {
        soil_con->AboveTreeLine[i] = 0;
        soil_con->BandElev[i] = 0.;
        soil_con->AreaFract[i] = 1.;
        soil_con->Pfactor[i] = 0.;
        soil_con->Tfactor[i] = 0.;
    }

    for (i = 0; i < MAX_LAYERS + 2; i++) {
        for (j = 0; j < MAX_ZWTVMOIST; j++) {
            soil_con->zwtvmoist_zwt[i][j] = 0.;
            soil_con->zwtvmoist_moist[i][j] = 0.;
        }
    }
}

/******************************************************************************
 * @brief    Initialize veg con sructure.
 *****************************************************************************/
void
initialize_veg_con(veg_con_struct *veg_con)
{
    extern option_struct options;
    size_t               i;

    veg_con->Cv = 0.;
    veg_con->Cv_sum = 0.;
    veg_con->veg_class = -1; // -1 to force a crash if inappropriate
    veg_con->vegetat_type_num = 0.;
    veg_con->sigma_slope = 0.;
    veg_con->lag_one = 0.;
    veg_con->fetch = 0.;
    veg_con->LAKE = 0;
    for (i = 0; i < MAX_LAYERS; i++) {
        veg_con->root[i] = 0.;
    }
    for (i = 0; i < options.ROOT_ZONES; i++) {
        veg_con->zone_depth[i] = 0.;
        veg_con->zone_fract[i] = 0.;
    }
    if (options.CARBON) {
        for (i = 0; i < options.Ncanopy; i++) {
            veg_con->CanopLayerBnd[i] = 0.;
        }
    }
}

/******************************************************************************
 * @brief    Initialize x2l_data_struct.
 *****************************************************************************/
void
initialize_x2l_data()
{
    extern x2l_data_struct *x2l_vic;
    extern domain_struct    local_domain;

    size_t                  i;

    log_info("Setting all x2l fields to %f", SHR_CONST_SPVAL);

    for (i = 0; i < local_domain.ncells; i++) {
        x2l_vic[i].x2l_Sa_z = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Sa_u = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Sa_v = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Sa_ptem = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Sa_shum = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Sa_pbot = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Sa_tbot = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Faxa_lwdn = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Faxa_rainc = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Faxa_rainl = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Faxa_snowc = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Faxa_snowl = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Faxa_swndr = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Faxa_swvdr = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Faxa_swndf = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Faxa_swvdf = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Sa_co2prog = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Sa_co2diag = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Faxa_bcphidry = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Faxa_bcphodry = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Faxa_bcphiwet = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Faxa_ocphidry = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Faxa_ocphodry = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Faxa_ocphiwet = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Faxa_dstwet1 = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Faxa_dstwet2 = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Faxa_dstwet3 = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Faxa_dstwet4 = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Faxa_dstdry1 = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Faxa_dstdry2 = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Faxa_dstdry3 = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Faxa_dstdry4 = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Flrr_flood = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_vars_set = false;
    }
}

/******************************************************************************
 * @brief    Initialize l2x_data_struct.
 *****************************************************************************/
void
initialize_l2x_data()
{
    extern l2x_data_struct *l2x_vic;
    extern domain_struct    local_domain;

    size_t                  i;

    log_info("Setting all l2x fields to %f", SHR_CONST_SPVAL);

    for (i = 0; i < local_domain.ncells; i++) {
        l2x_vic[i].l2x_Sl_t = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_Sl_tref = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_Sl_qref = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_Sl_avsdr = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_Sl_anidr = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_Sl_avsdf = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_Sl_anidf = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_Sl_snowh = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_Sl_u10 = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_Sl_ddvel = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_Sl_fv = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_Sl_ram1 = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_Sl_logz0 = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_Fall_taux = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_Fall_tauy = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_Fall_lat = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_Fall_sen = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_Fall_lwup = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_Fall_evap = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_Fall_swnet = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_Fall_fco2_lnd = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_Fall_flxdst1 = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_Fall_flxdst2 = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_Fall_flxdst3 = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_Fall_flxdst4 = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_Fall_flxvoc = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_Flrl_rofliq = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_Flrl_rofice = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_vars_set = true;
    }
}
