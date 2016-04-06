/******************************************************************************
 * @section DESCRIPTION
 *
 * Run function for image mode driver.
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
#include <vic_driver_image.h>
#include <synveg_cn.h>

// run vegetation

void synveg_run(void)
{
  extern option_struct       options;
  extern size_t              NR;
  extern size_t              current;
  extern all_vars_struct    *all_vars;
  extern atmos_data_struct  *atmos;
  extern dmy_struct         *dmy;
  extern domain_struct       local_domain;
  extern global_param_struct global_param;
  extern out_data_struct   **out_data;
  extern save_data_struct   *save_data;
  extern soil_con_struct    *soil_con;
  extern veg_con_struct    **veg_con;
  extern veg_hist_struct   **veg_hist;
  extern veg_lib_struct    **veg_lib;
  extern int                *Nveg;
  extern int                *Npfts;
  extern cn_data_struct     **cn;

  size_t                     i, iveg, band;
  int                        lbc, ubc, lbp, ubp;

  for (i = 0; i < local_domain.ncells; i++) {

    /* Convert LAI from global to local */
    if (current >= 0) {
        for (iveg = 0; iveg < Nveg[i]; iveg++) {
            for (band = 0; band < options.SNOW_BAND; band++) {
	      all_vars[i].veg_var[iveg][band].LAI = veg_hist[i][iveg].LAI[NR];
                all_vars[i].veg_var[iveg][band].LAI /= \
		  all_vars[i].veg_var[iveg][band].vegcover;
            }
        }
    }

    if(options.CARBON_MODEL == CN_NORMAL || options.CARBON_MODEL == CN_ADECOMP)
      {
	VICCNInterface(current, i, Nveg[i], Npfts[i], lbc, ubc, lbp, ubp, dmy, \
		   &global_param, &(atmos[i]), &(all_vars[i]), \
		   &(soil_con[i]), veg_con[i], veg_lib[i], &(cn[i]));
      }

    /* Convert LAI back to global */
    if (current >= 0) {
        for (iveg = 0; iveg < Nveg[i]; iveg++) {
            for (band = 0; band < options.SNOW_BAND; band++) {
	      veg_hist[i][iveg].LAI[NR] = all_vars[i].veg_var[iveg][band].LAI;
                veg_hist[i][iveg].LAI[NR] *= \
		  all_vars[i].veg_var[iveg][band].vegcover;
            }
        }
    }
  }

} 
