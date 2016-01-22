/******************************************************************************
 * @section DESCRIPTION
 *
 * Stand-alone image mode driver of the VIC model
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
#include <synveg_cn.h>

// initialize vegetation model
void synveg_init(void)
{
  extern domain_struct       local_domain;
  extern option_struct       options;
  extern int                *Nveg;
  extern int                *Npfts;
  extern int                 begg, endg, begc, endc, begp, endp;

  size_t                     j;
  int                        veg_class;
  int                       *vegclass;
  double                    *vegfract;
  double                     dt;

    /* Initialize CN state, MAB 8/29/15 */
  vegclass = calloc(local_domain.ncells * (MAX_VEG + 1), sizeof(int));
  vegfract = calloc(local_domain.ncells * (MAX_VEG + 1), sizeof(double));
  Nveg = calloc(local_domain.ncells, sizeof(int));
  Npfts = calloc(local_domain.ncells, sizeof(int));

  begg = 1;
  endg = local_domain.ncells;
  begc = 1;
  endc = options.SNOW_BAND;
  begp = 1;
  nlevgrnd = options.Nnode;
  for(i = 0; i < local_domain.ncells; i++) {
    endp += options.SNOW_BAND * veg_con[i][0].vegetat_type_num;
    Nveg[i] = veg_con[i][0].vegetat_type_num;
    Npfts[i] = 0;

    for(j = 0; j <= veg_con[i][0].vegetat_type_num; j++)
      {
	vegclass[i * (MAX_VEG + 1) + j] = veg_con[i][j].veg_class;
	vegfract[i * (MAX_VEG + 1) + j] = veg_con[i][j].Cv;
      }

  }

  dt = global_param.dt * 3600.0;

  clm_initialize2_(&dt, &nlevgrnd, &begg, &endg, &begc,			\
			 &endc, &begp, &endp, Nveg, vegclass,		\
			 vegfract, Npfts);
}
