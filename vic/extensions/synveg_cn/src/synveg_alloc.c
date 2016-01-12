/******************************************************************************
 * @section DESCRIPTION
 *
 * Allocate memory for VIC structures.
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

/******************************************************************************
 * @brief    Allocate memory for VIC structures.
 *****************************************************************************/
void alloc_cn(void)
{
    extern domain_struct       local_domain;
    extern option_struct       options;
    extern cn_data_struct    **cn;
    size_t                     i;
    size_t                     j;

    /* cn allocation, MAB 8/27/15 */
    cn = (cn_data_struct **) malloc((size_t) local_domain.ncells * \
				    sizeof(cn_data_struct *));
    if (cn == NULL) {
      log_err("Memory allocation error in vic_alloc().");
    }

    // allocate memory for individual grid cells
    for (i = 0; i < local_domain.ncells; i++) {
	/* allocate memory for cn, MAB 8/27/15 */
	alloc_cn(options.SNOW_BAND, options.Nnode, &(cn[i]));

    }
}
