#include <stdio.h>
#include <stdlib.h>
#include <vic_def.h>
#include <vic_run.h>
#include <string.h>

void get_CNatmos(double             lat,
		 atmos_data_struct *atmos,
		 cn_data_struct    *cn)

/**********************************************************************
	VICCNInterface	Michael Brunke		January 23, 2015

Puts VIC atmospheric data into cn structure.

**********************************************************************/

{

  double factor;

  /* Assign atmos structure variables to local variables, MAB 8/19/13 */

  cn[0].Tair = atmos->air_temp[NR];
  cn[0].vp = atmos->vp[NR];
  cn[0].vpd = atmos->vpd[NR];
  cn[0].psfc = atmos->pressure[NR];
  cn[0].lwrad = atmos->longwave[NR];
  cn[0].swrad = atmos->shortwave[NR];
  cn[0].precip = atmos->prec[NR];

  /* Calculation of direct and diffuse SW radiation elements, MAB 11/26/13 */
  /* Based off of CLM's data atmosphere model.  THIS SHOULD NOT BE ACTIVE  */
  /* WHEN COUPLED TO WRF!  These should be provided by WRF directly if     */
  /* possible. */

  cn[0].swrd[0] = 0.28 * cn[0].swrad;
  cn[0].swrd[1] = 0.31 * cn[0].swrad;
  cn[0].swri[0] = 0.24 * cn[0].swrad;
  cn[0].swri[1] = 0.17 * cn[0].swrad;

  factor = 1.0;
  if(lat > -60.0 && lat < -50.0)
    factor = 1.0 - (lat + 60.0) * (0.05 / 10.0);
  else if(lat >= -50.0 && lat <= 30.0)
    factor = 0.95;
  else if(lat > 30.0 && lat < 40.0)
    factor = 1.0 - (40.0 - lat) * (0.05 / 10.0);

  cn[0].swrd[0] *= factor;
  cn[0].swri[0] *= factor;
  cn[0].swrd[1] *= factor;
  cn[0].swri[1] *= factor;

}
