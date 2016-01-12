#include <stdio.h>
#include <stdlib.h>
#include <vic_def.h>
#include <vic_run.h>
#include <string.h>

void get_CNsoil(int Nbands, int Nnodes, int Nveg,
		soil_con_struct  *soil_con,
		all_vars_struct  *all_vars,
		cn_data_struct   *cn)

/**********************************************************************
	get_CNsoil	Michael Brunke		January 23, 2015

Puts VIC soil characteristics into cn structure.

**********************************************************************/

{

  int band, nidx, iveg;
  double factor;
  snow_data_struct  **snow;

  snow = all_vars->snow;

  /* Layer properties, MAB 10/10/13 */
  /* Can this be done earlier since this is constant? */

  for(band = 0; band < Nbands; band++)
    {

      /* Snow depth */
      cn[band].snowdep = 0.0;

      for(iveg = 0; iveg < Nveg; iveg++)
	{
	      cn[band].snowdep += snow[iveg][band].depth;
	}

      cn[band].snowdep /= Nveg;

      if(cn[band].snowdep > 1.e-3)
	{
	  cn[band].dz[0] = 1.0e-3;
	  cn[band].dz[1] = cn[band].snowdep - 1.0e-3;

	  cn[band].z[1] = -0.5 * (cn[band].snowdep - 1.0e-3);
	  cn[band].z[0] = -5.0e-4 - cn[band].dz[1];
	}
      else
	{
	  cn[band].z[0] = 0.0;
	  cn[band].z[1] = 0.0;

	  cn[band].dz[0] = 0.0;
	  cn[band].dz[1] = 0.0;
	}

  /* Get soil/snow roughness, MAB 10/11/13 */
      cn[band].z0 = soil_con->rough;
      cn[band].z0s = soil_con->snow_rough;

      for(nidx = 0; nidx < Nnodes; nidx++)
	{
	  cn[band].z[nidx+2] = soil_con->Zsum_node[nidx];
	  cn[band].dz[nidx+2] = soil_con->dz_node[nidx];
	}

    }

}
