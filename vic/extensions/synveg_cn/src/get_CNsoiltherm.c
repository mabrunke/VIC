#include <stdio.h>
#include <stdlib.h>
#include <vic_def.h>
#include <vic_run.h>
#include <string.h>

void get_CNsoiltherm(int Nbands, int Nnodes, int Nveg,
		     all_vars_struct  *all_vars,
		     cn_data_struct   *cn)

{

  int band, iveg, nidx;

  energy_bal_struct **energy;
  snow_data_struct  **snow;

  energy = all_vars->energy;
  snow = all_vars->snow;

  /* Combine snow and soil temperatures into one array */

  for(band = 0; band < Nbands; band++)
    {

      cn[band].alb = 0.0;
      for(nidx = 0; nidx < Nnodes; nidx++)
	{
	  cn[band].t_soisno[nidx] = 0.0;
	}

      for(iveg = 0; iveg < Nveg; iveg++)
	{
	      /* Add snow and pack temperatures */

	      cn[band].t_soisno[0] += snow[iveg][band].surf_temp;
	      cn[band].t_soisno[1] += snow[iveg][band].pack_temp;

              /* Add node temperatures */
	      for(nidx = 0; nidx < Nnodes; nidx++)
		{
		  cn[band].t_soisno[nidx + 2] += energy[iveg][band].T[nidx];
		}

	      cn[band].alb += energy[iveg][band].AlbedoUnder;

	}

      /* Average temperatures */

      for(nidx = 0; nidx < Nnodes + 2; nidx++)
	{
	  cn[band].t_soisno[nidx] /= Nveg;
	}

      /* Average albedo */
      cn[band].alb /= Nveg;

    }

}

