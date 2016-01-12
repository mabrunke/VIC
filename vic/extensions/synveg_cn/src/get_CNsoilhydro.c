#include <stdio.h>
#include <stdlib.h>
#include <vic_def.h>
#include <vic_run.h>
#include <string.h>

void get_CNsoilhydro(int Nbands, int Nnodes, int Nlayers, int Nveg, int Npft,
		     soil_con_struct  *soil_con,
		     veg_con_struct   *veg_con,
		     all_vars_struct  *all_vars,
		     cn_data_struct   *cn)

/**********************************************************************
	get_CNsoilhydro	Michael Brunke		January 25, 2015

Puts VIC soil hydrolog into cn structure.

**********************************************************************/

{

  int band, iveg, ipft, nidx, lidx, i;
  double Lsum;
  double rootf[MAX_PFT];
  double theta[MAX_LAYERS];
  double psand, pclay;
  double bsw[2][MAX_LAYERS];
  double watsat[MAX_LAYERS];
  double sucsat[2][MAX_LAYERS];
  double soisuc[2][MAX_LAYERS];
  double watfc[MAX_LAYERS];
  double ksat, fsat;
  
  char PAST_BOTTOM;

  cell_data_struct **cell;
  energy_bal_struct **energy;
  snow_data_struct  **snow;

  cell = all_vars->cell;
  energy = all_vars->energy;
  snow = all_vars->snow;

  /* Average baseflow and soil moisture over all bands for each 
     vegetated type */

  for(band = 0; band < Nbands; band++)
    {

      cn[band].baseflow = 0.0;
      for(nidx = 0; nidx < Nnodes + 2; nidx++)
	{
	  cn[band].moist[nidx] = 0.0;
          cn[band].ice[nidx] = 0.0;
	}

      for(iveg = 0; iveg <= Nveg; iveg++)
	{

	  if(iveg < Nveg)
	    {

		  cn[band].baseflow += cell[iveg][band].baseflow;

		  cn[band].moist[0] += snow[iveg][band].surf_water;
		  cn[band].moist[1] += snow[iveg][band].pack_water;
                  cn[band].ice[1] += snow[iveg][band].swq;
		  for(nidx = 0; nidx < Nnodes; nidx++)
		    {
		      cn[band].moist[nidx+2] += energy[iveg][band].moist[nidx];
		      cn[band].ice[nidx+2] += (energy[iveg][band].ice[nidx] * \
			snow[iveg][band].density);
		    }

	    }

	}

      cn[band].baseflow /= (Nveg - 1);

      for(nidx = 0; nidx < Nnodes + 2; nidx++)
	{
	  cn[band].moist[nidx] /= (Nveg - 1);
          cn[band].ice[nidx] /= (Nveg - 1);
	  if(nidx < 2)
	    {
	    cn[band].moist[nidx] *= CONST_RHOFW;
	    }
	  else
	    {
	    cn[band].moist[nidx] *= (soil_con->dz_node[nidx-2] * CONST_RHOFW);
	    cn[band].ice[nidx] *= (soil_con->dz_node[nidx-2] * CONST_RHOICE);
	    }
	}

    }

  /* Distribute root fraction across snow bands */
  for(band = 0; band < Nbands; band++)
    {

      Lsum = 0.0;
      lidx = 0;
      PAST_BOTTOM = false;

      for(nidx = 0; nidx < Nnodes; nidx++)
	{

	  for(ipft = 0; ipft < 21; ipft++)
	    {
	      rootf[ipft] = 0.0;
	    }

	  for(iveg = 0; iveg < Nveg; iveg++)
	    {

	      /* printf("%d %d %f %f\n", nidx, iveg, veg_con[iveg].root[lidx]); */
	      if(soil_con->Zsum_node[nidx] == Lsum + soil_con->depth[lidx] && \
		 nidx != 0 && lidx != Nlayers - 1)
		cov2pft((veg_con[iveg].root[lidx] +		\
			     veg_con[iveg].root[lidx+1]) / 2.0, \
			    veg_con[iveg].veg_class, veg_con[iveg].Cv, \
			    rootf);
	      else
		cov2pft(veg_con[iveg].root[lidx], veg_con[iveg].veg_class, \
			veg_con[iveg].Cv, rootf);

	      /* for(ipft = 0; ipft < 21; ipft++)
		{
		  printf("%d %f\n", ipft, rootf[ipft]);
		  } */
	    }

	  if(soil_con->Zsum_node[nidx] > Lsum + soil_con->depth[lidx] && \
		 !PAST_BOTTOM)
	    {
	      Lsum += soil_con->depth[lidx];
	      lidx++;
	      if(lidx == Nlayers)
		{
		  PAST_BOTTOM = true;
		  lidx = Nlayers - 1;
		}
	    }

	  for(ipft = 0; ipft < MAX_PFT; ipft++)
	    {
	      cn[band].rootfr[ipft][nidx] = rootf[ipft];
	      /* printf("%d %d %d %f %f\n", band, ipft, nidx, cn[band].rootfr[ipft][nidx], rootf[ipft]); */
	    }

	}

  /* Place saturated matric potential into array and calculate soil matric 
     potential for each layer (Need to convert from cm to MPa) */

      for(lidx = 0; lidx < Nlayers; lidx++)
	{
	  theta[lidx] = 0.0;
	}

      for(iveg = 0; iveg < Nveg; iveg++)
	{
	  for(lidx = 0; lidx < Nlayers; lidx++)
	    {
	      theta[lidx] += cell[iveg][band].layer[lidx].moist;
	    }
	}

      for(lidx = 0; lidx < Nlayers; lidx++)
	{
	  theta[lidx] /= (Nveg - 1);
	  theta[lidx] /= (soil_con->depth[lidx] * 1000.0);
	}

      for(lidx = 0; lidx < Nlayers; lidx++)
	{

	  psand = soil_con->quartz[lidx] * soil_con->bulk_dens_min[lidx] / \
	    soil_con->soil_dens_min[lidx] * 100.0;
          pclay = (1.0 - soil_con->quartz[lidx] - soil_con->porosity[lidx]) * \
	    soil_con->bulk_dens_min[lidx] / soil_con->soil_dens_min[lidx] * \
	    100.0;
	  if(pclay < 0.0)
	    pclay = 0.0;

	  bsw[0][lidx] = 2.91 + 0.159 * pclay;
          sucsat[0][lidx] = 10.0 * (10.0 * (1.88 - 0.0131 * psand));
	  watsat[lidx] = 0.489 - 0.00126 * psand;
          ksat = 0.0070556 * pow(10., -0.884 + 0.0153 * psand);
	  fsat = theta[band] / soil_con->porosity[lidx];
	  if(fsat < 0.001)
	    fsat = 0.001;
	  soisuc[0][lidx] = sucsat[0][lidx] * pow(fsat, bsw[0][lidx]);
	  if(soisuc[0][lidx] > 15.0)
	    soisuc[0][lidx] = 15.0;
	  if(soisuc[0][lidx] < 0.0)
	    soisuc[0][lidx] = 0.0;
	  watfc[lidx] = watsat[lidx] * pow(0.1 / ksat / 86400., 1.0 / \
					     (2.0 * bsw[0][lidx] + 3.0));

	  bsw[1][lidx] = -(3.10 + 0.157 * pclay - 0.003 * psand);
          sucsat[1][lidx] = -1.0 * (exp((1.54 - 0.0095 * psand + \
					       0.0063 *	(100.0	- psand \
							 - pclay )) *	\
				     log(10.0)) * 9.8e-5);
	  watsat[lidx] = (50.5 - 0.142 * psand - 0.037 * pclay) / 100.0;
	  fsat = theta[lidx] / soil_con->porosity[lidx];
	  if(fsat < 0.001)
	    fsat = 0.001;
	  soisuc[1][lidx] = sucsat[1][lidx] * pow(fsat, bsw[1][lidx]);
	  if(soisuc[1][lidx] < -15.0)
	    soisuc[1][lidx] = -15.0;
	  if(soisuc[1][lidx] > 0.0)
	    soisuc[1][lidx] = 0.0;

	}

      Lsum = 0.0;
      lidx = 0;
      PAST_BOTTOM = false;

      for(nidx = 0; nidx < Nnodes; nidx++)
	{

	  for(band = 0; band < Nbands; band++)
	    {

	      for(i = 0; i < 2; i++)

		{

		  if(soil_con->Zsum_node[nidx] == Lsum + \
		     soil_con->depth[lidx] && nidx != 0 && lidx != Nlayers - 1)
		    {
		      cn[band].bsw[nidx][i] = (bsw[i][lidx] + bsw[i][lidx+1]) \
			/ 2.0;
		      cn[band].sucsat[nidx][i] = (sucsat[i][lidx] +	\
					   sucsat[i][lidx+1]) / 2.0;
		      cn[band].soisuc[nidx][i] = (soisuc[i][lidx] +	\
					   soisuc[i][lidx+1]) / 2.0;
		    }
		  else
		    {
		      cn[band].bsw[nidx][i] = bsw[i][lidx];
		      cn[band].sucsat[nidx][i] = sucsat[i][lidx];
		      cn[band].soisuc[nidx][i] = soisuc[i][lidx];
		    }

		}

	      if(soil_con->Zsum_node[nidx] == Lsum + soil_con->depth[lidx] \
		 && nidx != 0 && lidx != Nlayers - 1)
		{
		  cn[band].watfc[nidx] = (watfc[lidx] + watfc[lidx+1]) / 2.0;
		}
	      else
		{
		  cn[band].watfc[nidx] = watfc[lidx];
		}

	      if(soil_con->Zsum_node[nidx] > Lsum + soil_con->depth[lidx] && \
		 !PAST_BOTTOM)
		{
		  Lsum += soil_con->depth[lidx];
		  lidx++;
		  if(lidx == Nlayers)
		    {
		      PAST_BOTTOM = true;
		      lidx = Nlayers - 1;
		    }
		}

	    }

	}
    }

}
