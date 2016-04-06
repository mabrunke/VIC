#include <vic_def.h>
#include <vic_run.h>
#include <synveg_cn.h>

void get_CNveg(int Nbands, int Nveg, double precip,
	       veg_con_struct   *veg_con,
	       veg_lib_struct   *veg_lib,
	       all_vars_struct  *all_vars,
	       cn_data_struct   *cn)

/**********************************************************************
	get_CNveg	Michael Brunke		January 26, 2015

Puts VIC vegetation characteristics into cn structure.

**********************************************************************/

{

  int band, iveg;
  double fsum[MAX_PFT+1];

  cell_data_struct **cell;
  energy_bal_struct **energy;
  veg_var_struct **veg;

  cell = all_vars->cell;
  energy = all_vars->energy;
  veg = all_vars->veg_var;

  /* Average canopy temperature and resistance */

  for(band = 0; band < Nbands; band++)
    {

      for(iveg = 0; iveg <= Nveg; iveg++)
	{
	  cn[band].t2m[iveg] = cn[0].Tair;
	  cn[band].Tveg[iveg] = 0.0;
	  cn[band].rveg[iveg] = 0.0;
          cn[band].zov[iveg] = 0.0;
          cn[band].displ[iveg] = 0.0;
	  cn[band].fpsn[iveg] = 0.0;
	  fsum[iveg] = 0.0;
	}

      for(iveg = 0; iveg < Nveg; iveg++)
	{

	  cn[band].Tveg[iveg] = energy[iveg][band].Tcanopy;
	  cn[band].rveg[iveg] = cell[iveg][band].aero_resist[1];
	  cn[band].zov[iveg] = veg_lib->roughness[veg_con[iveg].veg_class];
	  cn[band].displ[iveg] = veg_lib->displacement[veg_con[iveg].veg_class];
	  cn[band].fpsn[iveg] = veg[iveg][band].GPP * 1.e6 / \
	    veg[iveg][band].LAI;

            if(precip > 0.0)
	      cn[band].fwet[iveg] = 1.0;
	    else
	      cn[band].fwet[iveg] = 0.0;

	}

    }

}
