#include <stdio.h>
#include <stdlib.h>
#include <vic_def.h>
#include <vic_run.h>
#include <string.h>

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
  double fsum[MAX_CN_PFTS];

  cell_data_struct **cell;
  energy_bal_struct **energy;
  veg_var_struct **veg;

  cell = all_vars->cell;
  energy = all_vars->energy;
  veg = all_vars->veg_var;

  /* Average canopy temperature and resistance */

  for(band = 0; band < Nbands; band++)
    {

      for(iveg = 0; iveg <= MAX_CN_PFTS; iveg++)
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

	  switch(veg_con[iveg].veg_class)
	    {
	    case 0: cn[band].Tveg[2] += (energy[iveg][band].Tcanopy * veg_con[iveg].Cv);
	      cn[band].rveg[2] += (cell[iveg][band].aero_resist[1] * veg_con[iveg].Cv);
              cn[band].zov[2] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv);
              cn[band].displ[2] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] * \
		 veg_con[iveg].Cv);
	      cn[band].fpsn[2] += (veg[iveg][band].GPP * 1.e6 * \
				      veg_con[iveg].Cv / veg[iveg][band].LAI);
	      fsum[2] += veg_con[iveg].Cv;
	      break;
	    case 1: cn[band].Tveg[5] += (energy[iveg][band].Tcanopy * veg_con[iveg].Cv);
	      cn[band].rveg[5] += (cell[iveg][band].aero_resist[1] * veg_con[iveg].Cv);
              cn[band].zov[5] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv);
              cn[band].displ[5] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv);
	      cn[band].fpsn[5] += (veg[iveg][band].GPP * 1.e6 * \
				      veg_con[iveg].Cv / veg[iveg][band].LAI);
	      fsum[5] += veg_con[iveg].Cv;
	      break;
	    case 2: cn[band].Tveg[3] += (energy[iveg][band].Tcanopy * veg_con[iveg].Cv);
	      cn[band].rveg[3] += (cell[iveg][band].aero_resist[1] * veg_con[iveg].Cv);
              cn[band].zov[3] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv);
              cn[band].displ[3] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv);
	      cn[band].fpsn[3] += (veg[iveg][band].GPP * 1.e6 * \
				      veg_con[iveg].Cv / veg[iveg][band].LAI);
	      fsum[3] += veg_con[iveg].Cv;
	      break;
	    case 3: cn[band].Tveg[8] += (energy[iveg][band].Tcanopy * veg_con[iveg].Cv);
	      cn[band].rveg[8] += (cell[iveg][band].aero_resist[1] * veg_con[iveg].Cv);
              cn[band].zov[8] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv);
              cn[band].displ[8] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv);
	      cn[band].fpsn[8] += (veg[iveg][band].GPP * 1.e6 * \
				      veg_con[iveg].Cv / veg[iveg][band].LAI);
	      fsum[8] += veg_con[iveg].Cv;
	      break;
	    case 4: cn[band].Tveg[2] += (energy[iveg][band].Tcanopy * veg_con[iveg].Cv * 0.5);
	      cn[band].Tveg[8] += (energy[iveg][band].Tcanopy * veg_con[iveg].Cv * 0.5);
	      cn[band].rveg[2] += (cell[iveg][band].aero_resist[1] * veg_con[iveg].Cv * 0.5);
	      cn[band].rveg[8] += (cell[iveg][band].aero_resist[1] * veg_con[iveg].Cv * 0.5);
              cn[band].zov[8] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv * 0.5);
              cn[band].zov[2] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv * 0.5);
              cn[band].displ[8] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv * 0.5);
              cn[band].displ[2] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv * 0.5);
	      cn[band].fpsn[2] += (veg[iveg][band].GPP * 1.e6 * \
				      veg_con[iveg].Cv * 0.5 / veg[iveg][band].LAI);
	      cn[band].fpsn[8] += (veg[iveg][band].GPP * 1.e6 * \
				      veg_con[iveg].Cv * 0.5 / veg[iveg][band].LAI);
	      fsum[2] += veg_con[iveg].Cv * 0.5;
	      fsum[8] += veg_con[iveg].Cv * 0.5;
	      break;
	    case 5: cn[band].Tveg[2] += (energy[iveg][band].Tcanopy * veg_con[iveg].Cv * 0.4);
	      cn[band].Tveg[8] += (energy[iveg][band].Tcanopy * veg_con[iveg].Cv * 0.4);
	      cn[band].Tveg[9] += (energy[iveg][band].Tcanopy * veg_con[iveg].Cv * 0.1);
	      cn[band].Tveg[10] += (energy[iveg][band].Tcanopy * veg_con[iveg].Cv * 0.1);
	      cn[band].rveg[2] += (cell[iveg][band].aero_resist[1] * veg_con[iveg].Cv * 0.4);
	      cn[band].rveg[8] += (cell[iveg][band].aero_resist[1] * veg_con[iveg].Cv * 0.4);
	      cn[band].rveg[9] += (cell[iveg][band].aero_resist[1] * veg_con[iveg].Cv * 0.1);
	      cn[band].rveg[10] += (cell[iveg][band].aero_resist[1] * veg_con[iveg].Cv * 0.1);
              cn[band].zov[2] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv * 0.4);
              cn[band].zov[8] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv * 0.4);
              cn[band].zov[9] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv * 0.1);
              cn[band].zov[10] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv * 0.1);
              cn[band].displ[2] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv * 0.4);
              cn[band].displ[8] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv * 0.4);
              cn[band].displ[9] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv * 0.1);
              cn[band].displ[10] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv * 0.1);
	      cn[band].fpsn[2] += (veg[iveg][band].GPP * 1.e6 * \
				      veg_con[iveg].Cv * 0.4 / veg[iveg][band].LAI);
	      cn[band].fpsn[8] += (veg[iveg][band].GPP * 1.e6 * \
				      veg_con[iveg].Cv * 0.4 / veg[iveg][band].LAI);
	      cn[band].fpsn[9] += (veg[iveg][band].GPP * 1.e6 * \
				      veg_con[iveg].Cv * 0.1 / veg[iveg][band].LAI);
	      cn[band].fpsn[10] += (veg[iveg][band].GPP * 1.e6 * \
				      veg_con[iveg].Cv * 0.1 / veg[iveg][band].LAI);
	      fsum[2] += veg_con[iveg].Cv * 0.4;
	      fsum[8] += veg_con[iveg].Cv * 0.4;
	      fsum[9] += veg_con[iveg].Cv * 0.1;
	      fsum[10] += veg_con[iveg].Cv * 0.1;
	      break;
	    case 6: cn[band].Tveg[14] += (energy[iveg][band].Tcanopy * veg_con[iveg].Cv * 0.1);
	      cn[band].Tveg[13] += (energy[iveg][band].Tcanopy * veg_con[iveg].Cv * 0.3);
	      cn[band].Tveg[12] += (energy[iveg][band].Tcanopy * veg_con[iveg].Cv * 0.2);
	      cn[band].Tveg[2] += (energy[iveg][band].Tcanopy * veg_con[iveg].Cv * 0.2);
	      cn[band].Tveg[8] += (energy[iveg][band].Tcanopy * veg_con[iveg].Cv * 0.2);
	      cn[band].rveg[14] += (cell[iveg][band].aero_resist[1] * veg_con[iveg].Cv * 0.1);
	      cn[band].rveg[13] += (cell[iveg][band].aero_resist[1] * veg_con[iveg].Cv * 0.3);
	      cn[band].rveg[12] += (cell[iveg][band].aero_resist[1] * veg_con[iveg].Cv * 0.2);
	      cn[band].rveg[2] += (cell[iveg][band].aero_resist[1] * veg_con[iveg].Cv * 0.2);
	      cn[band].rveg[8] += (cell[iveg][band].aero_resist[1] * veg_con[iveg].Cv * 0.2);
              cn[band].zov[14] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv * 0.1);
              cn[band].zov[13] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv * 0.3);
              cn[band].zov[12] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv * 0.2);
              cn[band].zov[2] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv * 0.2);
              cn[band].zov[8] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv * 0.2);
              cn[band].displ[14] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv * 0.1);
              cn[band].displ[13] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv * 0.3);
              cn[band].displ[12] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv * 0.2);
              cn[band].displ[2] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv * 0.2);
              cn[band].displ[8] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv * 0.2);
	      cn[band].fpsn[14] += (veg[iveg][band].GPP * 1.e6 * \
				      veg_con[iveg].Cv * 0.1 / veg[iveg][band].LAI);
	      cn[band].fpsn[13] += (veg[iveg][band].GPP * 1.e6 * \
				      veg_con[iveg].Cv * 0.3 / veg[iveg][band].LAI);
	      cn[band].fpsn[12] += (veg[iveg][band].GPP * 1.e6 * \
				      veg_con[iveg].Cv * 0.2 / veg[iveg][band].LAI);
	      cn[band].fpsn[2] += (veg[iveg][band].GPP * 1.e6 * \
				      veg_con[iveg].Cv * 0.2 / veg[iveg][band].LAI);
	      cn[band].fpsn[8] += (veg[iveg][band].GPP * 1.e6 * \
				      veg_con[iveg].Cv * 0.2 / veg[iveg][band].LAI);
	      fsum[14] += veg_con[iveg].Cv * 0.1;
	      fsum[13] += veg_con[iveg].Cv * 0.3;
	      fsum[12] += veg_con[iveg].Cv * 0.2;
	      fsum[2] += veg_con[iveg].Cv * 0.2;
	      fsum[8] += veg_con[iveg].Cv * 0.2;
	      break;
	    case 7: cn[band].Tveg[9] += (energy[iveg][band].Tcanopy * veg_con[iveg].Cv * 0.3);
	      cn[band].Tveg[10] += (energy[iveg][band].Tcanopy * veg_con[iveg].Cv * 0.3);
	      cn[band].Tveg[11] += (energy[iveg][band].Tcanopy * veg_con[iveg].Cv * 0.2);
	      cn[band].Tveg[2] += (energy[iveg][band].Tcanopy * veg_con[iveg].Cv * 0.1);
	      cn[band].Tveg[8] += (energy[iveg][band].Tcanopy * veg_con[iveg].Cv * 0.1);
	      cn[band].rveg[9] += (cell[iveg][band].aero_resist[1] * veg_con[iveg].Cv * 0.3);
	      cn[band].rveg[10] += (cell[iveg][band].aero_resist[1] * veg_con[iveg].Cv * 0.3);
	      cn[band].rveg[11] += (cell[iveg][band].aero_resist[1] * veg_con[iveg].Cv * 0.2);
	      cn[band].rveg[2] += (cell[iveg][band].aero_resist[1] * veg_con[iveg].Cv * 0.1);
	      cn[band].rveg[8] += (cell[iveg][band].aero_resist[1] * veg_con[iveg].Cv * 0.1);
              cn[band].zov[9] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv * 0.3);
              cn[band].zov[10] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv * 0.3);
              cn[band].zov[11] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv * 0.2);
              cn[band].zov[8] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv * 0.1);
              cn[band].zov[2] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv * 0.1);
              cn[band].displ[9] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv * 0.3);
              cn[band].displ[10] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv * 0.3);
              cn[band].displ[11] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv * 0.2);
              cn[band].displ[8] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv * 0.1);
              cn[band].displ[2] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv * 0.1);
	      cn[band].fpsn[9] += (veg[iveg][band].GPP * 1.e6 * \
				      veg_con[iveg].Cv * 0.3 / veg[iveg][band].LAI);
	      cn[band].fpsn[10] += (veg[iveg][band].GPP * 1.e6 * \
				      veg_con[iveg].Cv * 0.3 / veg[iveg][band].LAI);
	      cn[band].fpsn[11] += (veg[iveg][band].GPP * 1.e6 * \
				      veg_con[iveg].Cv * 0.2 / veg[iveg][band].LAI);
	      cn[band].fpsn[8] += (veg[iveg][band].GPP * 1.e6 * \
				      veg_con[iveg].Cv * 0.1 / veg[iveg][band].LAI);
	      cn[band].fpsn[2] += (veg[iveg][band].GPP * 1.e6 * \
				      veg_con[iveg].Cv * 0.1 / veg[iveg][band].LAI);
	      fsum[9] += veg_con[iveg].Cv * 0.3;
	      fsum[10] += veg_con[iveg].Cv * 0.3;
	      fsum[11] += veg_con[iveg].Cv * 0.2;
	      fsum[8] += veg_con[iveg].Cv * 0.1;
	      fsum[2] += veg_con[iveg].Cv * 0.1;
	      break;
	    case 8: cn[band].Tveg[9] += (energy[iveg][band].Tcanopy * veg_con[iveg].Cv * 0.2);
	      cn[band].Tveg[10] += (energy[iveg][band].Tcanopy * veg_con[iveg].Cv * 0.2);
	      cn[band].Tveg[11] += (energy[iveg][band].Tcanopy * veg_con[iveg].Cv * 0.2);
	      cn[band].Tveg[14] += (energy[iveg][band].Tcanopy * veg_con[iveg].Cv * 0.05);
	      cn[band].Tveg[13] += (energy[iveg][band].Tcanopy * veg_con[iveg].Cv * 0.15);
	      cn[band].Tveg[12] += (energy[iveg][band].Tcanopy * veg_con[iveg].Cv * 0.1);
	      cn[band].rveg[9] += (cell[iveg][band].aero_resist[1] * veg_con[iveg].Cv * 0.2);
	      cn[band].rveg[10] += (cell[iveg][band].aero_resist[1] * veg_con[iveg].Cv * 0.2);
	      cn[band].rveg[11] += (cell[iveg][band].aero_resist[1] * veg_con[iveg].Cv * 0.2);
	      cn[band].rveg[14] += (cell[iveg][band].aero_resist[1] * veg_con[iveg].Cv * 0.05);
	      cn[band].rveg[13] += (cell[iveg][band].aero_resist[1] * veg_con[iveg].Cv * 0.15);
	      cn[band].rveg[12] += (cell[iveg][band].aero_resist[1] * veg_con[iveg].Cv * 0.1);
              cn[band].zov[9] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv * 0.2);
              cn[band].zov[10] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv * 0.2);
              cn[band].zov[11] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv * 0.2);
              cn[band].zov[14] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv * 0.05);
              cn[band].zov[13] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv * 0.15);
              cn[band].zov[12] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv * 0.1);
              cn[band].displ[9] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv * 0.2);
              cn[band].displ[10] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv * 0.2);
              cn[band].displ[11] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv * 0.2);
              cn[band].displ[14] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv * 0.05);
              cn[band].displ[13] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv * 0.15);
              cn[band].displ[12] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv * 0.1);
	      cn[band].fpsn[9] += (veg[iveg][band].GPP * 1.e6 * \
				      veg_con[iveg].Cv * 0.2 / veg[iveg][band].LAI);
	      cn[band].fpsn[10] += (veg[iveg][band].GPP * 1.e6 * \
				      veg_con[iveg].Cv * 0.2 / veg[iveg][band].LAI);
	      cn[band].fpsn[11] += (veg[iveg][band].GPP * 1.e6 * \
				      veg_con[iveg].Cv * 0.2 / veg[iveg][band].LAI);
	      cn[band].fpsn[14] += (veg[iveg][band].GPP * 1.e6 * \
				      veg_con[iveg].Cv * 0.05 / veg[iveg][band].LAI);
	      cn[band].fpsn[13] += (veg[iveg][band].GPP * 1.e6 * \
				      veg_con[iveg].Cv * 0.15 / veg[iveg][band].LAI);
	      cn[band].fpsn[12] += (veg[iveg][band].GPP * 1.e6 * \
				      veg_con[iveg].Cv * 0.1 / veg[iveg][band].LAI);
	      fsum[9] += veg_con[iveg].Cv * 0.2;
	      fsum[10] += veg_con[iveg].Cv * 0.2;
	      fsum[11] += veg_con[iveg].Cv * 0.2;
	      fsum[14] += veg_con[iveg].Cv * 0.05;
	      fsum[13] += veg_con[iveg].Cv * 0.15;
	      fsum[12] += veg_con[iveg].Cv * 0.1;
	      break;
	    case 9: cn[band].Tveg[14] += (energy[iveg][band].Tcanopy * veg_con[iveg].Cv * 0.1);
	      cn[band].Tveg[13] += (energy[iveg][band].Tcanopy * veg_con[iveg].Cv * 0.3);
	      cn[band].Tveg[12] += (energy[iveg][band].Tcanopy * veg_con[iveg].Cv * 0.2);
	      cn[band].rveg[14] += (cell[iveg][band].aero_resist[1] * veg_con[iveg].Cv * 0.1);
	      cn[band].rveg[13] += (cell[iveg][band].aero_resist[1] * veg_con[iveg].Cv * 0.3);
	      cn[band].rveg[12] += (cell[iveg][band].aero_resist[1] * veg_con[iveg].Cv * 0.2);
              cn[band].zov[14] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv * 0.1);
              cn[band].zov[13] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv * 0.3);
              cn[band].zov[12] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv * 0.2);
              cn[band].displ[14] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv * 0.1);
              cn[band].displ[13] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv * 0.3);
              cn[band].displ[12] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv * 0.2);
	      cn[band].fpsn[14] += (veg[iveg][band].GPP * 1.e6 * \
				      veg_con[iveg].Cv * 0.1 / veg[iveg][band].LAI);
	      cn[band].fpsn[13] += (veg[iveg][band].GPP * 1.e6 * \
				      veg_con[iveg].Cv * 0.3 / veg[iveg][band].LAI);
	      cn[band].fpsn[12] += (veg[iveg][band].GPP * 1.e6 * \
				      veg_con[iveg].Cv * 0.2 / veg[iveg][band].LAI);
	      fsum[14] += veg_con[iveg].Cv * 0.1;
	      fsum[13] += veg_con[iveg].Cv * 0.3;
	      fsum[12] += veg_con[iveg].Cv * 0.2;
	      break;
	    case 10: cn[band].Tveg[17] += (energy[iveg][band].Tcanopy * veg_con[iveg].Cv);
	      cn[band].rveg[17] += (cell[iveg][band].aero_resist[1] * veg_con[iveg].Cv);
              cn[band].zov[17] += (veg_lib->roughness[veg_con[iveg].veg_class] * \
			       veg_con[iveg].Cv);
              cn[band].displ[17] += \
		(veg_lib->displacement[veg_con[iveg].veg_class] *	\
		 veg_con[iveg].Cv);
	      cn[band].fpsn[17] += (veg[iveg][band].GPP * 1.e6 * \
				      veg_con[iveg].Cv / veg[iveg][band].LAI);
	      fsum[17] += veg_con[iveg].Cv;
	      break;
	    }

      cn[band].fwet[0] = 0.0;

      for(iveg = 1; iveg <= MAX_CN_PFTS; iveg++)

	{

	  if(fsum[iveg] != 0.0)
	    {
	      cn[band].Tveg[iveg] /= fsum[iveg];
	      cn[band].rveg[iveg] /= fsum[iveg];
	      cn[band].zov[iveg] /= fsum[iveg];
	      cn[band].displ[iveg] /= fsum[iveg];
              cn[band].fpsn[iveg] /= fsum[iveg];
	    }

	    if(precip > 0.0)
	      cn[band].fwet[iveg] = 1.0;
	    else
	      cn[band].fwet[iveg] = 0.0;
	      }

	}

    }

}
