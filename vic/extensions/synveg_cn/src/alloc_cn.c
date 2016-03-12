#include <vic_def.h>
#include <vic_driver_image.h>
#include <synveg_cn.h>

static char vcid[] = "$Id$";

/****************************************************************************/
/*			       alloc_cn()                                */
/****************************************************************************/
cn_data_struct *alloc_cn(int Nbands, int Nnode)
/*******************************************************************
  alloc_cn       Created by Michael Brunke

  Allocates data structure containing CN quantities in CN data
  structure.

  Modifications:

*******************************************************************/
{

  int iband, iveg, k, r, n;
  cn_data_struct *temp;

  temp = (cn_data_struct *) calloc(Nbands, sizeof(cn_data_struct *));
  if (temp == NULL)
      {
	log_err("Memory allocation error in alloc_cn().");
      }

      for(iband = 0; iband < Nbands; iband++)
	{
	  for(iveg = 0; iveg < MAX_PFT + 1; iveg++)
	    {
	      temp[iband].Tveg[iveg] = 0.0;
	      temp[iband].rveg[iveg] = 0.0;
	      temp[iband].zov[iveg] = 0.0;
	      temp[iband].displ[iveg] = 0.0;
	      temp[iband].fwet[iveg] = 0.0;
	      temp[iband].LAI[iveg] = 0.0;
	      temp[iband].dormant_flag[iveg] = 0.0;
	      temp[iband].days_active[iveg] = 0.0;
	      temp[iband].onset_flag[iveg] = 0.0;
	      temp[iband].onset_counter[iveg] = 0.0;
	      temp[iband].onset_gddflag[iveg] = 0.0;
	      temp[iband].onset_fdd[iveg] = 0.0;
	      temp[iband].onset_gdd[iveg] = 0.0;
	      temp[iband].onset_swi[iveg] = 0.0;
	      temp[iband].offset_flag[iveg] = 0.0;
	      temp[iband].offset_counter[iveg] = 0.0;
	      temp[iband].offset_fdd[iveg] = 0.0;
	      temp[iband].offset_swi[iveg] = 0.0;
	      temp[iband].lgsf[iveg] = 0.0;
	      temp[iband].bglfr[iveg] = 0.0;
	      temp[iband].bgtr[iveg] = 0.0;
	      temp[iband].dayl[iveg] = 0.0;
	      temp[iband].prev_dayl[iveg] = 0.0;
	      temp[iband].annavg_t2m[iveg] = 0.0;
	      temp[iband].tempavg_t2m[iveg] = 0.0;
	      temp[iband].gpp2[iveg] = 0.0;
	      temp[iband].availc[iveg] = 0.0;
	      temp[iband].xsmrpool_recover[iveg] = 0.0;
	      temp[iband].alloc_pnow[iveg] = 0.0;
	      temp[iband].c_allometry[iveg] = 0.0;
	      temp[iband].n_allometry[iveg] = 0.0;
	      temp[iband].tempsum_potential_gpp[iveg] = 0.0;
	      temp[iband].annsum_potential_gpp[iveg] = 0.0;
	      temp[iband].tempmax_retransn[iveg] = 0.0;
	      temp[iband].annmax_retransn[iveg] = 0.0;
	      temp[iband].avail_retransn[iveg] = 0.0;
	      temp[iband].plant_nalloc[iveg] = 0.0;
	      temp[iband].plant_calloc[iveg] = 0.0;
	      temp[iband].excess_cflux[iveg] = 0.0;
	      temp[iband].downreg[iveg] = 0.0;
	      temp[iband].prev_leafc_to_litter[iveg] = 0.0;
	      temp[iband].prev_frootc_to_litter[iveg] = 0.0;
	      temp[iband].tempsum_npp[iveg] = 0.0;
	      temp[iband].annsum_npp[iveg] = 0.0;
	      temp[iband].gpp[iveg] = 0.0;
	      temp[iband].npp[iveg] = 0.0;
	      temp[iband].ar[iveg] = 0.0;
	      temp[iband].leafc[iveg] = 0.0;
	      temp[iband].leafc_storage[iveg] = 0.0;
	      temp[iband].leafc_xfer[iveg] = 0.0;
	      temp[iband].frootc[iveg] = 0.0;
	      temp[iband].frootc_storage[iveg] = 0.0;
	      temp[iband].frootc_xfer[iveg] = 0.0;
	      temp[iband].livestemc[iveg] = 0.0;
	      temp[iband].livestemc_storage[iveg] = 0.0;
	      temp[iband].livestemc_xfer[iveg] = 0.0;
	      temp[iband].deadstemc[iveg] = 0.0;
	      temp[iband].deadstemc_storage[iveg] = 0.0;
	      temp[iband].deadstemc_xfer[iveg] = 0.0;
	      temp[iband].livecrootc[iveg] = 0.0;
	      temp[iband].livecrootc_storage[iveg] = 0.0;
	      temp[iband].livecrootc_xfer[iveg] = 0.0;
	      temp[iband].deadcrootc[iveg] = 0.0;
	      temp[iband].deadcrootc_storage[iveg] = 0.0;
	      temp[iband].deadcrootc_xfer[iveg] = 0.0;
	      temp[iband].gresp_storage[iveg] = 0.0;
	      temp[iband].gresp_xfer[iveg] = 0.0;
	      temp[iband].cpool[iveg] = 0.0;
	      temp[iband].xsmrpool[iveg] = 0.0;
	      temp[iband].pft_ctrunc[iveg] = 0.0;
	      temp[iband].totvegc[iveg] = 0.0;
	      temp[iband].woodc[iveg] = 0.0;
	      temp[iband].leafn[iveg] = 0.0;
	      temp[iband].leafn_storage[iveg] = 0.0;
	      temp[iband].leafn_xfer[iveg] = 0.0;
	      temp[iband].frootn[iveg] = 0.0;
	      temp[iband].frootn_storage[iveg] = 0.0;
	      temp[iband].frootn_xfer[iveg] = 0.0;
	      temp[iband].livestemn[iveg] = 0.0;
	      temp[iband].livestemn_storage[iveg] = 0.0;
	      temp[iband].livestemn_xfer[iveg] = 0.0;
	      temp[iband].deadstemn[iveg] = 0.0;
	      temp[iband].deadstemn_storage[iveg] = 0.0;
	      temp[iband].deadstemn_xfer[iveg] = 0.0;
	      temp[iband].livecrootn[iveg] = 0.0;
	      temp[iband].livecrootn_storage[iveg] = 0.0;
	      temp[iband].livecrootn_xfer[iveg] = 0.0;
	      temp[iband].deadcrootn[iveg] = 0.0;
	      temp[iband].deadcrootn_storage[iveg] = 0.0;
	      temp[iband].deadcrootn_xfer[iveg] = 0.0;
	      temp[iband].retransn[iveg] = 0.0;
	      temp[iband].npool[iveg] = 0.0;
	      temp[iband].pft_ntrunc[iveg] = 0.0;
	    }

	  temp[iband].Tair = 0.0;
	  temp[iband].vp = 0.0;
	  temp[iband].vpd = 0.0;
	  temp[iband].psfc = 0.0;
	  temp[iband].lwrad = 0.0;
	  temp[iband].swrad = 0.0;
	  temp[iband].precip = 0.0;
	  for(r = 0; r < 2; r++)
	    {
	      temp[iband].swrd[r] = 0.0;
	      temp[iband].swri[r] = 0.0;
	    }
	  temp[iband].alb = 0.0;
	  for(k = 0; k < Nnode + 2; k++)
	    {
	      temp[iband].t_soisno[k] = 0.0;
	      temp[iband].z[k] = 0.0;
	      temp[iband].dz[k] = 0.0;
	      temp[iband].ice[k] = 0.0;
	    }
	  temp[iband].z0 = 0.0;
	  temp[iband].z0s = 0.0;
	  temp[iband].baseflow = 0.0;
	  temp[iband].snowdep = 0.0;
	  temp[iband].decl = 0.0;
	  temp[iband].fpi = 0.0;
	  temp[iband].fpg = 0.0;
	  temp[iband].annsum_counter = 0.0;
	  temp[iband].cannsum_npp = 0.0;
	  temp[iband].cannavg_t2m = 0.0;
	  for(k = 0; k < Nnode; k++)
	    {
	      temp[iband].watfc[k] = 0.0;
	      temp[iband].moist[k] = 0.0;
	      for(iveg = 0; iveg < MAX_PFT; iveg++)
		  temp[iband].rootfr[iveg][k] = 0.0;
	      for(n = 0; n < 2; n++)
		{
		  temp[iband].bsw[k][n] = 0.0;
		  temp[iband].sucsat[k][n] = 0.0;
		  temp[iband].soisuc[k][n] = 0.0;
		}
	    }
	  temp[iband].me = 0.0;
	  temp[iband].fire_prob = 0.0;
	  temp[iband].mean_fire_prob = 0.0;
	  temp[iband].fireseasonl = 0.0;
	  temp[iband].ann_farea_burned = 0.0;
	  temp[iband].hr = 0.0;
	  temp[iband].lithr = 0.0;
	  temp[iband].nee = 0.0;
	  temp[iband].nep = 0.0;
	  temp[iband].cwdc = 0.0;
	  temp[iband].litr1c = 0.0;
	  temp[iband].litr2c = 0.0;
	  temp[iband].litr3c = 0.0;
	  temp[iband].soil1c = 0.0;
	  temp[iband].soil2c = 0.0;
	  temp[iband].soil3c = 0.0;
	  temp[iband].soil4c = 0.0;
	  temp[iband].seedc = 0.0;
	  temp[iband].col_ctrunc = 0.0;
	  temp[iband].totlitc = 0.0;
          temp[iband].totsomc = 0.0;
	  temp[iband].totcolc = 0.0;
	  temp[iband].prod10c = 0.0;
	  temp[iband].prod100c = 0.0;
	  temp[iband].cwdn = 0.0;
	  temp[iband].litr1n = 0.0;
	  temp[iband].litr2n = 0.0;
	  temp[iband].litr3n = 0.0;
	  temp[iband].soil1n = 0.0;
	  temp[iband].soil2n = 0.0;
	  temp[iband].soil3n = 0.0;
	  temp[iband].soil4n = 0.0;
	  temp[iband].sminn = 0.0;
	  temp[iband].seedn = 0.0;
	  temp[iband].col_ntrunc = 0.0;
	  temp[iband].totcoln = 0.0;
	  temp[iband].prod10n = 0.0;
	  temp[iband].prod100n = 0.0;

	  }

      /*    } */

  return temp;

}

/****************************************************************************/
/*	      		  free_cn()                                         */
/****************************************************************************/
void free_cn(cn_data_struct **cn)
/***************************************************************************
  Modifications:
***************************************************************************/
{

  if (*cn == NULL)
    return;

  free(*cn);
}
