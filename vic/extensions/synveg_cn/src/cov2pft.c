#include <vic_def.h>
#include <vic_run.h>
#include <synveg_cn.h>

void cov2pft(double cov_val, int veg_class, double veg_fract, double pft_vals[MAX_PFT+1])

/****************************************************************************
                                                                           
  pft2cov:  converts UMD vegetation cover values from CLM PFT values.

****************************************************************************/
{

  switch(veg_class)
    {
    case 0: pft_vals[2] += cov_val * veg_fract;
      break;
    case 1: pft_vals[5] += cov_val * veg_fract;
      break;
    case 2: pft_vals[3] += cov_val * veg_fract;
      break;
    case 3: pft_vals[8] += cov_val * veg_fract;
      break;
    case 4: pft_vals[2] += cov_val * 0.5 * veg_fract;
      pft_vals[8] += cov_val * 0.5 * veg_fract;
      break;
    case 5: pft_vals[2] += cov_val * 0.4 * veg_fract;
      pft_vals[8] += cov_val * 0.4 * veg_fract;
      pft_vals[9] += cov_val * 0.1 * veg_fract;
      pft_vals[10] += cov_val * 0.1 * veg_fract;
      break;
    case 6: pft_vals[14] += cov_val * 0.1 * veg_fract;
      pft_vals[13] += cov_val * 0.3 * veg_fract;
      pft_vals[12] += cov_val * 0.2 * veg_fract;
      pft_vals[2] += cov_val * 0.2 * veg_fract;
      pft_vals[8] += cov_val * 0.2 * veg_fract;
      break;
    case 7: pft_vals[9] += cov_val * 0.3 * veg_fract;
      pft_vals[10] += cov_val * 0.3 * veg_fract;
      pft_vals[11] += cov_val * 0.2 * veg_fract;
      pft_vals[2] += cov_val * 0.1 * veg_fract;
      pft_vals[8] += cov_val * 0.1 * veg_fract;
      break;
    case 8: pft_vals[9] += cov_val * 0.2 * veg_fract;
      pft_vals[10] += cov_val * 0.2 * veg_fract;
      pft_vals[11] += cov_val * 0.2 * veg_fract;
      pft_vals[14] += cov_val * 0.05 * veg_fract;
      pft_vals[13] += cov_val * 0.15 * veg_fract;
      pft_vals[14] += cov_val * 0.1 * veg_fract;
      pft_vals[0] += cov_val * 0.1 * veg_fract;
      break;
    case 9: pft_vals[14] += cov_val * 0.1 * veg_fract;
      pft_vals[13] += cov_val * 0.3 * veg_fract;
      pft_vals[12] += cov_val * 0.2 * veg_fract;
      pft_vals[0] += cov_val * 0.4 * veg_fract;
      break;
    case 10: pft_vals[17] += cov_val * veg_fract;
      break;
    default: pft_vals[0] += cov_val * veg_fract;
      break;
    }

}

