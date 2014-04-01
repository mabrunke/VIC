#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

double
linear_interp(double x,
              double lx,
              double ux,
              double ly,
              double uy)
{
    return((x - lx) / (ux - lx) * (uy - ly) + ly);
}

double
exp_interp(double x,
           double lx,
           double ux,
           double ly,
           double uy)
{
/**********************************************************************
   This subroutine interpolates the soil temperature at a given depth,
   under the assumption that the temperature decays exponentially with
   depth from a surface temperature of "ly" to an asymptotic limit of
   "uy".  "ux" here is the "damping" depth, at which difference between
   temperature at that depth and the asymptotic deep temperature is 1/e
   of the diffence at the surface.

   Modifications:

   2012-Feb-02 Fixed typo in the original formulation, which omitted the
              division by ux.						TJB
**********************************************************************/
    return(uy + (ly - uy) * exp(-(x - lx) / ux));
}

double
modify_Ksat(double Temp)
{
/**********************************************************************
        modify_Ksat	Keith Cherkauer		February 12, 1997

   This subroutine returns a parameter to multiply with Ksat to modify
   it for the effects of temperature on the viscosity and density of
   water.  It is assumed that the given Ksat value was measured at
   20C (68F).

   Viscosity and density taken from Linsley, "Hydrology for Engineers",
      A-10

        Temp	Rho	Mu	Factor
        C	kg/m^3	mPa-s
        0	999.84	1.79	0.560
        5	999.96	1.52	0.659
        10	999.70	1.31	0.770
        15	999.10	1.14	0.878
        20	998.21	1.00	1.00
        25	997.05	0.890	1.12
        30	995.65	0.798	1.25
        35	994.04	0.719	1.39
        40	992.22	0.653	1.52

**********************************************************************/

    extern option_struct options;

    double               Factor;

    /** formula generated by multiple regression against kinematic
        viscosity data from the Handbook of Chemistry and Physics **/
    if (options.FROZEN_SOIL) {
        Factor = 0.003557 / (0.006534 - 0.0002282 * Temp + 4.794e-6 * (Temp) *
                             (Temp) - 4.143e-8 * (Temp) * (Temp) * (Temp));
    }
    else {
        Factor = 1.;
    }

    if (Factor > 2.) {
        Factor = 2.;
    }

    /*return (Factor);*/
    return (1.0);
}
