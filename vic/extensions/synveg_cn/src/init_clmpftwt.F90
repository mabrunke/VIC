!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_clmpftwt
!
! !INTERFACE:
  subroutine init_clmpftwt(g, begc, endc, p, nveg, vegfrac)
!
! !DESCRIPTION:
! Initialize clmtype PFT weights
!
! !USES:
!
  use shr_kind_mod, only : r8 => shr_kind_r8
  use clmtype
!
! !ARGUMENTS:
    implicit none
    integer  :: begc, endc            ! clump beg and ending column indices
    integer  :: g, p                  ! current gridcell gridcell and PFT indices
    real(r8) :: vegfrac               ! VIC vegetation fractions
    integer  :: nveg                  ! # VIC vegetation types
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
! 10/23/12: Adapted for use in VIC by Michael Brunke
!
!EOP
!
! LOCAL VARAIBLES:
    integer :: c,m      ! indices
!------------------------------------------------------------------------

    ! Determine necessary indices

     do c = begc, endc
         clm3%g%c%p%wtgcell((g-1)*endc*nveg + (c-1)*nveg + p + 1) = &
		vegfrac / endc
          clm3%g%c%p%wtcol((g-1)*endc*nveg + (c-1)*nveg + p + 1) = &
		vegfrac / endc
     enddo

end subroutine init_clmpftwt

