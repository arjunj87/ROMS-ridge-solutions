#include "cppdefs.h"
 
 
      subroutine visc_rescale (tile)
      implicit none
      integer tile
#include "param.h"
#include "compute_tile_bounds.h"
      call visc_rescale_tile (istr,iend,jstr,jend)
      return
      end
 
 
 
      subroutine visc_rescale_tile (istr,iend,jstr,jend)
!
! Rescale variable horizontal viscosity and mixing coefficients
! according to the local grid box size. The maximum grid size
! "grdmax" is determined in metrics and it is assumed known in
! this subroutine.
!
      implicit none
      integer istr,iend,jstr,jend, i,j,itrc
      real cff
#include "param.h"
#include "grid.h"
#include "mixing.h"
#include "scalars.h"
!
#include "compute_auxiliary_bounds.h"
!
#if defined UV_VIS2 && defined VIS_GRID
      cff=visc2/grdmax
      do j=jstrR,jendR
        do i=istrR,iendR
          visc2_r(i,j)=cff*grdscl(i,j)
        enddo
      enddo
      cff=visc2*0.25/grdmax
      do j=jstr,jendR
        do i=istr,iendR
          visc2_p(i,j)=cff*( grdscl(i,j  )+grdscl(i-1,j  )
     &                      +grdscl(i,j-1)+grdscl(i-1,j-1))
        enddo
      enddo
#endif
#if defined TS_DIF2 && defined DIF_GRID
      do itrc=1,NT
        cff=tnu2(itrc)/grdmax
        do j=jstrR,jendR
          do i=istrR,iendR
            diff2(i,j,itrc)=cff*grdscl(i,j)
          enddo
        enddo
      enddo
#endif















      return
      end
 
