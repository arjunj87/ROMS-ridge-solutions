#include "cppdefs.h"
 
      subroutine t3dmix (tile)
      implicit none
      integer tile, itrc
#include "param.h"
#include "private_scratch.h"
#include "compute_tile_bounds.h"
      do itrc=1,NT
        call t3dmix_tile (istr,iend,jstr,jend, itrc, A2d(1,1),
     &                                     A2d(1,2), A2d(1,3))
      enddo
      return
      end
 
      subroutine t3dmix_tile (istr,iend,jstr,jend, itrc, FX,FE,LapT)
!
!  This routine computes horizontal mixing of tracers
!  along S-surfaces.
!
      implicit none
#   include "param.h"
      integer itrc, istr,iend,jstr,jend, i,j,k
      real FX(PRIVATE_2D_SCRATCH_ARRAY),     cff1,
     &     FE(PRIVATE_2D_SCRATCH_ARRAY),     cff2,
     &     LapT(PRIVATE_2D_SCRATCH_ARRAY)
#include "grid.h"
#include "ocean3d.h"
#include "mixing.h"
#include "climat.h"
#include "scalars.h"
!
#include "compute_auxiliary_bounds.h"
!
      do k=1,N
 
#ifdef TS_DIF2
!
!  Add in horizontal Laplacian diffusion along constant S-surfaces.
!--------------------------------------------------------------------
!  Compute XI- and ETA-components of diffusive tracer flux.
!

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++ Compute the fluxes for the second laplacian

        do j=jstr,jend
          do i=istr,iend+1
            FX(i,j)=
     &                   0.5*pmon_u(i,j)*(Hz(i,j,k)+Hz(i-1,j,k))*(
     &                     tstar(i,j,k,itrc)-tstar(i-1,j,k,itrc)
     &                                                              )
# ifdef MASKING
     &                                                   *umask(i,j)
# endif
          enddo
        enddo
        do j=jstr,jend+1
          do i=istr,iend
            FE(i,j)=
     &                    0.5*pnom_v(i,j)*(Hz(i,j,k)+Hz(i,j-1,k))*(
     &                      tstar(i,j,k,itrc)-tstar(i,j-1,k,itrc)
     &                                                              )
# ifdef MASKING
     &                                                   *vmask(i,j)
# endif
          enddo
        enddo
!+++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        do j=jstr,jend
          do i=istr,iend
          t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)+pm(i,j)*pn(i,j)*
     &    (FX(i+1,j)-FX(i,j)+FE(i,j+1)-FE(i,j))/Hz(i,j,k)                
          enddo
        enddo
!+++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#endif /* TS_DIF2 */
      enddo

# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
        call exchange_r3d_tile (istr,iend,jstr,jend,
     &                          t(START_2D_ARRAY,1,nnew,itrc))
# endif
      return
      end
 
