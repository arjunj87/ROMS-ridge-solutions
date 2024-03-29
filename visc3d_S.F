#include "cppdefs.h"
# if defined SOLVE3D && defined UV_VIS2
 
      subroutine visc3d (tile)
      implicit none
      integer tile
# include "param.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"
      call visc3d_tile (istr,iend,jstr,jend, A2d(1,1), A2d(1,2),
     &                                       A2d(1,3), A2d(1,4))
      return
      end
 
      subroutine visc3d_tile (istr,iend,jstr,jend, UFx,UFe,VFx,VFe)
!
! Computes harmonic mixing of momentum, along constant S-surfaces
! as horizontal divergence of the stress tensor. A transverse
! isotropy is assumed so the stress tensor is splitted into vertical
! and horizontal subtensors.
!
! Reference:
!
! Wajsowicz, R.C, 1993: A consistent formulation of the anisotropic
!     stress tensor for use in models of the large-scale ocean
!     circulation, JCP, 105, 333-338.
!
! Sadourny, R. and K. Maynard, 1997: Formulations of lateral
!     diffusion in geophysical fluid dynamics models, In "Numerical
!     Methods of Atmospheric and Oceanic Modelling". Lin, Laprise,
!     and Ritchie, Eds., NRC Research Press, 547-556.
!
! Griffies, S.M. and R.W. Hallberg, 2000: Biharmonic friction with
!     a Smagorinsky-like viscosity for use in large-scale eddy-
!     permitting ocean models, Mon. Wea. Rev., 128, 8, 2935-2946.
!
      implicit none
      integer istr,iend,jstr,jend, i,j,k, indx
      real UFe(PRIVATE_2D_SCRATCH_ARRAY),
     &     UFx(PRIVATE_2D_SCRATCH_ARRAY), cff,
     &     VFe(PRIVATE_2D_SCRATCH_ARRAY),
     &     VFx(PRIVATE_2D_SCRATCH_ARRAY)
# include "param.h"
# include "scalars.h"
# include "grid.h"
# include "ocean3d.h"
# include "coupling.h"
# include "mixing.h"
!
# include "compute_auxiliary_bounds.h"
!
      indx=3-nstp    !--> time index for target arrays;
!
!  Compute horizontal harmonic viscosity along constant S-surfaces.
!
!  Compute flux-components of the horizontal divergence of the stress
!  tensor (m5/s2) in XI- and ETA-directions.
!
      do k=1,N
        do j=jstrV-1,jend
          do i=istrU-1,iend
            cff=0.5*Hz(i,j,k)*visc2_r(i,j)*(
 
     &      dn_r(i,j)*pm(i,j)*( (pn(i  ,j)+pn(i+1,j))*u(i+1,j,k,nstp)
     &                         -(pn(i-1,j)+pn(i  ,j))*u(i  ,j,k,nstp)
     &                                                              )

     &     -dm_r(i,j)*pn(i,j)*( (pm(i,j  )+pm(i,j+1))*v(i,j+1,k,nstp)
     &                         -(pm(i,j-1)+pm(i,j  ))*v(i,j  ,k,nstp)
     &                                                             ))
            UFx(i,j)=  cff*dn_r(i,j)*dn_r(i,j)
            VFe(i,j)= -cff*dm_r(i,j)*dm_r(i,j)
          enddo
        enddo
        do j=jstr,jend+1
          do i=istr,iend+1
            cff=0.125*(Hz(i-1,j,k)+Hz(i,j,k)+Hz(i-1,j-1,k)
     &                         +Hz(i,j-1,k))*visc2_p(i,j)*(
 
     &      0.25*(pm(i-1,j)+pm(i,j)+pm(i-1,j-1)+pm(i,j-1))*dn_p(i,j)
     &                     *( (pn(i  ,j-1)+pn(i  ,j))*v(i  ,j,k,nstp)
     &                       -(pn(i-1,j-1)+pn(i-1,j))*v(i-1,j,k,nstp)
     &                                                              )

     &     +0.25*(pn(i-1,j)+pn(i,j)+pn(i-1,j-1)+pn(i,j-1))*dm_p(i,j)
     &                     *( (pm(i-1,j  )+pm(i,j  ))*u(i,j  ,k,nstp)
     &                       -(pm(i-1,j-1)+pm(i,j-1))*u(i,j-1,k,nstp)
     &                                                             ))
# ifdef MASKING
     &                                                   *pmask(i,j)
# endif
            UFe(i,j)= cff*dm_p(i,j)*dm_p(i,j)
            VFx(i,j)= cff*dn_p(i,j)*dn_p(i,j)
          enddo
        enddo
!
! Apply viscous terms. Note that at this stage arrays u,v(...,3-nstp)
! contain Hz*U and Hz*V with units of [m2/s].   Also compute vertical
! integral of viscous terms and add it into coupling terms for the
! barotropic mode
!
        do j=jstr,jend
          do i=istrU,iend
            cff=0.125*(pm(i-1,j)+pm(i,j))*(pn(i-1,j) +pn(i,j))
     &             *( (pn(i-1,j)+pn(i,j))*(UFx(i,j)-UFx(i-1,j))
     &               +(pm(i-1,j)+pm(i,j))*(UFe(i,j+1)-UFe(i,j))
     &                                                        )
            rufrc(i,j)=rufrc(i,j) + cff
            u(i,j,k,indx)=u(i,j,k,indx) + dt*cff
          enddo
        enddo
        do j=jstrV,jend
          do i=istr,iend
            cff=0.125*(pm(i,j)+pm(i,j-1))*(pn(i,j) +pn(i,j-1))
     &             *( (pn(i,j-1)+pn(i,j))*(VFx(i+1,j)-VFx(i,j))
     &               +(pm(i,j-1)+pm(i,j))*(VFe(i,j)-VFe(i,j-1))
     &                                                        )
            rvfrc(i,j)=rvfrc(i,j) + cff
            v(i,j,k,indx)=v(i,j,k,indx) + dt*cff
          enddo
        enddo
      enddo
      return
      end
#else
      subroutine visc3d_empty
      end
#endif
 
