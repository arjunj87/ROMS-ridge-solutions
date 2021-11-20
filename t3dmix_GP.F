#include "cppdefs.h"
#ifdef SOLVE3D
 
      subroutine t3dmix (tile)
      implicit none
      integer tile, itrc
# include "param.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"
      do itrc=1,NT
        call t3dmix_tile (istr,iend,jstr,jend, itrc, A3d(1,1),
     &                                    A2d(1, 1), A2d(1, 2),
     &                          A2d(1,3), A2d(1, 5), A2d(1, 7),
     &                          A2d(1,9), A2d(1,11), A2d(1,13))
      enddo
      return
      end
 
      subroutine t3dmix_tile (istr,iend,jstr,jend, itrc, LapT,
     &                    FX,FE, FC,dTdz, dTdx,dTde,dZdx,dZde)
      implicit none
# include "param.h"
      integer istr,iend,jstr,jend, itrc, i,j,k,k1,k2
      real LapT(PRIVATE_2D_SCRATCH_ARRAY,0:N),
     &       FX(PRIVATE_2D_SCRATCH_ARRAY),
     &       FE(PRIVATE_2D_SCRATCH_ARRAY),
     &       FC(PRIVATE_2D_SCRATCH_ARRAY,2),
     &     dTdz(PRIVATE_2D_SCRATCH_ARRAY,2),   cff,
     &     dTdx(PRIVATE_2D_SCRATCH_ARRAY,2),   cff1,
     &     dTde(PRIVATE_2D_SCRATCH_ARRAY,2),   cff2,
     &     dZdx(PRIVATE_2D_SCRATCH_ARRAY,2),   cff3,
     &     dZde(PRIVATE_2D_SCRATCH_ARRAY,2),   cff4
# include "grid.h"
# include "ocean3d.h"
# include "mixing.h"
# include "climat.h"
# include "scalars.h"
# ifdef DIAGNOSTICS_TS
# include "diags_ts.h"
# endif
# ifdef DIAGNOSTICS_PV
# include "diags_pv.h"
# endif
!
# include "compute_auxiliary_bounds.h"
!
# ifdef MASKING
#  define SWITCH *
# else
#  define SWITCH !
# endif

# ifdef TS_DIF2 
!                                            !!! WARNING: RECURSIVE
! Compute horizontal (along geopotential     !!!   BLOCKING SEQUENCE
! surfaces) Laplacian diffusion terms.       !!! 
!
      k2=1            ! vertical   dTdx,dTde(:,:,k2) k+1   rho-points 
      do k=0,N,+1     ! placement    dTdz,FC(:,:,k2) k+1/2   W-points 
        k1=k2         ! of scratch dTdx,dTde(:,:,k1) k     rho-points 
        k2=3-k1       ! variables    dTdz,FC(:,:,k1) k-1/2   W-points
        if (k.lt.N) then
          do j=jstr,jend
            do i=istr,iend+1
              cff=0.5*(pm(i,j)+pm(i-1,j)) SWITCH umask(i,j)
              dZdx(i,j,k2)=cff*(z_r(i,j,k+1) -z_r(i-1,j,k+1))
              dTdx(i,j,k2)=cff*(t(i,j,k+1,nrhs,itrc)
     &                               -t(i-1,j,k+1,nrhs,itrc)
#  ifdef CLIMAT_TS_MIXH
     &               -tclm(i,j,k+1,itrc)+tclm(i-1,j,k+1,itrc)
#  endif
     &                                                      )
            enddo
          enddo
          do j=jstr,jend+1
            do i=istr,iend
              cff=0.5*(pn(i,j)+pn(i,j-1)) SWITCH vmask(i,j)
              dZde(i,j,k2)=cff*(z_r(i,j,k+1) -z_r(i,j-1,k+1))
              dTde(i,j,k2)=cff*(t(i,j,k+1,nrhs,itrc)
     &                               -t(i,j-1,k+1,nrhs,itrc)
#  ifdef CLIMAT_TS_MIXH
     &               -tclm(i,j,k+1,itrc)+tclm(i,j-1,k+1,itrc)
#  endif
     &                                                      )
            enddo
          enddo
        endif
 
        if (k.eq.0 .or. k.eq.N) then
          do j=jstr-1,jend+1
            do i=istr-1,iend+1
              dTdz(i,j,k2)=0.0     ! Set no-flux bottom or
              FC(i,j,k2)=0.0       ! top boundary conditions.
            enddo
          enddo
        else
          do j=jstr-1,jend+1
            do i=istr-1,iend+1
              dTdz(i,j,k2)=(t(i,j,k+1,nrhs,itrc)-t(i,j,k,nrhs,itrc)
#  ifdef CLIMAT_TS_MIXH
     &                          -tclm(i,j,k+1,itrc)+tclm(i,j,k,itrc)
#  endif
     &                                   )/(z_r(i,j,k+1)-z_r(i,j,k))
            enddo
          enddo
        endif
 
        if (k.gt.0) then
          do j=jstr,jend
            do i=istr,iend+1
              FX(i,j)=0.25*(diff2(i,j,itrc)+diff2(i-1,j,itrc))
     &                     *(Hz(i,j,k)+Hz(i-1,j,k))*dn_u(i,j)*(
     &         dTdx(i  ,j,k1)-0.5*(
     &             min(dZdx(i,j,k1),0.)*(dTdz(i-1,j,k1)+dTdz(i,j,k2))
     &            +max(dZdx(i,j,k1),0.)*(dTdz(i-1,j,k2)+dTdz(i,j,k1))
     &                                                             ))
            enddo
          enddo
          do j=jstr,jend+1
            do i=istr,iend
              FE(i,j)=0.25*(diff2(i,j,itrc)+diff2(i,j-1,itrc))
     &                     *(Hz(i,j,k)+Hz(i,j-1,k))*dm_v(i,j)*(
     &         dTde(i,j,k1)-0.5*(
     &             min(dZde(i,j,k1),0.)*(dTdz(i,j-1,k1)+dTdz(i,j,k2))
     &            +max(dZde(i,j,k1),0.)*(dTdz(i,j-1,k2)+dTdz(i,j,k1))
     &                                                             ))
            enddo
          enddo
          if (k.lt.N) then
            do j=jstr,jend
              do i=istr,iend
 
                cff1=min(dZdx(i  ,j,k1),0.)
                cff2=min(dZdx(i+1,j,k2),0.)
 
                cff3=max(dZdx(i  ,j,k2),0.)
                cff4=max(dZdx(i+1,j,k1),0.)
 
                FC(i,j,k2)=0.5*diff2(i,j,itrc)*(
     &                       cff1*(cff1*dTdz(i,j,k2)-dTdx(i  ,j,k1))
     &                        +cff2*(cff2*dTdz(i,j,k2)-dTdx(i+1,j,k2))
 
     &                        +cff3*(cff3*dTdz(i,j,k2)-dTdx(i  ,j,k2))
     &                        +cff4*(cff4*dTdz(i,j,k2)-dTdx(i+1,j,k1))
     &                                                               )
                cff1=min(dZde(i,j  ,k1),0.)
                cff2=min(dZde(i,j+1,k2),0.)
 
                cff3=max(dZde(i,j  ,k2),0.)
                cff4=max(dZde(i,j+1,k1),0.)
 
                FC(i,j,k2)=FC(i,j,k2)+0.5*diff2(i,j,itrc)*(
     &                        cff1*(cff1*dTdz(i,j,k2)-dTde(i,j  ,k1))
     &                       +cff2*(cff2*dTdz(i,j,k2)-dTde(i,j+1,k2))
 
     &                       +cff3*(cff3*dTdz(i,j,k2)-dTde(i,j  ,k2))
     &                       +cff4*(cff4*dTdz(i,j,k2)-dTde(i,j+1,k1))
     &                                                              )
              enddo
            enddo
          endif
          do j=jstr,jend
            do i=istr,iend
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)+dt*( pm(i,j)*
     &                                   pn(i,j)*( FX(i+1,j)-FX(i,j)
     &                                            +FE(i,j+1)-FE(i,j))
     &                                         +FC(i,j,k2)-FC(i,j,k1)
     &                                                   )/Hz(i,j,k)

#if defined DIAGNOSTICS_PV && !defined DIAGNOSTICS_TS

              if (k.eq.N) then
                cff1=pm(i,j)*pn(i,j)
                pvDiab_sfc(i,j,itrc)=pvDiab_sfc(i,j,itrc)
     &                     + (FX(i+1,j)-FX(i,j)
     &                     + FE(i,j+1)-FE(i,j)
     &                     +(FC(i,j,k2)-FC(i,j,k1))/cff1)
     &                     *cff1/Hz(i,j,k)
              endif

              if (k.eq.1) then
                cff1=pm(i,j)*pn(i,j)
                pvDiab_bot(i,j,itrc)=pvDiab_bot(i,j,itrc)
     &                     + (FX(i+1,j)-FX(i,j)
     &                     + FE(i,j+1)-FE(i,j)
     &                     +(FC(i,j,k2)-FC(i,j,k1))/cff1)
     &                     *cff1/Hz(i,j,k)
              endif

#  endif /* DIAGNOSTICS */



#  ifdef DIAGNOSTICS_TS

              cff1=pm(i,j)*pn(i,j)
              THmix(i,j,k,itrc)=FX(i+1,j)-FX(i,j)
     &                         +FE(i,j+1)-FE(i,j)
     &                         +(FC(i,j,k2)-FC(i,j,k1))/cff1
#   ifdef MASKING
     &                                          * rmask(i,j)
#   endif

              Trate(i,j,k,itrc)=(Hz(i,j,k)*t(i,j,k,nnew,itrc)
     &                      -Hz_bak(i,j,k)*t(i,j,k,nstp,itrc))
     &                                              /(dt*cff1)
#   ifdef MASKING
     &                                          * rmask(i,j)
#   endif
!
! Divide all diagnostic terms by the cell volume
! (Hz(i,j,k,itrc)/(pm(i,j).*pn(i,j)). There after the unit
! of diag terms are: (unit of tracers)* s-1.
!
              THmix(i,j,k,itrc)=THmix(i,j,k,itrc)*cff1/Hz(i,j,k)
              Trate(i,j,k,itrc)=Trate(i,j,k,itrc)*cff1/Hz(i,j,k)

#  endif /* DIAGNOSTICS */

            enddo
          enddo
        endif  ! <-- k.gt.0
      enddo   ! --> k


#if defined DIAGNOSTICS_TS && defined DIAGNOSTICS_PV
         do j=jstr,jend
           do i=istr,iend
             pvDiab_sfc(i,j,itrc)=pvDiab_sfc(i,j,itrc) + THmix(i,j,N,itrc)
             pvDiab_bot(i,j,itrc)=pvDiab_bot(i,j,itrc) + THmix(i,j,1,itrc)
           enddo
         enddo
# endif /* DIAGNOSTICS_PV */










# else
      print*,'Please, choose TS_DIF2 instead of TS_DIF4'
      stop      
# endif /* T_DIF2 */

# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
        call exchange_r3d_tile (istr,iend,jstr,jend,
     &                          t(START_2D_ARRAY,1,nnew,itrc))
# endif

      return
      end
#else
      subroutine t3dmix_empty
      end
#endif /* SOLVE3D */
 
