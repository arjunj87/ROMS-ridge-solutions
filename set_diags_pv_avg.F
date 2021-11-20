#include "cppdefs.h"
#if (defined DIAGNOSTICS_PV && defined AVERAGES)
 
      subroutine set_diags_pv_avg (tile)
      implicit none
      integer tile
# include "param.h"
# include "compute_tile_bounds.h"
      call set_diags_pv_avg_tile (istr,iend,jstr,jend)
      return
      end
 
      subroutine set_diags_pv_avg_tile (istr,iend,jstr,jend)
!
! Compute time-averaged fields within a tile.
! ------- ------------- ------ ------ - -----
! Because of syncronization issues, the delayed mode averaging
! procedure is used. This procedure implies that all fields to be
! averaged are sampled during the next time step, rather than at
! the end of the time step when they were computed.
!
! Thought this algorithm results in somewhat ackwad controlling
! logic it has the advantage that that all fields to be sampled
! correspond to exactly the same time, which is time step "n".
! Particularly, this is done this way because vertical velocity
! corresponding to the newly computed horizontal velocities
! becomes available only during the following time step.
! The same applies to the density field.
!
! The algorithm consists of three logical blocks: (1) initialization
! of the averages arrays: when mod(iic-1,navg).eq.1 the target arrays
! are set to the first contribution; (2) accumulation of averaged
! data, when mod(iic-1,navg).gt.1; and (3) adding the last
! contribution and scaling.
!
      implicit none
      integer istr,iend,jstr,jend, i,j, ilc, iflux
# ifdef SOLVE3D
     &                       , itrc, k
# endif
      real cff
# include "param.h"
# include "scalars.h"
# include "ncvars.h"
# include "grid.h"
# include "ocean2d.h"
# include "ocean3d.h"

# include "averages.h"
#include "diags_pv.h"

!
#include "compute_auxiliary_bounds.h"
!
      ilc=1+iic-ntsdiags_pv_avg  ! number of time step since restart
!
! calculation of averaged fluxes will only be performed if ntsdiags_pv_avg
! is a positive number
!


      if (iic.gt.ntsdiags_pv_avg) then
        if (mod(iic-ntsdiags_pv_avg,nwrtdiags_pv_avg).eq.1) then

          if (ZEROTH_TILE) then
            timediags_pv_avg=time
      MPI_master_only write(*,*) 'started averaging',iic,
     &                       ntsdiags_pv_avg,nwrtdiags_pv_avg
          endif

          do iflux=1,2
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  pvForc_sfc_avg(i,j,iflux)=pvForc_sfc(i,j,iflux)  
           enddo
          enddo
          enddo
          do iflux=1,2
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  pvForc_bot_avg(i,j,iflux)=pvForc_bot(i,j,iflux)
           enddo
          enddo
          enddo
          do iflux=1,NT
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  pvDiab_sfc_avg(i,j,iflux)=pvDiab_sfc(i,j,iflux)
           enddo
          enddo
          enddo
          do iflux=1,NT
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  pvDiab_bot_avg(i,j,iflux)=pvDiab_bot(i,j,iflux)
           enddo
          enddo
          enddo



        elseif (mod(iic-ntsdiags_pv_avg,nwrtdiags_pv_avg).gt.1) then

          if (ZEROTH_TILE) timediags_pv_avg=timediags_pv_avg+time

          do iflux=1,2
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  pvForc_sfc_avg(i,j,iflux) = 
     &                       pvForc_sfc_avg(i,j,iflux) +
     &                       pvForc_sfc(i,j,iflux)
           enddo
          enddo
          enddo
          do iflux=1,2
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  pvForc_bot_avg(i,j,iflux) = 
     &                       pvForc_bot_avg(i,j,iflux) +
     &                       pvForc_bot(i,j,iflux)
           enddo
          enddo
          enddo
          do iflux=1,NT

             do j=JstrR,JendR
                do i=IstrR,IendR 
                  pvDiab_sfc_avg(i,j,iflux) = 
     &                       pvDiab_sfc_avg(i,j,iflux) +
     &                       pvDiab_sfc(i,j,iflux)
           enddo
          enddo

          enddo
          do iflux=1,NT
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  pvDiab_bot_avg(i,j,iflux) =
     &                       pvDiab_bot_avg(i,j,iflux) +
     &                       pvDiab_bot(i,j,iflux)
           enddo
          enddo
          enddo




        elseif (mod(iic-ntsdiags_pv_avg,nwrtdiags_pv_avg).eq.0) then

          cff=1./float(nwrtdiags_pv_avg)
          if (ZEROTH_TILE) then
                        timediags_pv_avg=cff*(timediags_pv_avg+time)
      MPI_master_only write(*,*) 'finish aver. pv',iic,
     &                       ntsdiags_pv_avg,nwrtdiags_pv_avg
          endif

          do iflux=1,2
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  pvForc_sfc_avg(i,j,iflux) = cff *
     &                       ( pvForc_sfc_avg(i,j,iflux) +
     &                         pvForc_sfc(i,j,iflux) )
             enddo
          enddo
          enddo
          do iflux=1,2
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  pvForc_bot_avg(i,j,iflux) = cff *
     &                       ( pvForc_bot_avg(i,j,iflux) +
     &                         pvForc_bot(i,j,iflux) )
           enddo
          enddo
          enddo
          do iflux=1,NT
             do j=JstrR,JendR
                do i=IstrR,IendR 
                  pvDiab_sfc_avg(i,j,iflux) = cff *
     &                       ( pvDiab_sfc_avg(i,j,iflux) +
     &                         pvDiab_sfc(i,j,iflux) )
           enddo
          enddo
          enddo
          do iflux=1,NT
              do j=JstrR,JendR
                do i=IstrR,IendR
                  pvDiab_bot_avg(i,j,iflux) = cff *
     &                       ( pvDiab_bot_avg(i,j,iflux) +
     &                         pvDiab_bot(i,j,iflux) )
           enddo
          enddo
          enddo


        endif
      endif       !<-- iic.gt.ntsavg

      return
      end
#else /* DIAGNOSTICS_PV && AVERAGES */
      subroutine set_diags_pv_avg_empty
      end
#endif /* DIAGNOSTICS_PV && AVERAGES */
