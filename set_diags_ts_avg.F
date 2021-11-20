#include "cppdefs.h"
#if (defined DIAGNOSTICS_TS && defined AVERAGES)
 
      subroutine set_diags_ts_avg (tile)
      implicit none
      integer tile
# include "param.h"
# include "compute_tile_bounds.h"
      call set_diags_ts_avg_tile (istr,iend,jstr,jend)
      return
      end
 
      subroutine set_diags_ts_avg_tile (istr,iend,jstr,jend)
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
#include "diags_ts.h"

!
#include "compute_auxiliary_bounds.h"
!
      ilc=1+iic-ntsdiags_ts_avg   ! number of time step since restart
!
! calculation of averaged fluxes will only be performed if ntsdiags_ts_avg
! is a positive number
!


      if (iic.gt.ntsdiags_ts_avg) then
        if (mod(iic-ntsdiags_ts_avg,nwrtdiags_ts_avg).eq.1) then

          if (ZEROTH_TILE) then
            timediags_ts_avg=time
      MPI_master_only write(*,*) 'started averaging',iic,
     &                       ntsdiags_ts_avg,nwrtdiags_ts_avg
          endif


          do iflux=1,NT
            do k=1,N
              do j=JstrR,JendR
                do i=IstrR,IendR 
                  TXadv_avg(i,j,k,iflux)=TXadv(i,j,k,iflux)  
                enddo
              enddo
            enddo
          enddo
          do iflux=1,NT
            do k=1,N
              do j=JstrR,JendR
                do i=IstrR,IendR 
                  TYadv_avg(i,j,k,iflux)=TYadv(i,j,k,iflux)
                enddo
              enddo
            enddo
          enddo
          do iflux=1,NT
            do k=1,N
              do j=JstrR,JendR
                do i=IstrR,IendR 
                  TVadv_avg(i,j,k,iflux)=TVadv(i,j,k,iflux)
                enddo
              enddo
            enddo
          enddo
          do iflux=1,NT
            do k=1,N
              do j=JstrR,JendR
                do i=IstrR,IendR
                  TForc_avg(i,j,k,iflux)=TForc(i,j,k,iflux)
                enddo
              enddo
            enddo
          enddo
          do iflux=1,NT
            do k=1,N
              do j=JstrR,JendR
                do i=IstrR,IendR
                  THmix_avg(i,j,k,iflux)=THmix(i,j,k,iflux)
                enddo
              enddo
            enddo
          enddo
          do iflux=1,NT
            do k=1,N
              do j=JstrR,JendR
                do i=IstrR,IendR 
                  TVmix_avg(i,j,k,iflux)=TVmix(i,j,k,iflux)
                end do
              end do
            end do
          end do


          do iflux=1,NT
            do k=1,N
              do j=JstrR,JendR
                do i=IstrR,IendR 
                  Trate_avg(i,j,k,iflux)=Trate(i,j,k,iflux)
                enddo
              enddo
            enddo
          enddo



        elseif (mod(iic-ntsdiags_ts_avg,nwrtdiags_ts_avg).gt.1) then

          if (ZEROTH_TILE) timediags_ts_avg=timediags_ts_avg+time

          do iflux=1,NT
            do k=1,N
              do j=JstrR,JendR
                do i=IstrR,IendR 
                  TXadv_avg(i,j,k,iflux) = 
     &                       TXadv_avg(i,j,k,iflux) +
     &                       TXadv(i,j,k,iflux)
                enddo
              enddo
            enddo
          enddo
          do iflux=1,NT
            do k=1,N
              do j=JstrR,JendR
                do i=IstrR,IendR 
                  TYadv_avg(i,j,k,iflux) = 
     &                       TYadv_avg(i,j,k,iflux) +
     &                       TYadv(i,j,k,iflux)
                enddo
              enddo
            enddo
          enddo
          do iflux=1,NT
            do k =1,N
              do j=JstrR,JendR
                do i=IstrR,IendR 
                  TVadv_avg(i,j,k,iflux) = 
     &                       TVadv_avg(i,j,k,iflux) +
     &                       TVadv(i,j,k,iflux)
                enddo
              enddo
            enddo
          enddo
          do iflux=1,NT
            do k=1,N
              do j=JstrR,JendR
                do i=IstrR,IendR
                  TForc_avg(i,j,k,iflux) =
     &                       TForc_avg(i,j,k,iflux) +
     &                       TForc(i,j,k,iflux)
                enddo
              enddo
            enddo
          enddo

          do iflux=1,NT
            do k=1,N
              do j=JstrR,JendR
                do i=IstrR,IendR
                  THmix_avg(i,j,k,iflux) =
     &                       THmix_avg(i,j,k,iflux) +
     &                       THmix(i,j,k,iflux)
                enddo
              enddo
            enddo
          enddo
          do iflux=1,NT
            do k=1,N
              do j=JstrR,JendR
                do i=IstrR,IendR 
                  TVmix_avg(i,j,k,iflux) = 
     &                       TVmix_avg(i,j,k,iflux) +
     &                       TVmix(i,j,k,iflux)
                enddo
              enddo
            enddo
          enddo


          do iflux=1,NT
            do k=1,N
              do j=JstrR,JendR
                do i=IstrR,IendR 
                  Trate_avg(i,j,k,iflux) = 
     &                    Trate_avg(i,j,k,iflux) + 
     &                    Trate(i,j,k,iflux)  
                enddo
              enddo
            enddo
          enddo


        elseif (mod(iic-ntsdiags_ts_avg,nwrtdiags_ts_avg).eq.0) then

          cff=1./float(nwrtdiags_ts_avg)
          if (ZEROTH_TILE) then
                   timediags_ts_avg=cff*(timediags_ts_avg+time)
      MPI_master_only write(*,*) 'finish avera. ts',iic,
     &                       ntsdiags_ts_avg,nwrtdiags_ts_avg

          endif


          do iflux=1,NT
            do k=1,N
              do j=JstrR,JendR
                do i=IstrR,IendR 
                  TXadv_avg(i,j,k,iflux) = cff *
     &                       ( TXadv_avg(i,j,k,iflux) +
     &                         TXadv(i,j,k,iflux) )
                   enddo
                enddo
             enddo
          enddo
          do iflux=1,NT
            do k=1,N
              do j=JstrR,JendR
                do i=IstrR,IendR 
                  TYadv_avg(i,j,k,iflux) = cff *
     &                       ( TYadv_avg(i,j,k,iflux) +
     &                         TYadv(i,j,k,iflux) )
                enddo
              enddo
            enddo
          enddo
          do iflux=1,NT
            do k=1,N
              do j=JstrR,JendR
                do i=IstrR,IendR 
                  TVadv_avg(i,j,k,iflux) = cff *
     &                       ( TVadv_avg(i,j,k,iflux) +
     &                         TVadv(i,j,k,iflux) )
                enddo
              enddo
            enddo
          enddo
          do iflux=1,NT
            do k =1,N
              do j=JstrR,JendR
                do i=IstrR,IendR
                  TForc_avg(i,j,k,iflux) = cff *
     &                       ( TForc_avg(i,j,k,iflux) +
     &                         TForc(i,j,k,iflux) )
                enddo
              enddo
            enddo
          enddo

          do iflux=1,NT
            do k =1,N
              do j=JstrR,JendR
                do i=IstrR,IendR
                  THmix_avg(i,j,k,iflux) = cff *
     &                       ( THmix_avg(i,j,k,iflux) +
     &                         THmix(i,j,k,iflux) )
                enddo
              enddo
            enddo
          enddo
          do iflux=1,NT
            do k=1,N
              do j=JstrR,JendR
                do i=IstrR,IendR 
                  TVmix_avg(i,j,k,iflux) = cff *
     &                       ( TVmix_avg(i,j,k,iflux) +
     &                         TVmix(i,j,k,iflux) )
                enddo
              enddo
            enddo
          enddo

          do iflux=1,NT
            do k=1,N
              do j=JstrR,JendR
                do i=IstrR,IendR 
                  Trate_avg(i,j,k,iflux) = cff *
     &                    ( Trate_avg(i,j,k,iflux) + 
     &                      Trate(i,j,k,iflux) )
                enddo
              enddo
            enddo
          enddo

        endif
      endif       !<-- iic.gt.ntsdiags_ts_avg

      return
      end
#else /* DIAGNOSTICS_TS && AVERAGES */
      subroutine set_diags_ts_avg_empty
      end
#endif /* DIAGNOSTICS_TS && AVERAGES */
