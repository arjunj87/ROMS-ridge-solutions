#include "cppdefs.h"
#if (defined DIAGNOSTICS_UV && defined AVERAGES)
 
      subroutine set_diags_uv_avg (tile)
      implicit none
      integer tile
# include "param.h"
# include "compute_tile_bounds.h"
      call set_diags_uv_avg_tile (istr,iend,jstr,jend)
      return
      end
 
      subroutine set_diags_uv_avg_tile (istr,iend,jstr,jend)
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
#include "diags_uv.h"

!
#include "compute_auxiliary_bounds.h"
!
      ilc=1+iic-ntsdiags_uv_avg   ! number of time step since restart
!
! calculation of averaged fluxes will only be performed if ntsdiags_uv_avg
! is a positive number
!


      if (iic.gt.ntsdiags_uv_avg) then
        if (mod(iic-ntsdiags_uv_avg,nwrtdiags_uv_avg).eq.1) then

          if (ZEROTH_TILE) then
            timediags_uv_avg=time
      MPI_master_only write(*,*) 'started averaging',iic,
     &                       ntsdiags_uv_avg,nwrtdiags_uv_avg
          endif


          do iflux=1,2
            do k=1,N
              do j=JstrR,JendR
                do i=IstrR,IendR 
                  MXadv_avg(i,j,k,iflux)=MXadv(i,j,k,iflux)  
                enddo
              enddo
            enddo
          enddo
          do iflux=1,2
            do k=1,N
              do j=JstrR,JendR
                do i=IstrR,IendR 
                  MYadv_avg(i,j,k,iflux)=MYadv(i,j,k,iflux)
                enddo
              enddo
            enddo
          enddo
          do iflux=1,2
            do k=1,N
              do j=JstrR,JendR
                do i=IstrR,IendR 
                  MVadv_avg(i,j,k,iflux)=MVadv(i,j,k,iflux)
                enddo
              enddo
            enddo
          enddo
          do iflux=1,2
            do k=1,N
              do j=JstrR,JendR
                do i=IstrR,IendR
                  MCor_avg(i,j,k,iflux)=MCor(i,j,k,iflux)
                enddo
              enddo
            enddo
          enddo
          do iflux=1,2
            do k=1,N
              do j=JstrR,JendR
                do i=IstrR,IendR
                  MPrsgrd_avg(i,j,k,iflux)=MPrsgrd(i,j,k,iflux)
                enddo
              enddo
            enddo
          enddo
          do iflux=1,2
            do k=1,N
              do j=JstrR,JendR
                do i=IstrR,IendR
                  MHmix_avg(i,j,k,iflux)=MHmix(i,j,k,iflux)
                enddo
              enddo
            enddo
          enddo
          do iflux=1,2
            do k=1,N
              do j=JstrR,JendR
                do i=IstrR,IendR 
                  MVmix_avg(i,j,k,iflux)=MVmix(i,j,k,iflux)
                end do
              end do
            end do
          end do
          do iflux=1,2
            do k=1,N
              do j=JstrR,JendR
                do i=IstrR,IendR 
                  MVmix2_avg(i,j,k,iflux)=MVmix2(i,j,k,iflux)
                end do
              end do
            end do
          end do

          do iflux=1,2
            do k=1,N
              do j=JstrR,JendR
                do i=IstrR,IendR 
                  Mrate_avg(i,j,k,iflux)=Mrate(i,j,k,iflux)
                enddo
              enddo
            enddo
          enddo



        elseif (mod(iic-ntsdiags_uv_avg,nwrtdiags_uv_avg).gt.1) then

          if (ZEROTH_TILE) timediags_uv_avg=timediags_uv_avg+time

          do iflux=1,2
            do k=1,N
              do j=JstrR,JendR
                do i=IstrR,IendR 
                  MXadv_avg(i,j,k,iflux) = 
     &                       MXadv_avg(i,j,k,iflux) +
     &                       MXadv(i,j,k,iflux)
                enddo
              enddo
            enddo
          enddo
          do iflux=1,2
            do k=1,N
              do j=JstrR,JendR
                do i=IstrR,IendR 
                  MYadv_avg(i,j,k,iflux) = 
     &                       MYadv_avg(i,j,k,iflux) +
     &                       MYadv(i,j,k,iflux)
                enddo
              enddo
            enddo
          enddo
          do iflux=1,2
            do k =1,N
              do j=JstrR,JendR
                do i=IstrR,IendR 
                  MVadv_avg(i,j,k,iflux) = 
     &                       MVadv_avg(i,j,k,iflux) +
     &                       MVadv(i,j,k,iflux)
                enddo
              enddo
            enddo
          enddo
          do iflux=1,2
            do k=1,N
              do j=JstrR,JendR
                do i=IstrR,IendR
                  MCor_avg(i,j,k,iflux) =
     &                       MCor_avg(i,j,k,iflux) +
     &                       MCor(i,j,k,iflux)
                enddo
              enddo
            enddo
          enddo
          do iflux=1,2
            do k=1,N
              do j=JstrR,JendR
                do i=IstrR,IendR
                  MPrsgrd_avg(i,j,k,iflux) =
     &                       MPrsgrd_avg(i,j,k,iflux) +
     &                       MPrsgrd(i,j,k,iflux)
                enddo
              enddo
            enddo
          enddo
          do iflux=1,2
            do k=1,N
              do j=JstrR,JendR
                do i=IstrR,IendR
                  MHmix_avg(i,j,k,iflux) =
     &                       MHmix_avg(i,j,k,iflux) +
     &                       MHmix(i,j,k,iflux)
                enddo
              enddo
            enddo
          enddo
          do iflux=1,2
            do k=1,N
              do j=JstrR,JendR
                do i=IstrR,IendR 
                  MVmix_avg(i,j,k,iflux) = 
     &                       MVmix_avg(i,j,k,iflux) +
     &                       MVmix(i,j,k,iflux)
                enddo
              enddo
            enddo
          enddo
          do iflux=1,2
            do k=1,N
              do j=JstrR,JendR
                do i=IstrR,IendR 
                  MVmix2_avg(i,j,k,iflux) = 
     &                       MVmix2_avg(i,j,k,iflux) +
     &                       MVmix2(i,j,k,iflux)
                enddo
              enddo
            enddo
          enddo

          do iflux=1,2
            do k=1,N
              do j=JstrR,JendR
                do i=IstrR,IendR 
                  Mrate_avg(i,j,k,iflux) = 
     &                    Mrate_avg(i,j,k,iflux) + 
     &                    Mrate(i,j,k,iflux)  
                enddo
              enddo
            enddo
          enddo


        elseif (mod(iic-ntsdiags_uv_avg,nwrtdiags_uv_avg).eq.0) then

          cff=1./float(nwrtdiags_uv_avg)
          if (ZEROTH_TILE) then
                   timediags_uv_avg=cff*(timediags_uv_avg+time)
      MPI_master_only write(*,*) 'finish avera. uv',iic,
     &                       ntsdiags_uv_avg,nwrtdiags_uv_avg

          endif


          do iflux=1,2
            do k=1,N
              do j=JstrR,JendR
                do i=IstrR,IendR 
                  MXadv_avg(i,j,k,iflux) = cff *
     &                       ( MXadv_avg(i,j,k,iflux) +
     &                         MXadv(i,j,k,iflux) )
                   enddo
                enddo
             enddo
          enddo
          do iflux=1,2
            do k=1,N
              do j=JstrR,JendR
                do i=IstrR,IendR 
                  MYadv_avg(i,j,k,iflux) = cff *
     &                       ( MYadv_avg(i,j,k,iflux) +
     &                         MYadv(i,j,k,iflux) )
                enddo
              enddo
            enddo
          enddo
          do iflux=1,2
            do k=1,N
              do j=JstrR,JendR
                do i=IstrR,IendR 
                  MVadv_avg(i,j,k,iflux) = cff *
     &                       ( MVadv_avg(i,j,k,iflux) +
     &                         MVadv(i,j,k,iflux) )
                enddo
              enddo
            enddo
          enddo
          do iflux=1,2
            do k =1,N
              do j=JstrR,JendR
                do i=IstrR,IendR
                  MCor_avg(i,j,k,iflux) = cff *
     &                       ( MCor_avg(i,j,k,iflux) +
     &                         MCor(i,j,k,iflux) )
                enddo
              enddo
            enddo
          enddo
          do iflux=1,2
            do k =1,N
              do j=JstrR,JendR
                do i=IstrR,IendR
                  MPrsgrd_avg(i,j,k,iflux) = cff *
     &                       ( MPrsgrd_avg(i,j,k,iflux) +
     &                         MPrsgrd(i,j,k,iflux) )
                enddo
              enddo
            enddo
          enddo
          do iflux=1,2
            do k =1,N
              do j=JstrR,JendR
                do i=IstrR,IendR
                  MHmix_avg(i,j,k,iflux) = cff *
     &                       ( MHmix_avg(i,j,k,iflux) +
     &                         MHmix(i,j,k,iflux) )
                enddo
              enddo
            enddo
          enddo
          do iflux=1,2
            do k=1,N
              do j=JstrR,JendR
                do i=IstrR,IendR 
                  MVmix_avg(i,j,k,iflux) = cff *
     &                       ( MVmix_avg(i,j,k,iflux) +
     &                         MVmix(i,j,k,iflux) )
                enddo
              enddo
            enddo
          enddo
          do iflux=1,2
            do k=1,N
              do j=JstrR,JendR
                do i=IstrR,IendR 
                  MVmix2_avg(i,j,k,iflux) = cff *
     &                       ( MVmix2_avg(i,j,k,iflux) +
     &                         MVmix2(i,j,k,iflux) )
                enddo
              enddo
            enddo
          enddo
          do iflux=1,2
            do k=1,N
              do j=JstrR,JendR
                do i=IstrR,IendR 
                  Mrate_avg(i,j,k,iflux) = cff *
     &                    ( Mrate_avg(i,j,k,iflux) + 
     &                      Mrate(i,j,k,iflux) )
                enddo
              enddo
            enddo
          enddo

        endif
      endif       !<-- iic.gt.ntsdiags_uv_avg

      return
      end
#else /* DIAGNOSTICS_UV && AVERAGES */
      subroutine set_diags_uv_avg_empty
      end
#endif /* DIAGNOSTICS_UV && AVERAGES */
