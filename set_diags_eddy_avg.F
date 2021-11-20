#include "cppdefs.h"
#if (defined DIAGNOSTICS_EDDY && defined AVERAGES)
 
      subroutine set_diags_eddy_avg (tile)
      implicit none
      integer tile
# include "param.h"
# include "compute_tile_bounds.h"
      call set_diags_eddy_avg_tile (istr,iend,jstr,jend)
      return
      end
 
      subroutine set_diags_eddy_avg_tile (istr,iend,jstr,jend)
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
     &                       ,imin,imax,jmin,jmax
# ifdef SOLVE3D
     &                       , itrc, k
# endif
      real cff, GRho
# include "param.h"
      real cff1,cff2,cff3
# include "work.h"
# include "scalars.h"
# include "ncvars.h"
# include "grid.h"
# include "ocean2d.h"
# include "ocean3d.h"
# include "eos_vars.h"
# include "averages.h"
#include "diags_eddy.h"

# include "private_scratch.h"
!

# include "compute_extended_bounds.h"
# ifdef EW_PERIODIC
      imin=max(istrR, istr-1)        ! The logic here is somewhat
      imax=min(iendR, iend+1)        ! opposite to the standard
#else
      imin=istr                      ! "auxiliary" bounds: extend 
      imax=iend                      ! the range of computation of  
# endif
# ifdef NS_PERIODIC
      jmin=max(jstrR, jstr-1)        ! "Wvlc", by one row of points
      jmax=min(jendR, jend+1)        ! on the side, but only in the
# else
      jmin=jstr                      ! case of periodic direction.
      jmax=jend
# endif
!
      ilc=1+iic-ntsdiags_eddy_avg   ! number of time step since restart
!
! calculation of averaged fluxes will only be performed if ntsdiags_eddy_avg
! is a positive number


      if (iic.gt.ntsdiags_eddy_avg) then

!###############################################

       ! Compute true vertical velocity (m/s).
         call Wvlcty_tile (istr,iend,jstr,jend, work, A2d(1,1),
     &                                   A2d(1,1), A2d(1,2))

!###############################################

        if (mod(iic-ntsdiags_eddy_avg,nwrtdiags_eddy_avg).eq.1) then

          if (ZEROTH_TILE) then
            timediags_eddy_avg=time
      MPI_master_only write(*,*) 'started averaging',iic,
     &                       ntsdiags_eddy_avg,nwrtdiags_eddy_avg
          endif

          do k=1,N
              do j=Jmin,jmax
                do i=imin,imax 

                  GRho=-1*g/rho0

                  cff1 = 0.5*(u(i+1,j,k,nstp)+u(i,j,k,nstp))
                  cff2 = 0.5*(v(i,j+1,k,nstp)+v(i,j,k,nstp))
#  ifdef SPLIT_EOS
                  cff3 =rho1(i,j,k)+qp1(i,j,k)
     &                          *(z_w(i,j,N)-z_r(i,j,k))
#  else
                  cff3 =rho(i,j,k)
#  endif
                  eddyuu_avg(i,j,k)=cff1**2
                  eddyvv_avg(i,j,k)=cff2**2
                  eddyuv_avg(i,j,k)=cff1*cff2
                  eddyub_avg(i,j,k)=GRho*cff1*cff3
                  eddyvb_avg(i,j,k)=GRho*cff2*cff3
                  eddywb_avg(i,j,k)=GRho*work(i,j,k)*cff3
                  eddynotused1_avg(i,j,k)=work(i,j,k)*cff1
                  eddynotused2_avg(i,j,k)=work(i,j,k)*cff2
                  eddynotused3_avg(i,j,k)=0.

                enddo
              enddo
            enddo



        elseif (mod(iic-ntsdiags_eddy_avg,nwrtdiags_eddy_avg).gt.1) then

          if (ZEROTH_TILE) timediags_eddy_avg=timediags_eddy_avg+time

          do k=1,N
              do j=Jmin,jmax
                do i=imin,imax 

                  cff1 = 0.5*(u(i+1,j,k,nstp)+u(i,j,k,nstp))
                  cff2 = 0.5*(v(i,j+1,k,nstp)+v(i,j,k,nstp))
#  ifdef SPLIT_EOS
                  cff3 =rho1(i,j,k)+qp1(i,j,k)
     &                          *(z_w(i,j,N)-z_r(i,j,k))
#  else
                  cff3 =rho(i,j,k)
#  endif
                  eddyuu_avg(i,j,k)=eddyuu_avg(i,j,k)
     &                                   +cff1**2
                  eddyvv_avg(i,j,k)=eddyvv_avg(i,j,k)
     &                                   +cff2**2
                  eddyuv_avg(i,j,k)=eddyuv_avg(i,j,k)
     &                                   +cff1*cff2
                  eddyub_avg(i,j,k)=eddyub_avg(i,j,k)
     &                          +GRho*cff1*cff3
                  eddyvb_avg(i,j,k)=eddyvb_avg(i,j,k)
     &                          +GRho*cff2*cff3
                  eddywb_avg(i,j,k)=eddywb_avg(i,j,k)
     &                          +GRho*work(i,j,k)*cff3
                  eddynotused1_avg(i,j,k)=eddynotused1_avg(i,j,k)
     &                          +work(i,j,k)*cff1
                  eddynotused2_avg(i,j,k)=eddynotused2_avg(i,j,k)
     &                          +work(i,j,k)*cff2
                  eddynotused3_avg(i,j,k)=0.

                enddo
              enddo
            enddo



        elseif (mod(iic-ntsdiags_eddy_avg,nwrtdiags_eddy_avg).eq.0) then

          cff=1./float(nwrtdiags_eddy_avg)
          if (ZEROTH_TILE) then
                   timediags_eddy_avg=cff*(timediags_eddy_avg+time)
      MPI_master_only write(*,*) 'finish avera. eddy',iic,
     &                       ntsdiags_eddy_avg,nwrtdiags_eddy_avg

          endif


          do k=1,N
              do j=Jmin,jmax
                do i=imin,imax 


                  cff1 = 0.5*(u(i+1,j,k,nstp)+u(i,j,k,nstp))
                  cff2 = 0.5*(v(i,j+1,k,nstp)+v(i,j,k,nstp))

#  ifdef SPLIT_EOS
                  cff3 =rho1(i,j,k)+qp1(i,j,k)
     &                          *(z_w(i,j,N)-z_r(i,j,k))
#  else
                  cff3 =rho(i,j,k)
#  endif
                  eddyuu_avg(i,j,k)= cff *
     &                       ( eddyuu_avg(i,j,k)
     &                         +cff1**2 )
                  eddyvv_avg(i,j,k)= cff *
     &                       ( eddyvv_avg(i,j,k)
     &                         +cff2**2 )
                  eddyuv_avg(i,j,k)= cff *
     &                       ( eddyuv_avg(i,j,k)
     &                         +cff1*cff2 )
                  eddyub_avg(i,j,k)= cff *
     &                       ( eddyub_avg(i,j,k)
     &                     +GRho*cff1*cff3 )
                  eddyvb_avg(i,j,k)= cff *
     &                       ( eddyvb_avg(i,j,k)
     &                     +GRho*cff2*cff3 )
                  eddywb_avg(i,j,k)= cff *
     &                       ( eddywb_avg(i,j,k)
     &                     +GRho*work(i,j,k)*cff3 )
                  eddynotused1_avg(i,j,k)= cff *
     &                       ( eddynotused1_avg(i,j,k)
     &                          +work(i,j,k)*cff1 )
                  eddynotused2_avg(i,j,k)= cff *
     &                       ( eddynotused2_avg(i,j,k)
     &                          +work(i,j,k)*cff2 )
                  eddynotused3_avg(i,j,k)=0.

                enddo
              enddo
            enddo

        endif
      endif       !<-- iic.gt.ntsdiags_eddy_avg

      return
      end
#else /* DIAGNOSTICS_EDDY && AVERAGES */
      subroutine set_diags_eddy_avg_empty
      end
#endif /* DIAGNOSTICS_EDDY && AVERAGES */
