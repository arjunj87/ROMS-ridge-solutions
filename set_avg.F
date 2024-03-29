#include "cppdefs.h"
#ifdef AVERAGES
 
      subroutine set_avg (tile)
      implicit none
      integer tile
# include "param.h"
# include "compute_tile_bounds.h"
      call set_avg_tile (istr,iend,jstr,jend)
      return
      end
 
      subroutine set_avg_tile (istr,iend,jstr,jend)
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
      integer istr,iend,jstr,jend, i,j
# ifdef SOLVE3D
     &                       , itrc, k
# endif
      real cff
# include "param.h"
# include "scalars.h"
# include "ncvars.h"
# include "grid.h"
# include "work.h"
# include "private_scratch.h"
# include "ocean2d.h"
# include "ocean3d.h"
# include "eos_vars.h"
# include "mixing.h"
# include "averages.h"
#ifdef BULK_FLUX
# include "forces.h"
#endif
!
#include "compute_auxiliary_bounds.h"
!

      if (iic.gt.ntsavg) then
        if (mod(iic-ntsavg,navg).eq.1) then
          if (ZEROTH_TILE) then
            time_avg=time
      MPI_master_only write(*,*) 'started averaging',iic,ntsavg,navg
          endif
          if (wrtavg(indxZ)) then
            do j=jstrR,jendR
              do i=istrR,iendR
                zeta_avg(i,j)=zeta(i,j,knew)
              enddo
            enddo
          endif
          if (wrtavg(indxUb)) then
            do j=jstrR,jendR
              do i=istrR,iendR
                ubar_avg(i,j)=ubar(i,j,knew)
              enddo
            enddo
          endif
          if (wrtavg(indxVb)) then
            do j=jstrR,jendR
              do i=istrR,iendR
                vbar_avg(i,j)=vbar(i,j,knew)
              enddo
            enddo
          endif
# ifdef SOLVE3D
          if (wrtavg(indxU)) then
            do k=1,N
              do j=jstrR,jendR
                do i=istrR,iendR
                  u_avg(i,j,k)=u(i,j,k,nstp)
#ifdef WRITE_SXSY
                  Sx_avg(i,j,k)=Sx_u(i,j,k)
#endif
                enddo
              enddo
            enddo
          endif
          if (wrtavg(indxV)) then
            do k=1,N
              do j=jstrR,jendR
                do i=istrR,iendR
                  v_avg(i,j,k)=v(i,j,k,nstp)
#ifdef WRITE_SXSY
                  Sy_avg(i,j,k)=Sy_v(i,j,k)  
#endif
                enddo
              enddo
            enddo
          endif
          do itrc=1,NT
           if (wrtavg(indxT+itrc-1)) then
             do k=1,N
               do j=jstrR,jendR
                 do i=istrR,iendR
                   t_avg(i,j,k,itrc)=t(i,j,k,nstp,itrc)
                 enddo
               enddo
             enddo
           endif
          enddo
          if (wrtavg(indxR)) then
            do k=1,N
              do j=jstrR,jendR
                do i=istrR,iendR
#  ifdef SPLIT_EOS
                  rho_avg(i,j,k)=rho1(i,j,k)+qp1(i,j,k)
     &                          *(z_w(i,j,N)-z_r(i,j,k))
#  else
                  rho_avg(i,j,k)=rho(i,j,k)
#  endif
                enddo
              enddo
            enddo
          endif
          if (wrtavg(indxO)) then
            do k=1,N
              do j=jstrR,jendR
                do i=istrR,iendR
                  w_avg(i,j,k)=w(i,j,k)
                enddo
              enddo
            enddo
          endif

          if (wrtavg(indxW)) then
!###############################################
       ! Compute true vertical velocity (m/s).
         call Wvlcty_tile (istr,iend,jstr,jend, work, A2d(1,1),
     &                                   A2d(1,1), A2d(1,2))
!###############################################
            do k=1,N
              do j=jstrR,jendR
                do i=istrR,iendR
                  wvlc_avg(i,j,k)=work(i,j,k)
                enddo
              enddo
            enddo
          endif
          if (wrtavg(indxAkv)) then
            do k=0,N
              do j=jstrR,jendR
                do i=istrR,iendR
                  akv_avg(i,j,k)=Akv(i,j,k)
                enddo
              enddo
            enddo
          endif
          if (wrtavg(indxAkt)) then
            do k=0,N
              do j=jstrR,jendR
                do i=istrR,iendR
                  akt_avg(i,j,k)=Akt(i,j,k,itemp)
                enddo
              enddo
            enddo
          endif
#  ifdef SALINITY
          if (wrtavg(indxAks)) then
            do k=0,N
              do j=jstrR,jendR
                do i=istrR,iendR
                  aks_avg(i,j,k)=Akt(i,j,k,isalt)
                enddo
              enddo
            enddo
          endif
#  endif
          
#  ifdef BULK_FLUX
            if (wrtavg(indxSustr_blk)) then
            do j=JstrR,JendR
              do i=IstrR,IendR
                sustr_blk_avg(i,j)=sustr_blk(i,j)
              enddo
            enddo
            endif
            if (wrtavg(indxSvstr_blk)) then
            do j=JstrR,JendR
              do i=IstrR,IendR
                svstr_blk_avg(i,j)=svstr_blk(i,j)
              enddo
            enddo
            endif
            if (wrtavg(indxShflx_rsw)) then
            do j=JstrR,JendR
              do i=IstrR,IendR
                shflx_rsw_avg(i,j)=shflx_rsw(i,j)
              enddo
            enddo
            endif
            if (wrtavg(indxShflx_rlw)) then
            do j=JstrR,JendR
              do i=IstrR,IendR
                shflx_rlw_avg(i,j)=shflx_rlw(i,j)
              enddo
            enddo
            endif
            if (wrtavg(indxShflx_lat)) then
            do j=JstrR,JendR
              do i=IstrR,IendR
                shflx_lat_avg(i,j)=shflx_lat(i,j)
              enddo
            enddo
            endif
            if (wrtavg(indxShflx_sen)) then
            do j=JstrR,JendR
              do i=IstrR,IendR
                shflx_sen_avg(i,j)=shflx_sen(i,j)
              enddo
            enddo
            endif
# endif                    
# endif
        elseif (mod(iic-ntsavg,navg).gt.1) then
          if (ZEROTH_TILE) time_avg=time_avg+time
          if (wrtavg(indxZ)) then
            do j=jstrR,jendR
              do i=istrR,iendR
                zeta_avg(i,j)=zeta_avg(i,j)+zeta(i,j,knew)
              enddo
            enddo
          endif
          if (wrtavg(indxUb)) then
            do j=jstrR,jendR
              do i=istrR,iendR
                ubar_avg(i,j)=ubar_avg(i,j)+ubar(i,j,knew)
              enddo
            enddo
          endif
          if (wrtavg(indxVb)) then
            do j=jstrR,jendR
              do i=istrR,iendR
                vbar_avg(i,j)=vbar_avg(i,j)+vbar(i,j,knew)
              enddo
            enddo
          endif
# ifdef SOLVE3D
          if (wrtavg(indxU)) then
            do k=1,N
              do j=jstrR,jendR
                do i=istrR,iendR
                  u_avg(i,j,k)=u_avg(i,j,k)+u(i,j,k,nstp)
#ifdef WRITE_SXSY
                  Sx_avg(i,j,k)=Sx_avg(i,j,k)
     &                           +Sx_u(i,j,k)
#endif
                enddo
              enddo
            enddo
          endif
          if (wrtavg(indxV)) then
            do k=1,N
              do j=jstrR,jendR
                do i=istrR,iendR
                  v_avg(i,j,k)=v_avg(i,j,k)+v(i,j,k,nstp)
#ifdef WRITE_SXSY
                  Sy_avg(i,j,k)=Sy_avg(i,j,k)
     &                           +Sy_v(i,j,k)  
#endif
                enddo
              enddo
            enddo
          endif
          do itrc=1,NT
            if (wrtavg(indxT+itrc-1)) then
              do k=1,N
                do j=jstrR,jendR
                  do i=istrR,iendR
                    t_avg(i,j,k,itrc)=t_avg(i,j,k,itrc)
     &                              +t(i,j,k,nstp,itrc)
                  enddo
                enddo
              enddo
            endif
          enddo         
          
          if (wrtavg(indxR)) then
            do k=1,N
              do j=jstrR,jendR
                do i=istrR,iendR
#  ifdef SPLIT_EOS
                  rho_avg(i,j,k)=rho_avg(i,j,k) + rho1(i,j,k)
     &                     +qp1(i,j,k)*(z_w(i,j,N)-z_r(i,j,k))
#  else
                  rho_avg(i,j,k)=rho_avg(i,j,k)+rho(i,j,k)
#  endif
                enddo
              enddo
            enddo
          endif
          if (wrtavg(indxO)) then
            do k=1,N
              do j=jstrR,jendR
                do i=istrR,iendR
                  w_avg(i,j,k)=w_avg(i,j,k)+w(i,j,k)
                enddo
              enddo
            enddo
          endif
          if (wrtavg(indxW)) then
!###############################################
       ! Compute true vertical velocity (m/s).
         call Wvlcty_tile (istr,iend,jstr,jend, work, A2d(1,1),
     &                                   A2d(1,1), A2d(1,2))
!###############################################
            do k=1,N
              do j=jstrR,jendR
                do i=istrR,iendR
                  wvlc_avg(i,j,k)=wvlc_avg(i,j,k)+work(i,j,k)
                enddo
              enddo
            enddo
          endif
          if (wrtavg(indxAkv)) then
            do k=0,N
              do j=jstrR,jendR
                do i=istrR,iendR
                  akv_avg(i,j,k)=akv_avg(i,j,k)+Akv(i,j,k)
                enddo
              enddo
            enddo
          endif
          if (wrtavg(indxAkt)) then
            do k=0,N
              do j=jstrR,jendR
                do i=istrR,iendR
                  akt_avg(i,j,k)=akt_avg(i,j,k)+Akt(i,j,k,itemp)
                enddo
              enddo
            enddo
          endif
#  ifdef SALINITY
          if (wrtavg(indxAks)) then
            do k=0,N
              do j=jstrR,jendR
                do i=istrR,iendR
                  aks_avg(i,j,k)=aks_avg(i,j,k)+Akt(i,j,k,isalt)
                enddo
              enddo
            enddo
          endif
#  endif
#  ifdef LMD_KPP
          if (wrtavg(indxHbl)) then
            do j=jstrR,jendR
              do i=istrR,iendR
                hbl_avg(i,j)=hbl_avg(i,j)+hbls(i,j,nstp)
# ifdef LMD_BKPP
                hbbl_avg(i,j)=hbbl_avg(i,j)+hbbls(i,j,nstp)
# endif                
              enddo
            enddo
          endif
#  endif
#  ifdef BULK_FLUX
            if (wrtavg(indxSustr_blk)) then
            do j=JstrR,JendR
              do i=IstrR,IendR
                sustr_blk_avg(i,j)=sustr_blk_avg(i,j)+sustr_blk(i,j)
              enddo
            enddo
            endif
            if (wrtavg(indxSvstr_blk)) then
            do j=JstrR,JendR
              do i=IstrR,IendR
                svstr_blk_avg(i,j)=svstr_blk_avg(i,j)+svstr_blk(i,j)
              enddo
            enddo
            endif
            if (wrtavg(indxShflx_rsw)) then
            do j=JstrR,JendR
              do i=IstrR,IendR
                shflx_rsw_avg(i,j)=shflx_rsw_avg(i,j)+shflx_rsw(i,j)
              enddo
            enddo
            endif
            if (wrtavg(indxShflx_rlw)) then
            do j=JstrR,JendR
              do i=IstrR,IendR
                shflx_rlw_avg(i,j)=shflx_rlw_avg(i,j)+shflx_rlw(i,j)
              enddo
            enddo
            endif
            if (wrtavg(indxShflx_lat)) then
            do j=JstrR,JendR
              do i=IstrR,IendR
                shflx_lat_avg(i,j)=shflx_lat_avg(i,j)+shflx_lat(i,j)
              enddo
            enddo
            endif
            if (wrtavg(indxShflx_sen)) then
            do j=JstrR,JendR
              do i=IstrR,IendR
                shflx_sen_avg(i,j)=shflx_sen_avg(i,j)+shflx_sen(i,j)
              enddo
            enddo
            endif
#  endif

# endif /* SOLVE3D */
        elseif (mod(iic-ntsavg,navg).eq.0) then
          cff=1./float(navg)
          if (ZEROTH_TILE) then
             time_avg=cff*(time_avg+time)
      MPI_master_only write(*,*) 'finish averaging',iic,ntsavg,navg
          endif
          if (wrtavg(indxZ)) then
            do j=jstrR,jendR
              do i=istrR,iendR
                zeta_avg(i,j)=cff*( zeta_avg(i,j)+zeta(i,j,knew))
              enddo
            enddo
          endif
          if (wrtavg(indxUb)) then
            do j=jstrR,jendR
              do i=istrR,iendR
                ubar_avg(i,j)=cff*( ubar_avg(i,j)+ubar(i,j,knew))
              enddo
            enddo
          endif
          if (wrtavg(indxVb)) then
            do j=jstrR,jendR
              do i=istrR,iendR
                vbar_avg(i,j)=cff*( vbar_avg(i,j)+vbar(i,j,knew))
              enddo
            enddo
          endif
# ifdef SOLVE3D
          if (wrtavg(indxU)) then
            do k=1,N
              do j=jstrR,jendR
                do i=istrR,iendR
                  u_avg(i,j,k)=cff*(u_avg(i,j,k)+u(i,j,k,nstp))
#ifdef WRITE_SXSY
                  Sx_avg(i,j,k)=cff*( Sx_avg(i,j,k)
     &                                    +  Sx_u(i,j,k))
#endif
                enddo
              enddo
            enddo
          endif
          if (wrtavg(indxV)) then
            do k=1,N
              do j=jstrR,jendR
                do i=istrR,iendR
                  v_avg(i,j,k)=cff*(v_avg(i,j,k)+v(i,j,k,nstp))
#ifdef WRITE_SXSY
                  Sy_avg(i,j,k)=cff*( Sy_avg(i,j,k)
     &                                    +  Sy_v(i,j,k)) 
#endif
                enddo
              enddo
            enddo
          endif
          do itrc=1,NT
            if (wrtavg(indxT+itrc-1)) then
              do k=1,N
                do j=jstrR,jendR
                  do i=istrR,iendR
                    t_avg(i,j,k,itrc)=cff*( t_avg(i,j,k,itrc)
     &                                    +t(i,j,k,nstp,itrc))
                  enddo
                enddo
              enddo
            endif
          enddo
                    
          
          if (wrtavg(indxR)) then
            do k=1,N
              do j=jstrR,jendR
                do i=istrR,iendR
#  ifdef SPLIT_EOS
                  rho_avg(i,j,k)=cff*( rho_avg(i,j,k) + rho1(i,j,k)
     &                          +qp1(i,j,k)*(z_w(i,j,N)-z_r(i,j,k)))
#  else
                  rho_avg(i,j,k)=cff*(rho_avg(i,j,k)+rho(i,j,k))
#  endif
                enddo
              enddo
            enddo
          endif
          if (wrtavg(indxO)) then
            do k=1,N
              do j=jstrR,jendR
                do i=istrR,iendR
                  w_avg(i,j,k)=cff*pm(i,j)*pn(i,j)*( w_avg(i,j,k)
     &                                                  +w(i,j,k))
                enddo
              enddo
            enddo
          endif
          if (wrtavg(indxW)) then
!###############################################
       ! Compute true vertical velocity (m/s).
         call Wvlcty_tile (istr,iend,jstr,jend, work, A2d(1,1),
     &                                   A2d(1,1), A2d(1,2))
!###############################################
            do k=1,N
              do j=jstrR,jendR
                do i=istrR,iendR
                  wvlc_avg(i,j,k)=cff*(wvlc_avg(i,j,k)+work(i,j,k))
                enddo
              enddo
            enddo
          endif
          if (wrtavg(indxAkv)) then
            do k=0,N
              do j=jstrR,jendR
                do i=istrR,iendR
                  akv_avg(i,j,k)=cff*(akv_avg(i,j,k)+Akv(i,j,k))
                enddo
              enddo
            enddo
          endif
          if (wrtavg(indxAkt)) then
            do k=0,N
              do j=jstrR,jendR
                do i=istrR,iendR
                  akt_avg(i,j,k)=cff*(akt_avg(i,j,k)+Akt(i,j,k,itemp))
                enddo
              enddo
            enddo
          endif
#  ifdef SALINITY
          if (wrtavg(indxAks)) then
            do k=0,N
              do j=jstrR,jendR
                do i=istrR,iendR
                  aks_avg(i,j,k)=cff*(aks_avg(i,j,k)+Akt(i,j,k,isalt))
                enddo
              enddo
            enddo
          endif
#  endif
#  ifdef LMD_KPP
          if (wrtavg(indxHbl)) then
            do j=jstrR,jendR
              do i=istrR,iendR
                hbl_avg(i,j)=cff*(hbl_avg(i,j)+hbls(i,j,nstp))
# ifdef LMD_BKPP
                hbbl_avg(i,j)=cff*(hbbl_avg(i,j)+hbbls(i,j,nstp))                
#endif                
              enddo
            enddo
          endif
#  endif


#  ifdef BULK_FLUX
           if (wrtavg(indxSustr_blk)) then
            do j=JstrR,JendR
              do i=IstrR,IendR
                sustr_blk_avg(i,j)=cff*( sustr_blk_avg(i,j)+
     &                                   sustr_blk(i,j) )
              enddo
            enddo
           endif
           if (wrtavg(indxSvstr_blk)) then
            do j=JstrR,JendR
              do i=IstrR,IendR
                svstr_blk_avg(i,j)=cff*( svstr_blk_avg(i,j)+
     &                                   svstr_blk(i,j) )
              enddo
            enddo
           endif
           if (wrtavg(indxShflx_rsw)) then
            do j=JstrR,JendR
              do i=IstrR,IendR
                shflx_rsw_avg(i,j)=cff*( shflx_rsw_avg(i,j)+
     &                                   shflx_rsw(i,j) )
              enddo
            enddo
           endif
           if (wrtavg(indxShflx_rlw)) then
            do j=JstrR,JendR
              do i=IstrR,IendR
                shflx_rlw_avg(i,j)=cff*( shflx_rlw_avg(i,j)+
     &                                   shflx_rlw(i,j) )
              enddo
            enddo
           endif
            if (wrtavg(indxShflx_lat)) then
            do j=JstrR,JendR
              do i=IstrR,IendR
                shflx_lat_avg(i,j)=cff*( shflx_lat_avg(i,j)+
     &                                   shflx_lat(i,j) )
              enddo
            enddo
            endif
            if (wrtavg(indxShflx_sen)) then
            do j=JstrR,JendR
              do i=IstrR,IendR
                shflx_sen_avg(i,j)=cff*( shflx_sen_avg(i,j)+
     &                                   shflx_sen(i,j) )
              enddo
            enddo
            endif
#  endif


# endif /* SOLVE3D */
        endif
      endif       !<-- iic.gt.ntsavg
      return
      end
#else
      subroutine set_avg_empty
      end
#endif /* AVERAGES */
