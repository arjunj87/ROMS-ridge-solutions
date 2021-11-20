#include "cppdefs.h"
#ifdef SOLVE3D
 
c--#define CONST_TRACERS

      subroutine step3d_t (tile)
      implicit none
      integer tile
# include "param.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"
      call step3d_t_tile (istr,iend,jstr,jend, A2d(1,1), A2d(1,2),
     &                               A2d(1,3), A2d(1,4), A2d(1,5),
     &                               A2d(1,1), A2d(1,2), A2d(1,3))
      return
      end
 
      subroutine step3d_t_tile (istr,iend,jstr,jend,
     &                   FX,FE,FF,WORK,dTdz,  FC,CF,DC)
      implicit none
# include "param.h"
      integer istr,iend,jstr,jend, itrc,
     &        imin,imax,jmin,jmax, i,j,k, indx
      real FX(PRIVATE_2D_SCRATCH_ARRAY),
     &     FE(PRIVATE_2D_SCRATCH_ARRAY),  cff,
     &     FF(PRIVATE_2D_SCRATCH_ARRAY),
     &   WORK(PRIVATE_2D_SCRATCH_ARRAY),  epsil,
     &   dTdz(PRIVATE_2D_SCRATCH_ARRAY,2),
     &     FC(PRIVATE_1D_SCRATCH_ARRAY,0:N),
     &     CF(PRIVATE_1D_SCRATCH_ARRAY,0:N),
     &     DC(PRIVATE_1D_SCRATCH_ARRAY,0:N)
      real gama     
      parameter (epsil=1.E-16)
# include "grid.h"
# include "ocean3d.h"
# include "forces.h"
# include "mixing.h"
# include "climat.h"
# include "scalars.h"
# include "sources.h"

# ifdef DIAGNOSTICS_TS
# include "diags_ts.h"
# elif defined DIAGNOSTICS_PV
      real  twrk_bot(PRIVATE_2D_SCRATCH_ARRAY,2),
     &      twrk_sfc(PRIVATE_2D_SCRATCH_ARRAY,2)
# endif

# ifdef DIAGNOSTICS_PV
# include "diags_pv.h"
# endif

!
# include "compute_auxiliary_bounds.h"
!
!
      do itrc=1,NT
        do k=1,N
 
# include "compute_horiz_tracer_fluxes.h"
 
          do j=jstr,jend
            do i=istr,iend
              t(i,j,k,nnew  ,itrc)=Hz_bak(i,j,k)*t(i,j,k,nstp,itrc)
     &                     -dt*pm(i,j)*pn(i,j)*( FX(i+1,j)-FX(i,j)
     &                                          +FE(i,j+1)-FE(i,j)
     &                                                           )

# ifdef DIAGNOSTICS_TS
              TXadv(i,j,k,itrc) = -(FX(i+1,j)-FX(i,j))
#  ifdef MASKING
     &              * rmask(i,j)
#  endif
              TYadv(i,j,k,itrc) = -(FE(i,j+1)-FE(i,j))
#  ifdef MASKING
     &              * rmask(i,j)
#  endif
# endif /* DIAGNOSTICS_TS */

            enddo
          enddo                 !-->  discard FX,FE
        enddo
      enddo
      
      




      do j=jstr,jend
        do itrc=1,NT
 
# include "compute_vert_tracer_fluxes.h"
 
          do k=1,N            ! Apply vertical advective fluxes.
            do i=istr,iend
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)-dt*pm(i,j)
     &                             *pn(i,j)*(FC(i,k)-FC(i,k-1))

# ifdef DIAGNOSTICS_TS
              TVadv(i,j,k,itrc) = -(FC(i,k)-FC(i,k-1))
#  ifdef MASKING
     &         * rmask(i,j)
#  endif
# endif /* DIAGNOSTICS_TS */  
            enddo
          enddo               !--> discard FC


!
! Add surface and bottom fluxes
!
#ifdef SEA_ICE_NOFLUX
          do i=istr,iend
            cff = 1.
            if (itrc.eq.itemp .and.
     &          (t(i,j,N,nnew,itrc) .le. -1.8)  ) cff = 0.
            t(i,j,N,nnew,itrc)=t(i,j,N,nnew,itrc)+dt*stflx(i,j,itrc)*cff
# ifdef DIAGNOSTICS_TS
            TForc(i,j,N,itrc)= stflx(i,j,itrc)*cff/(pm(i,j)*pn(i,j))
            do k=N-1,1,-1
                TForc(i,j,k,itrc)= 0.
            enddo
# elif defined DIAGNOSTICS_PV
            pvDiab_bot(i,j,itrc)=0.
            pvDiab_sfc(i,j,itrc)=stflx(i,j,itrc)*cff/(pm(i,j)*pn(i,j))
# endif  /* DIAGNOSTICS_TS */

# if defined DIAGNOSTICS_PV_SUP
            pvDiab_0(i,j,itrc)=0.
            pvDiab_1(i,j,itrc)=0.
            pvDiab_2(i,j,itrc)=0.
            pvDiab_3(i,j,itrc)=stflx(i,j,itrc)*cff
# endif  /* DIAGNOSTICS_PV_SUP */
          enddo
#else
!
! Add surface and bottom fluxes
!
          do i=istr,iend
            t(i,j,N,nnew,itrc)=t(i,j,N,nnew,itrc)+dt*stflx(i,j,itrc)
# ifdef DIAGNOSTICS_TS
            TForc(i,j,N,itrc)= stflx(i,j,itrc)/(pm(i,j)*pn(i,j))
            do k=N-1,1,-1
                TForc(i,j,k,itrc)= 0.
            enddo
# elif defined DIAGNOSTICS_PV
            pvDiab_bot(i,j,itrc)=0.
            pvDiab_sfc(i,j,itrc)=stflx(i,j,itrc)/(pm(i,j)*pn(i,j))
# endif  /* DIAGNOSTICS_TS */

# if defined DIAGNOSTICS_PV_SUP
            pvDiab_0(i,j,itrc)=0.
            pvDiab_1(i,j,itrc)=0.
            pvDiab_2(i,j,itrc)=0.
            pvDiab_3(i,j,itrc)=stflx(i,j,itrc)
# endif  /* DIAGNOSTICS_PV_SUP */

          enddo
#endif
  
# ifdef LMD_KPP
!
! Add the solar radiation flux in temperature equation. Also compute
! the nonlocal transport flux for unstable (convective) forcing
! conditions into matrix DC when using the Large et al. 1994 KPP
! scheme.
!
          if (itrc.eq.itemp) then
            do k=N-1,1,-1
              do i=istr,iend
                cff=srflx(i,j)*swr_frac(i,j,k)
#  ifdef LMD_NONLOCAL
     &                 -ghat(i,j,k)*(stflx(i,j,itemp)-srflx(i,j))
#  endif
                t(i,j,k+1,nnew,itemp)=t(i,j,k+1,nnew,itemp) -dt*cff
                t(i,j,k  ,nnew,itemp)=t(i,j,k  ,nnew,itemp) +dt*cff
# ifdef DIAGNOSTICS_TS
                TForc(i,j,k+1,itrc)= TForc(i,j,k+1,itrc)-cff/(pm(i,j)*pn(i,j))
                TForc(i,j,k,itrc)= TForc(i,j,k,itrc)+cff/(pm(i,j)*pn(i,j))
# elif defined DIAGNOSTICS_PV
               if (k.eq.N-1) then
                pvDiab_sfc(i,j,itrc)=pvDiab_sfc(i,j,itrc)-cff/(pm(i,j)*pn(i,j))
               endif

               if (k.eq.1) then
                pvDiab_bot(i,j,itrc)=cff/(pm(i,j)*pn(i,j))
               endif
# endif  /* DIAGNOSTICS_TS */

              enddo
            enddo

# if defined DIAGNOSTICS_PV_SUP
              do i=istr,iend
                pvDiab_2(i,j,itrc)=srflx(i,j)*swr_frac(i,j,N-1)
#  ifdef LMD_NONLOCAL
     &                 -ghat(i,j,N-1)*(stflx(i,j,itemp)-srflx(i,j))
#  endif
                pvDiab_1(i,j,itrc)=srflx(i,j)*swr_frac(i,j,1)
#  ifdef LMD_NONLOCAL
     &                 -ghat(i,j,1)*(stflx(i,j,itemp)-srflx(i,j))
#  endif
                pvDiab_0(i,j,itrc)=0.
              enddo
# endif  /* DIAGNOSTICS_PV_SUP */


c??
c?? WARNING: the above implies that light (if any) reaching all the
c?? way to the bottom is entirely absorbed within the botom-most grid
c?? box, rather that reflected/scatered back to the water column. It 
c?? is not clear, however, how to deal with this situation correctly
c??
c??         do i=istr,iend
c??           cff=srflx(i,j)*swr_frac(i,j,0)
c??#  ifdef LMD_NONLOCAL
c??  &               -ghat(i,j,0)*(stflx(i,j,itemp)-srflx(i,j))
c??#  endif
c??           t(i,j,1,nnew,itemp)=t(i,j,1,nnew,itemp) -dt*cff
c??         enddo


#  if defined LMD_NONLOCAL && defined SALINITY
          elseif (itrc.eq.isalt) then
            do k=N-1,1,-1
              do i=istr,iend
                cff=-dt*ghat(i,j,k)*stflx(i,j,isalt)
                t(i,j,k+1,nnew,isalt)=t(i,j,k+1,nnew,isalt) -cff
                t(i,j,k  ,nnew,isalt)=t(i,j,k  ,nnew,isalt) +cff
# ifdef DIAGNOSTICS_TS
                TForc(i,j,k+1,itrc)= TForc(i,j,k+1,itrc)-cff/(dt*pm(i,j)*pn(i,j))       
                TForc(i,j,k,itrc)= TForc(i,j,k,itrc)+cff/(dt*pm(i,j)*pn(i,j))
# elif defined DIAGNOSTICS_PV
               if (k.eq.N-1) then
                pvDiab_sfc(i,j,itrc)=pvDiab_sfc(i,j,itrc)-cff/(dt*pm(i,j)*pn(i,j))  
               endif
               if (k.eq.1) then
                pvDiab_bot(i,j,itrc)=cff/(dt*pm(i,j)*pn(i,j))
               endif
# endif  /* DIAGNOSTICS_TS */
              enddo
            enddo

# if defined DIAGNOSTICS_PV_SUP
              do i=istr,iend
                pvDiab_2(i,j,itrc)=-ghat(i,j,N-1)*stflx(i,j,isalt)
                pvDiab_1(i,j,itrc)=-ghat(i,j,1)*stflx(i,j,isalt)
                pvDiab_0(i,j,itrc)=0.
              enddo
# endif  /* DIAGNOSTICS_PV_SUP */



#  endif
          endif
# endif

!
! Perform implicit time step for vertical diffusion,
!
!   dq(k)     1     [         q(k+1)-q(k)             q(k)-q(k-1) ]
!  ------ = ----- * [ Akt(k)* ----------- - Akt(k-1)* ----------- ]
!    dt     Hz(k)   [            dz(k)                   dz(k-1)  ]
!
! where q(k) represents tracer field t(:,:,k,:,itrc). Doing so
! implies solution of a tri-diagonal system
!
!     -FC(k-1)*q_new(k-1) +[Hz(k)+FC(k-1)+FC(k)]*q_new(k)
!                       -FC(k)*q_new(k+1) = Hz(k)*q_old(k)
!
!                dt*Akt(k)
! where FC(k) = ----------- is normalized diffusivity coefficient 
!                  dz(k)
!
! defined at W-points; q_new(k) is the new-time-step (unknown) tracer
! field; q_old(k) is old-time-step tracer (known).  As long as
! vertical diffusivity Akt(k) is nonnegative, the tri-diagonal matrix
! is diagonally dominant which guarantees stability of a Gaussian
! elimination procedure, (e.g., Richtmeyer annd  Morton, 1967).
! Top and bottom boundary conditions are assumed to be no-flux,
! effectively Akt(N)=Akt(0)=0, hence FC(N)=FC(1)=0. This leads to
! equations for top and bottom grid boxes; 
!
!   -FC(N-1)*q_new(N-1) +[Hz(N)+FC(N-1)]*q_new(N) = Hz(N)*q_old(N)
!
!          [Hz(1)+FC(1)]*q_new(1) -FC(1)*q_new(2) = Hz(1)*q_old(1)
!
! The FC(N)=FC(0)=0 boundary conditions does not mean that physical
! boundary conditions are no flux: the forcing fluxes have been
! applied explicitly above.  Instead, the no-flux condition should
! be interpreted as that the implicit step merely redistributes the
! tracer concentration throughout the water column. At this moment
! the content of array t(:,:,:,nnew,itrc) has meaning of Hz*tracer.
! After the implicit step it becomes just tracer. 
!
# ifdef DIAGNOSTICS_TS
          do k=1,N
            do i=istr,iend
               TVmix(i, j, k, itrc)=t(i,j,k,nnew,itrc)
            enddo
         enddo
# elif defined DIAGNOSTICS_PV
            do i=istr,iend
               twrk_bot(i,j,itrc)=t(i,j,1,nnew,itrc)
               twrk_sfc(i,j,itrc)=t(i,j,N,nnew,itrc)
            enddo
# endif /* DIAGNOSTICS_TS */

          indx=itrc
          do i=istr,iend
            FC(i,1)=dt*Akt(i,j,1,indx)/(z_r(i,j,2)-z_r(i,j,1))
            cff=1./(Hz(i,j,1)+FC(i,1))
            CF(i,1)=cff*FC(i,1)
            DC(i,1)=cff*t(i,j,1,nnew,itrc)
          enddo
          do k=2,N-1,+1
            do i=istr,iend
              FC(i,k)=dt*Akt(i,j,k,indx)/(z_r(i,j,k+1)-z_r(i,j,k))
              cff=1./( Hz(i,j,k) +FC(i,k)+FC(i,k-1)*(1.-CF(i,k-1)) )
              CF(i,k)=cff*FC(i,k)
              DC(i,k)=cff*(t(i,j,k,nnew,itrc)+FC(i,k-1)*DC(i,k-1))
            enddo
          enddo
          do i=istr,iend
             t(i,j,N,nnew,itrc)=( t(i,j,N,nnew,itrc) +FC(i,N-1)
     &         *DC(i,N-1) )/(Hz(i,j,N)+FC(i,N-1)*(1.-CF(i,N-1)))
          enddo
          do k=N-1,1,-1
            do i=istr,iend
              t(i,j,k,nnew,itrc)=DC(i,k)+CF(i,k)*t(i,j,k+1,nnew,itrc)
            enddo
          enddo           !--> discard FC,CF,DC


# ifdef DIAGNOSTICS_TS
          do k=1,N
            do i=istr,iend
              TVmix(i,j,k,itrc) = 
     &            -(TVmix(i,j,k,itrc)-t(i,j,k,nnew,itrc)*Hz(i,j,k))
     &                                        /(dt*pm(i,j)*pn(i,j))
#  ifdef MASKING
     &                                                 * rmask(i,j)
#  endif
            enddo
          enddo     
# elif defined DIAGNOSTICS_PV
            do i=istr,iend
              pvDiab_bot(i,j,itrc) = pvDiab_bot(i,j,itrc)
     &             -(twrk_bot(i,j,itrc)-t(i,j,1,nnew,itrc)*Hz(i,j,1))
     &                                        /(dt*pm(i,j)*pn(i,j))  
#  ifdef MASKING
     &                                                 * rmask(i,j)
#  endif
              pvDiab_sfc(i,j,itrc) = pvDiab_sfc(i,j,itrc)
     &             -(twrk_sfc(i,j,itrc)-t(i,j,N,nnew,itrc)*Hz(i,j,N))
     &                                        /(dt*pm(i,j)*pn(i,j))
#  ifdef MASKING
     &                                                 * rmask(i,j)
#  endif
            enddo
# endif /* DIAGNOSTICS_TS */



#ifdef CONST_TRACERS
          do k=1,N
            do i=istr,iend
              t(i,j,k,nnew,itrc)=t(i,j,k,nstp,itrc)
#  ifdef DIAGNOSTICS_TS
              TVmix(i,j,k,itrc)=0.0
#   ifdef MASKING
     &              * rmask(i,j)
#   endif
#  endif /* DIAGNOSTICS_TS */
            enddo
          enddo
#endif
        enddo      ! <-- itrc
      enddo      ! <-- j
!
! Set lateral boundary conditions; nudge toward tracer climatology;
! apply land-sea mask and exchange periodic boundary conditions.
!
      do itrc=1,NT
        call t3dbc_tile (istr,iend,jstr,jend, itrc, WORK)
 
# if defined MASKING || defined TCLIMATOLOGY
#  ifdef EW_PERIODIC
#   define I_RANGE istr,iend
#  else
#   define I_RANGE istrR,iendR
#  endif
#  ifdef NS_PERIODIC
#   define J_RANGE jstr,jend
#  else
#   define J_RANGE jstrR,jendR
#  endif
        do k=1,N
          do j=J_RANGE
            do i=I_RANGE
#  ifdef TCLIMATOLOGY
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)
     &       +dt*Tnudgcof(i,j,itrc)*( tclm(i,j,k,itrc)
     &                             -t(i,j,k,nnew,itrc))
#  ifdef DIAGNOSTICS_TS
              TForc(i,j,k,itrc)=(TForc(i,j,k,itrc)
     &                    +Tnudgcof(i,j,itrc)
     &                    *(tclm(i,j,k,itrc)-t(i,j,k,nnew,itrc))
     &                    *(Hz(i,j,k)/(pm(i,j)*pn(i,j))))
#   ifdef MASKING
     &                    *rmask(i,j)
#   endif

# elif defined DIAGNOSTICS_PV

               if (k.eq.N) then
                pvDiab_sfc(i,j,itrc)=pvDiab_sfc(i,j,itrc)
     &                    +Tnudgcof(i,j,itrc)
     &                    *(tclm(i,j,k,itrc)-t(i,j,k,nnew,itrc))
     &                    *(Hz(i,j,k)/(pm(i,j)*pn(i,j)))
               endif

               if (k.eq.1) then
                pvDiab_bot(i,j,itrc)=pvDiab_bot(i,j,itrc)
     &                    +Tnudgcof(i,j,itrc)
     &                    *(tclm(i,j,k,itrc)-t(i,j,k,nnew,itrc))
     &                    *(Hz(i,j,k)/(pm(i,j)*pn(i,j)))
               endif

#  endif  /* DIAGNOSTICS_TS */


#  endif
#  ifdef MASKING
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)*rmask(i,j)
#  endif
            enddo
          enddo
        enddo



# ifdef DIAGNOSTICS_TS

!---------------------------------------------------------------
! Compute the tendency term of tracer diagnostics
! Divide all diagnostic terms by the cell volume
! (Hz(i,j,k,itrc)/(pm(i,j).*pn(i,j)). There after the unit
! of diagnostic terms will be: (unit of tracers)* s-1.
!
! Note: the Horizontal mixing term is computed in t3dmix
! where Trate is updated accordingly
!---------------------------------------------------------------

       do k=1,N
         do j=jstr,jend
           do i=istr,iend
              Trate(i,j,k,itrc)=(Hz(i,j,k)*t(i,j,k,nnew,itrc)
     &                          -Hz_bak(i,j,k)*t(i,j,k,nstp,itrc))
     &                                       /(dt*pm(i,j)*pn(i,j))
#  ifdef MASKING
     &                   *rmask(i,j)
#  endif

              cff=pm(i,j)*pn(i,j)/Hz(i,j,k)
              TXadv(i,j,k,itrc)=TXadv(i,j,k,itrc)*cff
              TYadv(i,j,k,itrc)=TYadv(i,j,k,itrc)*cff
              TVadv(i,j,k,itrc)=TVadv(i,j,k,itrc)*cff
              TVmix(i,j,k,itrc)=TVmix(i,j,k,itrc)*cff
              TForc(i,j,k,itrc)=TForc(i,j,k,itrc)*cff
              Trate(i,j,k,itrc)=Trate(i,j,k,itrc)*cff
              THmix(i,j,k,itrc)=0.
           enddo
         enddo
       enddo

# ifdef DIAGNOSTICS_PV
         do j=jstr,jend
           do i=istr,iend

              pvDiab_sfc(i,j,itrc)=TVmix(i,j,N,itrc) +
     &               TForc(i,j,N,itrc) 

              pvDiab_bot(i,j,itrc)=TVmix(i,j,1,itrc) +
     &               TForc(i,j,1,itrc) 


           enddo
         enddo
# endif /* DIAGNOSTICS_PV */
# endif /* DIAGNOSTICS_TS */



# if defined DIAGNOSTICS_PV && !defined DIAGNOSTICS_TS
         do j=jstr,jend
           do i=istr,iend

              cff=pm(i,j)*pn(i,j)/Hz(i,j,N)
              pvDiab_sfc(i,j,itrc)=pvDiab_sfc(i,j,itrc)*cff

              cff=pm(i,j)*pn(i,j)/Hz(i,j,1)
              pvDiab_bot(i,j,itrc)=pvDiab_bot(i,j,itrc)*cff

           enddo
         enddo

# endif










 


#  undef I_RANGE
#  undef J_RANGE
# endif
      enddo  ! <-- itrc
!
! Exchange computational margines and/or periodic boundaries:
!--------- ------------- -------- ------ -------- ----------- 
! This is done only here if the code is configured without explicit
! lateral diffusion. If there is one, then the exchange routine is
! called from the diffusive routine which immediately follows this
! step, therefore there is no need to call it here.
!
# if !defined TS_DIF2 && !defined TS_DIF2_SPONGE
#  if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
      do itrc=1,NT
        call exchange_r3d_tile (istr,iend,jstr,jend,
     &                          t(START_2D_ARRAY,1,nnew,itrc))
      enddo
#  endif
# endif

c+++++#if !defined EPVIC && defined MPI
c+++++      do itrc=1,NT
c+++++        call exchange_r3d_tile (istr,iend,jstr,jend,
c+++++     &                          t(START_2D_ARRAY,1,nnew,itrc))
c+++++      enddo
c++++#endif 
      return
      end




      subroutine check_trc_adv_switches (ierr)
!
! This routine keeps track of the status of CPP-switch settings in
! "compute_horiz_tracer_fluxes.h" and "compute_horiz_vert_fluxes.h".
! This code must be placed here rather than in a separate file in
! order to be exposed to the relevant CPP-settings.   It does not
! affect any model results, other than signature in global attribute
! "CPPS" in output netCDF files.
!
      implicit none
      integer ierr, is,ie, lenstr
# include "param.h"
# include "strings.h"
      ie=lenstr(cpps)
      is=ie+2
      ie=is+30
      if (ie.gt.max_opt_size) goto 99
      cpps(is:ie)='<compute_horiz_tracer_fluxes.h>'
      MPI_master_only write(stdout,'(1x,A)') cpps(is:ie) 
      is=ie+2

# ifdef UPSTREAM
      ie=is+7
      if (ie.gt.max_opt_size) goto 99
      cpps(is:ie)='UPSTREAM'
      MPI_master_only write(stdout,'(10x,A)') cpps(is:ie)
      is=ie+2
# endif
# ifdef AKIMA
      ie=is+4
      if (ie.gt.max_opt_size) goto 99
      cpps(is:ie)='AKIMA'
      MPI_master_only write(stdout,'(10x,A)') cpps(is:ie)
      is=ie+2
# endif
# ifdef CONST_TRACERS
      ie=is+7
      if (ie.gt.max_opt_size) goto 99
      cpps(is:ie)='CONST_TRACERS'
      MPI_master_only write(stdout,'(10x,A)') cpps(is:ie)
      is=ie+2
# endif

      ie=is+29
      if (ie.gt.max_opt_size) goto 99
      cpps(is:ie)='<compute_vert_tracer_fluxes.h>'
      MPI_master_only write(stdout,'(1x,A)') cpps(is:ie)
      is=ie+2
# ifdef SPLINES
      ie=is+6
      if (ie.gt.max_opt_size) goto 99
      cpps(is:ie)='SPLINES'
      MPI_master_only write(stdout,'(10x,A)') cpps(is:ie)
      is=ie+2
# endif
# ifdef AKIMA_V
      ie=is+6
      if (ie.gt.max_opt_size) goto 99
      cpps(is:ie)='AKIMA_V'
      MPI_master_only write(stdout,'(10x,A)') cpps(is:ie)
      is=ie+2
# endif
      return
  99  MPI_master_only write(stdout,'(/1x,2A/12x,A/)') '### ERROR: ',
     &  'Unsufficient lenght of string "cpps" in file "strings.h".',
     &        'Increase parameter "max_opt_size" it and recompile.'
      ierr=ierr+1
      return
      end
#else
      subroutine step3d_t_empty
      end
#endif  /* SOLVE3D */
 
