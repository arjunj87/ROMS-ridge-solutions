#include "cppdefs.h"
 
      subroutine init_arrays (tile)
      implicit none
      integer tile, i,j
#include "param.h"
#include "private_scratch.h"
#include "compute_tile_bounds.h"

#ifdef SOLVE3D
      do j=1,6           ! Initialize (first touch) private
        do i=1,N3d       ! scratch arrays in parallel by each
          A3d(i,j)=0.    ! thread.
        enddo
      enddo
      do i=1,N2d
        iA2d(i)=0
      enddo
#endif
      do j=1,32
        do i=1,N2d
          A2d(i,j)=0.
        enddo
      enddo

      call init_arrays_tile (istr,iend,jstr,jend)
      return
      end
 
      subroutine init_arrays_tile (istr,iend,jstr,jend)
!
! Initialize (first touch) globally accessable arrays. Most of them
! are assigned to zeros, vertical mixing coefficients are assinged
! to their background values. These will remain unchenged if no
! vertical mixing scheme is applied. Because of the "first touch"
! default data distribution policy, this operation actually performs
! distribution of the shared arrays accross the cluster, unless
! another distribution policy is specified to override the default.
!
      implicit none
      integer istr,iend,jstr,jend, i,j,k,itrc
      real init    !!!!  0xFFFA5A5A ==> NaN
      parameter (init=0.)
#define ALL_DATA
#include "param.h"
#include "scalars.h"
#include "grid.h"
#include "ocean2d.h"
#include "ocean3d.h"
#include "eos_vars.h"
#include "coupling.h"
#include "averages.h"
#include "mixing.h"
#include "forces.h"
#include "climat.h"

# ifdef DIAGNOSTICS_UV
#include "diags_uv.h"
# endif

# ifdef DIAGNOSTICS_VRT
#include "diags_vrt.h"
# endif

# ifdef DIAGNOSTICS_TS
#include "diags_ts.h"
# endif

# ifdef DIAGNOSTICS_PV
#include "diags_pv.h"
# endif

# ifdef DIAGNOSTICS_EDDY
#include "diags_eddy.h"
# endif

#undef ALL_DATA
#ifdef MPI
# include "mpif.h"
      integer status(MPI_STATUS_SIZE), blank, ierr
#endif
!
#include "compute_extended_bounds.h"
!
#ifdef NOT_
# ifdef MPI
      if (mynode.gt.0) then
        call MPI_Recv (blank, 1, MPI_INTEGER, mynode-1,
     &                 1, ocean_grid_comm, status, ierr)
      endif
      i=mynode
# else
      i=proc(2)
# endif
      write(*,'(I4/2(6x,A6,I3,3x,A6,I3))') i, 'istr =',istr,
     &        'iend =',iend,   'jstr =',jstr, 'jend =',jend
      write(*,'(4x,2(6x,A6,I3,3x,A6,I3)/)')   'istrR=',istrR,
     &        'iendR=',iendR, 'jstrR=',jstrR, 'jendR=',jendR
# ifdef MPI
      if (mynode .lt. NNODES) then
        call MPI_Send (blank, 1, MPI_INTEGER, mynode+1,
     &                        1, ocean_grid_comm,  ierr)
      endif
# endif
#endif
      do j=jstrR,jendR               ! Initialize
        do i=istrR,iendR             ! 2-D primitive
          zeta(i,j,1)=0. ! init      ! variables.
          zeta(i,j,2)=init
          zeta(i,j,3)=init
 
          ubar(i,j,1)=init
          ubar(i,j,2)=init
          ubar(i,j,3)=init
 
          vbar(i,j,1)=init
          vbar(i,j,2)=init
          vbar(i,j,3)=init
#ifdef AVERAGES
          zeta_avg(i,j)=init
          ubar_avg(i,j)=init
          vbar_avg(i,j)=init
#endif
#ifdef SOLVE3D
          rufrc(i,j)=init
          rufrc(i,j)=init
# ifdef VAR_RHO_2D
          rhoA(i,j)=0.
          rhoS(i,j)=0.
# endif
          Zt_avg1(i,j)=0. !init
          DU_avg1(i,j)=0. !init
          DV_avg1(i,j)=0. !init
          DU_avg2(i,j)=0. !init
          DV_avg2(i,j)=0. !init
#endif

#ifdef ADV_ISONEUTRAL
          usmax(i,j) = 1.
          vsmax(i,j) = 1.
#endif
!
! WARNING: "rmask" MUST BE initialized to all-one state in order
! to read grid variables (coordinates, metric, topography), which
! should not be masked. 
!
#ifdef MASKING
          rmask(i,j)=1.
#endif
        enddo
      enddo
#ifdef SOLVE3D
      do k=1,N                       ! Initialize
        do j=jstrR,jendR             ! 3-D primitive
          do i=istrR,iendR           ! variables.
            u(i,j,k,1)=init
            u(i,j,k,2)=init
 
            v(i,j,k,1)=init
            v(i,j,k,2)=init
 
            rho1(i,j,k)=init
# ifdef SPLIT_EOS
            qp1(i,j,k)=init
# else
            rho(i,j,k)=init
# endif
# ifdef AVERAGES
            rho_avg(i,j,k)=init
            u_avg(i,j,k)=init
            v_avg(i,j,k)=init
# endif

#ifdef ADV_ISONEUTRAL
            Sx_u  (i,j,k) = init
            Sy_v  (i,j,k) = init            
# ifdef WRITE_SXSY
            Sx_avg(i,j,k) = init 
            Sy_avg(i,j,k) = init 
# endif
            diff3u(i,j,k) = init
            diff3v(i,j,k) = init
            
            dRdX(i,j,k) = init
            dRdE(i,j,k) = init                      
#endif 




# ifdef DIAGNOSTICS_UV
            MXadv(i,j,k,1)=init ! indice 1 for u_momentum
            MXadv(i,j,k,2)=init ! indice 2 for v_momentum
            MYadv(i,j,k,1)=init ! ...
            MYadv(i,j,k,2)=init
            MVadv(i,j,k,1)=init
            MVadv(i,j,k,2)=init
            MCor(i,j,k,1)=init
            MCor(i,j,k,2)=init
            MPrsgrd(i,j,k,1)=init
            MPrsgrd(i,j,k,2)=init
            MHmix(i,j,k,1)=init
            MHmix(i,j,k,2)=init
            MVmix(i,j,k,1)=init
            MVmix(i,j,k,2)=init
            Mrate(i,j,k,1)=init
            Mrate(i,j,k,2)=init
            MVmix2(i,j,k,1)=init
            MVmix2(i,j,k,2)=init
#  ifdef AVERAGES
            MXadv_avg(i,j,k,1)=init
            MXadv_avg(i,j,k,2)=init
            MYadv_avg(i,j,k,1)=init
            MYadv_avg(i,j,k,2)=init
            MVadv_avg(i,j,k,1)=init
            MVadv_avg(i,j,k,2)=init
            MCor_avg(i,j,k,1)=init
            MCor_avg(i,j,k,2)=init
            MPrsgrd_avg(i,j,k,1)=init
            MPrsgrd_avg(i,j,k,2)=init
            MHmix_avg(i,j,k,1)=init
            MHmix_avg(i,j,k,2)=init
            MVmix_avg(i,j,k,1)=init
            MVmix_avg(i,j,k,2)=init
            Mrate_avg(i,j,k,1)=init
            Mrate_avg(i,j,k,2)=init
            MVmix2_avg(i,j,k,1)=init
            MVmix2_avg(i,j,k,2)=init
#  endif /* AVERAGES */
# endif /* DIAGNOSTICS_UV */

# ifdef DIAGNOSTICS_TS
            TXadv(i,j,k,1)=init ! indice 1 for u_momentum
            TXadv(i,j,k,2)=init ! indice 2 for v_momentum
            TYadv(i,j,k,1)=init ! ...
            TYadv(i,j,k,2)=init
            TVadv(i,j,k,1)=init
            TVadv(i,j,k,2)=init
            TForc(i,j,k,1)=init
            TForc(i,j,k,2)=init
            THmix(i,j,k,1)=init
            THmix(i,j,k,2)=init
            TVmix(i,j,k,1)=init
            TVmix(i,j,k,2)=init
            Trate(i,j,k,1)=init
            Trate(i,j,k,2)=init
#  ifdef AVERAGES
            TXadv_avg(i,j,k,1)=init
            TXadv_avg(i,j,k,2)=init
            TYadv_avg(i,j,k,1)=init
            TYadv_avg(i,j,k,2)=init
            TVadv_avg(i,j,k,1)=init
            TVadv_avg(i,j,k,2)=init
            TForc_avg(i,j,k,1)=init
            TForc_avg(i,j,k,2)=init
            THmix_avg(i,j,k,1)=init
            THmix_avg(i,j,k,2)=init
            TVmix_avg(i,j,k,1)=init
            TVmix_avg(i,j,k,2)=init
            Trate_avg(i,j,k,1)=init
            Trate_avg(i,j,k,2)=init

#  endif /* AVERAGES */
# endif /* DIAGNOSTICS_TS */

# ifdef DIAGNOSTICS_EDDY
#  ifdef AVERAGES
            eddyuu_avg(i,j,k)=init
            eddyvv_avg(i,j,k)=init
            eddyuv_avg(i,j,k)=init
            eddyub_avg(i,j,k)=init
            eddyvb_avg(i,j,k)=init
            eddywb_avg(i,j,k)=init
            eddynotused1_avg(i,j,k)=init
            eddynotused2_avg(i,j,k)=init
            eddynotused3_avg(i,j,k)=init
#  endif /* AVERAGES */
# endif /* DIAGNOSTICS_EDDY */

          enddo
        enddo
      enddo


# ifdef DIAGNOSTICS_VRT
        do j=jstrR,jendR             ! 3-D primitive
          do i=istrR,iendR           ! variables.

            vrtXadv(i,j)=init 
            vrtYadv(i,j)=init 
            vrtCor(i,j)=init
            vrtPrsgrd(i,j)=init
            vrtHmix(i,j)=init
            vrtVmix(i,j)=init
            vrtrate(i,j)=init
            vrtVmix2(i,j)=init
            vrtWind(i,j)=init
            vrtDrag(i,j)=init

#  ifdef AVERAGES
            vrtXadv_avg(i,j)=init
            vrtYadv_avg(i,j)=init
            vrtCor_avg(i,j)=init
            vrtPrsgrd_avg(i,j)=init
            vrtHmix_avg(i,j)=init
            vrtVmix_avg(i,j)=init
            vrtrate_avg(i,j)=init
            vrtVmix2_avg(i,j)=init
            vrtWind_avg(i,j)=init
            vrtDrag_avg(i,j)=init
#  endif /* AVERAGES */
          enddo
        enddo
# endif /* DIAGNOSTICS_VRT */





# ifdef DIAGNOSTICS_PV

        do j=jstrR,jendR    
          do i=istrR,iendR 

            pvForc_sfc(i,j,1)=init 
            pvForc_sfc(i,j,2)=init 
            pvForc_bot(i,j,1)=init 
            pvForc_bot(i,j,2)=init 

#  ifdef AVERAGES
            pvForc_sfc_avg(i,j,1)=init
            pvForc_sfc_avg(i,j,2)=init
            pvForc_bot_avg(i,j,1)=init
            pvForc_bot_avg(i,j,2)=init

#  endif /* AVERAGES */


          do itrc=1,NT

            pvDiab_sfc(i,j,itrc)=init 
            pvDiab_bot(i,j,itrc)=init 
#  ifdef AVERAGES
            pvDiab_sfc_avg(i,j,itrc)=init
            pvDiab_bot_avg(i,j,itrc)=init
#  endif /* AVERAGES */
           enddo

          enddo
         enddo
# endif /* DIAGNOSTICS_PV */







      do k=0,N
        do j=jstrR,jendR
          do i=istrR,iendR                       
            
            W(i,j,k)=init
# ifdef AVERAGES
            w_avg(i,j,k)=init
# endif
#ifdef ADV_ISONEUTRAL
            dZdR(i,j,k) = init
#endif  
          enddo
        enddo
      enddo
      
      do itrc=1,NT
        do k=1,N
          do j=jstrR,jendR
            do i=istrR,iendR
              t(i,j,k,1,itrc)=init
              t(i,j,k,2,itrc)=init
# ifdef AVERAGES
              t_avg(i,j,k,itrc)=init
# endif
            enddo
          enddo
        enddo
      enddo
#endif /* SOLVE3D */
!
!  Initialize forcing arrays (see "forces.h").
!
      do j=jstrR,jendR
        do i=istrR,iendR
          sustr(i,j)=init
          svstr(i,j)=init
#ifndef ANA_SMFLUX
          sustrg(i,j,1)=init
          svstrg(i,j,1)=init
          sustrg(i,j,2)=init
          svstrg(i,j,2)=init
#endif
        enddo
      enddo
#ifdef SOLVE3D
      do itrc=1,NT
        do j=jstrR,jendR
          do i=istrR,iendR
            stflx(i,j,itrc)=init
# if !defined ANA_STFLUX || !defined ANA_SSFLUX
            stflxg(i,j,1,itrc)=init
            stflxg(i,j,2,itrc)=init
# endif
          enddo
        enddo
      enddo
      do j=jstrR,jendR
        do i=istrR,iendR
# ifdef QCORRECTION
          dqdtg(i,j,1)=init
          sstg (i,j,1)=init
          dqdtg(i,j,2)=init
          sstg (i,j,2)=init
# endif
# if defined QCORRECTION || defined SALINITY_RESTORING
          sssg (i,j,1)=init
          sssg (i,j,2)=init
# endif
# ifdef BULK_FLUX
          tair (i,j)=init
          qair (i,j)=init
          rain (i,j)=init
          radlw (i,j)=init
          radsw (i,j)=init
          tau   (i,j)=init
#  ifdef BULK_SMFLUX
          uwnd(i,j)=init
          vwnd(i,j)=init
#  endif
          tairg(i,j,1)=init
          qairg (i,j,1)=init
          raing (i,j,1)=init
          radlwg (i,j,1)=init
          radswg (i,j,1)=init
          tairg(i,j,2)=init
          qairg (i,j,2)=init
          raing (i,j,2)=init
          radlwg (i,j,2)=init
          radswg (i,j,2)=init
#  ifdef BULK_SMFLUX
          uwndg(i,j,1)=init
          vwndg(i,j,1)=init
          uwndg(i,j,2)=init
          vwndg(i,j,2)=init
#  endif
          sustr_blk(i,j)=init
          svstr_blk(i,j)=init          
          shflx_rsw(i,j)=init
          shflx_rlw(i,j)=init
          shflx_lat(i,j)=init
          shflx_sen(i,j)=init
#  ifdef AVERAGES
          sustr_blk_avg(i,j)=init
          svstr_blk_avg(i,j)=init
          shflx_rsw_avg(i,j)=init
          shflx_rlw_avg(i,j)=init
          shflx_lat_avg(i,j)=init
          shflx_sen_avg(i,j)=init
#  endif
# endif /* BULK_FLUX */

          srflx(i,j)=init
# ifndef ANA_SRFLUX
          srflxg(i,j,1)=init
          srflxg(i,j,2)=init
# endif
        enddo
      enddo
#endif /* SOLVE3D */
!
! Initialize climatology arrays (see "climat.h").
!
#if defined M2NUDGING && !defined M2_FRC_BRY
      do j=jstrR,jendR
        do i=istrR,iendR
          ssh(i,j)=init
# ifndef ANA_SSH
          sshg(i,j,1)=init
          sshg(i,j,2)=init
# endif
        enddo
      enddo
#endif
#ifdef SOLVE3D
# ifdef TCLIMATOLOGY 
      do itrc=1,NT
        do j=jstrR,jendR
          do i=istrR,iendR
            Tnudgcof(i,j,itrc)=init
          enddo
        enddo
      enddo
# endif
# if (defined TCLIMATOLOGY || defined TNUDGING) && !defined ANA_TCLIMA
      do itrc=1,NT
        do k=1,N
          do j=jstrR,jendR
            do i=istrR,iendR
              tclm(i,j,k,itrc)=init
              tclima(i,j,k,1,itrc)=init
              tclima(i,j,k,2,itrc)=init
            enddo
          enddo
        enddo
      enddo
# endif
#endif

#ifdef UCLIMATOLOGY
# ifndef ANA_UCLIMA
      do j=jstrR,jendR
        do i=istrR,iendR
          ubclm(i,j)=init
          vbclm(i,j)=init
          ubclima(i,j,1)=init
          ubclima(i,j,2)=init
          vbclima(i,j,1)=init
          vbclima(i,j,2)=init
        enddo
      enddo
#  ifdef SOLVE3D
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
            uclm(i,j,k)=init
            vclm(i,j,k)=init
            uclima(i,j,k,1)=init
            uclima(i,j,k,2)=init
            vclima(i,j,k,1)=init
            vclima(i,j,k,2)=init
          enddo
        enddo
      enddo
#  endif
# endif
#endif
!
! Set variable horizontal viscosities and tracer diffusion
! coefficients (see "mixing.h") to their background values.
!
#ifdef UV_VIS2
        do j=jstrR,jendR
          do i=istrR,iendR
            visc2_r(i,j)=visc2
            visc2_p(i,j)=visc2
          enddo
        enddo
#endif
#ifdef SOLVE3D
# ifdef TS_DIF2
        do itrc=1,NT
          do j=jstrR,jendR
            do i=istrR,iendR
              diff2(i,j,itrc)=tnu2(itrc)
            enddo
          enddo
        enddo
# endif
!
! Initialize vertical mixing coefficients (see "mixing.h") to their
! background values. If no vertical closure scheme is selected, the
! vertical mixing coefficients are those specified by the background
! values.
!
      do k=0,N
        do j=jstrR,jendR
          do i=istrR,iendR
# if !defined LMD_MIXING && !defined BVF_MIXING\
  && !defined MY2_MIXING && !defined MY25_MIXING\
                         && !defined PP_MIXING
            Akv(i,j,k)=Akv_bak
# else
            Akv(i,j,k)=init
# endif
# if defined BVF_MIXING || defined LMD_MIXING  || defined LMD_KPP \
  || defined MY2_MIXING || defined MY25_MIXING || defined PP_MIXING
            bvf(i,j,k)=init
# endif
          enddo
        enddo
        do itrc=1,NT
          do j=jstrR,jendR
            do i=istrR,iendR
# if !defined LMD_MIXING && !defined BVF_MIXING\
  && !defined MY2_MIXING && !defined MY25_MIXING\
                         && !defined PP_MIXING
              Akt(i,j,k,itrc)=Akt_bak(itrc)
# else
              Akt(i,j,k,itrc)=init
# endif
            enddo
          enddo
        enddo
      enddo
 
# if defined LMD_KPP && defined LMD_NONLOCAL
      do k=1,N
        do j=jstrR,jendR
          do i=istrR,iendR
            ghat(i,j,k)=init
          enddo
        enddo
      enddo
# endif
#ifdef LMD_KPP
#ifdef DIAGS_CR
      do k=0,N
        do j=jstrR,jendR
          do i=istrR,iendR
             crzs(i,j,k) = 0.
             crzb(i,j,k) = 0.             
          enddo
        enddo
      enddo                    
#endif
!
!  Initialize depth of planetary boundary layer.
!
      do j=jstrR,jendR
        do i=istrR,iendR
          hbls(i,j,1)=0.  ! init
          hbls(i,j,2)=0.
#ifdef LMD_BKPP          
          hbbls(i,j,1)=0.  ! init
          hbbls(i,j,2)=0.
#endif          
        enddo
      enddo
# endif /* LMD_KPP */

#endif /* SOLVE3D */
      return
      end
 
