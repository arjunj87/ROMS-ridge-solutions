#include "cppdefs.h"
#if defined SOLVE3D && defined LMD_MIXING
# define SMOOTH_RIG  /* local Ri number smoothing */
 
      subroutine lmd_vmix (tile)
      implicit none
      integer tile
# include "param.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"
      call lmd_vmix_tile (istr,iend,jstr,jend, A3d(1,1), A3d(1,2),
     &                                         A3d(1,3), A3d(1,4))
# ifdef LMD_KPP
      call lmd_kpp_tile  (istr,iend,jstr,jend, A3d(1,1),  A3d(1, 2),
     &                                                     A3d(1,3), 

     &                   A2d(1, 1), A2d(1, 2), A2d(1, 3), A2d(1, 4),
     &                   A2d(1, 5), A2d(1, 6), A2d(1, 7), A2d(1, 8),
     &                   A2d(1, 9), A2d(1,10), A2d(1,11), A2d(1,12),
     &                   A2d(1,13), A2d(1,14), A2d(1,15), A2d(1,16),
     &                   A2d(1,17), A2d(1,18), A2d(1,19)          )
# endif
      return
      end
 
      subroutine lmd_vmix_tile (istr,iend,jstr,jend, Kv,Kt,Ks, Rig)
!
!  This subroutine computes vertical mixing coefficients for momentum
!  and tracers at the ocean interior using the Large, McWilliams and
!  Doney (1994) mixing scheme.
!
!  On Output:
!     Kv       vertical viscosity coefficient [m^2/s].
!     Kt       vertical diffusion coefficient for potential
!                 temperature [m^2/s].
!     Ks       vertical diffusion coefficient for salinity [m^2/s].
!
!  Reference:
!
!    Large, W.G., J.C. McWilliams, and S.C. Doney, 1994: A Review
!    and model with a nonlocal boundary layer parameterization,
!    Reviews of Geophysics, 32,363-403.
!
      implicit none
# include "param.h"
# include "grid.h"
# include "ocean3d.h"
# include "mixing.h"
# include "scalars.h"

      integer istr,iend,jstr,jend, i,j,k
      real Rig(PRIVATE_2D_SCRATCH_ARRAY,0:N),
     &      Kv(PRIVATE_2D_SCRATCH_ARRAY,0:N),
     &      Kt(PRIVATE_2D_SCRATCH_ARRAY,0:N),
     &      Ks(PRIVATE_2D_SCRATCH_ARRAY,0:N)

      real nu_sx, cff,dudz,dvdz
      real Ri0,  nuwm, nuws,  nu0m, nu0s, nu0c, lmd_nu, lmd_Rrho0,
     &     lmd_nuf,  lmd_fdd, lmd_tdd1, lmd_tdd2, lmd_tdd3, lmd_sdd1,
     &     lmd_sdd2, lmd_sdd3, eps

      parameter ( Ri0=0.7,  ! Critical gradient Richardson number
     &                      ! below which shear instabilty occurs.

     &     nu0m=50.e-4,     ! Maximum viscosity and diffusivity
     &     nu0s=50.e-4,     ! due to shear instability [m^2/s];

     &     nuwm=1e-4,     ! Interior viscosity and diffusivity
     &     nuws=0.1e-4,     ! due to wave breaking, [m^2/s]
!     &     nuwm=0.d0,     ! Interior viscosity and diffusivity
!     &     nuws=0.d0,     ! due to wave breaking, [m^2/s]
                            !
     &     nu0c=0.1,        ! Convective adjustment [m^2/s];
                            !
     &     lmd_nu=1.5e-6,   ! Molecular viscosity [m^2/s];
                            !
                            ! Value of double-diffusive density
     &     lmd_Rrho0=1.9,   ! ratio where diffusivity goes to
                            !              zero in salt fingering.
     &     lmd_nuf=10.0e-4, ! Scaling factors for double diffusion
     &     lmd_fdd=0.7,     ! coefficient in salt fingering.
                            !
     &     lmd_tdd1=0.909,  !
     &     lmd_tdd2=4.6,    ! Double diffusion constants
     &     lmd_tdd3=0.54,   ! for temperature (Marmorino and
     &     lmd_sdd1=0.15,   ! Caldwell, 1976) and salinity
     &     lmd_sdd2=1.85,   ! (Fedorov, 1988).
     &     lmd_sdd3=0.85)   !
 
# ifdef LMD_DDMIX
      real Rrho, ddDS, ddDT, nu_dds, nu_ddt, alfaobeta, Tt, Ts, Tp
 
      real A0,A1,A2,A3,A4, B0,B1, C0, D0,D1,D2, E0,F0,G0,H0, Smean
      parameter(A0=+0.665157E-01, A1=+0.170907E-01, A2=-0.203814E-03,
     &          A3=+0.298357E-05, A4=-0.255019E-07, B0=+0.378110E-02,
     &          B1=-0.846960E-04, C0=-0.678662E-05, D0=+0.380374E-04,
     &          D1=-0.933746E-06, D2=+0.791325E-08, E0=-0.164759E-06,
     &          F0=-0.251520E-11, G0=+0.512857E-12, H0=-0.302285E-13,
     &          Smean=35.0)
# endif
      parameter (eps=1.E-14)

# define tind nstp
# ifdef SMOOTH_RIG 
#  ifdef EW_PERIODIC
#   define I_EXT_RANGE istr-1,iend+1
#  else
#   define I_EXT_RANGE imin,imax
      integer imin,imax
#  endif
#  ifdef NS_PERIODIC
#   define J_EXT_RANGE jstr-1,jend+1
#  else
#   define J_EXT_RANGE jmin,jmax
      integer jmin,jmax
#  endif
#  ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
        imin=istr
      else
        imin=istr-1
      endif
      if (EASTERN_EDGE) then
        imax=iend
      else
        imax=iend+1
      endif
#  endif
#  ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
        jmin=jstr
      else
        jmin=jstr-1
      endif
      if (NORTHERN_EDGE) then
        jmax=jend
      else
        jmax=jend+1
      endif
#  endif
# else
#  define I_EXT_RANGE istr,iend
#  define J_EXT_RANGE jstr,jend
# endif
!                                    ! Compute horizontal velocity
! Compute local Richardson number:   ! shear  (du/dz)^2+(dv/dz)^2  at
!-------- ----- ---------- -------   ! horizontal RHO- and vertical
!                                    ! W-points. Then compute gradient
      do k=1,N-1                     ! Richardson number (already
        do j=J_EXT_RANGE             ! divided by its critical value.
          do i=I_EXT_RANGE
            cff=0.5/(z_r(i,j,k+1)-z_r(i,j,k))
            dudz=cff*( u(i  ,j,k+1,tind)-u(i  ,j,k,tind)
     &                +u(i+1,j,k+1,tind)-u(i+1,j,k,tind))
            dvdz=cff*( v(i,j  ,k+1,tind)-v(i,j  ,k,tind)
     &                +v(i,j+1,k+1,tind)-v(i,j+1,k,tind))

            Rig(i,j,k)=bvf(i,j,k)/( Ri0*max(
     &                    dudz*dudz+dvdz*dvdz, 1.E-10 ))
          enddo
        enddo

# ifdef SMOOTH_RIG
#  ifndef EW_PERIODIC
        if (WESTERN_EDGE) then
          do j=J_EXT_RANGE
            Rig(istr-1,j,k)=Rig(istr,j,k)
          enddo
        endif
        if (EASTERN_EDGE) then
          do j=J_EXT_RANGE
            Rig(iend+1,j,k)=Rig(iend,j,k)
          enddo
        endif
#  endif
#  ifndef NS_PERIODIC
        if (SOUTHERN_EDGE) then
          do i=I_EXT_RANGE
            Rig(i,jstr-1,k)=Rig(i,jstr,k)
          enddo
        endif
        if (NORTHERN_EDGE) then
          do i=I_EXT_RANGE
            Rig(i,jend+1,k)=Rig(i,jend,k)
          enddo
        endif
#   ifndef EW_PERIODIC
        if (WESTERN_EDGE .and.
     &      SOUTHERN_EDGE) then
          Rig(istr-1,jstr-1,k)=Rig(istr,jstr,k)
        endif
        if (WESTERN_EDGE .and.
     &      NORTHERN_EDGE) then
          Rig(istr-1,jend+1,k)=Rig(istr,jend,k)
        endif
        if (EASTERN_EDGE .and.
     &      SOUTHERN_EDGE) then
          Rig(iend+1,jstr-1,k)=Rig(iend,jstr,k)
        endif
        if (EASTERN_EDGE .and.
     &      NORTHERN_EDGE) then
          Rig(iend+1,jend+1,k)=Rig(iend,jend,k)
        endif
#   endif
#  endif
                                   ! Smooth Rig horizontally: use
        do j=jstr-1,jend           ! array Rig(:,:,0) as scratch. 
          do i=istr-1,iend
            Rig(i,j,0)=0.25*(Rig(i,j  ,k)+Rig(i+1,j  ,k)
     &                      +Rig(i,j+1,k)+Rig(i+1,j+1,k))
          enddo
        enddo
        do j=jstr,jend
          do i=istr,iend
            Rig(i,j,k)=0.25*(Rig(i,j  ,0)+Rig(i-1,j  ,0)
     &                      +Rig(i,j-1,0)+Rig(i-1,j-1,0))
          enddo
        enddo    !--> discard  Rig(:,:,0)
# endif /* SMOOTH_RIG */
# undef I_EXT_RANGE
# undef J_EXT_RANGE
!
! Compute "interior" viscosities and diffusivities everywhere
! as the superposition of three processes: local Richardson number
! instability due to resolved vertical shear, internal wave breaking,
! and double diffusion.
!
        do j=jstr,jend
          do i=istr,iend
# ifdef LMD_RIMIX
            cff=min(1., max(0., Rig(i,j,k)))  ! Compute mixing die
            nu_sx=1. - cff*cff                ! to shear instability
            nu_sx=nu_sx*nu_sx*nu_sx           ! and internal wave 
                                              ! breaking.
            Kv(i,j,k)=nuwm + nu0m*nu_sx
            Kt(i,j,k)=nuws + nu0s*nu_sx
# else
            Kv(i,j,k)=nuwm                    ! internal waves
            Kt(i,j,k)=nuws                    ! only
# endif
# if defined LMD_CONVEC && !defined LMD_KPP
            if (Rig(i,j,k).lt.0.) then        ! In the case of KPP
              Kv(i,j,k)=Kv(i,j,k) + nu0c      ! convective adjustment
              Kt(i,j,k)=Kt(i,j,k) + nu0c      ! will be added later
            endif                             ! and only below HBL
# endif
            Ks(i,j,k)=Kt(i,j,k)
          enddo
        enddo
# ifdef LMD_DDMIX
!
! Compute double-diffusive mixing.
!-----------------------------------------
! It can occur when vertical gradient of density is stable but the
! vertical gradient of salinity (salt figering) or temperature
! (diffusive convection) is unstable.
!
! Compute the ratio of thermal expansion and saline contraction
! coefficients at horizontal and vertical W-points.
!
! Compute double-diffusive density ratio, Rrho.
!
        do j=jstr,jend
          do i=istr,iend
            Tt=0.5*(t(i,j,k,tind,itemp)+t(i,j,k+1,tind,itemp))
            Ts=0.5*(t(i,j,k,tind,isalt)+t(i,j,k+1,tind,isalt))
     &                                                 -Smean
            Tp=-z_w(i,j,k)
 
            alfaobeta=A0+Tt*(A1+Tt*(A2+Tt*(A3+Tt*A4)))
     &                  +Ts*(B0+Tt*B1+Ts*C0)
     &                  +Tp*(D0+Tt*(D1+Tt*D2)+Ts*E0
     &                      +Tp*(Ts*F0+Tt*Tt*G0+Tp*H0))
 
            ddDT=t(i,j,k+1,tind,itemp)-t(i,j,k,tind,itemp)
            ddDS=t(i,j,k+1,tind,isalt)-t(i,j,k,tind,isalt)
            ddDS=sign(1.,ddDS)*max(abs(ddDS),eps)
            Rrho=alfaobeta*ddDT/ddDS
 
!
! Salt fingering case.
!--------------------------
!
            if (Rrho.gt.1. .and. ddDS.gt.0.) then
!
! Compute interior diffusivity for double diffusive mixing
! of salinity. Upper bound "Rrho" by "Rrho0"; (lmd_Rrho0=1.9,
! lmd_nuf=0.001).
!
              Rrho=min(Rrho,lmd_Rrho0)
              nu_dds=1.-((Rrho-1.)/(lmd_Rrho0-1.))**2
              nu_dds=lmd_nuf*nu_dds*nu_dds*nu_dds
!
! Compute interior diffusivity for double diffusive mixing
! of temperature (lmd_fdd=0.7).
!
              nu_ddt=lmd_fdd*nu_dds
!
! Diffusive convection case.
!-------------------------------
!
            elseif (Rrho.lt.1. .and. Rrho.gt.0. .and.
     &                               ddDS.lt.0.) then
!
! Compute interior diffusivity for double diffusive mixing of
! temperature (Marmorino and Caldwell, 1976); (lmd_nu=1.5e-6,
! lmd_tdd1=0.909, lmd_tdd2=4.6, lmd_tdd3=0.54).
!
              nu_ddt=lmd_nu*lmd_tdd1*
     &               exp(lmd_tdd2*exp(-lmd_tdd3*((1./Rrho)-1.)))
!
! Compute interior diffusivity for double diffusive mixing
! of salinity (lmd_sdd1=0.15, lmd_sdd2=1.85, lmd_sdd3=0.85).
!
              if (Rrho.lt.0.5) then
                nu_dds=nu_ddt*lmd_sdd1*Rrho
              else
                nu_dds=nu_ddt*(lmd_sdd2*Rrho-lmd_sdd3)
              endif
            else
              nu_ddt=0.
              nu_dds=0.
            endif
!
! Add double diffusion contribution to temperature and salinity
! mixing coefficients.
!
            Kt(i,j,k)=Kt(i,j,k)+nu_ddt
            Ks(i,j,k)=Ks(i,j,k)+nu_dds
          enddo
        enddo
# endif /* LMD_DDMIX */
      enddo     ! <-- k
# ifdef LMD_KPP
!
! Pad out surface and bottom values for lmd_blmix calculations.
! The interior values used here may not be the best values to
! use for the padding.
!
      do j=jstr,jend
        do i=istr,iend
          Kv(i,j,N)=Kv(i,j,N-1)
          Ks(i,j,N)=Ks(i,j,N-1)
          Kt(i,j,N)=Kt(i,j,N-1)
          Kv(i,j,0)=Kv(i,j,  1)
          Ks(i,j,0)=Ks(i,j,  1)
          Kt(i,j,0)=Kt(i,j,  1)
        enddo
      enddo
# else
!
! Finalize: Copy everything into shared arrays:
!
      do k=1,N-1
        do j=jstr,jend
          do i=istr,iend
            Akv(i,j,k)=Kv(i,j,k)
            Akt(i,j,k,itemp)=Kt(i,j,k)
# ifdef SALINITY
            Akt(i,j,k,isalt)=Ks(i,j,k)
# endif
          enddo
        enddo
      enddo
 
# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
      call exchange_w3d_tile (istr,iend,jstr,jend, Akv)
      call exchange_w3d_tile (istr,iend,jstr,jend,
     &                     Akt(START_2D_ARRAY,0,itemp))
#  ifdef SALINITY
      call exchange_w3d_tile (istr,iend,jstr,jend,
     &                     Akt(START_2D_ARRAY,0,isalt))
#  endif
# endif
 
# endif /* LMD_KPP */
      return
      end
#else
      subroutine lmd_vmix_empty
      end
#endif /* LMD_MIXING */
 
 
