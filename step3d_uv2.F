#include "cppdefs.h"
#ifdef SOLVE3D

# ifdef ADV_ISONEUTRAL
# define DIFF_COEFF
# define GAMMA 0.0833333333 
# endif

# define NO_SLIP_BOTTOM
 
      subroutine step3d_uv2 (tile)
      implicit none
      integer tile
# include "param.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"
      call step3d_uv2_tile (istr,iend,jstr,jend, A2d(1,1), A2d(1,2),
     &                                                     A2d(1,3))
      return
      end
 
      subroutine step3d_uv2_tile (istr,iend,jstr,jend, CF,FC,DC)
      implicit none
      integer istr,iend,jstr,jend, i,j,k
# include "param.h"
      real FC(PRIVATE_1D_SCRATCH_ARRAY,0:N), cff,
     &     CF(PRIVATE_1D_SCRATCH_ARRAY,0:N),
     &     DC(PRIVATE_1D_SCRATCH_ARRAY,0:N)
      real diff_clip
# include "grid.h"
# include "ocean2d.h"
# include "ocean3d.h"
# include "coupling.h"
# include "mixing.h"
# include "scalars.h"
# include "sources.h"
# ifdef DIAGNOSTICS_UV
#include "diags_uv.h"
# endif

# ifdef DIAGNOSTICS_VRT
#include "diags_vrt.h"
# endif

# ifdef DIAGNOSTICS_PV
# include "diags_pv.h"
# endif

!
# include "compute_auxiliary_bounds.h"
!

      do j=jstr,jend
        do i=istrU,iend                        !--> coupling
          CF(i,0)=0.5*(Hz(i,j,N)+Hz(i-1,j,N))
          DC(i,0)=u(i,j,N,nnew)
          u(i,j,N,nnew)=u(i,j,N,nnew)/CF(i,0)
        enddo
        do k=N-1,1,-1
          do i=istrU,iend
            cff=0.5*(Hz(i,j,k)+Hz(i-1,j,k))
            CF(i,0)=CF(i,0) +cff
            DC(i,0)=DC(i,0) +u(i,j,k,nnew)
            u(i,j,k,nnew)=u(i,j,k,nnew)/cff
          enddo
        enddo
        do i=istrU,iend
          DC(i,0)=( DC(i,0)*dn_u(i,j)-DU_avg1(i,j)
     &                        )/(CF(i,0)*dn_u(i,j))
        enddo
        do k=1,N
          do i=istrU,iend
            u(i,j,k,nnew)=(u(i,j,k,nnew)-DC(i,0))
#  ifdef MASKING
     &                                *umask(i,j)
#  endif
# ifdef DIAGNOSTICS_UV_old
            MVmix(i,j,k,1)= MVmix(i,j,k,1)- DC(i,0) *
     &                           0.5*(Hz(i,j,k)+Hz(i-1,j,k)) 
     &                                 *dm_u(i,j)*dn_u(i,j)/dt
#  ifdef MASKING
     &                            *umask(i,j)
#  endif 
# endif /* DIAGNOSTICS_UV */
# if defined DIAGNOSTICS_PV && !defined DIAGNOSTICS_UV
              if (k.eq.1) then
                pvForc_bot(i,j,1) = pvForc_bot(i,j,1) - DC(i,0) *
     &                           0.5*(Hz(i,j,k)+Hz(i-1,j,k)) 
     &                                 *dm_u(i,j)*dn_u(i,j)/dt
              endif

              if (k.eq.N) then
                pvForc_sfc(i,j,1) = pvForc_sfc(i,j,1) - DC(i,0) *
     &                           0.5*(Hz(i,j,k)+Hz(i-1,j,k)) 
     &                                 *dm_u(i,j)*dn_u(i,j)/dt
              endif
# endif /* DIAGNOSTICS_PV */



          enddo
        enddo

        if (j.ge.jstrV) then
          do i=istr,iend
            CF(i,0)=0.5*(Hz(i,j,N)+Hz(i,j-1,N))
            DC(i,0)=v(i,j,N,nnew)
            v(i,j,N,nnew)=v(i,j,N,nnew)/CF(i,0)
          enddo
          do k=N-1,1,-1
            do i=istr,iend
              cff=0.5*(Hz(i,j,k)+Hz(i,j-1,k))
              CF(i,0)=CF(i,0) +cff
              DC(i,0)=DC(i,0) +v(i,j,k,nnew)
              v(i,j,k,nnew)=v(i,j,k,nnew)/cff
            enddo
          enddo      !--> discard FC, keep DC,CF(:,0)
          do i=istr,iend
            DC(i,0)=( DC(i,0)*dm_v(i,j)-DV_avg1(i,j)
     &                          )/(CF(i,0)*dm_v(i,j))
          enddo
          do k=1,N
            do i=istr,iend
              v(i,j,k,nnew)=(v(i,j,k,nnew)-DC(i,0))
#  ifdef MASKING
     &                                  *vmask(i,j)
#  endif
# ifdef DIAGNOSTICS_UV_old
              MVmix(i,j,k,2)= MVmix(i,j,k,2) - DC(i,0)*
     &                             0.5*(Hz(i,j,k)+Hz(i,j-1,k)) 
     &                                   *dm_v(i,j)*dn_v(i,j)/dt
#  ifdef MASKING
     &                              *vmask(i,j)
#  endif
# endif /* DIAGNOSTICS_UV */

# if defined DIAGNOSTICS_PV && !defined DIAGNOSTICS_UV
              if (k.eq.1) then
                pvForc_bot(i,j,2) = pvForc_bot(i,j,2)- DC(i,0) *
     &                             0.5*(Hz(i,j,k)+Hz(i,j-1,k)) 
     &                                   *dm_v(i,j)*dn_v(i,j)/dt
              endif

              if (k.eq.N) then
                pvForc_sfc(i,j,2) = pvForc_sfc(i,j,2) - DC(i,0) *
     &                             0.5*(Hz(i,j,k)+Hz(i,j-1,k)) 
     &                                   *dm_v(i,j)*dn_v(i,j)/dt
              endif

# endif /* DIAGNOSTICS_PV */


            enddo
          enddo       !--> discard DC,CF(:,0)
        endif
      enddo      ! <-- j


!      write(*,*) 'step 3d: Istr,Iend,Jstr,Jend', Istr,Iend,Jstr,Jend

# ifdef DIAGNOSTICS_UV
      do k=1,N
        do j=Jstr,Jend
          do i=Istr,Iend

           Mrate(i,j,k,1)=( 0.5*(Hz(i,j,k)+Hz(i-1,j,k))*u(i,j,k,nnew)
     &           -0.5*(Hz_bak(i,j,k)+Hz_bak(i-1,j,k))*u(i,j,k,nstp))
     &                                           *dm_u(i,j)*dn_u(i,j)/dt
#  ifdef MASKING
     &                     *umask(i,j)
#  endif
           Mrate(i,j,k,2)=( 0.5*(Hz(i,j,k)+Hz(i,j-1,k))*v(i,j,k,nnew)
     &           -0.5*(Hz_bak(i,j,k)+Hz_bak(i,j-1,k))*v(i,j,k,nstp))
     &                                           *dm_v(i,j)*dn_v(i,j)/dt
#  ifdef MASKING
     &                     *vmask(i,j)
#  endif

! Divide all diagnostic terms by the cell volume Hz/(pm*pn).
! There after the unit of diag terms are :
! (unit of velocity) * s-1  =  m * s-2

            cff=0.5*(pm(i,j)+pm(i-1,j))
     &             *(pn(i,j)+pn(i-1,j))
     &             /(Hz(i,j,k)+Hz(i-1,j,k))
            
            MXadv(i,j,k,1)=MXadv(i,j,k,1)*cff
!            MXadv(i,j,k,1)=0.5*(Hz(i,j,k)+Hz(i-1,j,k))
            MYadv(i,j,k,1)=MYadv(i,j,k,1)*cff
            MVadv(i,j,k,1)=MVadv(i,j,k,1)*cff
            MCor(i,j,k,1)=MCor(i,j,k,1)*cff
            MPrsgrd(i,j,k,1)=MPrsgrd(i,j,k,1)*cff
# ifdef MASKING
     &                              *umask(i,j)
# endif  
            MHmix(i,j,k,1)=MHmix(i,j,k,1)*cff
            MVmix(i,j,k,1)=MVmix(i,j,k,1)*cff
            MVmix2(i,j,k,1)=MVmix2(i,j,k,1)*cff
            Mrate(i,j,k,1)=Mrate(i,j,k,1)*cff
!

            cff=0.5*(pm(i,j)+pm(i,j-1))
     &             *(pn(i,j)+pn(i,j-1))
     &             /(Hz(i,j,k)+Hz(i,j-1,k))

            MXadv(i,j,k,2)=MXadv(i,j,k,2)*cff
!            MXadv(i,j,k,2)=0.5*(Hz(i,j,k)+Hz(i,j-1,k))
            MYadv(i,j,k,2)=MYadv(i,j,k,2)*cff
            MVadv(i,j,k,2)=MVadv(i,j,k,2)*cff
            MCor(i,j,k,2)=MCor(i,j,k,2)*cff
            MPrsgrd(i,j,k,2)=MPrsgrd(i,j,k,2)*cff
# ifdef MASKING
     &                              *vmask(i,j)
# endif  
            MHmix(i,j,k,2)=MHmix(i,j,k,2)*cff
            MVmix(i,j,k,2)=MVmix(i,j,k,2)*cff
            MVmix2(i,j,k,2)=MVmix2(i,j,k,2)*cff
            Mrate(i,j,k,2)=Mrate(i,j,k,2)*cff
          enddo
        enddo
      enddo

! At this point we have Mrate = MXadv + MYadv + MVadv + MCor + MPrsgrd + MHmix + MVmix
! except at boundaries (1pt only _ I still need to check BC)


# elif defined DIAGNOSTICS_VRT

        do j=Jstr,Jend
          do i=Istr,Iend

            k=1
    
            wrkrate(i,j,1)=( 0.5*(Hz(i,j,k)+Hz(i-1,j,k))*u(i,j,k,nnew)
     &           -0.5*(Hz_bak(i,j,k)+Hz_bak(i-1,j,k))*u(i,j,k,nstp))
     &                                           *dm_u(i,j)*dn_u(i,j)/dt
#  ifdef MASKING
     &                     *umask(i,j)
#  endif
            wrkrate(i,j,2)=( 0.5*(Hz(i,j,k)+Hz(i,j-1,k))*v(i,j,k,nnew)
     &           -0.5*(Hz_bak(i,j,k)+Hz_bak(i,j-1,k))*v(i,j,k,nstp))
     &                                           *dm_v(i,j)*dn_v(i,j)/dt
#  ifdef MASKING
     &                     *vmask(i,j)
#  endif

            do k=2,N

                wrkrate(i,j,1) = wrkrate(i,j,1) + 
     &           ( 0.5*(Hz(i,j,k)+Hz(i-1,j,k))*u(i,j,k,nnew)
     &           -0.5*(Hz_bak(i,j,k)+Hz_bak(i-1,j,k))*u(i,j,k,nstp))
     &                                           *dm_u(i,j)*dn_u(i,j)/dt
#  ifdef MASKING
     &                     *umask(i,j)
#  endif
                wrkrate(i,j,2) = wrkrate(i,j,2) +
     &           ( 0.5*(Hz(i,j,k)+Hz(i,j-1,k))*v(i,j,k,nnew)
     &           -0.5*(Hz_bak(i,j,k)+Hz_bak(i,j-1,k))*v(i,j,k,nstp))
     &                                           *dm_v(i,j)*dn_v(i,j)/dt
#  ifdef MASKING
     &                     *vmask(i,j)
#  endif

            enddo

! Divide all diagnostic terms by (pm*pn).
! There after the unit of these terms are :
!  s-2

            cff=0.25*(pm(i,j)+pm(i-1,j))
     &             *(pn(i,j)+pn(i-1,j))

            
            wrkXadv(i,j,1)=wrkXadv(i,j,1)*cff
            wrkYadv(i,j,1)=wrkYadv(i,j,1)*cff
            wrkCor(i,j,1)=wrkCor(i,j,1)*cff
            wrkPrsgrd(i,j,1)=wrkPrsgrd(i,j,1)*cff
# ifdef MASKING
     &                              *umask(i,j)
# endif  
            wrkHmix(i,j,1)=wrkHmix(i,j,1)*cff
            wrkVmix(i,j,1)=wrkVmix(i,j,1)*cff
            wrkVmix2(i,j,1)=wrkVmix2(i,j,1)*cff
#  ifdef MASKING
     &                            *umask(i,j)
#  endif
             wrkrate(i,j,1)=wrkrate(i,j,1)*cff
!

            cff=0.25*(pm(i,j)+pm(i,j-1))
     &             *(pn(i,j)+pn(i,j-1))


            wrkXadv(i,j,2)=wrkXadv(i,j,2)*cff
            wrkYadv(i,j,2)=wrkYadv(i,j,2)*cff
            wrkCor(i,j,2)=wrkCor(i,j,2)*cff
            wrkPrsgrd(i,j,2)=wrkPrsgrd(i,j,2)*cff
# ifdef MASKING
     &                              *vmask(i,j)
# endif  
            wrkHmix(i,j,2)=wrkHmix(i,j,2)*cff
            wrkVmix(i,j,2)=wrkVmix(i,j,2)*cff
            wrkVmix2(i,j,2)=wrkVmix2(i,j,2)*cff
#  ifdef MASKING
     &                            *vmask(i,j)
#  endif
            wrkrate(i,j,2)=wrkrate(i,j,2)*cff
          enddo
        enddo


# endif /* DIAGNOSTICS_UV */





# if defined DIAGNOSTICS_PV && defined DIAGNOSTICS_UV

        do j=Jstr,Jend
          do i=Istr,Iend

            pvForc_sfc(i,j,1) = MHmix(i,j,N,1) 
     &                   + MVmix(i,j,N,1) + MVmix2(i,j,N,1)
            pvForc_sfc(i,j,2) = MHmix(i,j,N,2) 
     &                   + MVmix(i,j,N,2) + MVmix2(i,j,N,2)

            pvForc_bot(i,j,1) =MHmix(i,j,1,1) 
     &                   + MVmix(i,j,1,1) + MVmix2(i,j,1,1)
            pvForc_bot(i,j,2) = MHmix(i,j,1,2) 
     &                   + MVmix(i,j,1,2) + MVmix2(i,j,1,2)

          enddo
        enddo

# elif defined DIAGNOSTICS_PV
    
! Divide all diagnostic terms by (pm*pn).
! There after the unit of these terms are :
!  s-2
        do j=Jstr,Jend
          do i=Istr,Iend

            cff=0.5*(pm(i,j)+pm(i-1,j))
     &             *(pn(i,j)+pn(i-1,j))
     &             /(Hz(i,j,N)+Hz(i-1,j,N))

            pvForc_sfc(i,j,1) = pvForc_sfc(i,j,1)*cff

            cff=0.5*(pm(i,j)+pm(i-1,j))
     &             *(pn(i,j)+pn(i-1,j))
     &             /(Hz(i,j,1)+Hz(i-1,j,1))

            pvForc_bot(i,j,1) = pvForc_bot(i,j,1)*cff

            cff=0.5*(pm(i,j)+pm(i,j-1))
     &             *(pn(i,j)+pn(i,j-1))
     &             /(Hz(i,j,N)+Hz(i,j-1,N))

            pvForc_sfc(i,j,2) = pvForc_sfc(i,j,2)*cff

            cff=0.5*(pm(i,j)+pm(i,j-1))
     &             *(pn(i,j)+pn(i,j-1))
     &             /(Hz(i,j,1)+Hz(i,j-1,1))

            pvForc_bot(i,j,2) = pvForc_bot(i,j,2)*cff


          enddo
        enddo

# endif /* DIAGNOSTICS_PV */

# if defined DIAGNOSTICS_PV_SUP

! Divide all diagnostic terms by (pm*pn).
! There after the unit of these terms are :
!  s-2
        do j=Jstr,Jend
          do i=Istr,Iend

            cff=0.25*(pm(i,j)+pm(i-1,j))
     &             *(pn(i,j)+pn(i-1,j))

            pvForc_0(i,j,1) = pvForc_0(i,j,1)*cff
            pvForc_1(i,j,1) = pvForc_1(i,j,1)*cff
            pvForc_2(i,j,1) = pvForc_2(i,j,1)*cff
            pvForc_3(i,j,1) = pvForc_3(i,j,1)*cff

            cff=0.25*(pm(i,j)+pm(i,j-1))
     &             *(pn(i,j)+pn(i,j-1))

            pvForc_0(i,j,2) = pvForc_0(i,j,2)*cff
            pvForc_1(i,j,2) = pvForc_1(i,j,2)*cff
            pvForc_2(i,j,2) = pvForc_2(i,j,2)*cff
            pvForc_3(i,j,2) = pvForc_3(i,j,2)*cff

          enddo
        enddo

# endif /* DIAGNOSTICS_PV_SUP */

!
! Set PHYSICAL lateral boundary conditions.
!
      call u3dbc_tile (istr,iend,jstr,jend, DC)
      call v3dbc_tile (istr,iend,jstr,jend, DC) 
!
! Coupling 2D and 3D momentum equations:
!--------- -- --- -- -------- ----------
! Compute inaccurate vertical mean of the three-dimensional
! velocity field, then subtract it and replace it with the vertically
! integrated (barotropic) velocity field computed from the two-
! dimensional submodel. After that compute mass fluxes through grid
! box faces.
!
! Meaning of scratch variables below:
!
! DC(i,k) [where k=1,N]  height of grid box for U- or V-cell.
! DC(i,0) total depth of water column at horizontal U- or V-points.
! FC(i,0) vertically integrated mass flux, error, and correction for
!                       the new time step velocity u,v(:,:,:,nnew)
! CF(i,0) vertically integrated mass flux/error/correction for the
!             intermediate time step [n+1/2] mass fluxes FlxU,FlxV.
!
! This procedure also replaces instantaneous (with respect to the
! fast time step) vertically integrated (barotropic) velocities with
! their values based on fast-time-averaged mass fluxes. These are to
! be used as initial conditions for the barotropic mode at the new
! time step.
!
! Explanation of horizontal loop indices: in the case of periodic
! boundaries (in either direction) the coupling and computation of
! mass fluxes FlxU,FlxV is performed within the internal range of
! indices (excluding ghost zones), after that the ghost points for
! the newly computed arrays are exchanged; in the case of nonperiodic
! boundaries the coupling is done over the extended range of indices
! (that is including boundary points).
!
# ifdef EW_PERIODIC
#  define IU_RANGE istr,iend
#  define IV_RANGE istr,iend
# else
#  define IU_RANGE istr,iendR
#  define IV_RANGE istrR,iendR
# endif
 
# ifdef NS_PERIODIC
#  define JU_RANGE jstr,jend
# else
#  define JU_RANGE jstrR,jendR
# endif
 
      do j=JU_RANGE                  !--> XI-component
        do i=IU_RANGE
          DC(i,N)=0.5*(Hz(i,j,N)+Hz(i-1,j,N))*dn_u(i,j)
          DC(i,0)=DC(i,N)
          FC(i,0)=DC(i,N)*u(i,j,N,nnew)
        enddo
        do k=N-1,1,-1
          do i=IU_RANGE
            DC(i,k)=0.5*(Hz(i,j,k)+Hz(i-1,j,k))*dn_u(i,j)
            DC(i,0)=DC(i,0)+DC(i,k)
            FC(i,0)=FC(i,0)+DC(i,k)*u(i,j,k,nnew)
          enddo
        enddo
        do i=IU_RANGE
          DC(i,0)=1./DC(i,0)
          ubar(i,j,knew)=DC(i,0)*DU_avg1(i,j)
          FC(i,0)=DC(i,0)*(FC(i,0)-DU_avg1(i,j))
!>
          CF(i,0)=0.
        enddo          !--> keep DC(:,0) for later
        do k=N,1,-1
          do i=IU_RANGE
            u(i,j,k,nnew)=(u(i,j,k,nnew)-FC(i,0))
# ifdef MASKING
     &                                  *umask(i,j)
# endif

c--# define TR
# ifdef TR
            CF(i,k)=0.5*DC(i,k)*(u(i,j,k,nstp)+u(i,j,k,nnew))
#else
#define EPSIL 0.125
#define DELTA 0.75
 
            CF(i,k)=DELTA*FlxU(i,j,k) + EPSIL*DC(i,k)
     &                         *(u(i,j,k,nstp)+u(i,j,k,nnew))
#endif
            CF(i,0)=CF(i,0)+CF(i,k)
          enddo
        enddo
        do i=IU_RANGE
          CF(i,0)=DC(i,0)*(CF(i,0)-DU_avg2(i,j))
        enddo
        do k=1,N,+1
          do i=IU_RANGE
            FlxU(i,j,k)=CF(i,k)-DC(i,k)*CF(i,0)
#ifdef DIFF_COEFF
            diff3u(i,j,k) = sqrt( GAMMA*abs(u(i,j,k,nrhs))*dm_u(i,j)**3 )
#endif
          enddo
        enddo
 
        if (j.ge.jstr) then         !--> ETA-component
          do i=IV_RANGE
            DC(i,N)=0.5*(Hz(i,j,N)+Hz(i,j-1,N))*dm_v(i,j)
            DC(i,0)=DC(i,N)
            FC(i,0)=DC(i,N)*v(i,j,N,nnew)
          enddo
          do k=N-1,1,-1
            do i=IV_RANGE
              DC(i,k)=0.5*(Hz(i,j,k)+Hz(i,j-1,k))*dm_v(i,j)
              DC(i,0)=DC(i,0)+DC(i,k)
              FC(i,0)=FC(i,0)+DC(i,k)*v(i,j,k,nnew)
            enddo
          enddo
          do i=IV_RANGE
            DC(i,0)=1./DC(i,0)
            vbar(i,j,knew)=DC(i,0)*DV_avg1(i,j)
            FC(i,0)=DC(i,0)*(FC(i,0)-DV_avg1(i,j))
!>
            CF(i,0)=0.
          enddo        !--> keep DC(:,0) for later
          do k=N,1,-1
            do i=IV_RANGE
              v(i,j,k,nnew)=(v(i,j,k,nnew)-FC(i,0))
# ifdef MASKING
     &                                    *vmask(i,j)
# endif

# ifdef TR
              CF(i,k)=0.5*DC(i,k)*(v(i,j,k,nstp)+v(i,j,k,nnew))
# else
              CF(i,k)=DELTA*FlxV(i,j,k) + EPSIL*DC(i,k)
     &                           *(v(i,j,k,nstp)+v(i,j,k,nnew))
#  undef DELTA
#  undef EPSIL
# endif
              CF(i,0)=CF(i,0)+CF(i,k)
            enddo
          enddo
          do i=IV_RANGE
            CF(i,0)=DC(i,0)*(CF(i,0)-DV_avg2(i,j))
          enddo
          do k=1,N,+1
            do i=IV_RANGE
              FlxV(i,j,k)=CF(i,k)-DC(i,k)*CF(i,0)
#ifdef DIFF_COEFF
              diff3v(i,j,k) = sqrt( GAMMA*abs(v(i,j,k,nrhs))*dn_v(i,j)**3 )
#endif
            enddo
          enddo
        endif       !<-- j.ge.jstr
      enddo       ! <-- j
# undef IU_RANGE
# undef JU_RANGE
# undef IV_RANGE
!
! Exchange periodic boundaries and computational margins.
!
# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
      call exchange_u3d_tile (istr,iend,jstr,jend,
     &                                 u(START_2D_ARRAY,1,nnew))
      call exchange_v3d_tile (istr,iend,jstr,jend,
     &                                 v(START_2D_ARRAY,1,nnew))      
      
      call exchange_u3d_tile (istr,iend,jstr,jend, FlxU)
      call exchange_v3d_tile (istr,iend,jstr,jend, FlxV)

#ifdef DIFF_COEFF          
      call exchange_u3d_tile (istr,iend,jstr,jend, diff3U)
      call exchange_v3d_tile (istr,iend,jstr,jend, diff3V)
# endif 

! Need to exchange momentum terms for computation of vorticity terms
# if defined DIAGNOSTICS_UV && defined DIAGNOSTICS_VRT

      call exchange_u3d_tile (istr,iend,jstr,jend,
     &                                 MXadv(START_2D_ARRAY,1,1))
      call exchange_u3d_tile (istr,iend,jstr,jend,
     &                                 MYadv(START_2D_ARRAY,1,1))
      call exchange_u3d_tile (istr,iend,jstr,jend,
     &                                 MVadv(START_2D_ARRAY,1,1))
      call exchange_u3d_tile (istr,iend,jstr,jend,
     &                                 MPrsgrd(START_2D_ARRAY,1,1))
      call exchange_u3d_tile (istr,iend,jstr,jend,
     &                                 MCor(START_2D_ARRAY,1,1))
      call exchange_u3d_tile (istr,iend,jstr,jend,
     &                                 MHmix(START_2D_ARRAY,1,1))
      call exchange_u3d_tile (istr,iend,jstr,jend,
     &                                 MVmix(START_2D_ARRAY,1,1))
      call exchange_u3d_tile (istr,iend,jstr,jend,
     &                                 MVmix2(START_2D_ARRAY,1,1))
      call exchange_u3d_tile (istr,iend,jstr,jend,
     &                                 Mrate(START_2D_ARRAY,1,1))

      call exchange_u2d_tile (istr,iend,jstr,jend,
     &                                 wrkWind(START_2D_ARRAY,1))
      call exchange_u2d_tile (istr,iend,jstr,jend,
     &                                 wrkDrag(START_2D_ARRAY,1))



      call exchange_v3d_tile (istr,iend,jstr,jend,
     &                                 MXadv(START_2D_ARRAY,1,2))
      call exchange_v3d_tile (istr,iend,jstr,jend,
     &                                 MYadv(START_2D_ARRAY,1,2))
      call exchange_v3d_tile (istr,iend,jstr,jend,
     &                                 MVadv(START_2D_ARRAY,1,2))
      call exchange_v3d_tile (istr,iend,jstr,jend,
     &                                 MPrsgrd(START_2D_ARRAY,1,2))
      call exchange_v3d_tile (istr,iend,jstr,jend,
     &                                 MCor(START_2D_ARRAY,1,2))
      call exchange_v3d_tile (istr,iend,jstr,jend,
     &                                 MHmix(START_2D_ARRAY,1,2))
      call exchange_v3d_tile (istr,iend,jstr,jend,
     &                                 MVmix(START_2D_ARRAY,1,2))
      call exchange_v3d_tile (istr,iend,jstr,jend,
     &                                 MVmix2(START_2D_ARRAY,1,2))
      call exchange_v3d_tile (istr,iend,jstr,jend,
     &                                 Mrate(START_2D_ARRAY,1,2))

      call exchange_v2d_tile (istr,iend,jstr,jend,
     &                                 wrkWind(START_2D_ARRAY,2))
      call exchange_v2d_tile (istr,iend,jstr,jend,
     &                                 wrkDrag(START_2D_ARRAY,2))


# elif defined DIAGNOSTICS_VRT   


      call exchange_u2d_tile (istr,iend,jstr,jend,
     &                                 wrkXadv(START_2D_ARRAY,1))
      call exchange_u2d_tile (istr,iend,jstr,jend,
     &                                 wrkYadv(START_2D_ARRAY,1))
      call exchange_u2d_tile (istr,iend,jstr,jend,
     &                                 wrkPrsgrd(START_2D_ARRAY,1))
      call exchange_u2d_tile (istr,iend,jstr,jend,
     &                                 wrkCor(START_2D_ARRAY,1))
      call exchange_u2d_tile (istr,iend,jstr,jend,
     &                                 wrkHmix(START_2D_ARRAY,1))
      call exchange_u2d_tile (istr,iend,jstr,jend,
     &                                 wrkVmix(START_2D_ARRAY,1))
      call exchange_u2d_tile (istr,iend,jstr,jend,
     &                                 wrkVmix2(START_2D_ARRAY,1))
      call exchange_u2d_tile (istr,iend,jstr,jend,
     &                                 wrkrate(START_2D_ARRAY,1))

      call exchange_u2d_tile (istr,iend,jstr,jend,
     &                                 wrkWind(START_2D_ARRAY,1))
      call exchange_u2d_tile (istr,iend,jstr,jend,
     &                                 wrkDrag(START_2D_ARRAY,1))



      call exchange_v2d_tile (istr,iend,jstr,jend,
     &                                 wrkXadv(START_2D_ARRAY,2))
      call exchange_v2d_tile (istr,iend,jstr,jend,
     &                                 wrkYadv(START_2D_ARRAY,2))
      call exchange_v2d_tile (istr,iend,jstr,jend,
     &                                 wrkPrsgrd(START_2D_ARRAY,2))
      call exchange_v2d_tile (istr,iend,jstr,jend,
     &                                 wrkCor(START_2D_ARRAY,2))
      call exchange_v2d_tile (istr,iend,jstr,jend,
     &                                 wrkHmix(START_2D_ARRAY,2))
      call exchange_v2d_tile (istr,iend,jstr,jend,
     &                                 wrkVmix(START_2D_ARRAY,2))
      call exchange_v2d_tile (istr,iend,jstr,jend,
     &                                 wrkVmix2(START_2D_ARRAY,2))
      call exchange_v2d_tile (istr,iend,jstr,jend,
     &                                 wrkrate(START_2D_ARRAY,2))

      call exchange_v2d_tile (istr,iend,jstr,jend,
     &                                 wrkWind(START_2D_ARRAY,2))
      call exchange_v2d_tile (istr,iend,jstr,jend,
     &                                 wrkDrag(START_2D_ARRAY,2))


# endif 

      call exchange_u2d_tile (istr,iend,jstr,jend,
     &                                ubar(START_2D_ARRAY,knew))
      call exchange_v2d_tile (istr,iend,jstr,jend,
     &                                vbar(START_2D_ARRAY,knew))
# endif
      return
      end
#else
      subroutine step3d_uv_empty
      return
      end
#endif /* SOLVE3D */
 
