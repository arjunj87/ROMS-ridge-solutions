#include "cppdefs.h"
#if defined DIAGNOSTICS_VRT


      subroutine set_diags_vrt (tile)
      implicit none
      integer tile
# include "param.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"
      call set_diags_vrt_tile (istr,iend,jstr,jend,tile)
      return
      end
 
      subroutine set_diags_vrt_tile (istr,iend,jstr,jend,tile)
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
# include "param.h"
      integer istr,iend,jstr,jend, i,j, ilc, iflux,
     &     imin,imax,jmin,jmax,tile
# ifdef SOLVE3D
     &                       , itrc, k
# endif
      real cff,cff1,cff2,
#ifdef DIAGNOSTICS_UV
     &     wrkXadv(PRIVATE_2D_SCRATCH_ARRAY,2),
     &     wrkYadv(PRIVATE_2D_SCRATCH_ARRAY,2),
     &     wrkCor(PRIVATE_2D_SCRATCH_ARRAY,2),
     &     wrkPrsgrd(PRIVATE_2D_SCRATCH_ARRAY,2),
     &     wrkHmix(PRIVATE_2D_SCRATCH_ARRAY,2),
     &     wrkVmix(PRIVATE_2D_SCRATCH_ARRAY,2),
     &     wrkVmix2(PRIVATE_2D_SCRATCH_ARRAY,2),
     &     wrkrate(PRIVATE_2D_SCRATCH_ARRAY,2),
# endif
     &     dH(N), jstri(2), istri(2),
     &     jendi(2), iendi(2)
# include "scalars.h"
# include "ncvars.h"
# include "grid.h"
# include "ocean2d.h"
# include "ocean3d.h"

#if defined DIAGNOSTICS_UV
#include "diags_uv.h"
# endif

#include "diags_vrt.h"

#include "compute_auxiliary_bounds.h"




! ------- ------------- ------ ------ - -----

        jstri(1) = jstr-1
        jstri(2) = jstr

        istri(1) = istr
        istri(2) = istr-1

        jendi(1) = jend
        jendi(2) = jend

        iendi(1) = iend
        iendi(2) = iend

        !jstri(1) = jstrR
        !jstri(2) = jstrR

        !istri(1) = istrR
        !istri(2) = istrR

        !jendi(1) = jendR
        !jendi(2) = jendR

        !iendi(1) = iendR
        !iendi(2) = iendR


#ifdef DIAGNOSTICS_UV

        do itrc=1,2
         do j=jstri(itrc),jendi(itrc)
          do i=istri(itrc),iendi(itrc)


            if (itrc.eq.1) then
              do k=1,N
               dH(k) = (Hz(i,j,k)+Hz(i-1,j,k))/2.
              enddo
            else          
              do k=1,N
               dH(k) = (Hz(i,j,k)+Hz(i,j-1,k))/2.
              enddo
            endif

           ! ------- ------------- ------ ------ - -----
           wrkXadv(i,j,itrc)= MXadv(i,j,1,itrc)*dH(1)

           do k=2,N
             wrkXadv(i,j,itrc)= wrkXadv(i,j,itrc)+MXadv(i,j,k,itrc)*dH(k)
           enddo
           ! ------- ------------- ------ ------ - -----
           wrkYadv(i,j,itrc)= MYadv(i,j,1,itrc)*dH(1)

           do k=2,N
             wrkYadv(i,j,itrc)= wrkYadv(i,j,itrc)+MYadv(i,j,k,itrc)*dH(k)
           enddo
           ! ------- ------------- ------ ------ - -----
           wrkCor(i,j,itrc)= MCor(i,j,1,itrc)*dH(1)

           do k=2,N
             wrkCor(i,j,itrc)= wrkCor(i,j,itrc)+MCor(i,j,k,itrc)*dH(k)
           enddo
           ! ------- ------------- ------ ------ - -----
           wrkPrsgrd(i,j,itrc)= MPrsgrd(i,j,1,itrc)*dH(1)

           do k=2,N
             wrkPrsgrd(i,j,itrc)= wrkPrsgrd(i,j,itrc)+MPrsgrd(i,j,k,itrc)*dH(k)
           enddo
           ! ------- ------------- ------ ------ - -----
           wrkHmix(i,j,itrc)= MHmix(i,j,1,itrc)*dH(1)

           do k=2,N
             wrkHmix(i,j,itrc)= wrkHmix(i,j,itrc)+MHmix(i,j,k,itrc)*dH(k)
           enddo
           ! ------- ------------- ------ ------ - -----
           wrkVmix(i,j,itrc)= MVmix(i,j,1,itrc)*dH(1)

           do k=2,N
             wrkVmix(i,j,itrc)= wrkVmix(i,j,itrc)+MVmix(i,j,k,itrc)*dH(k)
           enddo
           ! ------- ------------- ------ ------ - -----
           wrkVmix2(i,j,itrc)= MVmix2(i,j,1,itrc)*dH(1)

           do k=2,N
             wrkVmix2(i,j,itrc)= wrkVmix2(i,j,itrc)+MVmix2(i,j,k,itrc)*dH(k)
           enddo
           ! ------- ------------- ------ ------ - -----
           wrkrate(i,j,itrc)= Mrate(i,j,1,itrc)*dH(1)

           do k=2,N
             wrkrate(i,j,itrc)= wrkrate(i,j,itrc)+Mrate(i,j,k,itrc)*dH(k)
           enddo


          enddo
         enddo
        enddo



#if defined BOB
        do itrc=1,2

        if (WESTERN_EDGE) then        
          do j=jstr,jend
             wrkXadv(imin-1,j,itrc) = wrkXadv(imin,j,itrc)
             wrkYadv(imin-1,j,itrc) = wrkYadv(imin,j,itrc)
             wrkCor(imin-1,j,itrc) = wrkCor(imin,j,itrc)
             wrkPrsgrd(imin-1,j,itrc) = wrkPrsgrd(imin,j,itrc)
             wrkHmix(imin-1,j,itrc) = wrkHmix(imin,j,itrc)
             wrkVmix(imin-1,j,itrc) = wrkVmix(imin,j,itrc)
             wrkVmix2(imin-1,j,itrc) = wrkVmix2(imin,j,itrc)
             wrkrate(imin-1,j,itrc) = wrkrate(imin,j,itrc)
          enddo
        endif
        if (EASTERN_EDGE) then
          do j=jstr,jend
             wrkXadv(imax+1,j,itrc) = wrkXadv(imax,j,itrc)
             wrkYadv(imax+1,j,itrc) = wrkYadv(imax,j,itrc)
             wrkCor(imax+1,j,itrc) = wrkCor(imax,j,itrc)
             wrkPrsgrd(imax+1,j,itrc) = wrkPrsgrd(imax,j,itrc)
             wrkHmix(imax+1,j,itrc) = wrkHmix(imax,j,itrc)
             wrkVmix(imax+1,j,itrc) = wrkVmix(imax,j,itrc)
             wrkVmix2(imax+1,j,itrc) = wrkVmix2(imax,j,itrc)
             wrkrate(imax+1,j,itrc) = wrkrate(imax,j,itrc)
          enddo
        endif

        if (SOUTHERN_EDGE) then        
          do i=istr,iend
             wrkXadv(i,jmin-1,itrc) = wrkXadv(i,jmin,itrc)
             wrkYadv(i,jmin-1,itrc) = wrkYadv(i,jmin,itrc)
             wrkCor(i,jmin-1,itrc) = wrkCor(i,jmin,itrc)
             wrkPrsgrd(i,jmin-1,itrc) = wrkPrsgrd(i,jmin,itrc)
             wrkHmix(i,jmin-1,itrc) = wrkHmix(i,jmin,itrc)
             wrkVmix(i,jmin-1,itrc) = wrkVmix(i,jmin,itrc)
             wrkVmix2(i,jmin-1,itrc) = wrkVmix2(i,jmin,itrc)
             wrkrate(i,jmin-1,itrc) = wrkrate(i,jmin,itrc)
          enddo
        endif
        if (NORTHERN_EDGE) then
          do i=istr,iend
             wrkXadv(i,jmax+1,itrc) = wrkXadv(i,jmax,itrc)
             wrkYadv(i,jmax+1,itrc) = wrkYadv(i,jmax,itrc)
             wrkCor(i,jmax+1,itrc) = wrkCor(i,jmax,itrc)
             wrkPrsgrd(i,jmax+1,itrc) = wrkPrsgrd(i,jmax,itrc)
             wrkHmix(i,jmax+1,itrc) = wrkHmix(i,jmax,itrc)
             wrkVmix(i,jmax+1,itrc) = wrkVmix(i,jmax,itrc)
             wrkVmix2(i,jmax+1,itrc) = wrkVmix2(i,jmax,itrc)
             wrkrate(i,jmax+1,itrc) = wrkrate(i,jmax,itrc)
          enddo
        endif

        enddo

#endif
#endif


      if (WESTERN_EDGE) then     ! Restrict extended ranges one
        imin=istr+1               ! point inward near the physical
      else                       ! boundary. Note that this version
        imin=istr            ! of code is suitable for MPI
      endif                      ! configuration too.

      if (SOUTHERN_EDGE) then
          jmin=jstr+1
      else
          jmin=jstr
      endif






         do j=jmin,jend
          do i=imin,iend      

           cff1 = 0.25*(pm(i,j) + pm(i-1,j) + pm(i,j-1) + pm(i-1,j-1))
           cff2 = 0.25*(pn(i,j) + pn(i-1,j) + pn(i,j-1) + pn(i-1,j-1))

            vrtXadv(i,j) = (wrkXadv(i,j,2) - wrkXadv(i-1,j,2)) * cff1
     &                     - (wrkXadv(i,j,1) - wrkXadv(i,j-1,1)) * cff2


            vrtYadv(i,j) = (wrkYadv(i,j,2) - wrkYadv(i-1,j,2)) * cff1
     &                     - (wrkYadv(i,j,1) - wrkYadv(i,j-1,1)) * cff2


            vrtCor(i,j) = (wrkCor(i,j,2) - wrkCor(i-1,j,2)) * cff1
     &                     - (wrkCor(i,j,1) - wrkCor(i,j-1,1)) * cff2


            vrtPrsgrd(i,j) = (wrkPrsgrd(i,j,2) - wrkPrsgrd(i-1,j,2)) * cff1
     &                     - (wrkPrsgrd(i,j,1) - wrkPrsgrd(i,j-1,1)) * cff2


            vrtHmix(i,j) = (wrkHmix(i,j,2) - wrkHmix(i-1,j,2)) * cff1
     &                     - (wrkHmix(i,j,1) - wrkHmix(i,j-1,1)) * cff2


            vrtVmix(i,j) = (wrkVmix(i,j,2) - wrkVmix(i-1,j,2)) * cff1
     &                     - (wrkVmix(i,j,1) - wrkVmix(i,j-1,1)) * cff2



            vrtVmix2(i,j) = (wrkVmix2(i,j,2) - wrkVmix2(i-1,j,2)) * cff1
     &                     - (wrkVmix2(i,j,1) - wrkVmix2(i,j-1,1)) * cff2


            vrtrate(i,j) = (wrkrate(i,j,2) - wrkrate(i-1,j,2)) * cff1
     &                     - (wrkrate(i,j,1) - wrkrate(i,j-1,1)) * cff2

            vrtWind(i,j) = (wrkWind(i,j,2) - wrkWind(i-1,j,2)) * cff1
     &                     - (wrkWind(i,j,1) - wrkWind(i,j-1,1)) * cff2

            vrtDrag(i,j) = (wrkDrag(i,j,2) - wrkDrag(i-1,j,2)) * cff1
     &                     - (wrkDrag(i,j,1) - wrkDrag(i,j-1,1)) * cff2


          enddo
         enddo

! ------- ------------- ------ ------ - -----
      return
      end
#else /* DIAGNOSTICS_VRT*/
      subroutine set_diags_vrt_empty
      end
#endif /* DIAGNOSTICS_VRT */










