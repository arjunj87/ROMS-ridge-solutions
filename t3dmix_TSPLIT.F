#include "cppdefs.h"

#define Z_CORR
 
      subroutine t3dmix_ISO (tile,niter)
      implicit none
      integer tile, itrc, niter
#include "param.h"
#include "private_scratch.h"
#include "compute_tile_bounds.h"
      do itrc=1,NT
        call t3dmix_ISO_tile (istr,iend,jstr,jend, itrc, A3d(1,1),
#ifdef Z_CORR
     &        A3d(1,2),
#endif          
     &                                    A2d(1, 1), A2d(1, 2),
     &                          A2d(1,3), A2d(1, 5), A2d(1, 7),
     &                          A2d(1,9)
#ifdef Z_CORR
     &        ,A2d(1,11), A2d(1,12), A2d(1,13), A2d(1,14)  
#endif           
     &  , niter)     
      enddo
      return
      end
 
      subroutine t3dmix_ISO_tile (istr,iend,jstr,jend, itrc, LapT, 
#ifdef Z_CORR
     &       Akz,
#endif        
     &       FX,FE,FC,dTdz, dTdx,dTde
#ifdef Z_CORR
     &       ,FFC,CF,CD,DC
#endif      
     & ,niter)
!
!--------------------------------------------------------------------
!  Compute isopycnic horizontal mixing terms for tracer equations.
!-------------------------------------------------------------------- 
!
      implicit none
#include "param.h"
      integer itrc, istr,iend,jstr,jend, i,j,k,k1,k2
      integer n1,n2,n3
      real LapT(PRIVATE_2D_SCRATCH_ARRAY,0:N),
# ifdef Z_CORR
     &      Akz (PRIVATE_2D_SCRATCH_ARRAY,0:N),   
     &       CF (PRIVATE_1D_SCRATCH_ARRAY,0:N),
     &       DC (PRIVATE_1D_SCRATCH_ARRAY,0:N),
     &       CD (PRIVATE_1D_SCRATCH_ARRAY,0:N),      
     &       FFC(PRIVATE_1D_SCRATCH_ARRAY,0:N),        
# endif
     &       FX(PRIVATE_2D_SCRATCH_ARRAY),
     &       FE(PRIVATE_2D_SCRATCH_ARRAY),
     &       FC(PRIVATE_2D_SCRATCH_ARRAY,2),
     &     dTdz(PRIVATE_2D_SCRATCH_ARRAY,2),   cff,
     &     dTdx(PRIVATE_2D_SCRATCH_ARRAY,2),   cff1,
     &     cff3,cff4,
     &     dTde(PRIVATE_2D_SCRATCH_ARRAY,2),   cff2
#if defined DIAGNOSTICS_TS
     &     ,work(PRIVATE_2D_SCRATCH_ARRAY,N)
# endif
#if defined DIAGNOSTICS_PV && !defined DIAGNOSTICS_TS
     &     ,work(PRIVATE_2D_SCRATCH_ARRAY,2)
# endif
       real dts,cff5
       integer imin,imax,jmin,jmax,niter
#include "grid.h"
#include "ocean3d.h"
#include "mixing.h"
#include "climat.h"
#include "scalars.h"
#include "eos_vars.h"

# ifdef DIAGNOSTICS_TS
# include "diags_ts.h"
# endif
# ifdef DIAGNOSTICS_PV
# include "diags_pv.h"
# endif

!
#include "compute_auxiliary_bounds.h"
!
#ifdef MASKING
# define SWITCH *
#else
# define SWITCH !
#endif
 

       dts = (1./NSPLITDIF)*dt
       
       if( mod(niter,2) == 0 ) then
          n1 = nrhs
          n2 = nnew
          n3 = n1
       else
          n1 = nnew
          n2 = nrhs
          n3 = n1
       endif       
       
       if( niter .eq. 1 ) then
          n1 = nstp
          n2 = nrhs
          n3 = nnew
       endif




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
#  else
      imin=istr-1
      imax=iend+1
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
#  else
      jmin=jstr-1
      jmax=jend+1
#  endif
     
!
!  Add in horizontal biharmonic diffusion along rotated surfaces.
!  The biharmonic operator is computed by applying the rotated
!  Laplacian operator twice.
!
!                                            !!! WARNING: RECURSIVE
!  The FIRST rotated Laplacian operator:           BLOCKING SEQUENCE
!--------------------------------------------------------------------
!                   !  Vertical    dTdx,dTde(:,:,k2) k+1   rho-points
      k2=1          !  placement     dTdr,FC(:,:,k2) k+1/2   W-points
      do k=0,N,+1   !  of fields:  dTdx,dTde(:,:,k1) k     pho-points
       k1=k2        !                dTdr,FC(:,:,k1) k-1/2   W-points
       k2=3-k1
        if (k.lt.N) then
          do j=jmin,jmax
            do i=imin,imax+1
              cff=0.5*(pm(i,j)+pm(i-1,j)) SWITCH umask(i,j)
              dTdx(i,j,k2)=cff*(t(i,j,k+1,n1,itrc)-t(i-1,j,k+1,n1,itrc) )
            enddo
          enddo
          do j=jmin,jmax+1
            do i=imin,imax
              cff=0.5*(pn(i,j)+pn(i,j-1)) SWITCH vmask(i,j)
              dTde(i,j,k2)=cff*( t(i,j,k+1,n1,itrc)-t(i,j-1,k+1,n1,itrc) )  
            enddo
          enddo
        endif
 
        if (k.eq.0 .or. k.eq.N) then
          do j=jmin-1,jmax+1
            do i=imin-1,imax+1
              FC(i,j,k2)=0.0
              dTdz(i,j,k2)=0.0
            enddo
          enddo
        else
          do j=jmin-1,jmax+1
            do i=imin-1,imax+1
              cff  =  1. / ( z_r(i,j,k+1)-z_r(i,j,k) )
              dTdz(i,j,k2)=cff*(  t(i,j,k+1,n1,itrc)-t(i,j,k  ,n1,itrc) ) SWITCH rmask(i,j)
            enddo
          enddo
        endif
 
        if (k.gt.0) then
          do j=jmin,jmax
            do i=imin,imax+1
              FX(i,j)=0.5*diff3u(i,j,k)*(Hz(i,j,k)+Hz(i-1,j,k))
     &         *dn_u(i,j)*(   dTdx(i  ,j,k1)+0.25*Sx_u(i,j,k)*(
     &                dTdz(i-1,j,k1)+dTdz(i-1,j,k2)     
     &         +      dTdz(i  ,j,k1)+dTdz(i  ,j,k2)       
     &                                              )   )
            enddo
          enddo
 
          do j=jmin,jmax+1
            do i=imin,imax
              FE(i,j)=0.5*diff3v(i,j,k)*(Hz(i,j,k)+Hz(i,j-1,k))
     &        *dm_v(i,j)*( dTde(i,j,k1)+0.25*Sy_v(i,j,k)*(
     &             dTdz(i,j-1,k1)+dTdz(i,j-1,k2)     
     &         +   dTdz(i,j  ,k1)+dTdz(i,j  ,k2)       
     &                                              )   )
            enddo
          enddo
 
          if (k.lt.N) then
            do j=jmin,jmax
              do i=imin,imax


             cff1=Sx_u(i  ,j,k  )
             cff2=Sx_u(i+1,j,k+1)
             cff3=Sx_u(i  ,j,k+1)
             cff4=Sx_u(i+1,j,k  )

             cff=
     &       diff3u(i  ,j,k  )*cff1*(cff1*dTdz(i,j,k2)+dTdx(i  ,j,k1))
     &      +diff3u(i+1,j,k+1)*cff2*(cff2*dTdz(i,j,k2)+dTdx(i+1,j,k2))
     &      +diff3u(i  ,j,k+1)*cff3*(cff3*dTdz(i,j,k2)+dTdx(i  ,j,k2))
     &      +diff3u(i+1,j,k  )*cff4*(cff4*dTdz(i,j,k2)+dTdx(i+1,j,k1))

             cff1=Sy_v(i,j  ,k  )
             cff2=Sy_v(i,j+1,k+1)
             cff3=Sy_v(i,j  ,k+1)
             cff4=Sy_v(i,j+1,k  )

             cff=cff+(
     &     diff3v(i,j  ,k  )*cff1*(cff1*dTdz(i,j,k2)+dTde(i,j  ,k1))
     &   + diff3v(i,j+1,k+1)*cff2*(cff2*dTdz(i,j,k2)+dTde(i,j+1,k2))
     &   + diff3v(i,j  ,k+1)*cff3*(cff3*dTdz(i,j,k2)+dTde(i,j  ,k2))
     &   + diff3v(i,j+1,k  )*cff4*(cff4*dTdz(i,j,k2)+dTde(i,j+1,k1))
     &               )
 
             FC(i,j,k2)=0.25*cff
             
              enddo
            enddo
          endif
          do j=jmin,jmax
            do i=imin,imax
              LapT(i,j,k)=( pm(i,j)*pn(i,j)*( FX(i+1,j)-FX(i,j)
     &                                       +FE(i,j+1)-FE(i,j))
     &                     +FC(i,j,k2)-FC(i,j,k1)    )/Hz(i,j,k)
            enddo
          enddo
        endif   ! <-- k.gt.0
      enddo   ! --> k
!
!  Apply boundary conditions to the Laplacian
!  (for the cases other than periodic: closed or gradient)
!
#  ifndef EW_PERIODIC
        if (WESTERN_EDGE) then            ! Apply lateral boundary
          do k=1,N                        ! conditions to the first
            do j=jmin,jmax                ! Laplacian (in cases other
#   ifndef OBC_WEST 
              LapT(istr-1,j,k)=0.         ! than periodic: closed
#   else
              LapT(istr-1,j,k)=LapT(istr,j,k)      !  or gradient)
#   endif
            enddo
          enddo
        endif
        if (EASTERN_EDGE) then
          do k=1,N
            do j=jmin,jmax
#   ifndef OBC_EAST
              LapT(iend+1,j,k)=0.
#   else
              LapT(iend+1,j,k)=LapT(iend,j,k)
#   endif
            enddo
          enddo
        endif
#  endif /* !EW_PERIODIC */
#  ifndef NS_PERIODIC
        if (SOUTHERN_EDGE) then
          do k=1,N
            do i=imin,imax
#   ifndef OBC_SOUTH
              LapT(i,jstr-1,k)=0.
#   else
              LapT(i,jstr-1,k)=LapT(i,jstr,k)
#   endif
            enddo
          enddo
        endif
        if (NORTHERN_EDGE) then
          do k=1,N
            do i=imin,imax
#   ifdef OBC_NORTH
              LapT(i,jend+1,k)=0.
#   else
              LapT(i,jend+1,k)=LapT(i,jend,k)
#   endif
            enddo
          enddo
        endif
#  endif /* !NS_PERIODIC */

!                                          !!! WARNING: RECURSIVE
!  The SECOND rotated Laplacian operator:        BLOCKING SEQUENCE
!--------------------------------------------------------------------
!                   !  Vertical    dTdx,dTde(:,:,k2) k+1   rho-points
      k2=1          !  placement     dTdr,FC(:,:,k2) k+1/2   W-points
      do k=0,N,+1   !  of fields:  dTdx,dTde(:,:,k1) k     pho-points
       k1=k2        !                dTdr,FC(:,:,k1) k-1/2   W-points
       k2=3-k1
 
        if (k.lt.N) then
          do j=jstr,jend
            do i=istr,iend+1
              cff=0.5*(pm(i,j)+pm(i-1,j)) SWITCH umask(i,j)
              dTdx(i,j,k2)=cff*(LapT(i,j,k+1)-LapT(i-1,j,k+1))
            enddo
          enddo
          do j=jstr,jend+1
            do i=istr,iend
              cff=0.5*(pn(i,j)+pn(i,j-1)) SWITCH vmask(i,j)
              dTde(i,j,k2)=cff*(LapT(i,j,k+1)-LapT(i,j-1,k+1))
            enddo
          enddo
        endif
 
        if (k.eq.0 .or. k.eq.N) then
          do j=jstr-1,jend+1
            do i=istr-1,iend+1
              FC(i,j,k2)=0.0
              dTdz(i,j,k2)=0.0
              Akz(i,j,k )= 0.0               
            enddo
          enddo
        else
          do j=jstr-1,jend+1
            do i=istr-1,iend+1
              cff  =  1. / (z_r(i,j,k+1)-z_r(i,j,k))
              dTdz(i,j,k2)=cff*(LapT(i,j,k+1)-LapT(i,j,k)) SWITCH rmask(i,j)
            enddo
          enddo
        endif
 
        if (k.gt.0) then
          do j=jstr,jend
            do i=istr,iend+1
              FX(i,j)=-0.5*diff3u(i,j,k)*(Hz(i,j,k)+Hz(i-1,j,k))
     &         *dn_u(i,j)*(   dTdx(i  ,j,k1)+0.25*Sx_u(i,j,k)*(
     &                dTdz(i-1,j,k1)+dTdz(i-1,j,k2)     
     &         +      dTdz(i  ,j,k1)+dTdz(i  ,j,k2)       
     &                                              )   )
            enddo
          enddo
          do j=jstr,jend+1
            do i=istr,iend
              FE(i,j)=-0.5*diff3v(i,j,k)*(Hz(i,j,k)+Hz(i,j-1,k))
     &        *dm_v(i,j)*( dTde(i,j,k1)+0.25*Sy_v(i,j,k)*(
     &             dTdz(i,j-1,k1)+dTdz(i,j-1,k2)     
     &         +   dTdz(i,j  ,k1)+dTdz(i,j  ,k2)       
     &                                              )   )
            enddo
          enddo
          if (k.lt.N) then
            do j=jstr,jend
              do i=istr,iend
             
             cff1=Sx_u(i  ,j,k  )
             cff2=Sx_u(i+1,j,k+1)
             cff3=Sx_u(i  ,j,k+1)
             cff4=Sx_u(i+1,j,k  )
 
             cff=    
     &       diff3u(i  ,j,k  )*cff1*(cff1*dTdz(i,j,k2)+dTdx(i  ,j,k1))
     &      +diff3u(i+1,j,k+1)*cff2*(cff2*dTdz(i,j,k2)+dTdx(i+1,j,k2)) 
     &      +diff3u(i  ,j,k+1)*cff3*(cff3*dTdz(i,j,k2)+dTdx(i  ,j,k2))
     &      +diff3u(i+1,j,k  )*cff4*(cff4*dTdz(i,j,k2)+dTdx(i+1,j,k1))

#ifdef Z_CORR
             cff5 = 1. / (z_r(i,j,k+1)-z_r(i,j,k))**2

             Akz(i,j,k) = 0.25*( 
     &       diff3u(i  ,j,k  )*(cff1*cff1*cff5+pm(i,j)*pm(i,j))
     &   +   diff3u(i+1,j,k+1)*(cff2*cff2*cff5+pm(i,j)*pm(i,j))
     &   +   diff3u(i  ,j,k+1)*(cff3*cff3*cff5+pm(i,j)*pm(i,j))
     &   +   diff3u(i+1,j,k  )*(cff4*cff4*cff5+pm(i,j)*pm(i,j))
     &                         )                !+++ sqrt(muX)*(sx^2/dz^2 + 1/dx^2) 
#endif 
 
             cff1=Sy_v(i,j  ,k  )
             cff2=Sy_v(i,j+1,k+1) 
             cff3=Sy_v(i,j  ,k+1)
             cff4=Sy_v(i,j+1,k  )

             cff=cff+(
     &     diff3v(i,j  ,k  )*cff1*(cff1*dTdz(i,j,k2)+dTde(i,j  ,k1))
     &   + diff3v(i,j+1,k+1)*cff2*(cff2*dTdz(i,j,k2)+dTde(i,j+1,k2))
     &   + diff3v(i,j  ,k+1)*cff3*(cff3*dTdz(i,j,k2)+dTde(i,j  ,k2))
     &   + diff3v(i,j+1,k  )*cff4*(cff4*dTdz(i,j,k2)+dTde(i,j+1,k1))  
     &               )
 
             FC(i,j,k2)=-0.25*cff

#ifdef Z_CORR
          Akz(i,j,k) =  Akz(i,j,k) + 0.25*(
     &        diff3v(i  ,j,k  )*(cff1*cff1*cff5+pn(i,j)*pn(i,j))
     &     +  diff3v(i,j+1,k+1)*(cff2*cff2*cff5+pn(i,j)*pn(i,j))
     &     +  diff3v(i,j  ,k+1)*(cff3*cff3*cff5+pn(i,j)*pn(i,j))
     &     +  diff3v(i,j+1,k  )*(cff4*cff4*cff5+pn(i,j)*pn(i,j)) 
     &                                    )                !+++ sqrt(muY)*(sy^2/dz^2 + 1/dy^2)  
#endif
              enddo
            enddo
          endif
          do j=jstr,jend
            do i=istr,iend



#  ifdef DIAGNOSTICS_TS
              work(i,j,k) = t(i,j,k,nnew,itrc)            
#  endif /* DIAGNOSTICS */

#if defined DIAGNOSTICS_PV && !defined DIAGNOSTICS_TS
              if (k.eq.1) work(i,j,1) = t(i,j,1,nnew,itrc) 
              if (k.eq.N) work(i,j,2) = t(i,j,N,nnew,itrc) 
#  endif /* DIAGNOSTICS */

#ifdef Z_CORR
              t(i,j,k,n2,itrc)=Hz(i,j,k)*t(i,j,k,n3,itrc) 
     &         + dts*(
     &                   pm(i,j)*pn(i,j)*( FX(i+1,j)-FX(i,j)
     &                                    +FE(i,j+1)-FE(i,j))
     &                  +FC(i,j,k2)-FC(i,j,k1)    )
#else
              t(i,j,k,n2,itrc)=t(i,j,k,n3,itrc) + dts*(
     &                   pm(i,j)*pn(i,j)*( FX(i+1,j)-FX(i,j)
     &                                    +FE(i,j+1)-FE(i,j))
     &                  +FC(i,j,k2)-FC(i,j,k1)    )/Hz(i,j,k)
#endif
            enddo
          enddo
        endif  ! <-- k.gt.0
      enddo   ! --> k


#ifdef Z_CORR
!+++++
      do k=1,N-1
       do j=jstr,jend
        do i=istr,iend
          cff        = dts*Akz(i,j,k)**2
          Akz(i,j,k) = max(0.25*(-1.+16.*cff),0.)      
          Akz(i,j,k) = Akz(i,j,k)*(z_r(i,j,k+1)-z_r(i,j,k))/dts   !-->   [Akz] = m  / s
        enddo
       enddo 
      enddo

      do j=jstr,jend 
!+++++             
            do i=istr,iend             
              do k=1,N-1
                CD(i,k) = Akz(i,j,k)*(t(i,j,k+1,n1,itrc)-t(i,j,k  ,n1,itrc))
              enddo 
              CD(i,0) = 0.
              CD(i,N) = 0.         
            enddo             
            do i=istr,iend  
              FFC(i,1)  = dts*Akz(i,j,1)        
              cff       = 1./(Hz(i,j,1)+FFC(i,1))
              CF(i,1)   = cff*FFC(i,1)
              DC(i,1)   = cff*( t(i,j,1,n2,itrc) 
     &                        - dts*(CD(i,1)-CD(i,0)) )
            enddo            
            do k=2,N-1,+1
             do i=istr,iend 
              FFC(i,k) = dts*Akz(i,j,k)
              cff      = 1./  
     &         ( Hz(i,j,k) + FFC(i,k) + FFC(i,k-1)*(1.-CF(i,k-1)) ) 
              CF(i,k)  = cff*FFC(i,k)
              DC(i,k)  = cff*( t(i,j,k,n2,itrc) 
     &                    -    dts*(CD(i,k)-CD(i,k-1)) 
     &                    +    FFC(i,k-1)*DC(i,k-1) )
             enddo
            enddo           
            do i=istr,iend 
              t(i,j,N,n2,itrc) = ( t(i,j,N,n2,itrc) 
     &        -  dts*(CD(i,N)-CD(i,N-1))   
     &        + FFC(i,N-1) *DC(i,N-1) ) /   
     &        ( Hz(i,j,N)+FFC(i,N-1)*(1.-CF(i,N-1)) )
            enddo                       
            do k=N-1,1,-1
             do i=istr,iend              
              t(i,j,k,n2,itrc) = DC(i,k)+CF(i,k)*t(i,j,k+1,n2,itrc)
              
              
              if( abs(t(i,j,k,n2,itrc)) .gt. 50. ) then
                  print*,'FUCK    = ',i,j,k,itrc
                  print*,'RE-FUCK = ',h(i,j),istr,jstr
                  print*,'HBL     = ',hbls(i,j,3-nstp)
                  print*,'HBBL    = ',hbbls(i,j,3-nstp)
              endif
              
             enddo              
            enddo                           
!+++++      
       enddo
#endif



#if defined DIAGNOSTICS_PV && !defined DIAGNOSTICS_TS
       do j=jstr,jend
        do i=istr,iend
            pvDiab_sfc(i,j,itrc)=pvDiab_sfc(i,j,itrc)
     &                     + ( t(i,j,N,nnew,itrc)
     &                     -  work(i,j,2) )/dt

            pvDiab_bot(i,j,itrc)=pvDiab_bot(i,j,itrc)
     &                     + ( t(i,j,1,nnew,itrc)
     &                     -  work(i,j,1) )/dt
        enddo                           
       enddo    
#  endif /* DIAGNOSTICS */




#  ifdef DIAGNOSTICS_TS
!
! All diagnostic terms have already divided by cell volume 
! (Hz(i,j,k,itrc)/(pm(i,j).*pn(i,j)) in t3dmix_GP.F. The unit
! of diag terms is already: (unit of tracers)* s-1.
!

       do k=1,N
        do j=jstr,jend
         do i=istr,iend

              THmix(i,j,k,itrc)=THmix(i,j,k,itrc)
     &                     + ( t(i,j,k,nnew,itrc)
     &                     -  work(i,j,k) )/dt
#   ifdef MASKING
     &                                          * rmask(i,j)
#   endif

              Trate(i,j,k,itrc)=Trate(i,j,k,itrc)
     &                     + THmix(i,j,k,itrc)
#   ifdef MASKING
     &                                          * rmask(i,j)
#   endif


         enddo                           
        enddo     
       enddo

#  endif /* DIAGNOSTICS */


#if defined DIAGNOSTICS_TS && defined DIAGNOSTICS_PV
         do j=jstr,jend
           do i=istr,iend
             pvDiab_sfc(i,j,itrc)=pvDiab_sfc(i,j,itrc) + THmix(i,j,N,itrc)
             pvDiab_bot(i,j,itrc)=pvDiab_bot(i,j,itrc) + THmix(i,j,1,itrc)
           enddo
         enddo
# endif /* DIAGNOSTICS_PV */






# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
        call exchange_r3d_tile (istr,iend,jstr,jend,
     &                          t(START_2D_ARRAY,1,n2,itrc))
# endif

      return
      end
 
