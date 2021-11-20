#include "cppdefs.h"

#define Z_CORR
#undef NEW_BC
#define DEBUG_EPVIC

#define   STAB1
#undef    STAB2
 
      subroutine t3dmix_ISO (tile)
      implicit none
      integer tile, itrc
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
     &    )     
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
     & )
!
!--------------------------------------------------------------------
!  Compute isopycnic horizontal mixing terms for tracer equations.
!-------------------------------------------------------------------- 
!
      implicit none
#include "param.h"
      integer itrc, istr,iend,jstr,jend, i,j,k,k1,k2
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
     &     cff3,cff4,cff5,
     &     dTde(PRIVATE_2D_SCRATCH_ARRAY,2),   cff2
       real cof1,cof2,cof3,cof4 
       integer imin,imax,jmin,jmax
#include "grid.h"
#include "ocean3d.h"
#include "mixing.h"
#include "climat.h"
#include "scalars.h"
#include "eos_vars.h"
!
#include "compute_auxiliary_bounds.h"
!
#ifdef MASKING
# define SWITCH *
#else
# define SWITCH !
#endif
 
#ifdef ADV_ISONEUTRAL

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
              dTdx(i,j,k2)=cff*(t(i,j,k+1,nstp,itrc)-t(i-1,j,k+1,nstp,itrc))
            enddo
          enddo
          do j=jmin,jmax+1
            do i=imin,imax
              cff=0.5*(pn(i,j)+pn(i,j-1)) SWITCH vmask(i,j)
              dTde(i,j,k2)=cff*(t(i,j,k+1,nstp,itrc)-t(i,j-1,k+1,nstp,itrc))
            enddo
          enddo
        endif
 
        if (k.eq.0 .or. k.eq.N) then
          do j=jmin-1,jmax+1
            do i=imin-1,imax+1
              FC(i,j,k2)=0.0
#ifndef NEW_BC
              dTdz(i,j,k2) = 0.0
#endif
             enddo
          enddo
#ifdef NEW_BC
         if( k.eq. 0 ) then
           do j=jstr-1,jend+1
           do i=istr-1,iend+1 
               dTdz(i,j,k2)=(t(i,j,2,nstp,itrc)-t(i,j,1,nstp,itrc)
     &                           )/(z_r(i,j,2)-z_r(i,j,1))               
           enddo
           enddo
        endif
       if(k.eq. N ) then
               do j=jstr-1,jend+1
               do i=istr-1,iend+1 
                   dTdz(i,j,k2)=dTdz(i,j,k1)
               enddo
              enddo
        endif
#endif 
         else
          do j=jmin-1,jmax+1
            do i=imin-1,imax+1
              cff  =  1. / ( z_r(i,j,k+1)-z_r(i,j,k) )
              dTdz(i,j,k2)=cff*(t(i,j,k+1,nstp,itrc)-t(i,j,k  ,nstp,itrc) 
     &                          ) SWITCH rmask(i,j)
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
#   ifndef OBC_NORTH
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
#ifndef NEW_BC
              dTdz(i,j,k2)=0.0
#endif
              Akz (i,j,k )= 0.0   
            enddo
          enddo
#ifdef NEW_BC
          if( k.eq. 0 ) then
             do j=jstr-1,jend+1
             do i=istr-1,iend+1                                       
               dTdz(i,j,k2)=(LapT(i,j,2)-LapT(i,j,1)
     &                      )/(z_r(i,j,2)-z_r(i,j,1))               
             enddo
             enddo
          endif
          if(k.eq.N) then
                  do j=jstr-1,jend+1
                  do i=istr-1,iend+1 
                      dTdz(i,j,k2)=dTdz(i,j,k1)
                  enddo
                  enddo
          endif
#endif 
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
     &     diff3u(i  ,j,k  )*cff1*(cff1*dTdz(i,j,k2)+dTdx(i  ,j,k1))
     &   + diff3u(i+1,j,k+1)*cff2*(cff2*dTdz(i,j,k2)+dTdx(i+1,j,k2))  
     &   + diff3u(i  ,j,k+1)*cff3*(cff3*dTdz(i,j,k2)+dTdx(i  ,j,k2))
     &   + diff3u(i+1,j,k  )*cff4*(cff4*dTdz(i,j,k2)+dTdx(i+1,j,k1)) 

#ifdef Z_CORR
             cff5 = 1. / (z_r(i,j,k+1)-z_r(i,j,k))**2

             cof2 = 0.25*( 
     &       diff3u(i  ,j,k  )*(cff1*cff1*cff5+pm(i,j)*pm(i,j))
     &   +   diff3u(i+1,j,k+1)*(cff2*cff2*cff5+pm(i+1,j)*pm(i+1,j))
     &   +   diff3u(i  ,j,k+1)*(cff3*cff3*cff5+pm(i,j)*pm(i,j))
     &   +   diff3u(i+1,j,k  )*(cff4*cff4*cff5+pm(i+1,j)*pm(i+1,j))
     &                         )                !+++ sqrt(muX)*(sx^2/dz^2 + 1/dx^2) 
#ifdef STAB2     
             cof1 = 0.25*( 
     &       diff3u(i  ,j,k  )*cff1*cff1
     &   +   diff3u(i+1,j,k+1)*cff2*cff2
     &   +   diff3u(i  ,j,k+1)*cff3*cff3
     &   +   diff3u(i+1,j,k  )*cff4*cff4
     &                           )  
#endif     
#endif  
             cff1=Sy_v(i,j  ,k  )
             cff2=Sy_v(i,j+1,k+1) 
             cff3=Sy_v(i,j  ,k+1)
             cff4=Sy_v(i,j+1,k  )
 
             cff=cff+(
     &        diff3v(i,j  ,k  )*cff1*(cff1*dTdz(i,j,k2)+dTde(i,j  ,k1))
     &      + diff3v(i,j+1,k+1)*cff2*(cff2*dTdz(i,j,k2)+dTde(i,j+1,k2))
     &      + diff3v(i,j  ,k+1)*cff3*(cff3*dTdz(i,j,k2)+dTde(i,j  ,k2))
     &      + diff3v(i,j+1,k  )*cff4*(cff4*dTdz(i,j,k2)+dTde(i,j+1,k1))
     &               )
 
             FC(i,j,k2)=-0.25*cff

#ifdef Z_CORR
          cof2 = cof2 + 0.25*(
     &        diff3v(i  ,j,k  )*(cff1*cff1*cff5+pn(i,j)*pn(i,j))
     &     +  diff3v(i,j+1,k+1)*(cff2*cff2*cff5+pn(i,j+1)*pn(i,j+1))
     &     +  diff3v(i,j  ,k+1)*(cff3*cff3*cff5+pn(i,j)*pn(i,j))
     &     +  diff3v(i,j+1,k  )*(cff4*cff4*cff5+pn(i,j+1)*pn(i,j+1)) 
     &                                    )   
#ifdef STAB2
          cof1 = cof1 + 0.25*( 
     &       diff3v(i  ,j,k  )*cff1*cff1
     &   +   diff3v(i,j+1,k+1)*cff2*cff2
     &   +   diff3v(i,j  ,k+1)*cff3*cff3
     &   +   diff3v(i,j+1,k  )*cff4*cff4
     &                           )  
          Akz(i,j,k) = 16.*cof1*cof2 / 
     &               (z_r(i,j,k+1)-z_r(i,j,k))
#else
          Akz(i,j,k) = cof2
#endif
#endif
              enddo
            enddo
          endif
          do j=jstr,jend
            do i=istr,iend
#ifdef Z_CORR
              t(i,j,k,nnew,itrc)=Hz(i,j,k)*t(i,j,k,nnew,itrc) 
     &         + dt*(
     &                   pm(i,j)*pn(i,j)*( FX(i+1,j)-FX(i,j)
     &                                    +FE(i,j+1)-FE(i,j))
     &                  +FC(i,j,k2)-FC(i,j,k1)    )
#else
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc) + dt*(
     &                   pm(i,j)*pn(i,j)*( FX(i+1,j)-FX(i,j)
     &                                    +FE(i,j+1)-FE(i,j))
     &                  +FC(i,j,k2)-FC(i,j,k1)    )/Hz(i,j,k)
#endif
            enddo
          enddo
        endif  ! <-- k.gt.0
      enddo   ! --> k


#ifdef Z_CORR
#ifdef STAB1
!+++++
      do k=1,N-1
       do j=jstr,jend
        do i=istr,iend        
          cff        = dt*Akz(i,j,k)**2
          Akz(i,j,k) = max(-1.+16.*cff,0.)      
          Akz(i,j,k) = Akz(i,j,k)*(z_r(i,   j,k+1)-z_r(i,j,k)) / dt   !-->   [Akz] = m  / s
        enddo
       enddo 
      enddo
#endif

      do j=jstr,jend 
!+++++             
            do i=istr,iend
              do k=1,N-1
                CD(i,k) = Akz(i,j,k)*(
     &             t(i,j,k+1,nstp,itrc)
     &      -      t(i,j,k  ,nstp,itrc)
     &                                    )
              enddo 
              CD(i,0) = 0.
              CD(i,N) = 0.         
            enddo             
            do i=istr,iend  
              FFC(i,1)  = dt*Akz(i,j,1)        
              cff       = 1./(Hz(i,j,1)+FFC(i,1))
              CF(i,1)   = cff*FFC(i,1)
              DC(i,1)   = cff*( t(i,j,1,nnew,itrc) 
     &                        - dt*(CD(i,1)-CD(i,0)) )
            enddo            
            do k=2,N-1,+1
             do i=istr,iend 
              FFC(i,k) = dt*Akz(i,j,k)
              cff      = 1./  
     &         ( Hz(i,j,k) + FFC(i,k) + FFC(i,k-1)*(1.-CF(i,k-1)) ) 
              CF(i,k)  = cff*FFC(i,k)
              DC(i,k)  = cff*( t(i,j,k,nnew,itrc) 
     &                    -    dt*(CD(i,k)-CD(i,k-1)) 
     &                    +    FFC(i,k-1)*DC(i,k-1) )
             enddo
            enddo           
            do i=istr,iend 
              t(i,j,N,nnew,itrc) = ( t(i,j,N,nnew,itrc) 
     &        -  dt*(CD(i,N)-CD(i,N-1))   
     &        + FFC(i,N-1) *DC(i,N-1) ) /   
     &        ( Hz(i,j,N)+FFC(i,N-1)*(1.-CF(i,N-1)) )
            enddo                       
            do k=N-1,1,-1
             do i=istr,iend              
              t(i,j,k,nnew,itrc) = DC(i,k)+CF(i,k)*t(i,j,k+1,nnew,itrc)
#ifdef DEBUG_EPVIC              
              if( abs(t(i,j,k,nnew,itrc)*rmask(i,j)) .gt. 50.) then
                  print*,'FUCK    = ',i,j,k,itrc
                  print*,'RE-FUCK = ',h(i,j),istr,jstr
             endif
#endif              
             enddo              
            enddo                           
!+++++      
       enddo
#endif

#  endif /* T_DIF4 */
     


# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
        call exchange_r3d_tile (istr,iend,jstr,jend,
     &                          t(START_2D_ARRAY,1,nnew,itrc))
# endif

      return
      end
 
