#include "cppdefs.h"
#if defined MIX_GP_TS || defined MIX_EN_TS
# undef GEOPOTENTIAL_DIFF
# define SMAX 0.01
# define SWITCH_TRIADS
# define HMIN_CONDITION
# define HHMIN 750.
c++++# define TAMERE
!++++++++++++++++++++++++++++++++++++++++++++ 
      subroutine Slopes (tile)
      implicit none
      integer tile
# include "param.h"
# include "private_scratch.h"
# include "compute_tile_bounds.h"
      call Slopes_tile ( istr,iend,jstr,jend, 
     &                   A3d(1,1),A3d(1,2),A2d(1,1),
     &                   A2d(1,3),A2d(1,5),A2d(1,7),
     &                   A2d(1,9),A2d(1,10))
      return
      end
 
      subroutine Slopes_tile (istr,iend,jstr,jend)
!====================================================================
      implicit none
# include "param.h"
      integer istr,iend,jstr,jend, i,j,k,k1,k2
      integer kstart,kend      
      real    cff,cff1,cff2,cff3,cff4,smax
      real    cffW1,cffW2
      real    hbbl,hbl,hmin,topo,eps
      real    z_to_rho,s_to_rho,s_to_z
      integer imin,imax,jmin,jmax,omp_get_thread_num
      parameter( hmin=50.   )
      parameter( eps =1.e-16)
!==================================================================== 
# include "grid.h"
# include "ocean3d.h"
# include "scalars.h"
# include "mixing.h"
# include "compute_auxiliary_bounds.h"
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ifndef EW_PERIODIC
      if (WESTERN_EDGE) then     ! Restrict extended ranges one
        imin=istrU               ! point inward near the physical
      else                       ! boundary. Note that this version
        imin=istrU-1             ! of code is suitable for MPI
      endif                      ! configuration too.
      if (EASTERN_EDGE) then
        imax=iend
      else
        imax=iend+1
      endif
# else
      imin=istr-1
      imax=iend+1
# endif
# ifndef NS_PERIODIC
        if (SOUTHERN_EDGE) then
          jmin=jstrV
        else
          jmin=jstrV-1
        endif
        if (NORTHERN_EDGE) then
          jmax=jend
        else
          jmax=jend+1
        endif
# else
        jmin=jstr-1
        jmax=jend+1
# endif
!+++
!+++ Chain rule gives :
!+++
!+++  sigma_to_rho = sigma_to_z + z_to_rho
!+++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      do k=1,N
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++       
        do j=jstr,jend
          do i=imin,imax
            
             if(umask(i,j) .gt. 0.5) then          
             cff1    = 0.5*(pm(i,j)+pm(i-1,j))
             s_to_z  = cff1*(z_r(i,j,k)-z_r(i-1,j,k))
             topo    = 0.5*(h(i,j)+h(i-1,j))                          
#ifdef SWITCH_TRIADS
             s_to_rho = -0.5*( 
     &             max(dRdx(i,j,k),0.)*(dzdr(i-1,j,k-1)+dzdr(i,j,k  ))
     &            +min(dRdx(i,j,k),0.)*(dzdr(i-1,j,k  )+dzdr(i,j,k-1))
     &                   ) 
#else
             s_to_rho = -0.25*( dzdr(i  ,j,k)+dzdr(i  ,j,k-1) 
     &            +             dzdr(i-1,j,k)+dzdr(i-1,j,k-1) )*drdx(i,j,k)
#endif  
#ifdef TAMERE 
             s_to_rho = -0.5*drdx(i,j,k)*(
     &       2.*dzdr(i  ,j,k)*dzdr(i  ,j,k-1)
     &    /    (dzdr(i  ,j,k)+dzdr(i  ,j,k-1))
     &    +  2.*dzdr(i-1,j,k)*dzdr(i-1,j,k-1)
     &    /    (dzdr(i-1,j,k)+dzdr(i-1,j,k-1))      
     &       )
#endif
             z_to_rho = s_to_rho + s_to_z     !++                  
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++                                      
#ifdef HMIN_CONDITION
             cff2 = (topo-hmin)/(HHMIN-hmin)                         
#endif
             cff3 = ( z_r(i,j,k)+h(i,j)+z_r(i-1,j,k)+h(i-1,j) )
     &            / ( z_w(i,j,5)+h(i,j)+z_w(i-1,j,5)+h(i-1,j) ) 
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
             smax = SMAX             
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++             
#ifdef HMIN_CONDITION
             if( cff2 .le. 1. ) then
                smax=(3.*cff2*cff2-2.*cff2*cff2*cff2)*SMAX
             endif
#endif

c++             if( cff3 .le. 1. ) then
c++                smax=(3.*cff3*cff3-2.*cff3*cff3*cff3)*SMAX           
c++             endif             
             
             cff1        =  0.5*
     &   (1.+tanh((0.4*smax-abs(z_to_rho))/(0.1*smax+eps)))              
           
               Sx_u(i,j,k) = cff1*z_to_rho - s_to_z

#ifdef GEOPOTENTIAL_DIFF
               Sx_u(i,j,k) = -s_to_z
#endif
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++                            
             else
               Sx_u(i,j,k) = 0.
             endif
          enddo
        enddo
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ifndef EW_PERIODIC
        if (WESTERN_EDGE) then         ! Extrapolate elementary
          do j=jstr,jend               ! differences near physical
              Sx_u(imin-1,j,k)=Sx_u(imin,j,k)            
          enddo
        endif
        if (EASTERN_EDGE) then
          do j=jstr,jend
              Sx_u(imax+1,j,k)=Sx_u(imax,j,k)            
          enddo
        endif
# endif
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        do j=jmin,jmax
          do i=istr,iend
             if(vmask(i,j) .gt. 0.5) then

             cff1  = 0.5*(pn(i,j)+pn(i,j-1))
             s_to_z   = cff1*(z_r(i,j,k)-z_r(i,j-1,k))            
             topo  = 0.5*(h(i,j)+h(i,j-1))             
#ifdef SWITCH_TRIADS
             s_to_rho = -0.5*( 
     &             max(dRde(i,j,k),0.)*(dzdr(i,j-1,k-1)+dzdr(i,j,k  ))
     &            +min(dRde(i,j,k),0.)*(dzdr(i,j-1,k  )+dzdr(i,j,k-1))
     &                 ) 
#else
             s_to_rho = -0.25*( dzdr(i,j  ,k)+dzdr(i,j  ,k-1) 
     &            +        dzdr(i,j-1,k)+dzdr(i,j-1,k-1) )*drde(i,j,k)
#endif                              
# ifdef TAMERE
             s_to_rho = -0.5*drde(i,j,k)*(
     &       2.*dzdr(i,j  ,k)*dzdr(i,j  ,k-1)
     &    /    (dzdr(i,j  ,k)+dzdr(i,j  ,k-1))
     &    +  2.*dzdr(i,j-1,k)*dzdr(i,j-1,k-1)
     &    /    (dzdr(i,j-1,k)+dzdr(i,j-1,k-1))      
     &       )
#endif
             z_to_rho = s_to_rho + s_to_z     !++ 
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++               
#ifdef HMIN_CONDITION
             cff2 = (topo-hmin)/(HHMIN-hmin) 
#endif           
             cff3 = ( z_r(i,j,k)+h(i,j)+z_r(i,j-1,k)+h(i,j-1) )
     &            / ( z_w(i,j,5)+h(i,j)+z_w(i,j-1,5)+h(i,j-1) )            
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
             smax = SMAX             
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++             
#ifdef HMIN_CONDITION
             if( cff2 .le. 1. ) then
                smax=(3.*cff2*cff2-2.*cff2*cff2*cff2)*SMAX
             endif
#endif
c++             if( cff3 .le. 1. ) then
c++                smax=(3.*cff3*cff3-2.*cff3*cff3*cff3)*SMAX           
c++             endif 
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
             cff1        =  0.5*
     &   (1.+tanh((0.4*smax-abs(z_to_rho))/(0.1*smax+eps)))     
                    
               Sy_v(i,j,k) = cff1*z_to_rho - s_to_z
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#ifdef GEOPOTENTIAL_DIFF
               Sy_v(i,j,k) = -s_to_z
#endif
             else
               Sy_v(i,j,k) = 0.
             endif

          enddo
        enddo   
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ifndef NS_PERIODIC
        if (SOUTHERN_EDGE) then
          do i=istr,iend
              Sy_v(i,jmin-1,k)=Sy_v(i,jmin,k)            
            enddo
        endif
        if (NORTHERN_EDGE) then
          do i=istr,iend
              Sy_v(i,jmax+1,k)=Sy_v(i,jmax,k)              
          enddo
        endif
# endif
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      enddo  ! <-- k
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI      
      call exchange_u3d_tile (istr,iend,jstr,jend, Sx_u)
      call exchange_v3d_tile (istr,iend,jstr,jend, Sy_v)      
# endif
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return 
      end
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#else
      subroutine Slopes_empty
      end
#endif /* MIX_TS_EN */
