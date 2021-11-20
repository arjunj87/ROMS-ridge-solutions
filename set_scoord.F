#include "cppdefs.h"
#ifdef SOLVE3D
                                     ! Setup S-coordinate system:
      subroutine set_scoord          ! input:  theta_s, theta_b, hc
      implicit none                  ! output: Cs_w(0:N), Cs_r(1:N) 
# include "param.h"
# include "scoord.h"
# include "scalars.h"
# include "comm_vars.h"
      integer k                      ! Compute vertical stretching
      real ds,sc, z1,zhc,z2,z3, CSF  ! curves at W- and RHO-points,
# ifndef NEW_S_COORD
      if (hc.le.hmin) then           !    -1 < Cs_r,Cs_w < 0, 
# endif
        ds=1.D0/dble(N)              ! then print z-coordinates 
        Cs_w(N)=0.D0                 ! of vertical S-levels over
        do k=N-1,1,-1                ! depths: minumum, maximum
          sc=ds*dble(k-N)            ! and half-way over slope.
          Cs_w(k)=CSF(sc, theta_s,theta_b)
        enddo
        Cs_w(0)=-1.D0
        do k=1,N
          sc=ds*(dble(k-N)-0.5D0)
          Cs_r(k)=CSF(sc, theta_s,theta_b)
        enddo

        MPI_master_only write(stdout,'(/1x,A/,/2x,A,7x,A/)')
     &        'Vertical S-coordinate system (z at W-points):',
     &             'level   S-coord    Cs-curve    Z at hmin',
     &                        'at hc    half way     at hmax'
        do k=N,0,-1
          sc=ds*dble(k-N)        
# ifdef NEW_S_COORD
          z1=hmin*(hc*sc + hmin*Cs_w(k))/(hc+hmin)
          zhc=0.5*hc*(sc + Cs_w(k))
          z2=0.5*hmax*(hc*sc + 0.5*hmax*Cs_w(k))/(hc+0.5*hmax)
          z3=hmax*(hc*sc + hmax*Cs_w(k))/(hc+hmax)
# else
          zhc=hc*sc
          z1=zhc +     (hmin-hc)*Cs_w(k)
          z2=zhc + (0.5*hmax-hc)*Cs_w(k)
          z3=zhc +     (hmax-hc)*Cs_w(k)
# endif
          MPI_master_only write(stdout,'(I7,F11.6,F12.7,4F12.3)')
     &                         k, ds*(k-N), Cs_w(k), z1,zhc,z2,z3
        enddo
# ifndef NEW_S_COORD 
      else
        write(stdout,'(/1x,2A,F7.2/12x,A,F7.2/)') '### ERROR: ',
     &    'Specified S-coordinate critical depth  hc   =', hc,
     &    'exceeds minimum unmasked topography.   hmin =', hmin
        may_day_flag=9
      endif
# endif
      return
      end


# ifdef NEW_S_COORD
                                            ! NOTE: Mathematical 
      function CSF (sc, theta_s,theta_b)    ! limits of CSF,csrf for
      implicit none                         ! theta_s, theta_b --> 0
      real*8 CSF, sc, theta_s,theta_b,csrf  ! match that under "else"
                                            ! logical branches.
      if (theta_s.gt.0.D0) then
        csrf=(1.D0-cosh(theta_s*sc))/(cosh(theta_s)-1.D0)
      else
        csrf=-sc**2
      endif
      if (theta_b.gt.0.D0) then
        CSF=(exp(theta_b*csrf)-1.D0)/(1.D0-exp(-theta_b))  
      else
        CSF=csrf
      endif
      return
      end

# else
      function CSF (sc, theta_s,theta_b)          ! Legacy
      real*8 CSF, sc, theta_s,theta_b             ! S-coordinate
                                                  ! transformation
      if (theta_s.gt.0.D0) then
        CSF=(1.D0-theta_b)*sinh(theta_s*sc)/sinh(theta_s)
     &           +theta_b*( 0.5D0*tanh(theta_s*(sc+0.5D0))
     &                       /tanh(0.5*theta_s)   -0.5D0 )
      else
        CSF=sc
      endif
      return
      end
# endif    /* NEW_S_COORD */



      subroutine check_scoord_switches (ierr)
      implicit none
      integer ierr, is,ie, lenstr
# include "param.h"
# include "strings.h"
      ie=lenstr(cpps)
      is=ie+2
      ie=is+10
      if (ie.gt.max_opt_size) goto 99
      cpps(is:ie)='<scoord.h>'
      MPI_master_only write(stdout,'(1x,A)') cpps(is:ie)
      is=ie+2
# ifdef NEW_S_COORD 
      ie=is+14
      if (ie.gt.max_opt_size) goto 99
      cpps(is:ie)='NEW_S_COORD'
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
      subroutine set_scoord_empty
      end
#endif /* SOLVE3D */
 
