#include "cppdefs.h"
# if defined SEAMOUNT || defined UPWELLING || defined PACIFIC \
                      || defined ISWAKE
#  define REPORT_MAX_UV
#  define REPORT_POTEN_ENRG
# elif defined USWC
#  define REPORT_MAX_UV
# endif

#define SUM_BY_PAIRS
#undef  ETALON_CHECK
#define MAX_ADV_CFL
#define REPORT_MAX_UV
 
      subroutine diag (tile)
      implicit none
      integer tile,   i,j
#include "param.h"
#include "scalars.h"
#include "private_scratch.h"
#include "compute_tile_bounds.h"

      j=min(iic,2*ninfo)-1          ! This introduces "logarithmic"
      i=2                           ! scaling of diag call frequency
      do while(i.lt.j)              ! at the beginning of model run.
        i=2*i
      enddo
      if (mod(iic-1,min(i,ninfo)) .eq. 0) then
        call diag_tile (istr,iend,jstr,jend, A2d(1,1),A2d(1,2),A2d(1,3)
#ifdef SOLVE3D
     &               , A2d(1,4), A2d(1,5), A2d(1,6), A2d(1,7), A2d(1,8)
#endif
     &                                                                )
      endif
      return
      end
 
      subroutine diag_tile (istr,iend,jstr,jend, dVol, ke,pe
#ifdef SOLVE3D
     &                            , ke2b, ke3bc, kesrf, ub,vb
#endif
     &                                                      )
      implicit none
      integer istr,iend,jstr,jend
#include "param.h"
      real, dimension(PRIVATE_2D_SCRATCH_ARRAY) :: dVol, ke, pe
#ifdef SOLVE3D
     &                             , ke2b, ke3bc, kesrf, ub, vb
#endif
#ifdef MPI
# include "mpif.h"
      integer size, step, status(MPI_STATUS_SIZE)
      real*QUAD buff(16)
      common /xyz/ buff
#endif
#include "grid.h"
#ifdef SOLVE3D
# include "ocean3d.h"
# include "eos_vars.h"
#else
# include "ocean2d.h"
#endif
#include "scalars.h"
#include "comm_vars.h"
#include "diag_vars.h"
#include "ncvars.h"

      real cff, dA, v2, my_v2d_max
      real*QUAD my_avzeta, my_ke,   my_pe
#ifdef SOLVE3D
     &        , my_ke3bc,  my_ke2b, my_kesrf
#endif
      integer i,j,k, nsubs, ierr, is,ie

#ifdef SOLVE3D
      real v2bc
# ifdef MAX_ADV_CFL
      real my_Cu_Adv, my_Cu_W, ciV, cx, cw
      integer my_i_cmax, my_j_cmax, my_k_cmax
# else
      real my_v3d_max, my_v3bc_max
# endif
#endif
#ifdef SUM_BY_PAIRS
      integer isize,itg, js,jsize,jtg
#endif

      integer, parameter :: max_check_line=128
      character check_line*(max_check_line), tstring*18

! Compute and report various diagnostics: volume-averaged kinetic and
! potential energy, horizontally averaged free-surface perturbation;
! maximum velocity, etc. Since this operation involves computation of
! global sums, it is done in three stages: at first, summation within
! the tile [subdomain of indices (istr:iend,jstr:jend)] independently
! by individual threads. In the case of three dimensions also perform
! verical summation at this stage; then summation accross the threads
! of the same MPI process (if any), and, finally, MPI reduction to
! compute global integrals.


      ierr=0
      my_v2d_max=0.
#ifdef SOLVE3D
# ifdef MAX_ADV_CFL
      my_Cu_Adv=0.
      my_Cu_W=0.
      my_i_cmax=0
      my_j_cmax=0
      my_k_cmax=0
# else
      my_v3d_max=0.
      my_v3bc_max=0.
#endif
      do j=jstr,jend+1
        do i=istr,iend+1
          ub(i,j)=(Hz(i,j,N)+Hz(i-1,j,N))*u(i,j,N,nstp)
          vb(i,j)=(Hz(i,j,N)+Hz(i,j-1,N))*v(i,j,N,nstp)
        enddo
        do k=N-1,2,-1
          do i=istr,iend+1
            ub(i,j)=ub(i,j)+(Hz(i,j,k)+Hz(i-1,j,k))*u(i,j,k,nstp)
            vb(i,j)=vb(i,j)+(Hz(i,j,k)+Hz(i,j-1,k))*v(i,j,k,nstp)
          enddo
        enddo
        do i=istr,iend+1
          ub(i,j)=(ub(i,j)+(Hz(i,j,1)+Hz(i-1,j,1))*u(i,j,1,nstp))
     &         /(z_w(i,j,N)+z_w(i-1,j,N)-z_w(i,j,0)-z_w(i-1,j,0))
          vb(i,j)=(vb(i,j)+(Hz(i,j,1)+Hz(i,j-1,1))*v(i,j,1,nstp))
     &         /(z_w(i,j,N)+z_w(i,j-1,N)-z_w(i,j,0)-z_w(i,j-1,0))
        enddo
      enddo
 
      cff=g/rho0
      do j=jstr,jend
        do i=istr,iend
          v2=0.5*(ub(i,j)**2+ub(i+1,j)**2 +vb(i,j)**2+vb(i,j+1)**2)
          my_v2d_max=max(my_v2d_max, v2)

          ke(i,j)=0.
          pe(i,j)=0.5*g*z_w(i,j,N)*z_w(i,j,N)

          ke2b(i,j)=0.5*(z_w(i,j,N)-z_w(i,j,0))*v2
          ke3bc(i,j)=0.
          kesrf(i,j)=0.25*( u(i,j,N,nstp)**2 + u(i+1,j,N,nstp)**2
     &                     +v(i,j,N,nstp)**2 + v(i,j+1,N,nstp)**2)
        enddo
        do k=N,1,-1
          do i=istr,iend
            v2=0.5*( u(i,j,k,nstp)**2 + u(i+1,j,k,nstp)**2
     &              +v(i,j,k,nstp)**2 + v(i,j+1,k,nstp)**2)
 
            v2bc=0.5*( (u(i  ,j,k,nstp)-ub(i  ,j))**2
     &                +(u(i+1,j,k,nstp)-ub(i+1,j))**2
     &                +(v(i,j  ,k,nstp)-vb(i,j  ))**2
     &                +(v(i,j+1,k,nstp)-vb(i,j+1))**2)
 
! Find maximum Advective Courant number: in the code segment below it
! is defined in the finite-volume sense as the sum of fluxes directed
! OUTWARD from grid box Hz divided by its control volume. Essentially
! it is the fraction of water replaced within the grid box during one
! time step.   Note that under this definition, Cu=1 is the condition
! when flux-split first-order upstream advection scheme looses its
! positive definiteness property (hence stability).  Reported values
! are: "Cu_Adv" is full (tri-dimensional) Courant number; "i,j,k" are 
! coordinates where its maximum occurs;  "Cu_W" is contribution into
! "Cu_Adv" due to vertical velocity.

# ifdef MAX_ADV_CFL
#  ifdef MASKING
            ciV=dt*rmask(i,j)*pm(i,j)*pn(i,j)/Hz(i,j,k)
#  else
            ciV=dt*pm(i,j)*pn(i,j)/Hz(i,j,k)
#  endif
            cw=ciV*( max(W(i,j,k), 0.) - min(W(i,j,k-1), 0.))

            cx=cw+ciV*( max(FlxU(i+1,j,k), 0.) -min(FlxU(i,j,k), 0.)
     &                   +max(FlxV(i,j+1,k), 0.) -min(FlxV(i,j,k), 0.))

            if (cx .gt. my_Cu_Adv) then
              my_Cu_Adv=cx
              my_Cu_W=cw
              my_i_cmax=i
              my_j_cmax=j
              my_k_cmax=k
            endif
# else
            my_v3d_max  = max(my_v3d_max,  v2)
            my_v3bc_max = max(my_v3bc_max, v2bc)
# endif
 
            ke(i,j)=ke(i,j) + 0.5*v2*Hz(i,j,k)
 
            pe(i,j)=pe(i,j) + cff*Hz(i,j,k)
# ifdef SPLIT_EOS
     &      *(rho1(i,j,k)+qp1(i,j,k)*(z_w(i,j,N)-z_r(i,j,k)))
# else
     &                                            *rho(i,j,k)
# endif
     &                               *(z_r(i,j,k)-z_w(i,j,0))
 
            ke3bc(i,j)=ke3bc(i,j) + 0.5*v2bc*Hz(i,j,k)
          enddo
        enddo
        do i=istr,iend
# ifdef MASKING
          dA=rmask(i,j)/(pm(i,j)*pn(i,j))
# else
          dA=1./(pm(i,j)*pn(i,j))
# endif
          dVol(i,j) = dA * z_w(i,j,N)
          ke(i,j)   = dA * ke(i,j)
          pe(i,j)   = dA * pe(i,j)
          ke2b(i,j) = dA * ke2b(i,j)
          ke3bc(i,j)= dA * ke3bc(i,j)
          kesrf(i,j)= dA * kesrf(i,j)
        enddo
      enddo
#else
      cff=0.5*g
      do j=jstr,jend
        do i=istr,iend
          v2=0.5*( ubar(i,j,kstp)**2+ubar(i+1,j,kstp)**2
     &            +vbar(i,j,kstp)**2+vbar(i,j+1,kstp)**2)
          my_v2d_max=max(my_v2d_max, v2)
# ifdef MASKING
          dA=rmask(i,j)/(pm(i,j)*pn(i,j))
# else
          dA=1./(pm(i,j)*pn(i,j))
# endif
          dVol(i,j)=dA*zeta(i,j,kstp)
          ke(i,j)=dA*0.5*(zeta(i,j,kstp)+h(i,j))*v2
          pe(i,j)=cff*dVol(i,j)*zeta(i,j,kstp)
        enddo
      enddo
#endif /* SOLVE3D */
                                               ! Horizontal summation
                                               ! within the subdomain:
#ifdef SUM_BY_PAIRS
      isize=iend-istr                  ! Horizontal summation within
      jsize=jend-jstr                  ! the subdomain: Compute partial
                                       ! sums within the tile using a
      do while (isize>0 .or. jsize>0)  ! split-directional reduction by
        if (jsize>0) then              ! pairs algorithm in order to
          js=(jsize+1)/2-1             ! avoid accumulation of roundoff
          do j=0,js                    ! errors.
            jtg=jstr+j
            do i=istr,istr+isize
               dVol(i,jtg) =  dVol(i,jtg+j) +  dVol(i,jtg+j+1)
                 ke(i,jtg) =    ke(i,jtg+j) +    ke(i,jtg+j+1)
                 pe(i,jtg) =    pe(i,jtg+j) +    pe(i,jtg+j+1)
# ifdef SOLVE3D
               ke2b(i,jtg) =  ke2b(i,jtg+j) +  ke2b(i,jtg+j+1)
              ke3bc(i,jtg) = ke3bc(i,jtg+j) + ke3bc(i,jtg+j+1)
              kesrf(i,jtg) = kesrf(i,jtg+j) + kesrf(i,jtg+j+1)
# endif
            enddo
          enddo
          if (2*js+1 < jsize) then     ! carry over points which do not
            js=js+1                    ! have pair to be added together
            jtg=jstr+js
            do i=istr,istr+isize
               dVol(i,jtg) =  dVol(i,jtg+js)
                 ke(i,jtg) =    ke(i,jtg+js)
                 pe(i,jtg) =    pe(i,jtg+js)
# ifdef SOLVE3D
               ke2b(i,jtg) =  ke2b(i,jtg+js)
              ke3bc(i,jtg) = ke3bc(i,jtg+js)
              kesrf(i,jtg) = kesrf(i,jtg+js)
# endif
            enddo
          endif
          jsize=js
        endif

        if (isize>0) then
          is=(isize+1)/2-1
          do j=jstr,jstr+jsize
            do i=0,is
              itg=istr+i
               dVol(itg,j) =  dVol(itg+i,j) +  dVol(itg+i+1,j)
                 ke(itg,j) =    ke(itg+i,j) +    ke(itg+i+1,j)
                 pe(itg,j) =    pe(itg+i,j) +    pe(itg+i+1,j)
# ifdef SOLVE3D
               ke2b(itg,j) =  ke2b(itg+i,j) +  ke2b(itg+i+1,j)
              ke3bc(itg,j) = ke3bc(itg+i,j) + ke3bc(itg+i+1,j)
              kesrf(itg,j) = kesrf(itg+i,j) + kesrf(itg+i+1,j)
# endif
            enddo
          enddo
          if (2*is+1 < isize) then
            is=is+1
            itg=istr+is
            do j=jstr,jstr+jsize
               dVol(itg,j) =  dVol(itg+is,j)
                 ke(itg,j) =    ke(itg+is,j)
                 pe(itg,j) =    pe(itg+is,j)
# ifdef SOLVE3D
               ke2b(itg,j) =  ke2b(itg+is,j)
              ke3bc(itg,j) = ke3bc(itg+is,j)
              kesrf(itg,j) = kesrf(itg+is,j)
# endif
            enddo
          endif
          isize=is
        endif
      enddo  !<-- while()
#else
! Alternative: a two-stage (one for each direction) straightforward
! summation algorithm as an effort to mitigate roundoff errors.

      do j=jstr+1,jend
        do i=istr,iend
          dVol(i,jstr)  = dVol(i,jstr)  + dVol(i,j)
          ke(i,jstr)    = ke(i,jstr)    + ke(i,j)
          pe(i,jstr)    = pe(i,jstr)    + pe(i,j)
# ifdef SOLVE3D
          ke2b(i,jstr)  = ke2b(i,jstr)  + ke2b(i,j)
          ke3bc(i,jstr) = ke3bc(i,jstr) + ke3bc(i,j)
          kesrf(i,jstr) = kesrf(i,jstr) + kesrf(i,j)
# endif
        enddo
      enddo
      do i=istr+1,iend
        dVol(istr,jstr)  = dVol(istr,jstr)  + dVol(i,jstr)
        ke(istr,jstr)    = ke(istr,jstr)    + ke(i,jstr)
        pe(istr,jstr)    = pe(istr,jstr)    + pe(i,jstr)
# ifdef SOLVE3D
        ke2b(istr,jstr)  = ke2b(istr,jstr)  + ke2b(i,jstr)
        ke3bc(istr,jstr) = ke3bc(istr,jstr) + ke3bc(i,jstr)
        kesrf(istr,jstr) = kesrf(istr,jstr) + kesrf(i,jstr)
# endif
      enddo
#endif /* SUM_BY_PAIRS */
 
      my_avzeta=dVol(istr,jstr)
      my_ke=ke(istr,jstr)
      my_pe=pe(istr,jstr)
#ifdef SOLVE3D
      my_ke2b=ke2b(istr,jstr)
      my_ke3bc=ke3bc(istr,jstr)
      my_kesrf=kesrf(istr,jstr)
#endif
      if (SINGLE_TILE_MODE) then
        nsubs=1
      else
        nsubs=NSUB_X*NSUB_E
      endif

! Perform global summation: whoever gets first to the critical region
! resets global sums before global summation starts; after the global
! summation is completed, thread, which is the last one to enter the
! critical region, finalizes the computation of diagnostics and
! prints them out.

C$OMP CRITICAL (diag_cr_rgn)
      if (tile_count.eq.0) then                ! Initialize global
        avzeta=my_avzeta                       ! sums for summation
        avke=my_ke                             ! among the threads
        avpe=my_pe                             ! within shared-memory 
        v2d_max=my_v2d_max
#ifdef SOLVE3D
        avke2b=my_ke2b
        avke3bc=my_ke3bc
        avkesrf=my_kesrf
# ifdef MAX_ADV_CFL
        Cu_Adv=my_Cu_Adv
        Cu_W=my_Cu_W
        i_cx_max=my_i_cmax
        j_cx_max=my_j_cmax
        k_cx_max=my_k_cmax
# else
        v3d_max=my_v3d_max
        v3bc_max=my_v3bc_max
# endif
#endif
      else                                     ! Perform global
        avzeta=avzeta+my_avzeta                ! summation among
        avke=avke+my_ke                        ! the threads within
        avpe=avpe+my_pe                        ! MPI process.
        v2d_max=max(v2d_max, my_v2d_max)
#ifdef SOLVE3D
        avke2b=avke2b+my_ke2b
        avke3bc=avke3bc+my_ke3bc
        avkesrf=avkesrf+my_kesrf
# ifdef MAX_ADV_CFL
        if (my_Cu_Adv.gt.Cu_Adv) then
          Cu_Adv=my_Cu_Adv
          Cu_W=my_Cu_W
          i_cx_max=my_i_cmax
          j_cx_max=my_j_cmax
          k_cx_max=my_k_cmax
        endif
# else
        v3d_max=max(v3d_max, my_v3d_max)
        v3bc_max=max(v3bc_max,my_v3bc_max)
# endif
#endif
      endif
                                               ! This counter
      tile_count=tile_count+1                  ! identifies the last
                                               ! thread, whoever it
      if (tile_count.eq.nsubs) then            ! happens to be, not
        tile_count=0                           ! necessarily master.
#ifdef MPI
# ifdef SOLVE3D
#   ifdef MAX_ADV_CFL
        i_cx_max=i_cx_max + iSW_corn
        j_cx_max=j_cx_max + jSW_corn
#  endif
# endif
        if (NNODES.gt.1) then                  ! Perform global
          size=NNODES                          ! summation among
          do while (size.gt.1)                 ! MPI processes
            step=(size+1)/2
            if (mynode.ge.step .and. mynode.lt.size) then
              buff(1)=avzeta
              buff(2)=avke                     ! This is essentially
              buff(3)=avpe                     ! MPI_Reduce operation.
              buff(4)=avke2b                   ! Note that because
              buff(5)=avke3bc                  ! buff() may be real*16
              buff(6)=avkesrf                  ! the size of message
              buff(7)=v2d_max                  ! in MPI_Send/Recv is 
# ifdef SOLVE3D
#   ifdef MAX_ADV_CFL
              buff(8)=Cu_Adv                   ! doubled.
              buff(9)=Cu_W
              buff(10)=i_cx_max
              buff(11)=j_cx_max
              buff(12)=k_cx_max
#   else
              buff(8)=v3d_max
              buff(9)=v3bc_max
#  endif
# endif
 
!>            write(*,*) 'sending ', mynode, '-->', mynode-step
 
              call MPI_Send (buff,  24, MPI_DOUBLE_PRECISION,
     &               mynode-step, 17, ocean_grid_comm,      ierr)
            elseif (mynode .lt. size-step) then
              call MPI_Recv (buff,  24, MPI_DOUBLE_PRECISION,
     &            mynode+step, 17, ocean_grid_comm, status, ierr)
 
!<            write(*,*) 'received ',  mynode, '<--', mynode+step
 
              avzeta=avzeta+buff(1)
              avke=avke+buff(2)                ! here "v2" is used as 
              avpe=avpe+buff(3)                ! a scratch variable to 
              avke2b=avke2b+buff(4)            ! avoid type conflict 
              avke3bc=avke3bc+buff(5)          ! in the case if "buff"
              avkesrf=avkesrf+buff(6)          ! is 16-Byte, while
              v2=buff(7)                       ! "v2", "v2d_max", etc..
              v2d_max=max(v2d_max, v2)         ! are only 8-byte long.

# ifdef SOLVE3D
#   ifdef MAX_ADV_CFL

              if (buff(8).gt.Cu_Adv) then
                Cu_Adv=buff(8)
                Cu_W=buff(9)
                i_cx_max=buff(10)
                j_cx_max=buff(11)
                k_cx_max=buff(12)
              endif
#   else
              v2=buff(8)
              v3d_max=max(v3d_max, v2)
              v2=buff(9)
              v3bc_max=max(v3bc_max, v2)
#  endif
# endif
            endif
            size=step
          enddo
        endif
        if (mynode.eq.0) then
#endif
          avke=avke/(volume+avzeta)            ! Compute and print
          avpe=avpe/(volume+avzeta)            ! global diagnostics
          avke2b=avke2b/(volume+avzeta)        ! (the last thread of 
          avke3bc=avke3bc/(volume+avzeta)      ! master MPI process
          avkesrf=avkesrf/area                 ! only).
          v2d_max=sqrt(v2d_max)

# ifdef SOLVE3D
#   ifdef MAX_ADV_CFL

#   else
          v3d_max=sqrt(v3d_max)
          v3bc_max=sqrt(v3bc_max)
#  endif
# endif


#ifdef SOLVE3D
# ifdef MAX_ADV_CFL
          if (first_time.eq.0) then
            first_time=1
            write(*,1)   'STEP', 'time[DAYS]',  'KINETIC_ENRG',
     &                'BAROTR_KE', 'MAX_ADV_CFL', 'MAX_VERT_CFL',
     &                                     'i_cx', 'j_cx', 'k_c'
C$   &                                                   , 'trd'
          endif
          write(check_line,2)  avke,  avke2b,  Cu_Adv,  Cu_W,
     &                            i_cx_max,  j_cx_max,  k_cx_max
C$   &                                                 , proc(2)

  1       format(1x,A,2x,A,1x,A,5x,A,8x,A,2(5x,A),3x,A,2(1x,A))
  2       format(ES18.11, ES17.10, 2F16.12, 2I7,I5, I3)

# elif defined REPORT_MAX_UV && defined REPORT_POTEN_ENRG

          if (first_time.eq.0) then
            first_time=1
            write(*,1) 'STEP', 'time[DAYS]', 'KINETIC_ENRG',
     &                   'BAROTR_KE', 'MAX_U3D', 'MAX_UBAR',
     &                   'MAX_UBCL', 'POTEN_ENRG'
C$   &                                                 , 'trd'
          endif
          write(check_line,2)  avke,    avke2b,   v3d_max,
     &                           v2d_max, v3bc_max, avpe
C$   &                                            , proc(2)
  1       format(1x,A,2x,A,1x,A,5x,A,6x,A,7x,A,6x,A,7x,A,3x,A)
  2       format(ES18.11, ES16.9, 3(ES15.8), ES16.8,I3)

# elif defined REPORT_MAX_UV
          if (first_time.eq.0) then
            first_time=1
            write(*,1) 'STEP', 'time[DAYS]', 'KINETIC_ENRG',
     &                 'BAROTR_ENRG', 'MAX_U3D', 'MAX_UBAR',
     &                 'MAX_UBCL', 'FREE_SURFACE'
C$   &                                                 , 'trd'
          endif
          write(check_line,2)  avke, avke2b, v3d_max, v2d_max,
     &                                    v3bc_max, avzeta/area
C$   &                                                , proc(2)
  1       format(1x,A,2x,A,1x,A,5x,A,4x,A,7x,A,6x,A,7x,A,1x,A)
  2       format(ES18.11,ES16.9,3(ES15.8),ES16.8,I3)
# else
          if (first_time.eq.0) then
            first_time=1
            write(*,'(1x,A,2x,A,1x,A,5x,A,4x,A,4x,A,1x,A)')
     &              'STEP',       'time[DAYS]', 'KINETIC_ENRG',
     &              'SURFACE_KE', 'POTEN_ENRG', 'FREE_SURFACE'
C$   &                                                  ,'trd'
          endif
          write(check_line,'(ES18.11,3(ES15.8),I3)') avke, avkesrf,
     &                                               avpe, avzeta/area
C$   &                                                   , proc(2)
# endif
#else
          if (first_time.eq.0) then
            first_time=1
            write(*,'(1x,A,2x,A,1x,A,5x,A,4x,A,9x,A,1x,A)')
     &                   'STEP',    'time[DAYS]', 'KINETIC_ENRG',
     &                   'POTEN_ENRG', 'KE+PE',   'FREE_SURFACE'
C$   &                                                 , 'trd'
          endif
          write(check_line,'(ES18.11,3(ES15.8),I3)')
     &                   avke, avpe, avke+avpe, avzeta/area
C$   &                                                 , proc(2)
#endif
          ie=max_check_line
          do while (check_line(ie:ie).eq.' ' .and. ie.gt.0)
            ie=ie-1
          enddo
          is=1
          do while (check_line(is:is).eq.' ' .and. is.lt.ie)
            is=is+1
          enddo                            ! Suppress floating-point
          i=is-1                           ! "E"s in order to shorten
          do while (i.lt.ie)               ! printed diagnostic line
            i=i+1
            if (check_line(i:i).eq.'E' .or.check_line(i:i).eq.'e') then
              check_line(i:ie-1)=check_line(i+1:ie)
              check_line(ie:ie)=' '
              ie=ie-1

! Appearance of any non-digit character other than blank,"+", "-", or
! "." indicates floating point exception condition (usually caused by
! blowup).    Force writing history file at every time step and set
! may_day_flag to signal termination.  The code behavior in this
! situation is to write the last state of all variables into history
! file and quit. 

            elseif (ichar(check_line(i:i)).lt.48 .or.
     &              ichar(check_line(i:i)).gt.57) then
              if (check_line(i:i).ne.' ' .and.
     &            check_line(i:i).ne.'+'  .and.
     &            check_line(i:i).ne.'-'  .and.
     &            check_line(i:i).ne.'.') then
                nwrt=1
                if (may_day_flag.eq.0) may_day_flag=1  !--> BLOW UP
              endif
            endif
          enddo  ! <-- while(..)

          write(tstring,'(F18.8)') tdays
          i=1
          do while (i.lt.18 .and. tstring(i:i).eq.' ')
            i=i+1
          enddo
          write(*,'(I7,2(1x,A))') iic-1, tstring(i:i+8),
     &                                  check_line(is:ie)
          call flush(6)

#ifdef  ETALON_CHECK
# define WRITE_ETALON

! The following segment performs comparison of model output from
! current run (stored and formatted as "check_line") with pre-stored
! data ("etalon_line"). This version is designed to handle situation
! where the number of significant digits for each parameter to be
! compared is unknown (determined dynamically) and it is only assumed
! that format contains mantissa and power, which are recognized as
! A1,P1 and A0,P0 for current and etalon data respectively. This is
! necessary to rescale the mantissa of the difference (if occurs) in
! such a way that its power is the same as for the current output.
! This code organization also has the property that no adjustment is
! necessary, if output format (see above) is modified (i.e., number
! of controlled parameters or number of significant digits in each
! parameter is changed).

C$        ie=ie-3             ! <--  disregard thread number
          do icheck=1,ncheck
            if (iic-1 .eq. check_point(icheck)) then
# ifdef WRITE_ETALON
              if (mcheck.eq.0) open(unit=19, file='etalon_data.NEW')
              mcheck=mcheck+1
              write(19,'(6x,A,I2,A,I6/6x,A,I2,A/5x,A,1x,3A)')
     &             'check_point(', mcheck, ')=', check_point(icheck),
     &             'etalon_line(', mcheck, ')=', '&',  '''',
     &                                        check_line(is:ie), ''''
              if (icheck.lt.ncheck) then
                call flush(19)
              else
                close(19)
                write(*,*) 'finished etalon_data.NEW.'
              endif
# endif

! Extract arrays of mantissa and exponent A1,P1 and A0,P0 from check
! and etalon lines respectively. Perform the comparison, choosing the
! number of controlled parameters to be the smaller of one from check
! and etalon line. Rescale "A0" of etalon line to match scale of
! mantissa of check line.

              call read_check_line (check_line,ie, A1,P1,j, npars)
              call read_check_line (etalon_line(icheck),
     &                                           ie, A0,P0,k, npars)
              k=min(j,k)
              if (k.gt.0) then
                ierr=0
                do i=1,k
                  if (P0(i).eq.P1(i) .and. A0(i).eq.A1(i)) then
                    A1(i)=0.
                  else
                    A1(i)=A1(i)-A0(i)*10.**(P0(i)-P1(i))
                    ierr=ierr+1
                  endif
                enddo            !--> discard A0,P0,P1, keep A1

! Print out comparison results. If difference is discovered, print
! out only its mantissa, retaining the same number of digits after
! decimal point as in the original format. To do so, analyse the
! format of check line and save the location of decimal point as
! P0(j) and last digit of mantissa P1(j) for each controlled
! parameter indexed j=1,..,k.

                if (ierr.eq.0) then
                  write(*,*) 'PASSED_ETALON_CHECK'
                else
                  do j=1,k             ! P0(j) is location of decimal
                    P0(j)=0            ! point of parameter "j";
                    P1(j)=0            ! P1(j) is location of the last
                  enddo                ! digit of its mantissa
                  j=0
                  do i=1,ie
                    if (check_line(i:i).eq.'.'.and.j.lt.npars) then
                      j=j+1
                      P0(j)=i
                    elseif (check_line(i:i).eq.' '.and.j.gt.0) then
                      if (P0(j).gt.0 .and. P1(j).eq.0) P1(j)=i-1
                    endif
                  enddo                ! unless j.eq.0, "i" is
                                       ! location the first symbol
                  do j=1,k             ! after the last digit of
                    if (j.gt.1) then   ! mantissa of the previous
                      i=P1(j-1)+1      ! parameter. Note that the
                    else               ! whole string is shifted by
                      i=-5             ! 6 characters when printing.
                    endif
                    write(fmt,'(2(A,I2),A)') '(F', P1(j)-i+1, '.',
     &                                             P1(j)-P0(j), ')'
                    write(check_line(i+6:P1(j)+6),fmt) A1(j)
                  enddo
                                       ! For visial purposes replace 0s
                  j=0                  ! after decimal point with dots. 
                  do i=2,P1(k)+6
                    if (check_line(i-1:i).eq.'0.') then
                      check_line(i:i)=':'
                      j=1
                    elseif (check_line(i:i).eq.'0' .and.
     &                                     j.eq.1) then
                      check_line(i:i)='.'
                    else
                      j=0
                    endif
                  enddo
                  write(*,'(1x,2A/1x,A)')    'difference',
     &            check_line(1:P1(k)+6), 'DIFFERENCE DISCOVERED'
                endif      !<-- ierr.eq.0
              else
                write(*,*) 'WARNING: No parameters can ',
     &                     'be compared by ETALON CHECK.'
              endif     !<-- k.gt.0
            endif     !<-- iic-1 .eq. check_point(icheck)
          enddo     !<-- icheck
#endif /* ETALON_CHECK */


c          if (prev_ke .gt. 1.D0-20) then
c            if (avke .gt. 1.002*prev_ke) then
c              nwrt=1
c              write(*,*) 'Abrupt change in KE: setting nwrt to 1.'
c            endif
c          endif
c          prev_ke=avke


#ifdef MPI
        endif    ! <-- mynode.eq.0

        buff(1)=may_day_flag
        buff(2)=nwrt

        call MPI_Bcast(buff, 4, MPI_DOUBLE_PRECISION,
     &                      0, ocean_grid_comm, ierr)
        may_day_flag=buff(1)
        nwrt=buff(2)
#endif
      endif  ! <-- tile_count.eq.nsubs
C$OMP END CRITICAL (diag_cr_rgn)
      return
      end

#ifdef ETALON_CHECK
      subroutine read_check_line (string,lstr, A,P,np, npars)

! For a given character string of length lstr, which is assumed to
! contain a sequence of floating point numbers written in Fw.d, Ew.d
! or ESw.d format unknown in advance, determine how many numbers "np"
! are written and extract mantissa A and exponent P separately for
! each number. If exponent is absent (Fw.d format), the corresponding
! P is set to zero.

      implicit none
      character(len=*) string
      integer  lstr, npars, P(npars), np,  i,k, ierr
      real A(npars)

      np=0
      i=0
      do while (i.lt.lstr .and. np.lt.npars)
        i=i+1
        if (string(i:i).eq.'.') then
          k=i-1
          do while (k.gt.1 .and. ichar(string(k:k)).ge.48
     &                     .and. ichar(string(k:k)).le.57)
            k=k-1
          enddo
          if (string(k:k).eq.' ') k=k+1
          i=i+1
          do while (ichar(string(i:i)).ge.48
     &        .and. ichar(string(i:i)).le.57)
            i=i+1
          enddo
          np=np+1
          A(np)=0.
          P(np)=0
          read(string(k:i-1),*, iostat=ierr) A(np)
          if (string(i:i).eq.'E' .or.
     &        string(i:i).eq.'e') then
            string(i:i)=' '
            i=i+1
          endif
          if (string(i:i).eq.'+' .or. string(i:i).eq.'-') then
            k=i
            i=i+1
            do while (ichar(string(i:i)).ge.48
     &          .and. ichar(string(i:i)).le.57)
              i=i+1
            enddo
            if (string(i:i).ne.'.') then
              read(string(k:i-1),*,iostat=ierr) P(np)
              string(k:i-1)=' '
            endif
          endif
c**       write(*,*) np, A(np),P(np)
        endif
      enddo
      return
      end
#endif
