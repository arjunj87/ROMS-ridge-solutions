#include "cppdefs.h"
#ifdef SOLVE3D
 
      subroutine u3dbc_tile (istr,iend,jstr,jend, grad)
!
! Set lateral boundary conditions for XI-component velocity
! u(:,:,:,nnew).
!
      implicit none
      integer istr,iend,jstr,jend, i,j,k
      real grad(PRIVATE_2D_SCRATCH_ARRAY), eps, dtfwd,
     &     cff, cx,cy, dft,dfx,dfy, tau,tau_in,tau_out
      parameter (eps=1.E-20)
# include "param.h"
# include "scalars.h"
# include "grid.h"
# include "ocean3d.h"
# ifdef M3_FRC_BRY
#  include "boundary.h"
# else
#  include "climat.h"
# endif

!
# include "compute_auxiliary_bounds.h"
!
      if (PRE_STEP) then      ! because predictor sub-step advances
        dtfwd=0.5*dt          ! u, v, t to "n+1/2", the forward step
      else                    ! employed here for upstream advection
        dtfwd=dt              ! in the vicinity of open boundary is
      endif                   ! actually a halfstep.

# define dt illegal
# if defined M3_FRC_BRY || defined M3NUDGING
!      tau_in=dtfwd*tauM3_in
!      tau_out=dtfwd*tauM3_out
       tau_in=0.d0
       tau_out=0.d0
# endif

# ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
#  ifdef OBC_WEST
#   if defined OBC_M3ORLANSKI || defined OBC_M3ORLANSKI_WESTONLY
        do k=1,N                             ! Western edge radiation
          do j=jstr,jend+1                   ! ======= ==== =========
            grad(istr  ,j)=(u(istr  ,j,k,nstp)-u(istr  ,j-1,k,nstp))
#    ifdef MASKING
     &                                                *pmask(istr,j)
#    endif
            grad(istr+1,j)=(u(istr+1,j,k,nstp)-u(istr+1,j-1,k,nstp))
#    ifdef MASKING
     &                                              *pmask(istr+1,j)
#    endif
          enddo
          do j=jstr,jend
            dft=u(istr+1,j,k,nstp)-u(istr+1,j,k,nnew)
            dfx=u(istr+1,j,k,nnew)-u(istr+2,j,k,nnew)
 
            if (dfx*dft .lt. 0.) then
              dft=0.                       ! <-- cancel cx, if inflow
#    if defined M3_FRC_BRY || defined M3NUDGING
              tau=tau_in
            else
              tau=tau_out
#    endif
            endif
 
            if (dft*(grad(istr+1,j)+grad(istr+1,j+1)) .gt. 0.) then
              dfy=grad(istr+1,j)
            else
              dfy=grad(istr+1,j+1)
            endif
 
#    ifdef OBC_RAD_NORMAL
            dfy=0.
#    endif
            cff=max(dfx*dfx+dfy*dfy, eps)
            cx=dft*dfx
#    ifdef OBC_RAD_NPO
            cy=0.
#    else
            cy=min(cff,max(dft*dfy,-cff))
#    endif
 
            u(istr,j,k,nnew)=( cff*u(istr,j,k,nstp)
     &                        +cx*u(istr+1,j,k,nnew)
     &                    -max(cy,0.)*grad(istr,j  )
     &                    -min(cy,0.)*grad(istr,j+1)
     &                                   )/(cff+cx)
#    if defined M3_FRC_BRY  || defined M3NUDGING 
            u(istr,j,k,nnew)=(1.-tau)*u(istr,j,k,nnew)
#     ifdef M3_FRC_BRY
     &                                +tau*u_west(j,k)
#     else     
     &                             +tau*uclm(istr,j,k)
#     endif
#    endif
#    ifdef MASKING
            u(istr,j,k,nnew)=u(istr,j,k,nnew)*umask(istr,j)
#    endif
          enddo
        enddo
#   else            /* alternative open */
        do k=1,N
          do j=jstr,jend
#    if defined OBC_M3SPECIFIED && !defined OBC_M3ORLANSKI_WESTONLY
#     ifdef M3_FRC_BRY
            u(istr,j,k,nnew)=u_west(j,k)         ! specified
#     else
            u(istr,j,k,nnew)=uclm(istr,j,k)
#     endif
#    else
            u(istr,j,k,nnew)=u(istr+1,j,k,nnew)  ! Gradient: default
#    endif
#    ifdef MASKING
     &                           *umask(istr,j)
#    endif
          enddo
        enddo
#   endif
#  else           /* alternative to open */
        do k=1,N                               ! Western edge closed
          do j=jstr,jend                       ! ======= ==== ======
            u(istr,j,k,nnew)=0.                !  (no-flux: default)
          enddo
        enddo
#  endif              /* OBC_WEST */
      endif         !<-- WESTERN_EDGE

#ifdef ISWAKE
# define OBC_M3ORLANSKI
#endif
 
      if (EASTERN_EDGE) then
#  ifdef OBC_EAST
#   if defined OBC_M3ORLANSKI || defined OBC_M3ORLANSKI_EASTONLY
        do k=1,N                             ! Eastern edge radiation
          do j=jstr,jend+1                   ! ======= ==== =========
            grad(iend  ,j)=(u(iend  ,j,k,nstp)-u(iend  ,j-1,k,nstp))
#    ifdef MASKING
     &                                                *pmask(iend,j)
#    endif
            grad(iend+1,j)=(u(iend+1,j,k,nstp)-u(iend+1,j-1,k,nstp))
#    ifdef MASKING
     &                                              *pmask(iend+1,j)
#    endif
          enddo
          do j=jstr,jend
            dft=u(iend,j,k,nstp)-u(iend  ,j,k,nnew)
            dfx=u(iend,j,k,nnew)-u(iend-1,j,k,nnew)
 
            if (dfx*dft .lt. 0.) then
              dft=0.                      ! <-- cancel cx, if inflow
#    if defined M3_FRC_BRY || defined M3NUDGING
              tau=tau_in
            else
              tau=tau_out
#    endif
            endif
 
            if (dft*(grad(iend,j)+grad(iend,j+1)) .gt. 0.) then
              dfy=grad(iend,j)
            else
              dfy=grad(iend,j+1)
            endif
 
#    ifdef OBC_RAD_NORMAL
            dfy=0.
#    endif
            cff=max(dfx*dfx+dfy*dfy, eps)
            cx=dft*dfx
#    ifdef OBC_RAD_NPO
            cy=0.
#    else
            cy=min(cff,max(dft*dfy,-cff))
#    endif
 
            u(iend+1,j,k,nnew)=( cff*u(iend+1,j,k,nstp)
     &                              +cx*u(iend,j,k,nnew)
     &                      -max(cy,0.)*grad(iend+1,j  )
     &                      -min(cy,0.)*grad(iend+1,j+1)
     &                                       )/(cff+cx)
#    if defined M3_FRC_BRY  || defined M3NUDGING 
            u(iend+1,j,k,nnew)=(1.-tau)*u(iend+1,j,k,nnew)
#     ifdef M3_FRC_BRY
     &                                    +tau*u_east(j,k)
#     else     
     &                               +tau*uclm(iend+1,j,k)
#     endif
#    endif
#    ifdef MASKING
            u(iend+1,j,k,nnew)=u(iend+1,j,k,nnew)*umask(iend+1,j)
#    endif
          enddo
        enddo
#   else                /* alternative open */
        do k=1,N
          do j=jstr,jend
#    if defined OBC_M3SPECIFIED && !defined OBC_M3ORLANSKI_EASTONLY
#     ifdef M3_FRC_BRY
            u(iend+1,j,k,nnew)=u_east(j,k)       ! specified
#     else
            u(iend+1,j,k,nnew)=uclm(iend+1,j,k)
#     endif
#    else
            u(iend+1,j,k,nnew)=u(iend,j,k,nnew)  ! gradient (default)
#    endif
#    ifdef MASKING
     &                         *umask(iend+1,j)
#    endif
          enddo
        enddo
#   endif
#  else
        do k=1,N                                ! Eastern edge closed
          do j=jstr,jend                        ! ======= ==== ======
            u(iend+1,j,k,nnew)=0.               !  (no-flux: default)
          enddo
        enddo
#  endif
      endif         !<-- EASTERN_EDGE
# endif             /* !EW_PERIODIC */
 
 
 
# ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
#  ifdef OBC_SOUTH
#   ifdef OBC_M3ORLANSKI
        do k=1,N                            ! Southern edge radiation
          do i=istrU-1,iend                 ! ======== ==== =========
            grad(i,jstr-1)=u(i+1,jstr-1,k,nstp)-u(i,jstr-1,k,nstp)
            grad(i,jstr  )=u(i+1,jstr  ,k,nstp)-u(i,jstr  ,k,nstp)
          enddo
          do i=istrU,iend
/* --->
            dft=u(i,jstr,k,nstp)-u(i,jstr  ,k,nnew)
            dfx=u(i,jstr,k,nnew)-u(i,jstr+1,k,nnew)
 
            if (dfx*dft .lt. 0.) then
              dft=0.                       ! <-- cancel cx, if inflow
#    if defined M3_FRC_BRY || defined M3NUDGING
              tau=tau_in
            else
              tau=tau_out
#    endif
            endif
 
            if (dft*(grad(i-1,jstr)+grad(i,jstr)) .gt. 0.) then
              dfy=grad(i-1,jstr)
            else
              dfy=grad(i  ,jstr)
            endif
 
#    ifdef OBC_RAD_NORMAL
            dfy=0.
#    endif
            cff=max(dfx*dfx+dfy*dfy, eps)
            cx=dft*dfx
#    ifdef OBC_RAD_NPO
            cy=0.
#    else
            cy=min(cff,max(dft*dfy,-cff))
#    endif
            u(i,jstr-1,k,nnew)=( cff*u(i,jstr-1,k,nstp)
     &                              +cx*u(i,jstr,k,nnew)
     &                      -max(cy,0.)*grad(i-1,jstr-1)
     &                      -min(cy,0.)*grad(i  ,jstr-1)
     &                                       )/(cff+cx)

---> */


          cx=-0.125*dtfwd*(v(i,jstr,k,nrhs)+v(i-1,jstr,k,nrhs))
     &                            *( pn(i,jstr-1)+pn(i-1,jstr-1)
     &                                +pn(i,jstr)+pn(i-1,jstr) )

          cy= 0.125*dtfwd*(u(i,jstr-1,k,nrhs)+u(i,jstr,k,nrhs))
     &                            *( pm(i,jstr-1)+pm(i-1,jstr-1)
     &                                +pm(i,jstr)+pm(i-1,jstr) )

          if (cx.gt.0.) then
            tau=0.
          else
            tau=-cx
            cx=0.
          endif

          u(i,jstr-1,k,nnew)=(1.-cx)*(   u(i,jstr-1,k,nstp)
     &                          -max(cy,0.)*grad(i-1,jstr-1)
     &                          -min(cy,0.)*grad(i  ,jstr-1)
     &                                                     )
     &                       +cx*(         u(i,jstr,k,nstp)
     &                            -max(cy,0.)*grad(i-1,jstr)
     &                            -min(cy,0.)*grad(i  ,jstr)
     &                                                     )







#    if defined M3_FRC_BRY  || defined M3NUDGING 
           u(i,jstr-1,k,nnew)=(1.-tau)*u(i,jstr-1,k,nnew)
#     ifdef M3_FRC_BRY
     &                                  +tau*u_south(i,k)
#     else     
     &                              +tau*uclm(i,jstr-1,k)
#     endif
#    endif
#    ifdef MASKING
            u(i,jstr-1,k,nnew)=u(i,jstr-1,k,nnew)*umask(i,jstr-1)
#    endif
          enddo
        enddo
#   else
        do k=1,N
          do i=istrU,iend
#    ifdef OBC_M3SPECIFIED
#     ifdef M3_FRC_BRY
            u(i,jstr-1,k,nnew)=u_south(i,k)      ! specified
#     else
            u(i,jstr-1,k,nnew)=uclm(i,jstr-1,k)
#     endif
#    else
            u(i,jstr-1,k,nnew)=u(i,jstr,k,nnew)  ! gradient (default)
#    endif
#    ifdef MASKING
     &                         *umask(i,jstr-1)
#    endif
          enddo
        enddo
#   endif
#  else
#   ifdef EW_PERIODIC
#    define I_RANGE istrU,iend
#   else
#    define I_RANGE istr,iendR
#   endif
        do k=1,N                        ! Wall: free-slip (gamma2=+1)
          do i=I_RANGE                  ! =====   no-slip (gamma2=-1)
            u(i,jstr-1,k,nnew)=gamma2*u(i,jstr,k,nnew)
#   ifdef MASKING
     &                                *umask(i,jstr-1)
#   endif
          enddo
        enddo
#   undef I_RANGE
#  endif
      endif              !<-- SOUTHERN_EDGE
 
 
 
      if (NORTHERN_EDGE) then
#  ifdef OBC_NORTH
#   if defined OBC_M3ORLANSKI || defined OBC_M3ORLANSKI_NORTHONLY
        do k=1,N                            ! Northern edge radiation
          do i=istrU-1,iend                 ! ======== ==== =========
            grad(i,jend  )=u(i+1,jend  ,k,nstp)-u(i,jend  ,k,nstp)
            grad(i,jend+1)=u(i+1,jend+1,k,nstp)-u(i,jend+1,k,nstp)
          enddo
          do i=istrU,iend
/* --->
            dft=u(i,jend,k,nstp)-u(i,jend  ,k,nnew)
            dfx=u(i,jend,k,nnew)-u(i,jend-1,k,nnew)
 
            if (dfx*dft .lt. 0.) then
              dft=0.                                 ! <-- INFLOW
#    if defined M3_FRC_BRY || defined M3NUDGING
              tau=tau_in
            else
              tau=tau_out
#    endif
            endif
 
            if (dft*(grad(i-1,jend)+grad(i,jend)) .gt. 0.) then
              dfy=grad(i-1,jend)
            else
              dfy=grad(i  ,jend)
            endif
 
#    ifdef OBC_RAD_NORMAL
            dfy=0.
#    endif
            cff=max(dfx*dfx+dfy*dfy, eps)
            cx=dft*dfx
#    ifdef OBC_RAD_NPO
            cy=0.
#    else
            cy=min(cff,max(dft*dfy,-cff))
#    endif
 
            u(i,jend+1,k,nnew)=( cff*u(i,jend+1,k,nstp)
     &                              +cx*u(i,jend,k,nnew)
     &                      -max(cy,0.)*grad(i-1,jend+1)
     &                      -min(cy,0.)*grad(i  ,jend+1)
     &                                       )/(cff+cx)
---> */


          cx=0.125*dtfwd*(v(i,jend+1,k,nrhs)+v(i-1,jend+1,k,nrhs))
     &                               *( pn(i,jend+1)+pn(i-1,jend+1)
     &                                   +pn(i,jend)+pn(i-1,jend) )

          cy=0.125*dtfwd*(u(i,jend,k,nrhs)+u(i,jend+1,k,nrhs))
     &                               *( pm(i,jend+1)+pm(i-1,jend+1)
     &                                   +pm(i,jend)+pm(i-1,jend) )

          if (cx.gt.0.) then
            tau=0.
          else
            tau=-cx
            cx=0.
          endif

          u(i,jend+1,k,nnew)=(1.-cx)*(  u(i,jend+1,k,nstp)
     &                          -max(cy,0.)*grad(i-1,jend+1)
     &                          -min(cy,0.)*grad(i  ,jend+1)
     &                                                     )
     &                       +cx*(         u(i,jend,k,nstp)
     &                            -max(cy,0.)*grad(i-1,jend)
     &                            -min(cy,0.)*grad(i  ,jend)
     &                                                     )



#    if defined M3_FRC_BRY  || defined M3NUDGING 
            u(i,jend+1,k,nnew)=(1.-tau)*u(i,jend+1,k,nnew)
#     ifdef M3_FRC_BRY
     &                                   +tau*u_north(i,k)
#     else     
     &                               +tau*uclm(i,jend+1,k)
#     endif
#    endif
#    ifdef MASKING
            u(i,jend+1,k,nnew)=u(i,jend+1,k,nnew)*umask(i,jend+1)
#    endif
          enddo
        enddo
#   else               /* alternative open */
        do k=1,N
          do i=istrU,iend
#    if defined OBC_M3SPECIFIED && !defined OBC_M3ORLANSKI_NORTHONLY && !defined OBC_M3GRAD_NORTHONLY
#     ifdef M3_FRC_BRY
            u(i,jend+1,k,nnew)=u_north(i,k)      ! specified
#     else
            u(i,jend+1,k,nnew)=uclm(i,jend+1,k)
#     endif
#    else
            u(i,jend+1,k,nnew)=u(i,jend,k,nnew)  ! gradient (default)
#    endif
#    ifdef MASKING
     &                         *umask(i,jend+1)
#    endif
          enddo
        enddo
#   endif
#  else
#   ifdef EW_PERIODIC
#    define I_RANGE istrU,iend
#   else
#    define I_RANGE istr,iendR
#   endif
        do k=1,N                        ! Wall: free-slip (gamma2=+1)
          do i=I_RANGE                  ! =====   no-slip (gamma2=-1)
            u(i,jend+1,k,nnew)=gamma2*u(i,jend,k,nnew)
#   ifdef MASKING
     &                                *umask(i,jend+1)
#   endif
          enddo
        enddo
#   undef I_RANGE
#  endif
      endif   !<-- NORTHERN_EDGE
# endif        /* !NS_PERIODIC */
 
                           ! Corners between adjacent open boundaries
                           ! ======= ======= ======== ==== ==========
 
# if defined OBC_SOUTH && defined OBC_WEST
      if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
        do k=1,N
          u(istr,jstr-1,k,nnew)=0.5*( u(istr+1,jstr-1,k,nnew)
     &                               +u(istr  ,jstr  ,k,nnew))
        enddo
      endif
# endif
# if defined OBC_SOUTH && defined OBC_EAST
      if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
        do k=1,N
          u(iend+1,jstr-1,k,nnew)=0.5*( u(iend,jstr-1,k,nnew)
     &                                 +u(iend+1,jstr,k,nnew))
        enddo
      endif
# endif
# if defined OBC_NORTH && defined OBC_WEST
      if (WESTERN_EDGE .and. NORTHERN_EDGE) then
        do k=1,N
          u(istr,jend+1,k,nnew)=0.5*( u(istr+1,jend+1,k,nnew)
     &                               +u(istr  ,jend  ,k,nnew))
        enddo
      endif
# endif
# if defined OBC_NORTH && defined OBC_EAST
      if (EASTERN_EDGE .and. NORTHERN_EDGE) then
        do k=1,N
          u(iend+1,jend+1,k,nnew)=0.5*( u(iend,jend+1,k,nnew)
     &                                 +u(iend+1,jend,k,nnew))
        enddo
      endif
# endif
      return
      end
#else
      subroutine u3dbc_empty
      end
#endif /* SOLVE3D */
