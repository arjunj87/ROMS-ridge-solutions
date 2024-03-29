#include "cppdefs.h"
#ifdef SOLVE3D
 
      subroutine v3dbc_tile (istr,iend,jstr,jend, grad)
!
! Set lateral boundary conditions for ETA-component velocity
! v(:,:,:,nnew) 
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

# ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
#  ifdef OBC_SOUTH
#   ifdef OBC_M3ORLANSKI
        do k=1,N                            ! Southern edge radiation
          do i=istr,iend+1                  ! ======== ==== =========
            grad(i,jstr  )=(v(i,jstr  ,k,nstp)-v(i-1,jstr  ,k,nstp))
#    ifdef MASKING
     &                                                *pmask(i,jstr)
#    endif
            grad(i,jstr+1)=(v(i,jstr+1,k,nstp)-v(i-1,jstr+1,k,nstp))
#    ifdef MASKING
     &                                              *pmask(i,jstr+1)
#    endif
          enddo
          do i=istr,iend
            dft=v(i,jstr+1,k,nstp)-v(i,jstr+1,k,nnew)
            dfx=v(i,jstr+1,k,nnew)-v(i,jstr+2,k,nnew)
 
            if (dfx*dft .lt. 0.) then
              dft=0.                      ! <-- cancel cx, if inflow
#    if defined M3_FRC_BRY || defined M3NUDGING
              tau=tau_in
            else
              tau=tau_out
#    endif
            endif
 
            if (dft*(grad(i,jstr+1)+grad(i+1,jstr+1)) .gt. 0.) then
              dfy=grad(i,jstr+1)
            else
              dfy=grad(i+1,jstr+1)
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
 
            v(i,jstr,k,nnew)=( cff*v(i,jstr,k,nstp)
     &                        +cx*v(i,jstr+1,k,nnew)
     &                    -max(cy,0.)*grad(i  ,jstr)
     &                    -min(cy,0.)*grad(i+1,jstr)
     &                                   )/(cff+cx)
#    if defined M3_FRC_BRY  || defined M3NUDGING 
            v(i,jstr,k,nnew)=(1.-tau)*v(i,jstr,k,nnew)
#     ifdef M3_FRC_BRY
     &                               +tau*v_south(i,k)
#     else     
     &                             +tau*vclm(i,jstr,k)
#     endif
#    endif
#    ifdef MASKING
            v(i,jstr,k,nnew)=v(i,jstr,k,nnew)*vmask(i,jstr)
#    endif
          enddo
        enddo
#   else           /* alternative open */
        do k=1,N
          do i=istr,iend
#    ifdef OBC_M3SPECIFIED
#     ifdef M3_FRC_BRY
            v(i,jstr,k,nnew)=v_south(i,k)        ! specified
#     else
            v(i,jstr,k,nnew)=vclm(i,jstr,k)
#     endif
#    else
            v(i,jstr,k,nnew)=v(i,jstr+1,k,nnew)  ! gradient (default)
#    endif
#    ifdef MASKING
     &                           *vmask(i,jstr)
#    endif
          enddo
        enddo
#   endif
#  else
        do k=1,N                               ! Southern edge closed
          do i=istr,iend                       ! ======== ==== ======
            v(i,jstr,k,nnew)=0.                !  (no-flux: default)
          enddo
        enddo
#  endif              /* OBC_SOUTH */
      endif         !<-- SOUTHERN_EDGE
 
 
 
      if (NORTHERN_EDGE) then
#  ifdef OBC_NORTH
#   if defined OBC_M3ORLANSKI || defined OBC_M3ORLANSKI_NORTHONLY
        do k=1,N                            ! Northern edge radiation
          do i=istr,iend+1                  ! ======== ==== =========
            grad(i,jend  )=(v(i,jend  ,k,nstp)-v(i-1,jend  ,k,nstp))
#    ifdef MASKING
     &                                                *pmask(i,jend)
#    endif
            grad(i,jend+1)=(v(i,jend+1,k,nstp)-v(i-1,jend+1,k,nstp))
#    ifdef MASKING
     &                                              *pmask(i,jend+1)
#    endif
          enddo
          do i=istr,iend
            dft=v(i,jend,k,nstp)-v(i,jend  ,k,nnew)
            dfx=v(i,jend,k,nnew)-v(i,jend-1,k,nnew)
 
            if (dfx*dft .lt. 0.) then
              dft=0.                       ! <-- cancel cx, if inflow
#    if defined M3_FRC_BRY || defined M3NUDGING
              tau=tau_in
            else
              tau=tau_out
#    endif
            endif
 
            if (dft*(grad(i,jend)+grad(i+1,jend)) .gt. 0.) then
              dfy=grad(i,jend)
            else
              dfy=grad(i+1,jend)
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
 
            v(i,jend+1,k,nnew)=( cff*v(i,jend+1,k,nstp)
     &                              +cx*v(i,jend,k,nnew)
     &                      -max(cy,0.)*grad(i  ,jend+1)
     &                      -min(cy,0.)*grad(i+1,jend+1)
     &                                      )/(cff+cx)
#    if defined M3_FRC_BRY  || defined M3NUDGING 
            v(i,jend+1,k,nnew)=(1.-tau)*v(i,jend+1,k,nnew)
#     ifdef M3_FRC_BRY
     &                                   +tau*v_north(i,k)
#     else     
     &                               +tau*vclm(i,jend+1,k)
#     endif
#    endif
#    ifdef MASKING
            v(i,jend+1,k,nnew)=v(i,jend+1,k,nnew)*vmask(i,jend+1)
#    endif
          enddo
        enddo
#   else                /* alternative open */
        do k=1,N
          do i=istr,iend
#    if defined OBC_M3SPECIFIED && !defined OBC_M3ORLANSKI_NORTHONLY && !defined OBC_M3GRAD_NORTHONLY
#     ifdef M3_FRC_BRY
            v(i,jend+1,k,nnew)=v_north(i,k)      ! specified
#     else
            v(i,jend+1,k,nnew)=vclm(i,jend+1,k)
#     endif
#    else
            v(i,jend+1,k,nnew)=v(i,jend,k,nnew)  ! gradient (default)
#    endif
#    ifdef MASKING
     &                         *vmask(i,jend+1)
#    endif
          enddo
        enddo
#   endif
#  else
        do k=1,N                               ! Northern edge closed
          do i=istr,iend                       ! ======== ==== ======
            v(i,jend+1,k,nnew)=0.              !   (no-flux: default)
          enddo
        enddo
#  endif
      endif     !<--  NORTHERN_EDGE
# endif          /* !NS_PERIODIC */
 
 
 
# ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
#  ifdef OBC_WEST
#   if defined OBC_M3ORLANSKI || defined OBC_M3ORLANSKI_WESTONLY
        do k=1,N                             ! Western edge radiation
          do j=jstrV-1,jend                  ! ======= ==== =========
            grad(istr-1,j)=v(istr-1,j+1,k,nstp)-v(istr-1,j,k,nstp)
            grad(istr  ,j)=v(istr  ,j+1,k,nstp)-v(istr  ,j,k,nstp)
          enddo
          do j=jstrV,jend
/* --->
            dft=v(istr,j,k,nstp)-v(istr  ,j,k,nnew)
            dfx=v(istr,j,k,nnew)-v(istr+1,j,k,nnew)
 
            if (dfx*dft .lt. 0.) then
              dft=0.                       ! <-- cancel cx, if inflow
#    if defined M3_FRC_BRY || defined M3NUDGING
              tau=tau_in
            else
              tau=tau_out
#    endif
            endif
 
            if (dft*(grad(istr,j-1)+grad(istr,j)) .gt. 0.) then
              dfy=grad(istr,j-1)
            else
              dfy=grad(istr,j  )
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
 
            v(istr-1,j,k,nnew)=( cff*v(istr-1,j,k,nstp)
     &                              +cx*v(istr,j,k,nnew)
     &                      -max(cy,0.)*grad(istr-1,j-1)
     &                      -min(cy,0.)*grad(istr-1,j  )
     &                                       )/(cff+cx)
---> */


          cx=-0.125*dtfwd*(u(istr,j,k,nrhs)+u(istr,j-1,k,nrhs))
     &                            *( pm(istr-1,j)+pm(istr-1,j-1)
     &                                +pm(istr,j)+pm(istr,j-1) )

          cy= 0.125*dtfwd*(v(istr-1,j,k,nrhs)+v(istr,j,k,nrhs))
     &                            *( pn(istr-1,j)+pn(istr-1,j-1)
     &                                +pn(istr,j)+pn(istr,j-1) )

          if (cx.gt.0.) then
            tau=0.
          else
            tau=-cx
            cx=0.
          endif

          v(istr-1,j,k,nnew)=(1.-cx)*( v(istr-1,j,k,nstp)
     &                          -max(cy,0.)*grad(istr-1,j-1)
     &                          -min(cy,0.)*grad(istr-1,j  )
     &                                                     )
     &                       +cx*(        v(istr,j,k,nstp)
     &                            -max(cy,0.)*grad(istr,j-1)
     &                            -min(cy,0.)*grad(istr,j  )
     &                                                     )




#    if defined M3_FRC_BRY  || defined M3NUDGING 
            v(istr-1,j,k,nnew)=(1.-tau)*v(istr-1,j,k,nnew)
#     ifdef M3_FRC_BRY
     &                                    +tau*v_west(j,k)
#     else     
     &                               +tau*vclm(istr-1,j,k)
#     endif
#    endif
#    ifdef MASKING
            v(istr-1,j,k,nnew)=v(istr-1,j,k,nnew)*vmask(istr-1,j)
#    endif
          enddo
        enddo
#   else
        do k=1,N
          do j=jstrV,jend
#    if defined OBC_M3SPECIFIED && !defined OBC_M3ORLANSKI_WESTONLY 
#     ifdef M3_FRC_BRY
            v(istr-1,j,k,nnew)=v_west(j,k)       ! specified
#     else
            v(istr-1,j,k,nnew)=vclm(istr-1,j,k)
#     endif
#    else
            v(istr-1,j,k,nnew)=v(istr,j,k,nnew)  ! gradient (default)
#    endif
#    ifdef MASKING
     &                         *vmask(istr-1,j)
#    endif
          enddo
        enddo
#   endif
#  else
#   ifdef NS_PERIODIC
#    define J_RANGE jstrV,jend
#   else
#    define J_RANGE jstr,jendR
#   endif
        do k=1,N                        ! Wall: free-slip (gamma2=+1)
          do j=J_RANGE                  ! =====   no-slip (gamma2=-1)
            v(istr-1,j,k,nnew)=gamma2*v(istr,j,k,nnew)
#   ifdef MASKING
     &                                *vmask(istr-1,j)
#   endif
          enddo
        enddo
#   undef J_RANGE
#  endif
      endif          !<-- WESTERN_EDGE
 
#ifdef ISWAKE
# define OBC_M3ORLANSKI
#endif
 
      if (EASTERN_EDGE) then
#  ifdef OBC_EAST
#   if defined OBC_M3ORLANSKI || defined OBC_M3ORLANSKI_EASTONLY
        do k=1,N                             ! Eastern edge radiation
          do j=jstrV-1,jend                  ! ======= ==== =========
            grad(iend  ,j)=v(iend  ,j+1,k,nstp)-v(iend  ,j,k,nstp)
            grad(iend+1,j)=v(iend+1,j+1,k,nstp)-v(iend+1,j,k,nstp)
          enddo
          do j=jstrV,jend
/* --->
            dft=v(iend,j,k,nstp)-v(iend  ,j,k,nnew)
            dfx=v(iend,j,k,nnew)-v(iend-1,j,k,nnew)
 
            if (dfx*dft .lt. 0.) then
              dft=0.                       ! <-- cancel cx, if inflow
#    if defined M3_FRC_BRY || defined M3NUDGING
              tau=tau_in
            else
              tau=tau_out
#    endif
            endif
 
            if (dft*(grad(iend,j-1)+grad(iend,j)) .gt. 0.) then
              dfy=grad(iend,j-1)
            else
              dfy=grad(iend,j  )
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
 
            v(iend+1,j,k,nnew)=( cff*v(iend+1,j,k,nstp)
     &                              +cx*v(iend,j,k,nnew)
     &                      -max(cy,0.)*grad(iend+1,j-1)
     &                      -min(cy,0.)*grad(iend+1,j  )
     &                                       )/(cff+cx)

---> */

          cx=0.125*dtfwd*(u(iend+1,j,k,nrhs)+u(iend+1,j-1,k,nrhs))
     &                               *( pm(iend+1,j)+pm(iend+1,j-1)
     &                                   +pm(iend,j)+pm(iend,j-1) )

          cy=0.125*dtfwd*(v(iend,j,k,nrhs)+v(iend+1,j,k,nrhs))
     &                               *( pn(iend+1,j)+pn(iend+1,j-1)
     &                                   +pn(iend,j)+pn(iend,j-1) )

          if (cx.gt.0.) then
            tau=0.
          else
            tau=-cx
            cx=0.
          endif

          v(iend+1,j,k,nnew)=(1.-cx)*(   v(iend+1,j,k,nstp)
     &                          -max(cy,0.)*grad(iend+1,j-1)
     &                          -min(cy,0.)*grad(iend+1,j  )
     &                                                     )
     &                       +cx*(         v(iend,j,k,nstp)
     &                            -max(cy,0.)*grad(iend,j-1)
     &                            -min(cy,0.)*grad(iend,j  )
     &                                                     )



#    if defined M3_FRC_BRY  || defined M3NUDGING 
            v(iend+1,j,k,nnew)=(1.-tau)*v(iend+1,j,k,nnew)
#     ifdef M3_FRC_BRY
     &                                    +tau*v_east(j,k)
#     else     
     &                               +tau*vclm(iend+1,j,k)
#     endif
#    endif
#    ifdef MASKING
            v(iend+1,j,k,nnew)=v(iend+1,j,k,nnew)*vmask(iend+1,j)
#    endif
          enddo
        enddo
#   else               /* alternative open */
        do k=1,N
          do j=jstr,jend
#    if defined OBC_M3SPECIFIED && !defined OBC_M3ORLANSKI_EASTONLY
#     ifdef M3_FRC_BRY
            v(iend+1,j,k,nnew)=v_east(j,k)       ! specified
#     else
            v(iend+1,j,k,nnew)=vclm(iend+1,j,k)
#     endif
#    else
            v(iend+1,j,k,nnew)=v(iend,j,k,nnew)  ! gradient (default)
#    endif
#    ifdef MASKING
     &                         *vmask(iend+1,j)
#    endif
          enddo
        enddo
#   endif
#  else
#   ifdef NS_PERIODIC
#    define J_RANGE jstrV,jend
#   else
#    define J_RANGE jstr,jendR
#   endif
        do k=1,N                        ! Wall: free-slip (gamma2=+1)
          do j=J_RANGE                  ! ====    no-slip (gamma2=-1)
            v(iend+1,j,k,nnew)=gamma2*v(iend,j,k,nnew)
#   ifdef MASKING
     &                                *vmask(iend+1,j)
#   endif
          enddo
        enddo
#   undef J_RANGE
#  endif
      endif     !<-- EASTERN_EDGE
# endif          /* !EW_PERIODIC */
 
                           ! Corners between adjacent open boundaries
                           ! ======= ======= ======== ==== ==========
 
# if defined OBC_SOUTH && defined OBC_WEST
      if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
        do k=1,N
          v(istr-1,jstr,k,nnew)=0.5*( v(istr-1,jstr+1,k,nnew)
     &                               +v(istr  ,jstr  ,k,nnew))
        enddo
      endif
# endif
# if defined OBC_SOUTH && defined OBC_EAST
      if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
        do k=1,N
          v(iend+1,jstr,k,nnew)=0.5*( v(iend+1,jstr+1,k,nnew)
     &                               +v(iend  ,jstr  ,k,nnew))
        enddo
      endif
# endif
# if defined OBC_NORTH && defined OBC_WEST
      if (WESTERN_EDGE .and. NORTHERN_EDGE) then
        do k=1,N
          v(istr-1,jend+1,k,nnew)=0.5*( v(istr-1,jend,k,nnew)
     &                                 +v(istr,jend+1,k,nnew))
        enddo
      endif
# endif
# if defined OBC_NORTH && defined OBC_EAST
      if (EASTERN_EDGE .and. NORTHERN_EDGE) then
        do k=1,N
          v(iend+1,jend+1,k,nnew)=0.5*( v(iend+1,jend,k,nnew)
     &                                 +v(iend,jend+1,k,nnew))
        enddo
      endif
# endif
      return
      end
#else
      subroutine v3dbc_empty
      end
#endif /* SOLVE3D */
 
