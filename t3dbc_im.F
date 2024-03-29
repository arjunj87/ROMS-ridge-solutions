#include "cppdefs.h"
#ifdef SOLVE3D
 
      subroutine t3dbc_tile (istr,iend,jstr,jend, itrc, grad)
!
! Set lateral boundary conditions for tracer field t(:,:,:,itrc)
!
      implicit none
      integer istr,iend,jstr,jend, itrc, i,j,k
      real grad(PRIVATE_2D_SCRATCH_ARRAY), eps, dtfwd,
     &     cff, cx,cy, dft,dfx,dfy, tau,tau_in,tau_out
      parameter (eps=1.E-20)
# include "param.h"
# include "grid.h"
# include "ocean3d.h"
# include "climat.h"
# include "scalars.h"
# include "boundary.h"
!
# include "compute_auxiliary_bounds.h"
!
      if (PRE_STEP) then      ! because predictor sub-step advances
        dtfwd=0.5*dt          ! u, v, t to "n+1/2", the forward step
      else                    ! employed here for upstream advection
        dtfwd=dt              ! in the vicinity of open boundary is
      endif                   ! actually a halfstep.

# define dt illegal
# if defined T_FRC_BRY || defined TNUDGING
      tau_in=dtfwd*tauT_in
      tau_out=dtfwd*tauT_out
# endif

# ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
#  if defined OBC_WEST && defined OBC_TORLANSKI || defined OBC_TORLANSKI_WESTONLY
        do k=1,N
          do j=jstr,jend+1
            grad(istr-1,j)=( t(istr-1,j  ,k,nstp,itrc)
     &                      -t(istr-1,j-1,k,nstp,itrc))
#   ifdef MASKING
     &                                 *vmask(istr-1,j)
#   endif
            grad(istr  ,j)=( t(istr  ,j  ,k,nstp,itrc)
     &                      -t(istr  ,j-1,k,nstp,itrc))
#   ifdef MASKING
     &                                   *vmask(istr,j)
#   endif
          enddo
          do j=jstr,jend
/* --->
            dft=t(istr,j,k,nstp,itrc)-t(istr  ,j,k,nnew,itrc)
            dfx=t(istr,j,k,nnew,itrc)-t(istr+1,j,k,nnew,itrc)
 
            if (dfx*dft .lt. 0.) then
              dft=0.                  ! <-- cancel, if inflow 
#   if defined T_FRC_BRY || defined TNUDGING
              tau=tau_in
            else
              tau=tau_out
#   endif
            endif
 
            if (dft*(grad(istr,j)+grad(istr,j+1)) .gt. 0.) then
              dfy=grad(istr,j)
            else
              dfy=grad(istr,j+1)
            endif
 
#   ifdef OBC_RAD_NORMAL
            dfy=0.
#   endif
            cff=max(dfx*dfx+dfy*dfy, eps)
            cx=dft*dfx
#   ifdef OBC_RAD_NPO
            cy=0.
#   else
            cy=min(cff,max(dft*dfy,-cff))
#   endif
            t(istr-1,j,k,nnew,itrc)=( cff*t(istr-1,j,k,nstp,itrc)
     &                                   +cx*t(istr,j,k,nnew,itrc)
     &                                -max(cy,0.)*grad(istr-1,j  )
     &                                -min(cy,0.)*grad(istr-1,j+1)
     &                                                 )/(cff+cx)

---> */



            cx=-dtfwd*u(istr,j,k,nrhs)*pm(istr-1,j)
            cy=0.5*dtfwd*(v(istr-1,j,k,nrhs)+v(istr-1,j+1,k,nrhs))
     &                                               *pn(istr-1,j)

            if (cx.gt.0.) then
              tau=0.
            else
              tau=-cx
              cx=0.
            endif

            t(istr-1,j,k,nnew,itrc)=(1.-cx)*( t(istr-1,j,k,nstp,itrc)
     &                                   -max(cy,0.)*grad(istr-1,j  )
     &                                   -min(cy,0.)*grad(istr-1,j+1)
     &                                                              )
     &                                  +cx*(  t(istr,j,k,nstp,itrc)
     &                                     -max(cy,0.)*grad(istr,j  )
     &                                     -min(cy,0.)*grad(istr,j+1)
     &                                                              )



#   if defined T_FRC_BRY || defined TNUDGING
            t(istr-1,j,k,nnew,itrc)=(1.-tau)*t(istr-1,j,k,nnew,itrc)
#    ifdef T_FRC_BRY
     &                                    +tau*t_west(j,k,itrc)     
#    else
     &                                    +tau*tclm(istr-1,j,k,itrc)
#    endif
#   endif
#   ifdef MASKING
            t(istr-1,j,k,nnew,itrc)=t(istr-1,j,k,nnew,itrc)
     &                                      *rmask(istr-1,j)
#   endif
          enddo
        enddo
#  else
        do k=1,N
          do j=jstr,jend
#   if defined OBC_WEST && defined OBC_TSPECIFIED && !defined OBC_TORLANSKI_WESTONLY
#    ifdef T_FRC_BRY
            t(istr-1,j,k,nnew,itrc)=t_west(j,k,itrc)
#    else
            t(istr-1,j,k,nnew,itrc)=tclm(istr-1,j,k,itrc)
#    endif
#   else
            t(istr-1,j,k,nnew,itrc)=t(istr,j,k,nnew,itrc)
#   endif
#   ifdef MASKING
     &                                   *rmask(istr-1,j)
#   endif
          enddo
        enddo
#  endif
      endif     ! <-- WESTERN_EDGE

#ifdef ISWAKE
# define OBC_TORLANSKI
#endif

      if (EASTERN_EDGE) then
#  if defined OBC_EAST && defined OBC_TORLANSKI || defined OBC_TORLANSKI_EASTONLY
!
!                                        !  Eastern edge radiation BC
        do k=1,N                         !  ======= ==== ========= ==
          do j=jstr,jend+1
           grad(iend  ,j)=( t(iend  ,j  ,k,nstp,itrc)
     &                     -t(iend  ,j-1,k,nstp,itrc))
#   ifdef MASKING
     &                                  *vmask(iend,j)
#   endif
           grad(iend+1,j)=( t(iend+1,j  ,k,nstp,itrc)
     &                     -t(iend+1,j-1,k,nstp,itrc))
#   ifdef MASKING
     &                                *vmask(iend+1,j)
#   endif
          enddo
          do j=jstr,jend
/* --->
            dft=t(iend,j,k,nstp,itrc)-t(iend  ,j,k,nnew,itrc)
            dfx=t(iend,j,k,nnew,itrc)-t(iend-1,j,k,nnew,itrc)
 
            if (dfx*dft .lt. 0.) then
              dft=0.                  ! <-- cancel, if inflow 
#   if defined T_FRC_BRY || defined TNUDGING
              tau=tau_in
            else
              tau=tau_out
#   endif
            endif
 
            if (dft*(grad(iend,j)+grad(iend,j+1)) .gt. 0.) then
              dfy=grad(iend,j)
            else
              dfy=grad(iend,j+1)
            endif
 
#   ifdef OBC_RAD_NORMAL
            dfy=0.
#   endif
            cff=max(dfx*dfx+dfy*dfy, eps)
            cx=dft*dfx
#   ifdef OBC_RAD_NPO
            cy=0.
#   else
            cy=min(cff,max(dft*dfy,-cff))
#   endif
 
            t(iend+1,j,k,nnew,itrc)=( cff*t(iend+1,j,k,nstp,itrc)
     &                                   +cx*t(iend,j,k,nnew,itrc)
     &                                -max(cy,0.)*grad(iend+1,j  )
     &                                -min(cy,0.)*grad(iend+1,j+1)
     &                                                 )/(cff+cx)

---> */


            cx=dtfwd*u(iend+1,j,k,nrhs)*pm(iend+1,j)
            cy=0.5*dtfwd*(v(iend+1,j,k,nrhs)+v(iend+1,j+1,k,nrhs))
     &                                               *pn(iend+1,j)
          
            if (cx.gt.0.) then
              tau=0.
            else
              tau=-cx
              cx=0.
            endif

            t(iend+1,j,k,nnew,itrc)=(1.-cx)*( t(iend+1,j,k,nstp,itrc)
     &                                   -max(cy,0.)*grad(iend+1,j  )
     &                                   -min(cy,0.)*grad(iend+1,j+1)
     &                                                              )
     &                                  +cx*(  t(iend,j,k,nnew,itrc)
     &                                     -max(cy,0.)*grad(iend,j  )
     &                                     -min(cy,0.)*grad(iend,j+1)
     &                                                              )


#   if defined T_FRC_BRY || defined TNUDGING
            t(iend+1,j,k,nnew,itrc)=(1.-tau)*t(iend+1,j,k,nnew,itrc)
#    ifdef T_FRC_BRY
     &                                         +tau*t_east(j,k,itrc)
#    else
     &                                    +tau*tclm(iend+1,j,k,itrc)
#    endif
#   endif
#   ifdef MASKING
            t(iend+1,j,k,nnew,itrc)=t(iend+1,j,k,nnew,itrc)
     &                                     *rmask(iend+1,j)
#   endif
          enddo
        enddo
#  else
        do k=1,N
          do j=jstr,jend
#   if defined OBC_EAST && defined OBC_TSPECIFIED && !defined OBC_TORLANSKI_EASTONLY
#    ifdef T_FRC_BRY
            t(iend+1,j,k,nnew,itrc)=t_east(j,k,itrc) 
#    else
            t(iend+1,j,k,nnew,itrc)=tclm(iend+1,j,k,itrc)
#    endif
#   else
            t(iend+1,j,k,nnew,itrc)=t(iend,j,k,nnew,itrc)
#   endif
#   ifdef MASKING
     &                                   *rmask(iend+1,j)
#   endif
          enddo
        enddo
#  endif
      endif    ! <-- EASTERN_EDGE
# endif        /* !EW_PERIODIC */






 
# ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
#  if defined OBC_SOUTH && defined OBC_TORLANSKI
        do k=1,N
          do i=istr,iend+1
            grad(i,jstr  )=( t(i  ,jstr  ,k,nstp,itrc)
     &                      -t(i-1,jstr  ,k,nstp,itrc))
#   ifdef MASKING
     &                                   *umask(i,jstr)
#   endif
            grad(i,jstr-1)=( t(i  ,jstr-1,k,nstp,itrc)
     &                      -t(i-1,jstr-1,k,nstp,itrc))
#   ifdef MASKING
     &                                *umask(i,jstr-1)
#   endif
          enddo
          do i=istr,iend
/* --->
            dft=t(i,jstr,k,nstp,itrc)-t(i,jstr  ,k,nnew,itrc)
            dfx=t(i,jstr,k,nnew,itrc)-t(i,jstr+1,k,nnew,itrc)
 
            if (dfx*dft .lt. 0.) then
              dft=0.                   ! <-- cancel, if inflow 
#   if defined T_FRC_BRY || defined TNUDGING
              tau=tau_in
            else
              tau=tau_out
#   endif
            endif
 
            if (dft*(grad(i,jstr)+grad(i+1,jstr)) .gt. 0.) then
              dfy=grad(i,jstr)
            else
              dfy=grad(i+1,jstr)
            endif
 
#   ifdef OBC_RAD_NORMAL
            dfy=0.
#   endif
            cff=max(dfx*dfx+dfy*dfy, eps)
            cx=dft*dfx
#   ifdef OBC_RAD_NPO
            cy=0.
#   else
            cy=min(cff,max(dft*dfy,-cff))
#   endif
 
            t(i,jstr-1,k,nnew,itrc)=( cff*t(i,jstr-1,k,nstp,itrc)
     &                                   +cx*t(i,jstr,k,nnew,itrc)
     &                                -max(cy,0.)*grad(i  ,jstr-1)
     &                                -min(cy,0.)*grad(i+1,jstr-1)
     &                                                 )/(cff+cx)

---> */


            cx=-dtfwd*v(i,jstr,k,nrhs)*pn(i,jstr-1)
            cy=0.5*dtfwd*(u(i,jstr-1,k,nrhs)+u(i+1,jstr-1,k,nrhs))
     &                                            *pm(i,jstr-1)

            if (cx.gt.0.) then
              tau=0.
            else
              tau=-cx
              cx=0.
            endif

            t(i,jstr-1,k,nnew,itrc)=(1.-cx)*( t(i,jstr-1,k,nstp,itrc)
     &                                   -max(cy,0.)*grad(i  ,jstr-1)
     &                                   -min(cy,0.)*grad(i+1,jstr-1)
     &                                                              )
     &                                  +cx*(  t(i,jstr,k,nstp,itrc)
     &                                     -max(cy,0.)*grad(i  ,jstr)
     &                                     -min(cy,0.)*grad(i+1,jstr)
     &                                                              )


#   if defined T_FRC_BRY || defined TNUDGING
            t(i,jstr-1,k,nnew,itrc)=(1.-tau)*t(i,jstr-1,k,nnew,itrc)
#    ifdef T_FRC_BRY
     &                                        +tau*t_south(i,k,itrc)
#    else
     &                                    +tau*tclm(i,jstr-1,k,itrc)
#    endif
#   endif
#   ifdef MASKING
            t(i,jstr-1,k,nnew,itrc)=t(i,jstr-1,k,nnew,itrc)
     &                                     *rmask(i,jstr-1)
#   endif
          enddo
        enddo
#  else
        do k=1,N
          do i=istr,iend
#   if defined OBC_SOUTH && defined OBC_TSPECIFIED
#    ifdef T_FRC_BRY
            t(i,jstr-1,k,nnew,itrc)=t_south(i,k,itrc)
#    else
            t(i,jstr-1,k,nnew,itrc)=tclm(i,jstr-1,k,itrc)
#    endif
#   else
            t(i,jstr-1,k,nnew,itrc)=t(i,jstr,k,nnew,itrc)
#   endif
#   ifdef MASKING
     &                                   *rmask(i,jstr-1)
#   endif
          enddo
        enddo
#  endif
      endif    ! <-- SOUTHERN_EDGE


      if (NORTHERN_EDGE) then
#  if defined OBC_NORTH && defined OBC_TORLANSKI || defined OBC_TORLANSKI_NORTHONLY
        do k=1,N
          do i=istr,iend+1
            grad(i,jend  )=( t(i  ,jend  ,k,nstp,itrc)
     &                      -t(i-1,jend  ,k,nstp,itrc))
#   ifdef MASKING
     &                                   *umask(i,jend)
#   endif
            grad(i,jend+1)=( t(i  ,jend+1,k,nstp,itrc)
     &                      -t(i-1,jend+1,k,nstp,itrc))
#   ifdef MASKING
     &                                 *umask(i,jend+1)
#   endif
          enddo
          do i=istr,iend
/* --->
            dft=t(i,jend,k,nstp,itrc)-t(i,jend  ,k,nnew,itrc)
            dfx=t(i,jend,k,nnew,itrc)-t(i,jend-1,k,nnew,itrc)
 
            if (dfx*dft .lt. 0.) then
              dft=0.                          ! <-- cancel, if inflow 
#   if defined T_FRC_BRY || defined TNUDGING
              tau=tau_in
            else
              tau=tau_out
#   endif
            endif
 
            if (dft*(grad(i,jend)+grad(i+1,jend)) .gt. 0.) then
              dfy=grad(i,jend)
            else
              dfy=grad(i+1,jend)
            endif
 
#   ifdef OBC_RAD_NORMAL
            dfy=0.
#   endif
            cff=max(dfx*dfx+dfy*dfy, eps)
            cx=dft*dfx
#   ifdef OBC_RAD_NPO
            cy=0.
#   else
            cy=min(cff,max(dft*dfy,-cff))
#   endif
            t(i,jend+1,k,nnew,itrc)=( cff*t(i,jend+1,k,nstp,itrc)
     &                                 +cx*t(i,jend  ,k,nnew,itrc)
     &                                -max(cy,0.)*grad(i  ,jend+1)
     &                                -min(cy,0.)*grad(i+1,jend+1)
     &                                                 )/(cff+cx)

---> */
            cx=dtfwd*v(i,jend+1,k,nrhs)*pn(i,jend+1)
            cy=0.5*dtfwd*(u(i,jend+1,k,nrhs)+u(i+1,jend+1,k,nrhs))
     &                                            *pm(i,jend+1)

            if (cx.gt.0.) then
              tau=0.
            else
              tau=-cx
              cx=0.
            endif

            t(i,jend+1,k,nnew,itrc)=(1.-cx)*( t(i,jend+1,k,nstp,itrc)
     &                                   -max(cy,0.)*grad(i  ,jend+1)
     &                                   -min(cy,0.)*grad(i+1,jend+1)
     &                                                              )
     &                                  +cx*(  t(i,jend,k,nnew,itrc)
     &                                     -max(cy,0.)*grad(i  ,jend)
     &                                     -min(cy,0.)*grad(i+1,jend)
     &                                                              )





#   if defined T_FRC_BRY || defined TNUDGING
            t(i,jend+1,k,nnew,itrc)=(1.-tau)*t(i,jend+1,k,nnew,itrc)
#    ifdef T_FRC_BRY
     &                                        +tau*t_north(i,k,itrc)
#    else
     &                                    +tau*tclm(i,jend+1,k,itrc)
#    endif
#   endif
#   ifdef MASKING
            t(i,jend+1,k,nnew,itrc)=t(i,jend+1,k,nnew,itrc)
     &                                     *rmask(i,jend+1)
#   endif
          enddo
        enddo
#  else
        do k=1,N
          do i=istr,iend
#   if defined OBC_NORTH && defined OBC_TSPECIFIED && !defined OBC_TORLANSKI_NORTHONLY && !defined OBC_M3GRAD_NORTHONLY
#    ifdef T_FRC_BRY
            t(i,jend+1,k,nnew,itrc)=t_north(i,k,itrc)
#    else
            t(i,jend+1,k,nnew,itrc)=tclm(i,jend+1,k,itrc)
#    endif
#   else
            t(i,jend+1,k,nnew,itrc)=t(i,jend,k,nnew,itrc)
#   endif
#   ifdef MASKING
     &                                   *rmask(i,jend+1)
#   endif
          enddo
        enddo
#  endif
      endif    ! <-- NORTHERN_EDGE
# endif /* ! NS_PERIODIC */
!
! Corner points between adjacent boundaries. Note that because
! boundary conditions for tracers are of gradient nature -- either
! Neumann, or radiation, the corner points have to be set in all
! cases, even if the adjacent boundaries are closed. This differs
! from setting boundaries for velocities, where there is anisotropy
! associated with normal and tangential components. In the case when 
! one/or both points near the corner is/are masked, the use of masked
! points is avoided. 

# ifndef EW_PERIODIC
#  ifndef NS_PERIODIC
      if (SOUTHERN_EDGE .and. WESTERN_EDGE) then
#   ifdef MASKING
        cff=rmask(istr,jstr-1)+rmask(istr-1,jstr)
        if (cff.gt.0.) then
          cff=1./cff
          do k=1,N
            t(istr-1,jstr-1,k,nnew,itrc)=cff*(
     &              rmask(istr,jstr-1)*t(istr,jstr-1,k,nnew,itrc)
     &             +rmask(istr-1,jstr)*t(istr-1,jstr,k,nnew,itrc))
          enddo
        else
          do k=1,N
            t(istr-1,jstr-1,k,nnew,itrc)=0.
          enddo
        endif
#   else
        do k=1,N
          t(istr-1,jstr-1,k,nnew,itrc)=0.5*( t(istr,jstr-1,k,nnew,
     &                           itrc)+t(istr-1,jstr,k,nnew,itrc))
        enddo
#   endif
      endif

      if (SOUTHERN_EDGE .and. EASTERN_EDGE) then
#   ifdef MASKING
        cff=rmask(iend,jstr-1)+rmask(iend+1,jstr)
        if (cff.gt.0.) then
          cff=1./cff
          do k=1,N
            t(iend+1,jstr-1,k,nnew,itrc)=cff*(
     &              rmask(iend,jstr-1)*t(iend,jstr-1,k,nnew,itrc)
     &             +rmask(iend+1,jstr)*t(iend+1,jstr,k,nnew,itrc))
          enddo
        else
          do k=1,N
            t(iend+1,jstr-1,k,nnew,itrc)=0.
          enddo
        endif
#   else
        do k=1,N
          t(iend+1,jstr-1,k,nnew,itrc)=0.5*(t(iend,jstr-1,k,nnew,
     &                           itrc)+t(iend+1,jstr,k,nnew,itrc))
        enddo
#   endif
      endif

      if (NORTHERN_EDGE .and. WESTERN_EDGE) then
#   ifdef MASKING
        cff=rmask(istr,jend+1)+rmask(istr-1,jend)
        if (cff.gt.0.) then
          cff=1./cff
          do k=1,N
            t(istr-1,jend+1,k,nnew,itrc)=cff*(
     &              rmask(istr,jend+1)*t(istr,jend+1,k,nnew,itrc)
     &             +rmask(istr-1,jend)*t(istr-1,jend,k,nnew,itrc))
          enddo
        else
          do k=1,N
            t(istr-1,jend+1,k,nnew,itrc)=0.
          enddo
        endif
#   else
        do k=1,N
          t(istr-1,jend+1,k,nnew,itrc)=0.5*( t(istr,jend+1,k,nnew,
     &                           itrc)+t(istr-1,jend,k,nnew,itrc))
        enddo
#   endif
      endif

      if (NORTHERN_EDGE .and. EASTERN_EDGE) then
#   ifdef MASKING
        cff=rmask(iend,jend+1)+rmask(iend+1,jend)
        if (cff.gt.0.) then
          cff=1./cff
          do k=1,N
            t(iend+1,jend+1,k,nnew,itrc)=cff*(
     &              rmask(iend,jend+1)*t(iend,jend+1,k,nnew,itrc)
     &             +rmask(iend+1,jend)*t(iend+1,jend,k,nnew,itrc))
          enddo
        else
          do k=1,N
            t(iend+1,jend+1,k,nnew,itrc)=0.
          enddo
        endif
#   else
        do k=1,N
          t(iend+1,jend+1,k,nnew,itrc)=0.5*( t(iend,jend+1,k,nnew,
     &                           itrc)+t(iend+1,jend,k,nnew,itrc))
        enddo
#   endif
      endif
#  endif
# endif
      return
      end
#else
      subroutine t3dbc_empty
      end
#endif /* SOLVE3D */
