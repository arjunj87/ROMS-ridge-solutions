! This code segment computes horizontal fluxes for tracer variables.
! Basically it interpolates tracer values from their native locations
! on C grid to horizontal velocity points. Curently three options are
! supported: 4-point symmetric fourth-order method (default); 3-point
! upstream-biased parabolic interpolation (UPSTREAM); and 4-point
! scheme where arithmetic averaging of elementary differences is
! replaced by harmonic averaging (AKIMA), resulting in mid-point
! values bounded by two nearest values at native location, regardless
! of grid-scale roughness of the interpolated field, while still
! retaining asymptotic fourth-order behavior for smooth fields.
! This code is extracted into a special module ibecause it is used
! twice, in predictor and corrector substeps for tracer variables. 
! 

#ifdef ADV_ISONEUTRAL
# undef  UPSTREAM
#else
# define UPSTREAM
#endif

# undef AKIMA

c--# define CONST_TRACERS

# ifdef UPSTREAM
#  define curv WORK
# else
#  define grad WORK
# endif
# ifdef EW_PERIODIC
          imin=istr-1                  ! Determine extended index
          imax=iend+2                  ! range for computation of
# else
          if (WESTERN_EDGE) then       ! elementary differences:
            imin=istr                  ! it needs to be restricted 
          else                         ! because in the vicinity of  
            imin=istr-1                ! physical boundary the extra
          endif                        ! point may be not available,
          if (EASTERN_EDGE) then       ! and extrapolation of slope
            imax=iend+1                ! is used instead.
          else
            imax=iend+2
          endif
# endif
          do j=jstr,jend
            do i=imin,imax
              FX(i,j)=(t(i,j,k,nrhs,itrc)-t(i-1,j,k,nrhs,itrc))
# ifdef MASKING
     &                                               *umask(i,j)
# endif
            enddo
          enddo            !--> discard imin,imax
# ifndef EW_PERIODIC
          if (WESTERN_EDGE) then
            do j=jstr,jend
              FX(istr-1,j)=FX(istr,j)
            enddo
          endif
          if (EASTERN_EDGE) then
            do j=jstr,jend
              FX(iend+2,j)=FX(iend+1,j)
            enddo
          endif
# endif
          do j=jstr,jend
            do i=istr-1,iend+1
# if defined UPSTREAM
              curv(i,j)=FX(i+1,j)-FX(i,j)
# elif defined AKIMA
              cff=2.*FX(i+1,j)*FX(i,j)
              if (cff.gt.epsil) then
                grad(i,j)=cff/(FX(i+1,j)+FX(i,j))
              else
                grad(i,j)=0.
              endif
# else
              grad(i,j)=0.5*(FX(i+1,j)+FX(i,j))
# endif
            enddo
          enddo             !--> discard FX
          do j=jstr,jend
            do i=istr,iend+1
# ifdef UPSTREAM
              FX(i,j)=0.5*(t(i,j,k,nrhs,itrc)+t(i-1,j,k,nrhs,itrc))
     &                                                  *FlxU(i,j,k)
     &            -0.166666666666*( curv(i-1,j)*max(FlxU(i,j,k),0.)
     &                             +curv(i  ,j)*min(FlxU(i,j,k),0.))
# else
              FX(i,j)=0.5*( t(i,j,k,nrhs,itrc)+t(i-1,j,k,nrhs,itrc)
     &                      -0.333333333333*( grad(i,j)-grad(i-1,j))
     
     &                                                 )*FlxU(i,j,k)
# endif
            enddo           !--> discard curv,grad, keep FX
          enddo
 
# ifdef NS_PERIODIC
          jmin=jstr-1
          jmax=jend+2
# else
          if (SOUTHERN_EDGE) then
            jmin=jstr
          else
            jmin=jstr-1
          endif
          if (NORTHERN_EDGE) then
            jmax=jend+1
          else
            jmax=jend+2
          endif
# endif
          do j=jmin,jmax
            do i=istr,iend
              FE(i,j)=(t(i,j,k,nrhs,itrc)-t(i,j-1,k,nrhs,itrc))
# ifdef MASKING
     &                                               *vmask(i,j)
# endif
            enddo
          enddo         !--> discard jmin,jmax
# ifndef NS_PERIODIC
          if (SOUTHERN_EDGE) then
            do i=istr,iend
              FE(i,jstr-1)=FE(i,jstr)
            enddo
          endif
          if (NORTHERN_EDGE) then
            do i=istr,iend
              FE(i,jend+2)=FE(i,jend+1)
            enddo
          endif
# endif
          do j=jstr-1,jend+1
            do i=istr,iend
# if defined UPSTREAM
              curv(i,j)=FE(i,j+1)-FE(i,j)
# elif defined AKIMA
              cff=2.*FE(i,j+1)*FE(i,j)
              if (cff.gt.epsil) then
                grad(i,j)=cff/(FE(i,j+1)+FE(i,j))
              else
                grad(i,j)=0.
              endif
# else
              grad(i,j)=0.5*(FE(i,j+1)+FE(i,j))
# endif
            enddo
          enddo            !--> discard FE
 
          do j=jstr,jend+1
            do i=istr,iend
# ifdef UPSTREAM
              FE(i,j)=0.5*(t(i,j,k,nrhs,itrc)+t(i,j-1,k,nrhs,itrc))
     &                                                  *FlxV(i,j,k)
     &            -0.166666666666*( curv(i,j-1)*max(FlxV(i,j,k),0.)
     &                             +curv(i,j  )*min(FlxV(i,j,k),0.))
# else
              FE(i,j)=0.5*( t(i,j,k,nrhs,itrc)+t(i,j-1,k,nrhs,itrc)
     &                       -0.333333333333*(grad(i,j)-grad(i,j-1))
     &                                                 )*FlxV(i,j,k)
# endif
            enddo
          enddo            !--> discard curv,grad, keep FE
