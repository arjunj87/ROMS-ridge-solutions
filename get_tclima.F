#include "cppdefs.h"
#if defined SOLVE3D &&( defined TCLIMATOLOGY ||( defined TNUDGING &&\
              defined T_FRC_BRY )) && !defined ANA_TCLIMA

 
                                       ! Read climatology fields for
      subroutine get_tclima (ierr)     ! tracers from climatological
      implicit none                    ! file at appropriate time.
# define TCLIMA_DATA
# include "param.h"
# include "scalars.h"
# include "climat.h"
# include "ncvars.h"
# include "netcdf.inc"
      real cff
      integer i,itrc, lstr,lvar,lenstr, ierr, nf_fread
!
! Initialization: Check, whether climatological data file is already
!================ opened, and if not, open it.  Find and save netCDF
! IDs for relevant variables.  Determine whether there is cycling to
! reuse the input data and find cycling period "tclm_cycle", set
! initial cycling index "tclm_ncycle" and record index "tclm_rec",
! time index "itclm" and both time record bounds to large negative
! artificial values, so that it will trigger the logic in reading
! part below.
!
      ierr=nf_noerr
      lstr=lenstr(clm_file)
      if (iic.eq.ntstart) then
        if (ncidclm.eq.-1) then
          ierr=nf_open (clm_file(1:lstr), nf_nowrite, ncidclm)
          if (ierr .ne. nf_noerr) write(stdout,'(/1x,4A/)') 'ERROR ',
     &                  'in get_tclima: can not open netCDF file ''',
     &                                      clm_file(1:lstr),   '''.'
        endif
        do itrc=1,NT
          if (ierr.eq.nf_noerr) then
            if (itrc.eq.itemp) then
              ierr=nf_inq_varid (ncidclm,'tclm_time',tclm_tid(itrc))
              if (ierr .ne. nf_noerr) write(stdout,1) 'tclm_time',
     &                                         itrc, clm_file(1:lstr)
# ifdef SALINITY
            elseif (itrc.eq.isalt) then
              ierr=nf_inq_varid (ncidclm, 'sclm_time',tclm_tid(itrc))
              if (ierr .ne. nf_noerr) write(stdout,1) 'sclm_time',
     &                                         itrc, clm_file(1:lstr)
# endif
            endif
            if (ierr.eq.nf_noerr .and. tclm_tid(itrc).gt.0) then
              lvar=lenstr(vname(1,indxT+itrc-1))
              ierr=nf_inq_varid (ncidclm, vname(1,indxT+itrc-1)
     &                                  (1:lvar), tclm_id(itrc))
              if (ierr .eq. nf_noerr) then
                call set_cycle (ncidclm,tclm_tid(itrc), nttclm(itrc),
     &                          tclm_cycle(itrc),  tclm_ncycle(itrc),
     &                                         tclm_rec(itrc), ierr)
                ittclm(itrc)=1
                tclm_time(1,itrc)=-1.E+20
                tclm_time(2,itrc)=-1.E+20
              else
                write(stdout,1) vname(1,indxT+itrc-1)(1:lvar),
     &                                   itrc, clm_file(1:lstr)
              endif
            endif
          endif
        enddo
      endif
  1   format(/1x,'ERROR in get_tclima: cannot find variable ''',
     &    A, ''', itrc =',I3/7x, 'in netCDF file ''', A, '''.' /)
!
! Read data from the file:  Check if model time is bounded by past
!===== ==== ==== === =====  and future data times: if not, increment
! record and cycling indices, flip time index and read a new portion
! of data. Repeat until model time falls between the two data times.
!
      do itrc=1,NT
        if (tclm_tid(itrc).gt.0) then
          do while (tclm_time(ittclm(itrc),itrc).lt.time+0.5*dt .and.
     &                                              ierr.eq.nf_noerr)
            call advance_cycle (tclm_cycle(itrc),  nttclm(itrc),
     &                 tclm_ncycle(itrc), tclm_rec(itrc), ierr)
            if (ierr.eq.nf_noerr) then
              ierr=nf_get_var1_FTYPE(ncidclm, tclm_tid(itrc),
     &                                  tclm_rec(itrc), cff)
              if (ierr.eq.nf_noerr) then
                ittclm(itrc)=min(3-ittclm(itrc),nttclm(itrc))
                tclm_time(ittclm(itrc),itrc)=cff*day2sec
     &               + tclm_cycle(itrc)*tclm_ncycle(itrc)
 
                ierr=nf_fread (tclima(START_2D_ARRAY,1,ittclm(itrc),
     &                                itrc), ncidclm, tclm_id(itrc),
     &                                      tclm_rec(itrc), r3dvar)
                if (ierr.eq.nf_noerr) then
                  write(stdout,'(6x,A,I2,1x,A,F12.4,2(1x,A,I4))')
     &            'get_tclima -- read climatology for tracer', itrc,
     &            'tclm_time =', cff,  'rec =', tclm_rec(itrc) MYID
                  if (nttclm(itrc).eq.1) goto 3
                else
                  lvar=lenstr(vname(1,indxT+itrc-1))
                  write(stdout,4) vname(1,indxT+itrc-1)(1:lvar),
     &                                     itrc, tclm_rec(itrc)
                endif
              else
                write(stdout,4) 'Xclm_time', itrc, tclm_rec(itrc)
              endif
            else
              write(stdout,'(/1x,A,I3,1x,A,I3/7x,3A/7x,2(A,G12.4)/)')
     &       'ERROR in get_tclima: requested time record ',
     &        tclm_rec(itrc), 'exeeds the last record', nttclm(itrc),
     &       'available in netCDF file ''',  clm_file(1:lstr),  '''',
     &       'tdays = ', tdays,' but the last available tclm_time =',
     &                          tclm_time(ittclm(itrc),itrc)*sec2day
            endif
          enddo
        endif
  3    continue
      enddo
  4   format(/1x,'ERROR in get_tclima: cannot read variable ''',
     &                  A, ''', itrc =', I3, ' for record', I3/)
      return
      end
 
 
      subroutine set_tclima_tile (istr,iend,jstr,jend, ierr)
!
! Set-up tracer climatology for current tile.
!
      implicit none
      integer istr,iend,jstr,jend, ierr, itrc,i,j,k, it1,it2
      real cff, cff1, cff2
# define TCLIMA_DATA
# include "param.h"
# include "climat.h"
# include "scalars.h"
# include "ocean3d.h"
!
# include "compute_extended_bounds.h"
!
! Set coefficients for interpolation. Check that for the next time
! step [when time=time+dt] both weights will still be positive, and
! if not, set synchro_flag to signal that new data should be read
! from an appropriate netCDF input file (master thread only).
! After that either load time-invariant data, or interpolate in time
! or complain about error and signal to quit, if interpolation is
! needed, but not possible.
!
      do itrc=1,NT
        if (tclm_tid(itrc).gt.0) then
          it1=3-ittclm(itrc)
          it2=ittclm(itrc)
          cff1=tclm_time(it2,itrc)-time
          cff2=time-tclm_time(it1,itrc)

          if (cff1.lt.dt) synchro_flag=.true.
 
          if (nttclm(itrc).eq.1) then         ! Load time-invariant
            if (iic.eq.0) then                ! tracer climatology.
              do k=1,N
                do j=jstrR,jendR
                  do i=istrR,iendR
                    tclm(i,j,k,itrc)=tclima(i,j,k,ittclm(itrc),itrc)
                  enddo
                enddo
              enddo
            endif
          elseif (cff1.ge.0. .and. cff2.ge.0.) then
            cff=1./(cff1+cff2)                        ! Interpolate
            cff1=cff1*cff                             ! tracer
            cff2=cff2*cff                             ! climatology
            do k=1,N                                  ! climatology
              do j=jstrR,jendR                        ! in time.
                do i=istrR,iendR
                  tclm(i,j,k,itrc)=cff1*tclima(i,j,k,it1,itrc)
     &                            +cff2*tclima(i,j,k,it2,itrc)
                enddo
              enddo
            enddo
          elseif (ZEROTH_TILE) then
            write(stdout,'(/2A/3(1x,A,F16.10)/)') ' SET_TCLIMA_TILE',
     &         ' --  model time is outside bounds of ''tclm_time''.',
     &         'TCLM_TSTART=', tclm_time(it1,itrc)*sec2day, 'TDAYS=',
     &          tdays, 'TCLM_TEND=',     tclm_time(it2,itrc)*sec2day
            ierr=ierr+1
          endif
        endif   ! <-- tclm_tid(itrc).gt.0
      enddo   ! <-- itrc
      return
      end
#else
      subroutine get_tclima_empty
      end
#endif
 
