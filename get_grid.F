#include "cppdefs.h"
#ifndef ANA_GRID
                             ! Read grid information
      subroutine get_grid    ! from GRID NetCDF file.
 
      implicit none
      character*1 char1
      integer ierr, ncid, varid, lstr, lenstr, checkdims, nf_fread
# include "param.h"
# include "scalars.h"
# include "grid.h"
# include "ncvars.h"
# include "netcdf.inc"
!
! Open grid netCDF file for reading. Check that dimensions in that
! file are consistent with the model, then read all necessary
! variables.
!
      lstr=lenstr(grdname)
      ierr=nf_open(grdname(1:lstr), nf_nowrite, ncid)
      if (ierr .eq. nf_noerr) then
        ierr=checkdims (ncid, grdname(1:lstr), lstr, varid)
        if (ierr. ne. nf_noerr) goto 99
      else
        write(stdout,'(/3(1x,A)/)') 'GET_GRID ERROR: Cannot open',
     &                       'input NetCDF file:', grdname(1:lstr)
        goto 99                                           !--> ERROR
      endif
!
! Logical switch for spherical grid configuration:
!
      ierr=nf_inq_varid (ncid, 'spherical', varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_get_var1_text (ncid, varid, 1, char1)
        if (ierr .ne. nf_noerr) then
          write(stdout,2) 'spherical', grdname(1:lstr)
          goto 99                                         !--> ERROR
        endif
      else
        write(stdout,1) 'spherical', grdname(1:lstr)
        goto 99                                           !--> ERROR
      endif
 
      if (char1.eq.'t' .or. char1.eq.'T') then
# ifdef SPHERICAL
        MPI_master_only write(stdout,'(/1x,A/)')
     &                          'Spherical grid detected.'
# else
        write(stdout,'(/1x,A/12x,A/)')
     &      'GET_GRID - ERROR: Spherical grid detected, but',
     &               'model SPHERICAL CPP-switch is not set.'
        goto 99                                            !--> ERROR
# endif
      endif 
!
! Read in grid arrays.
!===== == ==== =======
! Bathymetry:
!
      ierr=nf_inq_varid (ncid, 'h', varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread(h (START_2D_ARRAY), ncid, varid, 0, r2dvar)
        if (ierr .ne. nf_noerr) then
          write(stdout,2) 'h', grdname(1:lstr)
          goto 99                                         !--> ERROR
        endif
      else
        write(stdout,1) 'h', grdname(1:lstr)
        goto 99                                           !--> ERROR
      endif
!
! Coriolis parameter.
!
      ierr=nf_inq_varid (ncid, 'f', varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (f(START_2D_ARRAY), ncid, varid, 0, r2dvar)
        if (ierr .ne. nf_noerr) then
          write(stdout,2) 'f', grdname(1:lstr)
          goto 99                                         !--> ERROR
        endif
      else
        write(stdout,1) 'f', grdname(1:lstr)
        goto 99                                           !--> ERROR
      endif
!
! Coordinate transfomation metrics (m,n) associated with the
! differential distances in XI and ETA.
!
      ierr=nf_inq_varid (ncid, 'pm', varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (pm(START_2D_ARRAY), ncid, varid, 0, r2dvar)
        if (ierr .ne. nf_noerr) then
          write(stdout,2) 'pm', grdname(1:lstr)
          goto 99                                         !--> ERROR
        endif
      else
        write(stdout,1) 'pm', grdname(1:lstr)
        goto 99                                           !--> ERROR
      endif
 
      ierr=nf_inq_varid (ncid, 'pn', varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (pn(START_2D_ARRAY), ncid, varid, 0, r2dvar)
        if (ierr .ne. nf_noerr) then
          write(stdout,2) 'pn', grdname(1:lstr)
          goto 99                                         !--> ERROR
        endif
      else
        write(stdout,1) 'pn', grdname(1:lstr)
        goto 99                                           !--> ERROR
      endif
!
! Coordinates (lon,lat [degrees]) or (x,y [meters]) at RHO-points.
!
# ifdef SPHERICAL
      ierr=nf_inq_varid (ncid, 'lon_rho', varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (lonr(START_2D_ARRAY), ncid, varid, 0, r2dvar)
        if (ierr .ne. nf_noerr) then
          write(stdout,2) 'lon_rho', grdname(1:lstr)
          goto 99                                         !--> ERROR
        endif
      else
        write(stdout,1) 'lon_rho', grdname(1:lstr)
        goto 99                                           !--> ERROR
      endif
 
      ierr=nf_inq_varid (ncid, 'lat_rho', varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (latr(START_2D_ARRAY), ncid, varid, 0, r2dvar)
        if (ierr .ne. nf_noerr) then
          write(stdout,2) 'lat_rho', grdname(1:lstr)
          goto 99                                         !--> ERROR
        endif
      else
        write(stdout,1) 'lat_rho', grdname(1:lstr)
        goto 99                                           !--> ERROR
      endif
# else
      ierr=nf_inq_varid (ncid, 'x_rho', varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (xr(START_2D_ARRAY), ncid, varid, 0, r2dvar)
        if (ierr .ne. nf_noerr) then
          write(stdout,2) 'x_rho', grdname(1:lstr)
          goto 99                                         !--> ERROR
        endif
      else
        write(stdout,1) 'x_rho', grdname(1:lstr)
        goto 99                                           !--> ERROR
      endif
 
      ierr=nf_inq_varid (ncid, 'y_rho', varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (yr(START_2D_ARRAY), ncid, varid, 0, r2dvar)
        if (ierr .ne. nf_noerr) then
          write(stdout,2) 'y_rho', grdname(1:lstr)
          goto 99                                           !--> ERROR
        endif
      else
        write(stdout,1) 'y_rho', grdname(1:lstr)
        goto 99                                           !--> ERROR
      endif
# endif
# ifdef CURVGRID
!
! Angle (radians) between XI-axis and EAST at RHO-points.
!
      ierr=nf_inq_varid (ncid, 'angle', varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread (angler(START_2D_ARRAY), ncid, varid, 0, r2dvar)
        if (ierr .ne. nf_noerr) then
          write(stdout,2) 'angle', grdname(1:lstr)
          goto 99                                         !--> ERROR
        endif
      else
        write(stdout,1) 'angle', grdname(1:lstr)
        goto 99                                           !--> ERROR
      endif
# endif
# ifdef MASKING
!
! Mask at RHO-points.
!
      ierr=nf_inq_varid (ncid, 'mask_rho', varid)
      if (ierr .eq. nf_noerr) then
        ierr=nf_fread(rmask(START_2D_ARRAY), ncid, varid, 0, r2dvar)
        if (ierr .ne. nf_noerr) then
          write(stdout,2) 'mask_rho', grdname(1:lstr)
          goto 99                                         !--> ERROR
        endif
      else
        write(stdout,1) 'mask_rho', grdname(1:lstr)
        goto 99                                           !--> ERROR
      endif
# endif
      ierr=nf_close(ncid)
      return                                      !--> NORMAL RETURN
 
   1  format(/1x, 'GET_GRID - unable to find grid variable:',
     &           1x, A / 12x, 'in grid netCDF file:', 1x, A/)
   2  format(/1x, 'GET_GRID - error while reading variable:',
     &           1x, A / 12x, 'in grid netCDF file:', 1x, A/)
  99  may_day_flag=2
      return                                              !--> ERROR
      end
#else
      subroutine get_grid_empty
      end
#endif /* !ANA_GRID */
 
