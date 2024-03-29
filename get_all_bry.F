#include "cppdefs.h"
#if !defined ANA_BRY && (defined T_FRC_BRY || defined M2_FRC_BRY || \
                         defined M3_FRC_BRY || defined Z_FRC_BRY )

      subroutine get_all_bry (ierr)     ! Read side boundary forcing 
      implicit none                     ! fields from boundary file
# include "param.h"
# include "scalars.h"
# include "ncvars.h"
# include "netcdf.inc"
# include "boundary.h"
      real cff
      integer lstr, lenstr, ierr, ierr_all, itrc
      integer nf_read_bry_EW, nf_read_bry_NS 
# ifdef OBC_WEST
      logical west_bry_active
# endif
# ifdef OBC_EAST
      logical east_bry_active
# endif
# ifdef OBC_SOUTH
      logical south_bry_active
# endif
# ifdef OBC_NORTH
      logical north_bry_active
# endif
!
! Determine which side(s) require boundary conditions and check
! whether it is needed at all.  Note that "itrc" is used here as
! a trigger variable: it it remains 0, then it means that this
! MPI subdomain does not have edge adjacent to an open boundary. 
!
      itrc=0 
# ifdef OBC_WEST
#  ifdef MPI
      west_bry_active=.not.west_inter
#  else
      west_bry_active=.true.
#  endif
      if (west_bry_active) itrc=1
# endif
# ifdef OBC_EAST
#  ifdef MPI
      east_bry_active=.not.east_inter
#  else
      east_bry_active=.true.
#  endif
      if (east_bry_active) itrc=1
# endif
# ifdef OBC_SOUTH
#  ifdef MPI
      south_bry_active=.not.south_inter
#  else
      south_bry_active=.true.
#  endif
      if (south_bry_active) itrc=1
# endif
# ifdef OBC_NORTH
#  ifdef MPI
      north_bry_active=.not.north_inter
#  else
      north_bry_active=.true.
#  endif
      if (north_bry_active) itrc=1
# endif
      if (itrc.eq.0) return  !--> no boundary forcing is needed.


!
! Initialization: Check, whether boundary forcing file is already
!---------------- opened, and if not, open it.  Find and save netCDF
! IDs for relevant variables.  Determine whether there is cycling to
! reuse the input data and find cycling period "bry_cycle", set
! initial cycling index "bry_ncycle" and record index "bry_rec",
! time index "itbry" and both time record bounds to large negative
! artificial values, so that it will trigger the logic in reading
! part below.
!
      ierr=nf_noerr
      lstr=lenstr(bry_file)
      if (iic.eq.ntstart) then
        if (bry_id .eq. -1) then
          ierr=nf_open (bry_file(1:lstr), nf_nowrite, bry_id)
          if (ierr.ne.nf_noerr) write(stdout,'(/1x,4A/)') 'ERROR ',
     &               'in get_all_bry: can not open netCDF file ''',
     &                                   bry_file(1:lstr),   '''.'
        endif
        if (ierr.eq.nf_noerr) then
          ierr_all=0 
          ierr=nf_inq_varid (bry_id, 'bry_time',  bry_time_id)
          if (ierr.ne.nf_noerr) write(stdout,1) 'bry_time'
          ierr_all=ierr_all+ierr

# ifdef OBC_WEST
          if (west_bry_active) then
#  ifdef Z_FRC_BRY
            ierr=nf_inq_varid (bry_id, 'zeta_west', zeta_west_id)
            if (ierr.ne.nf_noerr) write(stdout,1) 'zeta_west'
            ierr_all=ierr_all+ierr
#  endif
#  ifdef M2_FRC_BRY
            ierr=nf_inq_varid (bry_id, 'ubar_west', ubar_west_id)
            if (ierr.ne.nf_noerr) write(stdout,1) 'ubar_west'
            ierr_all=ierr_all+ierr
            ierr=nf_inq_varid (bry_id, 'vbar_west', vbar_west_id)
            if (ierr.ne.nf_noerr) write(stdout,1) 'vbar_west'
            ierr_all=ierr_all+ierr
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
            ierr=nf_inq_varid (bry_id,    'u_west', u_west_id)
            if (ierr.ne.nf_noerr)  write(stdout,1) 'u_west'
            ierr_all=ierr_all+ierr
            ierr=nf_inq_varid (bry_id,    'v_west', v_west_id)
            if (ierr.ne.nf_noerr)  write(stdout,1) 'v_west'
            ierr_all=ierr_all+ierr
#   endif
#   ifdef T_FRC_BRY
            ierr=nf_inq_varid (bry_id, 'temp_west', t_west_id(itemp))
            if (ierr.ne.nf_noerr) write(stdout,1) 'temp_west'
            ierr_all=ierr_all+ierr
#    ifdef SALINITY
            ierr=nf_inq_varid (bry_id, 'salt_west', t_west_id(isalt))
            if (ierr.ne.nf_noerr) write(stdout,1) 'salt_west'
            ierr_all=ierr_all+ierr
#    endif
#   endif
#  endif
          endif
# endif

# ifdef OBC_EAST
          if (east_bry_active) then 
#  ifdef Z_FRC_BRY
            ierr=nf_inq_varid (bry_id, 'zeta_east', zeta_east_id)
            if (ierr.ne.nf_noerr) write(stdout,1) 'zeta_east'
            ierr_all=ierr_all+ierr
#  endif
#  ifdef M2_FRC_BRY
            ierr=nf_inq_varid (bry_id, 'ubar_east', ubar_east_id)
            if (ierr.ne.nf_noerr) write(stdout,1) 'ubar_east'
            ierr_all=ierr_all+ierr
            ierr=nf_inq_varid (bry_id, 'vbar_east', vbar_east_id)
            if (ierr.ne.nf_noerr) write(stdout,1) 'vbar_east'
            ierr_all=ierr_all+ierr
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
            ierr=nf_inq_varid (bry_id,    'u_east',   u_east_id)
            if (ierr.ne.nf_noerr)  write(stdout,1)   'u_east'
            ierr_all=ierr_all+ierr
            ierr=nf_inq_varid (bry_id,    'v_east',   v_east_id)
            if (ierr.ne.nf_noerr)  write(stdout,1)   'v_east'
            ierr_all=ierr_all+ierr
#   endif
#   ifdef T_FRC_BRY
            ierr=nf_inq_varid (bry_id, 'temp_east', t_east_id(itemp))
            if (ierr.ne.nf_noerr) write(stdout,1) 'temp_east'
            ierr_all=ierr_all+ierr
#    ifdef SALINITY
            ierr=nf_inq_varid (bry_id, 'salt_east', t_east_id(isalt))
            if (ierr.ne.nf_noerr) write(stdout,1) 'salt_east'
            ierr_all=ierr_all+ierr
#    endif
#   endif
#  endif
          endif
# endif

# ifdef OBC_SOUTH
          if (south_bry_active) then
#  ifdef Z_FRC_BRY
            ierr=nf_inq_varid (bry_id, 'zeta_south', zeta_south_id)
            if (ierr.ne.nf_noerr) write(stdout,1) 'zeta_south'
            ierr_all=ierr_all+ierr
#  endif
#  ifdef M2_FRC_BRY
            ierr=nf_inq_varid (bry_id, 'ubar_south', ubar_south_id)
            if (ierr.ne.nf_noerr) write(stdout,1) 'ubar_south'
            ierr_all=ierr_all+ierr
            ierr=nf_inq_varid (bry_id, 'vbar_south', vbar_south_id)
            if (ierr.ne.nf_noerr) write(stdout,1) 'vbar_south'
            ierr_all=ierr_all+ierr
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
            ierr=nf_inq_varid (bry_id,    'u_south',  u_south_id)
            if (ierr.ne.nf_noerr) write(stdout,1)    'u_south'
            ierr_all=ierr_all+ierr
            ierr=nf_inq_varid (bry_id,    'v_south',  v_south_id)
            if (ierr.ne.nf_noerr) write(stdout,1)    'v_south'
            ierr_all=ierr_all+ierr
#   endif
#   ifdef T_FRC_BRY
            ierr=nf_inq_varid (bry_id,'temp_south',t_south_id(itemp))
            if (ierr.ne.nf_noerr) write(stdout,1) 'temp_south'
            ierr_all=ierr_all+ierr
#    ifdef SALINITY
            ierr=nf_inq_varid (bry_id,'salt_south',t_south_id(isalt))
            if (ierr.ne.nf_noerr) write(stdout,1) 'salt_south'
            ierr_all=ierr_all+ierr
#    endif
#   endif
#  endif
          endif
# endif

# ifdef OBC_NORTH
          if (north_bry_active) then
#  ifdef Z_FRC_BRY
            ierr=nf_inq_varid (bry_id, 'zeta_north', zeta_north_id)
            if (ierr.ne.nf_noerr) write(stdout,1) 'zeta_north'
            ierr_all=ierr_all+ierr
#  endif
#  ifdef M2_FRC_BRY
            ierr=nf_inq_varid (bry_id, 'ubar_north', ubar_north_id)
            if (ierr.ne.nf_noerr) write(stdout,1) 'ubar_north'
            ierr_all=ierr_all+ierr
            ierr=nf_inq_varid (bry_id, 'vbar_north', vbar_north_id)
            if (ierr.ne.nf_noerr) write(stdout,1) 'vbar_north'
            ierr_all=ierr_all+ierr
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
            ierr=nf_inq_varid (bry_id,    'u_north',  u_north_id)
            if (ierr.ne.nf_noerr) write(stdout,1)    'u_north'
            ierr_all=ierr_all+ierr
            ierr=nf_inq_varid (bry_id,    'v_north',  v_north_id)
            if (ierr.ne.nf_noerr)   write(stdout,1)  'v_north'
            ierr_all=ierr_all+ierr
#   endif
#   ifdef T_FRC_BRY
            ierr=nf_inq_varid (bry_id,'temp_north',t_north_id(itemp))
            if (ierr.ne.nf_noerr) write(stdout,1) 'temp_north'
            ierr_all=ierr_all+ierr
#    ifdef SALINITY
            ierr=nf_inq_varid (bry_id,'salt_north',t_north_id(isalt))
            if (ierr.ne.nf_noerr) write(stdout,1) 'salt_north'
            ierr_all=ierr_all+ierr
#    endif
#   endif
#  endif
          endif
# endif
          ierr=ierr_all
          if (ierr.eq.nf_noerr) then
            call set_cycle (bry_id, bry_time_id, ntbry,
     &           bry_cycle, bry_ncycle, bry_rec, ierr)
            itbry=1
            bry_time(1)=-1.E+20
            bry_time(2)=-1.E+20
          else
            write(stdout,'(8x,4A)') '### ERROR(s) occur while ',
     &                   'examining content of netCDF file ''',
     &                    bry_file(1:lstr), '''.'
          endif
        endif
      endif
  1   format(' ERROR in get_all_bry: cannot find variable ''',A,'''')
!
! Read data from the file:  Check if model time is bounded by past
!----- ---- ---- --- -----  and future data times: if not, increment
! record and cycling indices, flip time index and read a new portion
! of data. Repeat until model time falls between the two data times.
!
      do while (bry_time(itbry).lt.time+0.5*dt .and.ierr.eq.nf_noerr)
        call advance_cycle (bry_cycle,ntbry,bry_ncycle,bry_rec,ierr)
        if (ierr.eq.nf_noerr) then
          ierr=nf_get_var1_FTYPE (bry_id, bry_time_id, bry_rec, cff)
          if (ierr.eq.nf_noerr) then
            itbry=min(3-itbry,ntbry)
            bry_time(itbry)=cff*day2sec + bry_cycle*bry_ncycle
            ierr_all=0
#ifdef OBC_WEST
            if (west_bry_active) then
#  ifdef Z_FRC_BRY
              ierr=nf_read_bry_EW (zeta_west_dt(0,itbry), bry_id,
     &                             zeta_west_id, bry_rec, r2dvar)
              if (ierr.ne.nf_noerr) write(stdout,2) 'zeta_west'
              ierr_all=ierr_all+ierr
#  endif
#  ifdef M2_FRC_BRY
              ierr=nf_read_bry_EW (ubar_west_dt(0,itbry), bry_id,
     &                           ubar_west_id, bry_rec, u2dvar)
              if (ierr.ne.nf_noerr) write(stdout,2) 'ubar_west'
              ierr_all=ierr_all+ierr
              ierr=nf_read_bry_EW (vbar_west_dt(0,itbry), bry_id,
     &                           vbar_west_id, bry_rec, v2dvar)
              if (ierr.ne.nf_noerr) write(stdout,2) 'vbar_west'
              ierr_all=ierr_all+ierr
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
              ierr=nf_read_bry_EW (u_west_dt(0,1,itbry), bry_id,
     &                           u_west_id,  bry_rec,  u3dvar)
              if (ierr.ne.nf_noerr)  write(stdout,2)  'u_west'
              ierr_all=ierr_all+ierr
              ierr=nf_read_bry_EW (v_west_dt(0,1,itbry), bry_id,
     &                           v_west_id,  bry_rec,  v3dvar)
              if (ierr.ne.nf_noerr)  write(stdout,2)  'v_west'
              ierr_all=ierr_all+ierr
#   endif
#   ifdef T_FRC_BRY
              do itrc=1,NT
                ierr=nf_read_bry_EW (t_west_dt(0,1,itbry,itrc),
     &                 bry_id, t_west_id(itrc),  bry_rec, r3dvar)
                if (ierr.ne.nf_noerr)  write(stdout,2)  'trc_west'
                ierr_all=ierr_all+ierr
              enddo
#   endif
#  endif
            endif
# endif

# ifdef OBC_EAST
            if (east_bry_active) then
#  ifdef Z_FRC_BRY
              ierr=nf_read_bry_EW (zeta_east_dt(0,itbry), bry_id,
     &                             zeta_east_id, bry_rec, r2dvar)
              if (ierr.ne.nf_noerr) write(stdout,2) 'zeta_east'
              ierr_all=ierr_all+ierr
#  endif
#  ifdef M2_FRC_BRY
              ierr=nf_read_bry_EW (ubar_east_dt(0,itbry), bry_id,
     &                             ubar_east_id, bry_rec, u2dvar)
              if (ierr.ne.nf_noerr) write(stdout,2) 'ubar_east'
              ierr_all=ierr_all+ierr
              ierr=nf_read_bry_EW (vbar_east_dt(0,itbry), bry_id,
     &                             vbar_east_id, bry_rec, v2dvar)
              if (ierr.ne.nf_noerr) write(stdout,2) 'vbar_east'
              ierr_all=ierr_all+ierr
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
              ierr=nf_read_bry_EW (u_east_dt(0,1,itbry), bry_id,
     &                             u_east_id,  bry_rec,  u3dvar)
              if (ierr.ne.nf_noerr)  write(stdout,2)  'u_east'
              ierr_all=ierr_all+ierr
              ierr=nf_read_bry_EW (v_east_dt(0,1,itbry), bry_id, 
     &                             v_east_id,  bry_rec,  v3dvar)
              if (ierr.ne.nf_noerr)  write(stdout,2)  'v_east'
              ierr_all=ierr_all+ierr
#   endif
#   ifdef T_FRC_BRY
              do itrc=1,NT
                ierr=nf_read_bry_EW (t_east_dt(0,1,itbry,itrc),
     &                bry_id, t_east_id(itrc),  bry_rec, r3dvar)
                if (ierr.ne.nf_noerr)  write(stdout,2)  'trc_east'
                ierr_all=ierr_all+ierr
              enddo
#   endif
#  endif
            endif
# endif

# ifdef OBC_SOUTH
            if (south_bry_active) then
#  ifdef Z_FRC_BRY
              ierr=nf_read_bry_NS (zeta_south_dt(0,itbry), bry_id,
     &                            zeta_south_id, bry_rec,  r2dvar)
              if (ierr.ne.nf_noerr) write(stdout,2) 'zeta_south'
              ierr_all=ierr_all+ierr
#  endif
#  ifdef M2_FRC_BRY
              ierr=nf_read_bry_NS (ubar_south_dt(0,itbry), bry_id,
     &                             ubar_south_id, bry_rec, u2dvar)
              if (ierr.ne.nf_noerr) write(stdout,2) 'ubar_south'
              ierr_all=ierr_all+ierr
              ierr=nf_read_bry_NS (vbar_south_dt(0,itbry), bry_id,
     &                             vbar_south_id, bry_rec, v2dvar)
              if (ierr.ne.nf_noerr) write(stdout,2) 'vbar_south'
              ierr_all=ierr_all+ierr
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
              ierr=nf_read_bry_NS (u_south_dt(0,1,itbry), bry_id,
     &                             u_south_id,  bry_rec,  u3dvar)
              if (ierr.ne.nf_noerr)  write(stdout,2)  'u_south'
              ierr_all=ierr_all+ierr

              ierr=nf_read_bry_NS (v_south_dt(0,1,itbry), bry_id,
     &                             v_south_id,  bry_rec,  v3dvar)
              if (ierr.ne.nf_noerr)  write(stdout,2)  'v_south'
              ierr_all=ierr_all+ierr
#   endif
#   ifdef T_FRC_BRY
              do itrc=1,NT
                ierr=nf_read_bry_NS (t_south_dt(0,1,itbry,itrc),
     &                bry_id, t_south_id(itrc),  bry_rec, r3dvar)
                if (ierr.ne.nf_noerr) write(stdout,2) 'trc_south'
                ierr_all=ierr_all+ierr
              enddo
#   endif
#  endif
            endif
# endif

# ifdef OBC_NORTH
            if (north_bry_active) then
#  ifdef Z_FRC_BRY
              ierr=nf_read_bry_NS (zeta_north_dt(0,itbry), bry_id,
     &                             zeta_north_id, bry_rec, r2dvar)
              if (ierr.ne.nf_noerr) write(stdout,2) 'zeta_north'
              ierr_all=ierr_all+ierr
#  endif
#  ifdef M2_FRC_BRY
              ierr=nf_read_bry_NS (ubar_north_dt(0,itbry), bry_id,
     &                             ubar_north_id, bry_rec, u2dvar)
              if (ierr.ne.nf_noerr) write(stdout,2) 'ubar_north'
              ierr_all=ierr_all+ierr
              ierr=nf_read_bry_NS (vbar_north_dt(0,itbry), bry_id, 
     &                             vbar_north_id, bry_rec, v2dvar)
              if (ierr.ne.nf_noerr) write(stdout,2) 'vbar_north'
              ierr_all=ierr_all+ierr
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
              ierr=nf_read_bry_NS (u_north_dt(0,1,itbry), bry_id,
     &                             u_north_id,  bry_rec,  u3dvar)
              if (ierr.ne.nf_noerr)  write(stdout,2)  'u_north'
              ierr_all=ierr_all+ierr
              ierr=nf_read_bry_NS (v_north_dt(0,1,itbry), bry_id,
     &                             v_north_id,  bry_rec,  v3dvar)
              if (ierr.ne.nf_noerr)  write(stdout,2)  'v_north'
              ierr_all=ierr_all+ierr
#   endif
#   ifdef T_FRC_BRY
              do itrc=1,NT
                ierr=nf_read_bry_NS (t_north_dt(0,1,itbry,itrc),
     &                 bry_id, t_north_id(itrc),  bry_rec, r3dvar)
                if (ierr.ne.nf_noerr) write(stdout,2) 'temp_north'
                ierr_all=ierr_all+ierr
              enddo
#   endif
#  endif
            endif
# endif
            ierr=ierr_all
            if (ierr.eq.0) then
              write(stdout,'(6x,A,2x,A,F12.4,2(1x,A,I4))')
     &            'get_all_bry -- read boundary for all fields',
     &            'bry_time =', cff, 'rec =', bry_rec MYID
            else
              write(stdout,'(1x,2A,I4/8x,3A)')
     &           '### ERROR(s) occur while trying to read record ',
     &             bry_rec,  'in file ''',bry_file(1:lstr),'''.'
            endif 
          else
            write(stdout,2) 'bry_time'
          endif
        else
           write(stdout,'(/1x,A,I4,1x,A,I4/7x,4A/7x,2(A,G12.4)/)')
     &    'ERROR in get_bry_all: requested time record ',  bry_rec,
     &    'exeeds the last record', ntbry,  'available in netCDF ',
     &    'file ''',   bry_file(1:lstr),  '''',  'tdays = ', tdays,
     &    '  but the last available  bry_time =',
     &                                   bry_time(itbry)*sec2day
        endif
      enddo
  2   format(/1x,'### ERROR in get_bry_all: Cannot read variable ''',
     &                                                     A, '''.')
      return
      end 




      subroutine set_all_bry_tile (istr,iend,jstr,jend, ierr)
      implicit none
      integer istr,iend,jstr,jend, ierr, imin,imax,jmin,jmax,
     &                                   i,j,k, it1,it2, itrc
      real cff, cff1,cff2
# include "param.h"
# include "scalars.h"
# include  "boundary.h"
                                        ! NOTE: Time interpolation
      if (WESTERN_EDGE) then            ! of BRY forcing variables is 
        imin=istr-1                     ! done over non-standard 
# ifdef MPI
      elseif (WEST_INTER) then          ! ranges of spatial indices
        imin=istr-1                     ! which is designed to avoid
# endif
      else                              ! the need of MPI-excange of
        imin=istr                       ! ghost points of BRY arrays 
      endif                             ! and, at the same time, to
                                        ! have only one row of ghost
      if (EASTERN_EDGE) then            ! points (hence, standard
        imax=iend+1                     ! "extended" bounds) are not
# ifdef MPI
      elseif (EAST_INTER) then          ! applicable here.
        imax=iend+1 
# endif
      else
        imax=iend
      endif

      if (SOUTHERN_EDGE) then
        jmin=jstr-1
# ifdef MPI
      elseif (SOUTH_INTER) then
        jmin=jstr-1
# endif
      else
        jmin=jstr
      endif

      if (NORTHERN_EDGE) then
        jmax=jend+1
# ifdef MPI
      elseif (NORTH_INTER) then
        jmax=jend+1
# endif
      else
        jmax=jend
      endif


      it1=0
# ifdef OBC_WEST
      if (WESTERN_EDGE) it1=1 
# endif
# ifdef OBC_EAST
      if (EASTERN_EDGE) it1=1
# endif
# ifdef OBC_SOUTH
      if (SOUTHERN_EDGE) it1=1
# endif
# ifdef OBC_NORTH
      if (NORTHERN_EDGE) it1=1 
# endif
      if (it1.eq.0) return



      it1=3-itbry
      it2=itbry
      cff1=bry_time(it2)-time
      cff2=time-bry_time(it1)

      if (cff1.lt.dt) synchro_flag=.true.

      if (cff1.ge.0. .and. cff2.ge.0.) then
        cff=1./(cff1+cff2)                      ! interpolate
        cff1=cff1*cff                           ! boundary values 
        cff2=cff2*cff                           ! in time 

# ifdef OBC_WEST
        if (WESTERN_EDGE) then
#  ifdef Z_FRC_BRY
          do j=jmin,jmax
            zeta_west(j)=cff1*zeta_west_dt(j,it1)
     &                  +cff2*zeta_west_dt(j,it2)
          enddo
#  endif
#  ifdef M2_FRC_BRY
          do j=jmin,jmax
            ubar_west(j)=cff1*ubar_west_dt(j,it1)
     &                  +cff2*ubar_west_dt(j,it2)
            vbar_west(j)=cff1*vbar_west_dt(j,it1)
     &                  +cff2*vbar_west_dt(j,it2)
          enddo
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
          do k=1,N
            do j=jmin,jmax
              u_west(j,k)=cff1*u_west_dt(j,k,it1)
     &                   +cff2*u_west_dt(j,k,it2)
              v_west(j,k)=cff1*v_west_dt(j,k,it1)
     &                   +cff2*v_west_dt(j,k,it2)
            enddo
          enddo
#   endif
#   ifdef T_FRC_BRY
          do itrc=1,NT
            do k=1,N
              do j=jmin,jmax
                t_west(j,k,itrc)=cff1*t_west_dt(j,k,it1,itrc)
     &                          +cff2*t_west_dt(j,k,it2,itrc)
              enddo
            enddo
          enddo
#   endif
#  endif
        endif
# endif

# ifdef OBC_EAST
        if (EASTERN_EDGE) then
#  ifdef Z_FRC_BRY
          do j=jmin,jmax
            zeta_east(j)=cff1*zeta_east_dt(j,it1)
     &                  +cff2*zeta_east_dt(j,it2)
          enddo
#  endif
#  ifdef M2_FRC_BRY
          do j=jmin,jmax
            ubar_east(j)=cff1*ubar_east_dt(j,it1)
     &                  +cff2*ubar_east_dt(j,it2)
            vbar_east(j)=cff1*vbar_east_dt(j,it1)
     &                  +cff2*vbar_east_dt(j,it2)
          enddo
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
          do k=1,N
            do j=jmin,jmax
              u_east(j,k)=cff1*u_east_dt(j,k,it1)
     &                   +cff2*u_east_dt(j,k,it2)
              v_east(j,k)=cff1*v_east_dt(j,k,it1)
     &                   +cff2*v_east_dt(j,k,it2)
            enddo
          enddo
#   endif
#   ifdef T_FRC_BRY
          do itrc=1,NT
            do k=1,N
              do j=jmin,jmax
                t_east(j,k,itrc)=cff1*t_east_dt(j,k,it1,itrc)
     &                          +cff2*t_east_dt(j,k,it2,itrc)
              enddo
            enddo
          enddo
#   endif
#  endif
        endif
# endif

# ifdef OBC_SOUTH
        if (SOUTHERN_EDGE) then
#  ifdef Z_FRC_BRY
          do i=imin,imax
            zeta_south(i)=cff1*zeta_south_dt(i,it1)
     &                   +cff2*zeta_south_dt(i,it2)
          enddo
#  endif
#  ifdef M2_FRC_BRY
          do i=imin,imax
            ubar_south(i)=cff1*ubar_south_dt(i,it1)
     &                   +cff2*ubar_south_dt(i,it2)
            vbar_south(i)=cff1*vbar_south_dt(i,it1)
     &                   +cff2*vbar_south_dt(i,it2)
          enddo
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
          do k=1,N
            do i=imin,imax
              u_south(i,k)=cff1*u_south_dt(i,k,it1)
     &                    +cff2*u_south_dt(i,k,it2)
              v_south(i,k)=cff1*v_south_dt(i,k,it1)
     &                    +cff2*v_south_dt(i,k,it2)
            enddo
          enddo
#   endif
#   ifdef T_FRC_BRY
          do itrc=1,NT
            do k=1,N
              do i=imin,imax
                t_south(i,k,itrc)=cff1*t_south_dt(i,k,it1,itrc)
     &                           +cff2*t_south_dt(i,k,it2,itrc)
              enddo
            enddo
          enddo
#   endif
#  endif
        endif
# endif

# ifdef OBC_NORTH
        if (NORTHERN_EDGE) then
#  ifdef Z_FRC_BRY
          do i=imin,imax
            zeta_north(i)=cff1*zeta_north_dt(i,it1)
     &                   +cff2*zeta_north_dt(i,it2)
          enddo
#  endif
#  ifdef M2_FRC_BRY
          do i=imin,imax
            ubar_north(i)=cff1*ubar_north_dt(i,it1)
     &                   +cff2*ubar_north_dt(i,it2)
            vbar_north(i)=cff1*vbar_north_dt(i,it1)
     &                   +cff2*vbar_north_dt(i,it2)
          enddo
#  endif
#  ifdef SOLVE3D
#   ifdef M3_FRC_BRY
          do k=1,N
            do i=imin,imax
              u_north(i,k)=cff1*u_north_dt(i,k,it1)
     &                    +cff2*u_north_dt(i,k,it2)
              v_north(i,k)=cff1*v_north_dt(i,k,it1)
     &                    +cff2*v_north_dt(i,k,it2)
            enddo
          enddo
#   endif
#   ifdef T_FRC_BRY
          do itrc=1,NT
            do k=1,N
              do i=imin,imax
                t_north(i,k,itrc)=cff1*t_north_dt(i,k,it1,itrc)
     &                           +cff2*t_north_dt(i,k,it2,itrc)
              enddo
            enddo
          enddo
#   endif
#  endif
        endif
# endif
      else
        write(stdout,'(/2(1x,A)/3(1x,A,F16.10)/)')
     &   'set_all_bry_tile - current model time is out of bounds of',
     &   '''bry_time''.',     'BRY_TSTART=',   bry_time(it1)*sec2day,
     &   'TDAYS=',  tdays,    'BRY_TEND=',     bry_time(it2)*sec2day
        ierr=ierr+1
      endif
      return
      end
#else
      subroutine get_bry_empty
      end
#endif

