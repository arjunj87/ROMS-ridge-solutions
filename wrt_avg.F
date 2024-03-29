#include "cppdefs.h"
#ifdef AVERAGES
                                     ! Write time-averaged
      subroutine wrt_avg             ! fields into averages
                                     ! netCDF file.
      implicit none
      integer ierr, record, lstr, lvar, lenstr
     &  , start(2), count(2), ibuff(4), nf_fwrite
# ifdef SOLVE3D
     &            , itrc
# endif
# include "param.h"
# include "scalars.h"
# include "averages.h"
# include "ncvars.h"
# include "netcdf.inc"
# if defined MPI && !defined PARALLEL_FILES
#  include "mpif.h"
      integer status(MPI_STATUS_SIZE), blank
      if (mynode.gt.0) then
        call MPI_Recv (blank, 1, MPI_INTEGER, mynode-1,
     &                 3, ocean_grid_comm, status, ierr)
      endif
# endif
!
! Create/open averages file; write grid arrays, if so needed,
! check for errors, then set record within the file "record".
!
      call def_avg (ncidavg, nrecavg, ierr)
      lstr=lenstr(avgname)
      if (ierr .ne. nf_noerr) goto 99       !!! WARNING: Here it is
      nrecavg=max(nrecavg,1)                !!! assumed that global
      if (nrpfavg.eq.0) then                !!! averages record index
        record=nrecavg                      !!! "nrecavg" is already
      else                                  !!! advanced by main.
        record=1+mod(nrecavg-1, nrpfavg)
      endif
!
! Write out time-averaged variables:
! ----- --- ------------- ----------
!
! Time step and record indices.
!
      ibuff(1)=iic
      ibuff(2)=nrecrst
      ibuff(3)=nrechis
      ibuff(4)=nrecavg
 
      start(1)=1
      start(2)=record
      count(1)=4
      count(2)=1
      ierr=nf_put_vara_int (ncidavg, avgTstep, start, count, ibuff)
      if (ierr .ne. nf_noerr) then
        write(stdout,1) 'time_step', record,ierr MYID
        goto 99                                           !--> ERROR
      endif
!
! Averaged time
!
      ierr=nf_put_var1_FTYPE (ncidavg, avgTime, record, time_avg)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxTime))
        write(stdout,1) vname(1,indxTime)(1:lvar), record, ierr
     &                  MYID
        goto 99                                           !--> ERROR
      endif
!
! Barotropic mode variables: free-surface and 2D momentum
! components in XI-,ETA-directions.
!
      if (wrtavg(indxZ)) then
        ierr=nf_fwrite (zeta_avg(START_2D_ARRAY), ncidavg, avgZ,
     &                                            record, r2dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxZ))
          write(stdout,1) vname(1,indxTime)(1:lvar), record, ierr
     &                    MYID
          goto 99                                         !--> ERROR
        endif
      endif
 
      if (wrtavg(indxUb)) then
        ierr=nf_fwrite (ubar_avg(START_2D_ARRAY), ncidavg, avgUb,
     &                                            record, u2dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxUb))
          write(stdout,1) vname(1,indxUb)(1:lvar), record, ierr
     &                    MYID
          goto 99                                         !--> ERROR
        endif
      endif
 
      if (wrtavg(indxVb)) then
        ierr=nf_fwrite (vbar_avg(START_2D_ARRAY), ncidavg, avgVb,
     &                                            record, v2dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxVb))
          write(stdout,1) vname(1,indxVb)(1:lvar), record, ierr
     &                    MYID
          goto 99                                         !--> ERROR
        endif
      endif
# ifdef SOLVE3D
!
! 3D momentum components in XI- and ETA-directions.
!
      if (wrtavg(indxU)) then
        ierr=nf_fwrite (u_avg(START_2D_ARRAY,1), ncidavg, avgU,
     &                                           record, u3dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxU))
          write(stdout,1) vname(1,indxU)(1:lvar), record, ierr
     &                    MYID
          goto 99                                         !--> ERROR
        endif
      endif
 
      if (wrtavg(indxV)) then
        ierr=nf_fwrite (v_avg(START_2D_ARRAY,1), ncidavg, avgV,
     &                                           record, v3dvar)
        if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxV))
          write(stdout,1) vname(1,indxV)(1:lvar), record, ierr
     &                    MYID
          goto 99                                         !--> ERROR
        endif
      endif
!
! Tracer variables.
!
      do itrc=1,NT
        if (wrtavg(indxT+itrc-1)) then
          ierr=nf_fwrite (t_avg(START_2D_ARRAY,1,itrc), ncidavg,
     &                              avgT(itrc), record, r3dvar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxT+itrc-1))
            write(stdout,1) vname(1,indxT+itrc-1)(1:lvar),
     &                      record, ierr MYID
            goto 99                                       !--> ERROR
          endif
        endif
      enddo
!
! Density anomaly.
!
      if (wrtavg(indxR)) then
        ierr=nf_fwrite (rho_avg(START_2D_ARRAY,1), ncidavg, avgR,
     &                                             record, r3dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxR))
          write(stdout,1) vname(1,indxR)(1:lvar), record, ierr
     &                    MYID
          goto 99                                         !--> ERROR
        endif
      endif
!
!  Write out S-coordinate omega vertical velocity (m/s).
!
      if (wrtavg(indxO)) then
        ierr=nf_fwrite (w_avg, ncidavg, avgO, record, w3dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxO))
          write(stdout,1) vname(1,indxO)(1:lvar), record, ierr
     &                    MYID
          goto 99                                         !--> ERROR
        endif
      endif
!
!  Write out  vertical velocity (m/s).
!
      if (wrtavg(indxW)) then
        ierr=nf_fwrite (wvlc_avg, ncidavg, avgW, record, r3dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxW))
          write(stdout,1) vname(1,indxW)(1:lvar), record, ierr
     &                    MYID
          goto 99                                         !--> ERROR
        endif
      endif

!
!  Write out vertical diffusion coefficient for momentum.
!
      if (wrtavg(indxAkv)) then
        ierr=nf_fwrite (akv_avg(START_2D_ARRAY,0), ncidavg,
     &                                avgAkv, record, w3dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxAkv))
          write(stdout,1) vname(1,indxAkv)(1:lvar), record, ierr,
     &                    nf_strerror(ierr) MYID
          goto 99                                         !--> ERROR
        endif
      endif
!
!  Write out vertical diffusion coefficient for potential temperature.
!
      if (wrtavg(indxAkt)) then
        ierr=nf_fwrite (akt_avg(START_2D_ARRAY,0), ncidavg,
     &                                avgAkt, record, w3dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxAkt))
          write(stdout,1) vname(1,indxAkt)(1:lvar), record, ierr,
     &                    nf_strerror(ierr) MYID
          goto 99                                         !--> ERROR
        endif
      endif
# ifdef SALINITY
!
!  Write out vertical diffusion coefficient for salinity.
!
      if (wrtavg(indxAks)) then
        ierr=nf_fwrite (aks_avg(START_2D_ARRAY,0), ncidavg,
     &                                avgAks, record, w3dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxAks))
          write(stdout,1) vname(1,indxAks)(1:lvar), record, ierr,
     &                    nf_strerror(ierr) MYID
          goto 99                                         !--> ERROR
        endif
      endif
# endif
# ifdef LMD_KPP
!
!  Write out depth of planetary boundary layer (m).
!
      if (wrtavg(indxHbl)) then
        ierr=nf_fwrite (hbl_avg(START_2D_ARRAY), ncidavg,
     &                              avgHbl, record, r2dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxHbl))
          write(stdout,1) vname(1,indxHbl)(1:lvar), record, ierr,
     &                    nf_strerror(ierr) MYID
          goto 99                                         !--> ERROR
        endif
# ifdef LMD_BKPP        
        ierr=nf_fwrite (hbbl_avg(START_2D_ARRAY), ncidavg,
     &                              avgHbbl, record, r2dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxHbbl))
          write(stdout,1) vname(1,indxHbbl)(1:lvar), record, ierr,
     &                    nf_strerror(ierr) MYID
          goto 99                                         !--> ERROR
        endif        
# endif                
      endif
# endif



#if defined ADV_ISONEUTRAL & defined WRITE_SXSY
        ierr=nf_fwrite (Sx_avg(START_2D_ARRAY,1), ncidavg,
     &                                avgSX, record, u3dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxSX))
          write(stdout,1) vname(1,indxSX)(1:lvar), record, ierr,
     &                    nf_strerror(ierr) MYID
          goto 99                                         !--> ERROR
        endif

        ierr=nf_fwrite (Sy_avg(START_2D_ARRAY,1), ncidavg,
     &                                avgSY, record, v3dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxSY))
          write(stdout,1) vname(1,indxSY)(1:lvar), record, ierr,
     &                    nf_strerror(ierr) MYID
          goto 99                                         !--> ERROR
        endif
#endif  /* WRITE_SXSY */


#  ifdef BULK_FLUX
!
!  Write out surface heat flux component (degC/s)
!
      if (wrtavg(indxSustr_blk)) then
       ierr=nf_fwrite(sustr_blk_avg(START_2D_ARRAY), ncidavg, 
     &                                       avgSustr_blk, record, 
# ifdef WND_AT_RHO_POINTS
     &                                                     r2dvar)
# else
     &                                                     u2dvar)
# endif
       if (ierr .ne. nf_noerr) then
         lvar=lenstr(vname(1,indxSustr_blk))
         write(stdout,1) vname(1,indxSustr_blk)(1:lvar), record, ierr
     &                  MYID
         goto 99                                         !--> ERROR
       endif
      endif

      if (wrtavg(indxSvstr_blk)) then
        ierr=nf_fwrite(svstr_blk_avg(START_2D_ARRAY), ncidavg, 
     &                                        avgSvstr_blk, record,
# ifdef WND_AT_RHO_POINTS
     &                                                     r2dvar)
# else
     &                                                     v2dvar)
# endif

       if (ierr .ne. nf_noerr) then
           lvar=lenstr(vname(1,indxSvstr_blk))
           write(stdout,1) vname(1,indxSvstr_blk)(1:lvar), record, ierr
     &                  MYID
            goto 99                                         !--> ERROR
       endif
      endif

      if (wrtavg(indxShflx_rsw)) then
       ierr=nf_fwrite(shflx_rsw_avg(START_2D_ARRAY), ncidavg, 
     &                                       avgShflx_rsw, record, 
     &                                                     r2dvar)
       if (ierr .ne. nf_noerr) then
         lvar=lenstr(vname(1,indxShflx_rsw))
         write(stdout,1) vname(1,indxShflx_rsw)(1:lvar), record, ierr
     &                  MYID
         goto 99                                         !--> ERROR
       endif
      endif

      if (wrtavg(indxShflx_rlw)) then
        ierr=nf_fwrite(shflx_rlw_avg(START_2D_ARRAY), ncidavg, 
     &                                        avgShflx_rlw, record,
     &                                                     r2dvar)
       if (ierr .ne. nf_noerr) then
           lvar=lenstr(vname(1,indxShflx_rlw))
           write(stdout,1) vname(1,indxShflx_rlw)(1:lvar), record, ierr
     &                  MYID
            goto 99                                         !--> ERROR
       endif
      endif

      if (wrtavg(indxShflx_lat)) then
        ierr=nf_fwrite(shflx_lat_avg(START_2D_ARRAY), ncidavg, 
     &                                        avgShflx_lat, record,
     &                                                     r2dvar)
       if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxShflx_lat))
          write(stdout,1) vname(1,indxShflx_lat)(1:lvar), record, ierr
     &                  MYID
          goto 99                                         !--> ERROR
       endif
      endif

      if (wrtavg(indxShflx_sen)) then
         ierr=nf_fwrite(shflx_sen_avg(START_2D_ARRAY), ncidavg, 
     &    avgShflx_sen, record,
     &                                                     r2dvar)
       if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxShflx_sen))
          write(stdout,1) vname(1,indxShflx_sen)(1:lvar), record, ierr
     &                  MYID
          goto 99                                         !--> ERROR
       endif
      endif
#  endif /* BULK_FLUX */


# endif /* SOLVE3D */
  1   format(/' WRT_AVG - ERROR while writing variable(',1x,a,1x,
     &               ')into averages file.',/,11x,'Time record:',
     &                      i6,3x,'netCDF error code',i4,3x,a,i4)
      goto 100
  99  may_day_flag=3
 100  continue
!
! Close or synchronize netCDF file to disk to allow other
! processes to access data immediately after it is written.
!
# if defined MPI && !defined PARALLEL_FILES
      ierr=nf_close(ncidavg)
      if (nrpfavg.gt.0 .and. record.ge.nrpfavg) ncidavg=-1
# else
      if (nrpfavg.gt.0 .and. record.ge.nrpfavg) then
        ierr=nf_close(ncidavg)
        ncidavg=-1
      else
        ierr=nf_sync(ncidavg)
      endif
# endif
      if (ierr .eq. nf_noerr) then
        write(stdout,'(6x,A,2(A,I4,1x),A,I3)') 'WRT_AVG -- wrote ',
     &            'averaged fields into time record =', record, '/',
     &             nrecavg  MYID
      else
        write(stdout,'(/1x,2A/)') 'WRT_AVG ERROR: Cannot ',
     &             'synchronize/close averages netCDF file.'
        may_day_flag=3
      endif
 
# if defined MPI && !defined PARALLEL_FILES
      if (mynode .lt. NNODES-1) call MPI_Send (blank, 1, MPI_INTEGER, 
     &                            mynode+1, 3, ocean_grid_comm, ierr)
# endif
      return
      end
#else
      subroutine wrt_avg_empty
      end
#endif /* AVERAGES */
 
