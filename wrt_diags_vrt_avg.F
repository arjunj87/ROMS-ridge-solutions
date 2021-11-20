#include "cppdefs.h"
#if defined DIAGNOSTICS_VRT && defined AVERAGES

!
!---------------------------------------------------------------
! Write time-averaged diags flux fields into averages netCDF file
! Writes requested model fields at requested levels into
! diagnostics netCDF file.
!---------------------------------------------------------------
!
      subroutine wrt_diags_vrt_avg     

      implicit none
      integer ierr,  record, lstr, lvar, lenstr
     &   , start(2), count(2), ibuff(4), nf_fwrite
#include "param.h"
#include "scalars.h"
#include "grid.h"
#include "ocean2d.h"
#include "ocean3d.h"
#include "eos_vars.h"
#include "mixing.h"
#include "ncvars.h"
#include "netcdf.inc"

#include "diags_vrt.h"

#ifdef SOLVE3D
      integer tile, itrc,i,j,k,ivar
# include "work.h"
#endif
#if defined MPI && !defined PARALLEL_FILES
# include "mpif.h"
      integer status(MPI_STATUS_SIZE), blank
      if (mynode.gt.0) call MPI_Recv (blank, 1, MPI_INTEGER,
     &            mynode-1, 2, ocean_grid_comm, status, ierr)
#endif
!
! Create/open history file; write grid arrays, if so needed.
!


      call def_diags_vrt_avg (nciddiags_vrt_avg, nrecdiags_vrt_avg, ierr)
      if (ierr .ne. nf_noerr) goto 99
      lstr=lenstr(diags_vrt_avgname)




!                                          !!! WARNING: Once time  
! Set record within the file.              !!! stepping has been 
!                                          !!! started, it is assumed
      nrecdiags_vrt_avg=max(nrecdiags_vrt_avg,1) !!! that global hist.
      if (nrpfdiags_vrt_avg.eq.0) then          !!! record index is
        record=nrecdiags_vrt_avg                !!! advanced by main.  
      else
        record=1+mod(nrecdiags_vrt_avg-1, nrpfdiags_vrt_avg)
      endif




 
 
c-#define CR
CR      write(*,*) 'wrt_his: Entry  ' MYID
!
! Write out evolving model variables:
! ----- --- -------- ----- ----------
!
! Time step number and record numbers.
!
      ibuff(1)=iic
      ibuff(2)=nrecrst
      ibuff(3)=nrechis
      ibuff(4)=nrecdiags_vrt_avg
      start(1)=1
      start(2)=record
      count(1)=4
      count(2)=1
      ierr=nf_put_vara_int (nciddiags_vrt_avg, diags_vrt_avgTstep,
     &                                             start, count, ibuff)
      if (ierr .ne. nf_noerr) then
        write(stdout,1) 'time_step_avg',record, ierr, nf_strerror(ierr)
     &                   MYID
        goto 99                                           !--> ERROR
      endif
!
! Averaged diag Time
!
      ierr=nf_put_var1_FTYPE (nciddiags_vrt_avg, diags_vrt_avgTime,
     &                                record,timediags_vrt_avg)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxTime))
        write(stdout,1) vname(1,indxTime)(1:lvar), record, ierr,
     &                  nf_strerror(ierr) MYID
        goto 99                                           !--> ERROR
      endif
 



!---------------------------------------------------------------
! Momentum diagnostic variables.
!---------------------------------------------------------------
!
!

        if (wrtdiags_vrt_avg(2)) then

            itrc = 1
            ivar=p2dvar


!------------------------------
!
! indxvrtXadv
!

          ierr=nf_fwrite(vrtXadv_avg(START_2D_ARRAY), nciddiags_vrt_avg,
     &                               diags_vrtXadv_avg(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxvrtXadv+itrc-1))
            write(stdout,1) vname(1,indxvrtXadv+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif   
!
! indxvrtYadv
!
          ierr=nf_fwrite(vrtYadv_avg(START_2D_ARRAY), nciddiags_vrt_avg,
     &                               diags_vrtYadv_avg(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxvrtYadv+itrc-1))
            write(stdout,1) vname(1,indxvrtYadv+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif

!
! indxvrtCor
!
          ierr=nf_fwrite(vrtCor_avg(START_2D_ARRAY), nciddiags_vrt_avg,
     &                               diags_vrtCor_avg(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxvrtCor+itrc-1))
            write(stdout,1) vname(1,indxvrtCor+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif
!
! indxvrtPrsgrd
!
          ierr=nf_fwrite(vrtPrsgrd_avg(START_2D_ARRAY), nciddiags_vrt_avg,
     &                               diags_vrtPrsgrd_avg(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxvrtPrsgrd+itrc-1))
            write(stdout,1) vname(1,indxvrtPrsgrd+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif
!
! indxvrtHmix
!
          ierr=nf_fwrite(vrtHmix_avg(START_2D_ARRAY), nciddiags_vrt_avg,
     &                               diags_vrtHmix_avg(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxvrtHmix+itrc-1))
            write(stdout,1) vname(1,indxvrtHmix+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif
!
! indxvrtVmix
!
          ierr=nf_fwrite(vrtVmix_avg(START_2D_ARRAY), nciddiags_vrt_avg,
     &                               diags_vrtVmix_avg(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxvrtVmix+itrc-1))
            write(stdout,1) vname(1,indxvrtVmix+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif
!
! indxvrtrate
!
          ierr = nf_fwrite(vrtrate_avg(START_2D_ARRAY), 
     &        nciddiags_vrt_avg, diags_vrtrate_avg(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxvrtrate+itrc-1))
            write(stdout,1) vname(1,indxvrtrate+itrc-1)(1:lvar),
     &           record, ierr, nf_strerror(ierr) MYID 
            goto 99
          endif


!
! indxvrtVmix2
!
          ierr=nf_fwrite(vrtVmix2_avg(START_2D_ARRAY), nciddiags_vrt_avg,
     &                               diags_vrtVmix2_avg(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxvrtVmix2+itrc-1))
            write(stdout,1) vname(1,indxvrtVmix2+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif


!
! indxvrtWind
!
          ierr=nf_fwrite(vrtWind_avg(START_2D_ARRAY), nciddiags_vrt_avg,
     &                               diags_vrtWind_avg(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxvrtWind+itrc-1))
            write(stdout,1) vname(1,indxvrtWind+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif


!
! indxvrtDrag
!
          ierr=nf_fwrite(vrtDrag_avg(START_2D_ARRAY), nciddiags_vrt_avg,
     &                               diags_vrtDrag_avg(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxvrtDrag+itrc-1))
            write(stdout,1) vname(1,indxvrtDrag+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif

      endif

!


!---------------------------------------------------------------


  1   format(/1x, 'WRT_DIAGS_VRT_AVG while writing variable ''', A,
     &         ''' into diags_vrt_avg file.'  /11x, 'Time record:', I6,
     &         3x,'netCDF error code',i4 /11x,'Cause of error: ',
     &                                            A, 3x, A, i4)
      goto 100
  99  may_day_flag=3
 100  continue 
 
!
! Synchronize netCDF file to disk to allow other processes
! to access data immediately after it is written.
!
#if defined MPI && !defined PARALLEL_FILES
      ierr=nf_close (nciddiags_vrt_avg)
      if (nrpfdiags_vrt_avg.gt.0 .and. record.ge.nrpfdiags_vrt_avg) nciddiags_vrt_avg=-1
#else
      if (nrpfdiags_vrt_avg.gt.0 .and. record.ge.nrpfdiags_vrt_avg) then
        ierr=nf_close (nciddiags_vrt_avg)
        nciddiags_vrt_avg=-1
      else
        ierr=nf_sync(nciddiags_vrt_avg)
      endif
#endif
      if (ierr .eq. nf_noerr) then
        write(stdout,'(6x,2A,I4,1x,A,I4,1x,A,I3)')  'WRT_DIAGS_VRT_AVG -- ',
     &               'wrote fields into record =', record,
     &               '/', nrechis MYID
      else
        write(stdout,'(/1x,2A/)') 'ERROR in wrt_DIAGS_VRT_AVG: Can not ',
     &                 'synchronize/close netCDF file.'
        may_day_flag=3
      endif
#if defined MPI && !defined PARALLEL_FILES
      if (mynode .lt. NNODES-1) call MPI_Send (blank, 1, MPI_INTEGER,
     &                            mynode+1, 2, ocean_grid_comm, ierr)
#endif
      return
      end
#endif
