#include "cppdefs.h"
#if defined DIAGNOSTICS_TS && defined AVERAGES

!
!---------------------------------------------------------------
! Write time-averaged diags flux fields into averages netCDF file
! Writes requested model fields at requested levels into
! diagnostics netCDF file.
!---------------------------------------------------------------
!
      subroutine wrt_diags_ts_avg     

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

#include "diags_ts.h"

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


      call def_diags_ts_avg (nciddiags_ts_avg, nrecdiags_ts_avg, ierr)
      if (ierr .ne. nf_noerr) goto 99
      lstr=lenstr(diags_ts_avgname)




!                                          !!! WARNING: Once time  
! Set record within the file.              !!! stepping has been 
!                                          !!! started, it is assumed
      nrecdiags_ts_avg=max(nrecdiags_ts_avg,1) !!! that global hist.
      if (nrpfdiags_ts_avg.eq.0) then          !!! record index is
        record=nrecdiags_ts_avg                !!! advanced by main.  
      else
        record=1+mod(nrecdiags_ts_avg-1, nrpfdiags_ts_avg)
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
      ibuff(4)=nrecdiags_ts_avg
      start(1)=1
      start(2)=record
      count(1)=4
      count(2)=1
      ierr=nf_put_vara_int (nciddiags_ts_avg, diags_ts_avgTstep,
     &                                             start, count, ibuff)
      if (ierr .ne. nf_noerr) then
        write(stdout,1) 'time_step_avg',record, ierr, nf_strerror(ierr)
     &                   MYID
        goto 99                                           !--> ERROR
      endif
!
! Averaged diag Time
!

      write(*,*) 'ts diags_ts_avgTime, timediags_ts_avg',
     &                   diags_ts_avgTime, timediags_ts_avg

      ierr=nf_put_var1_FTYPE (nciddiags_ts_avg, diags_ts_avgTime,
     &                                record,timediags_ts_avg)
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
      do itrc=1,NT
        if (wrtdiags_ts_avg(2)) then



            ivar=r3dvar


!------------------------------
!
! indxTXadv
!

          ierr=nf_fwrite(TXadv_avg(START_2D_ARRAY,1,itrc), nciddiags_ts_avg,
     &                               diags_tsXadv_avg(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxTXadv+itrc-1))
            write(stdout,1) vname(1,indxTXadv+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif   
!
! indxTYadv
!
          ierr=nf_fwrite(TYadv_avg(START_2D_ARRAY,1,itrc), nciddiags_ts_avg,
     &                               diags_tsYadv_avg(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxTYadv+itrc-1))
            write(stdout,1) vname(1,indxTYadv+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif
!
! indxTVadv
!
          ierr=nf_fwrite(TVadv_avg(START_2D_ARRAY,1,itrc), nciddiags_ts_avg,
     &                               diags_tsVadv_avg(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxTVadv+itrc-1))
            write(stdout,1) vname(1,indxTVadv+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif
!
! indxTForc
!
          ierr=nf_fwrite(TForc_avg(START_2D_ARRAY,1,itrc), nciddiags_ts_avg,
     &                               diags_tsForc_avg(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxTForc+itrc-1))
            write(stdout,1) vname(1,indxTForc+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif

!
! indxTHmix
!
          ierr=nf_fwrite(THmix_avg(START_2D_ARRAY,1,itrc), nciddiags_ts_avg,
     &                               diags_tsHmix_avg(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxTHmix+itrc-1))
            write(stdout,1) vname(1,indxTHmix+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif
!
! indxTVmix
!
          ierr=nf_fwrite(TVmix_avg(START_2D_ARRAY,1,itrc), nciddiags_ts_avg,
     &                               diags_tsVmix_avg(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxTVmix+itrc-1))
            write(stdout,1) vname(1,indxTVmix+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif
!
! indxTrate
!
          ierr = nf_fwrite(Trate_avg(START_2D_ARRAY,1,itrc), 
     &        nciddiags_ts_avg, diags_tsrate_avg(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxTrate+itrc-1))
            write(stdout,1) vname(1,indxTrate+itrc-1)(1:lvar),
     &           record, ierr, nf_strerror(ierr) MYID 
            goto 99
          endif



      endif
      enddo
!


!---------------------------------------------------------------


  1   format(/1x, 'WRT_DIAGS_TS_AVG while writing variable ''', A,
     &         ''' into diags_ts_avg file.'  /11x, 'Time record:', I6,
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
      ierr=nf_close (nciddiags_ts_avg)
      if (nrpfdiags_ts_avg.gt.0 .and. record.ge.nrpfdiags_ts_avg) nciddiags_ts_avg=-1
#else
      if (nrpfdiags_ts_avg.gt.0 .and. record.ge.nrpfdiags_ts_avg) then
        ierr=nf_close (nciddiags_ts_avg)
        nciddiags_ts_avg=-1
      else
        ierr=nf_sync(nciddiags_ts_avg)
      endif
#endif
      if (ierr .eq. nf_noerr) then
        write(stdout,'(6x,2A,I4,1x,A,I4,1x,A,I3)')  'WRT_DIAGS_TS_AVG -- ',
     &  'wrote fields into record =', record,
     &  '/', nrechis MYID
      else
        write(stdout,'(/1x,2A/)') 'ERROR in wrt_DIAGS_TS_AVG: Can not ',
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
