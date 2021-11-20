#include "cppdefs.h"
#if defined DIAGNOSTICS_EDDY && defined AVERAGES

!
!---------------------------------------------------------------
! Write time-averaged diags flux fields into averages netCDF file
! Writes requested model fields at requested levels into
! diagnostics netCDF file.
!---------------------------------------------------------------
!
      subroutine wrt_diags_eddy_avg     

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

#include "diags_eddy.h"

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


      call def_diags_eddy_avg (nciddiags_eddy_avg, nrecdiags_eddy_avg, ierr)
      if (ierr .ne. nf_noerr) goto 99
      lstr=lenstr(diags_eddy_avgname)




!                                          !!! WARNING: Once time  
! Set record within the file.              !!! stepping has been 
!                                          !!! started, it is assumed
      nrecdiags_eddy_avg=max(nrecdiags_eddy_avg,1) !!! that global hist.
      if (nrpfdiags_eddy_avg.eq.0) then          !!! record index is
        record=nrecdiags_eddy_avg                !!! advanced by main.  
      else
        record=1+mod(nrecdiags_eddy_avg-1, nrpfdiags_eddy_avg)
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
      ibuff(4)=nrecdiags_eddy_avg
      start(1)=1
      start(2)=record
      count(1)=4
      count(2)=1
      ierr=nf_put_vara_int (nciddiags_eddy_avg, diags_eddy_avgTstep,
     &                                             start, count, ibuff)
      if (ierr .ne. nf_noerr) then
        write(stdout,1) 'time_step_avg',record, ierr, nf_strerror(ierr)
     &                   MYID
        goto 99                                           !--> ERROR
      endif
!
! Averaged diag Time
!



      ierr=nf_put_var1_FTYPE (nciddiags_eddy_avg, diags_eddy_avgTime,
     &                                record,timediags_eddy_avg)
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
        itrc=1
        if (wrtdiags_eddy_avg(itrc+1)) then



        ivar=r3dvar


!------------------------------
!
! indxeddyuu
!

          ierr=nf_fwrite(eddyuu_avg(START_2D_ARRAY,1), nciddiags_eddy_avg,
     &                               diags_eddyuu_avg(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxeddyuu+itrc-1))
            write(stdout,1) vname(1,indxeddyuu+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif   
!
! indxeddyvv
!
          ierr=nf_fwrite(eddyvv_avg(START_2D_ARRAY,1), nciddiags_eddy_avg,
     &                               diags_eddyvv_avg(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxeddyvv+itrc-1))
            write(stdout,1) vname(1,indxeddyvv+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif
!
! indxeddyuv
!
          ierr=nf_fwrite(eddyuv_avg(START_2D_ARRAY,1), nciddiags_eddy_avg,
     &                               diags_eddyuv_avg(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxeddyuv+itrc-1))
            write(stdout,1) vname(1,indxeddyuv+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif
!
! indxeddyub
!
          ierr=nf_fwrite(eddyub_avg(START_2D_ARRAY,1), nciddiags_eddy_avg,
     &                               diags_eddyub_avg(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxeddyub+itrc-1))
            write(stdout,1) vname(1,indxeddyub+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif
!
! indxeddyvb
!
          ierr=nf_fwrite(eddyvb_avg(START_2D_ARRAY,1), nciddiags_eddy_avg,
     &                               diags_eddyvb_avg(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxeddyvb+itrc-1))
            write(stdout,1) vname(1,indxeddyvb+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif
!
! indxeddywb
!
          ierr=nf_fwrite(eddywb_avg(START_2D_ARRAY,1), nciddiags_eddy_avg,
     &                               diags_eddywb_avg(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxeddywb+itrc-1))
            write(stdout,1) vname(1,indxeddywb+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif
!
! indxeddynotused1
!
          ierr=nf_fwrite(eddynotused1_avg(START_2D_ARRAY,1), nciddiags_eddy_avg,
     &                               diags_eddynotused1_avg(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxeddynotused1+itrc-1))
            write(stdout,1) vname(1,indxeddynotused1+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif
!
! indxeddynotused2
!
          ierr = nf_fwrite(eddynotused2_avg(START_2D_ARRAY,1), 
     &        nciddiags_eddy_avg, diags_eddynotused2_avg(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxeddynotused2+itrc-1))
            write(stdout,1) vname(1,indxeddynotused2+itrc-1)(1:lvar),
     &           record, ierr, nf_strerror(ierr) MYID 
            goto 99
          endif


!
! indxeddynotused3
!
          ierr=nf_fwrite(eddynotused3_avg(START_2D_ARRAY,1), nciddiags_eddy_avg,
     &                               diags_eddynotused3_avg(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxeddynotused3+itrc-1))
            write(stdout,1) vname(1,indxeddynotused3+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif



      endif

!


!---------------------------------------------------------------


  1   format(/1x, 'WRT_DIAGS_EDDY_AVG while writing variable ''', A,
     &         ''' into diags_eddy_avg file.'  /11x, 'Time record:', I6,
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
      ierr=nf_close (nciddiags_eddy_avg)
      if (nrpfdiags_eddy_avg.gt.0 .and. record.ge.nrpfdiags_eddy_avg) nciddiags_eddy_avg=-1
#else
      if (nrpfdiags_eddy_avg.gt.0 .and. record.ge.nrpfdiags_eddy_avg) then
        ierr=nf_close (nciddiags_eddy_avg)
        nciddiags_eddy_avg=-1
      else
        ierr=nf_sync(nciddiags_eddy_avg)
      endif
#endif
      if (ierr .eq. nf_noerr) then
        write(stdout,'(6x,2A,I4,1x,A,I4,1x,A,I3)')  'WRT_DIAGS_EDDY_AVG -- ',
     &  'wrote fields into record =', record,
     &  '/', nrechis MYID
      else
        write(stdout,'(/1x,2A/)') 'ERROR in wrt_DIAGS_EDDY_AVG: Can not ',
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
