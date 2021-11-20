#include "cppdefs.h"
#if defined DIAGNOSTICS_PV && defined AVERAGES

!
!---------------------------------------------------------------
! Write time-averaged diags flux fields into averages netCDF file
! Writes requested model fields at requested levels into
! diagnostics netCDF file.
!---------------------------------------------------------------
!
      subroutine wrt_diags_pv_avg     

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

#include "diags_pv.h"

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


      call def_diags_pv_avg (nciddiags_pv_avg, nrecdiags_pv_avg, ierr)
      if (ierr .ne. nf_noerr) goto 99
      lstr=lenstr(diags_pv_avgname)




!                                          !!! WARNING: Once time  
! Set record within the file.              !!! stepping has been 
!                                          !!! started, it is assumed
      nrecdiags_pv_avg=max(nrecdiags_pv_avg,1) !!! that global hist.
      if (nrpfdiags_pv_avg.eq.0) then          !!! record index is
        record=nrecdiags_pv_avg                !!! advanced by main.  
      else
        record=1+mod(nrecdiags_pv_avg-1, nrpfdiags_pv_avg)
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
      ibuff(4)=nrecdiags_pv_avg
      start(1)=1
      start(2)=record
      count(1)=4
      count(2)=1
      ierr=nf_put_vara_int (nciddiags_pv_avg, diags_pv_avgTstep,
     &                                             start, count, ibuff)
      if (ierr .ne. nf_noerr) then
        write(stdout,1) 'time_step_avg',record, ierr, nf_strerror(ierr)
     &                   MYID
        goto 99                                           !--> ERROR
      endif
!
! Averaged diag Time
!
      ierr=nf_put_var1_FTYPE (nciddiags_pv_avg, diags_pv_avgTime,
     &                                record,timediags_pv_avg)
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

        if (wrtdiags_pv_avg(2)) then

          do itrc=1,2


          if (itrc.eq.1) then
            ivar=u2dvar
          else
            ivar=v2dvar
          endif



!------------------------------
!
! indxpvForc_sfc
!

          ierr=nf_fwrite(pvForc_sfc_avg(START_2D_ARRAY,itrc), nciddiags_pv_avg,
     &                               diags_pvForc_sfc_avg(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxpvForc_sfc+itrc-1))
            write(stdout,1) vname(1,indxpvForc_sfc+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif   
!
! indxpvForc_bot
!
          ierr=nf_fwrite(pvForc_bot_avg(START_2D_ARRAY,itrc), nciddiags_pv_avg,
     &                               diags_pvForc_bot_avg(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxpvForc_bot+itrc-1))
            write(stdout,1) vname(1,indxpvForc_bot+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif



      enddo

      do itrc=1,NT

            ivar=r2dvar



!
! indxpvDiab_sfc
!
          ierr=nf_fwrite(pvDiab_sfc_avg(START_2D_ARRAY,itrc), nciddiags_pv_avg,
     &                               diags_pvDiab_sfc_avg(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxpvDiab_sfc+itrc-1))
            write(stdout,1) vname(1,indxpvDiab_sfc+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif
!
! indxpvDiab_bot
!
          ierr=nf_fwrite(pvDiab_bot_avg(START_2D_ARRAY,itrc), nciddiags_pv_avg,
     &                               diags_pvDiab_bot_avg(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxpvDiab_bot+itrc-1))
            write(stdout,1) vname(1,indxpvDiab_bot+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif



      enddo
      endif

!


!---------------------------------------------------------------


  1   format(/1x, 'WRT_DIAGS_PV_AVG while writing variable ''', A,
     &         ''' into diags_pv_avg file.'  /11x, 'Time record:', I6,
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
      ierr=nf_close (nciddiags_pv_avg)
      if (nrpfdiags_pv_avg.gt.0 .and. record.ge.nrpfdiags_pv_avg) nciddiags_pv_avg=-1
#else
      if (nrpfdiags_pv_avg.gt.0 .and. record.ge.nrpfdiags_pv_avg) then
        ierr=nf_close (nciddiags_pv_avg)
        nciddiags_pv_avg=-1
      else
        ierr=nf_sync(nciddiags_pv_avg)
      endif
#endif
      if (ierr .eq. nf_noerr) then
        write(stdout,'(6x,2A,I4,1x,A,I4,1x,A,I3)')  'WRT_DIAGS_PV_AVG -- ',
     &               'wrote fields into record =', record,
     &               '/', nrechis MYID
      else
        write(stdout,'(/1x,2A/)') 'ERROR in wrt_DIAGS_PV_AVG: Can not ',
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
