#include "cppdefs.h"
#if defined DIAGNOSTICS_UV && defined AVERAGES

!
!---------------------------------------------------------------
! Write time-averaged diags flux fields into averages netCDF file
! Writes requested model fields at requested levels into
! diagnostics netCDF file.
!---------------------------------------------------------------
!
      subroutine wrt_diags_uv_avg     

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

#include "diags_uv.h"

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


      call def_diags_uv_avg (nciddiags_uv_avg, nrecdiags_uv_avg, ierr)
      if (ierr .ne. nf_noerr) goto 99
      lstr=lenstr(diags_uv_avgname)




!                                          !!! WARNING: Once time  
! Set record within the file.              !!! stepping has been 
!                                          !!! started, it is assumed
      nrecdiags_uv_avg=max(nrecdiags_uv_avg,1) !!! that global hist.
      if (nrpfdiags_uv_avg.eq.0) then          !!! record index is
        record=nrecdiags_uv_avg                !!! advanced by main.  
      else
        record=1+mod(nrecdiags_uv_avg-1, nrpfdiags_uv_avg)
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
      ibuff(4)=nrecdiags_uv_avg
      start(1)=1
      start(2)=record
      count(1)=4
      count(2)=1
      ierr=nf_put_vara_int (nciddiags_uv_avg, diags_uv_avgTstep,
     &                                             start, count, ibuff)
      if (ierr .ne. nf_noerr) then
        write(stdout,1) 'time_step_avg',record, ierr, nf_strerror(ierr)
     &                   MYID
        goto 99                                           !--> ERROR
      endif
!
! Averaged diag Time
!

      write(*,*) 'uv diags_uv_avgTime, timediags_uv_avg',
     &                   diags_uv_avgTime, timediags_uv_avg

      ierr=nf_put_var1_FTYPE (nciddiags_uv_avg, diags_uv_avgTime,
     &                                record,timediags_uv_avg)
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
      do itrc=1,2
        if (wrtdiags_uv_avg(itrc+1)) then


          if (itrc.eq.1) then
            ivar=u3dvar
          else
            ivar=v3dvar
          endif

!------------------------------
!
! indxMXadv
!

          ierr=nf_fwrite(MXadv_avg(START_2D_ARRAY,1,itrc), nciddiags_uv_avg,
     &                               diags_uvXadv_avg(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxMXadv+itrc-1))
            write(stdout,1) vname(1,indxMXadv+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif   
!
! indxMYadv
!
          ierr=nf_fwrite(MYadv_avg(START_2D_ARRAY,1,itrc), nciddiags_uv_avg,
     &                               diags_uvYadv_avg(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxMYadv+itrc-1))
            write(stdout,1) vname(1,indxMYadv+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif
!
! indxMVadv
!
          ierr=nf_fwrite(MVadv_avg(START_2D_ARRAY,1,itrc), nciddiags_uv_avg,
     &                               diags_uvVadv_avg(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxMVadv+itrc-1))
            write(stdout,1) vname(1,indxMVadv+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif
!
! indxMCor
!
          ierr=nf_fwrite(MCor_avg(START_2D_ARRAY,1,itrc), nciddiags_uv_avg,
     &                               diags_uvCor_avg(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxMCor+itrc-1))
            write(stdout,1) vname(1,indxMCor+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif
!
! indxMPrsgrd
!
          ierr=nf_fwrite(MPrsgrd_avg(START_2D_ARRAY,1,itrc), nciddiags_uv_avg,
     &                               diags_uvPrsgrd_avg(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxMPrsgrd+itrc-1))
            write(stdout,1) vname(1,indxMPrsgrd+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif
!
! indxMHmix
!
          ierr=nf_fwrite(MHmix_avg(START_2D_ARRAY,1,itrc), nciddiags_uv_avg,
     &                               diags_uvHmix_avg(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxMHmix+itrc-1))
            write(stdout,1) vname(1,indxMHmix+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif
!
! indxMVmix
!
          ierr=nf_fwrite(MVmix_avg(START_2D_ARRAY,1,itrc), nciddiags_uv_avg,
     &                               diags_uvVmix_avg(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxMVmix+itrc-1))
            write(stdout,1) vname(1,indxMVmix+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif
!
! indxMrate
!
          ierr = nf_fwrite(Mrate_avg(START_2D_ARRAY,1,itrc), 
     &        nciddiags_uv_avg, diags_uvrate_avg(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxMrate+itrc-1))
            write(stdout,1) vname(1,indxMrate+itrc-1)(1:lvar),
     &           record, ierr, nf_strerror(ierr) MYID 
            goto 99
          endif


!
! indxMVmix2
!
          ierr=nf_fwrite(MVmix2_avg(START_2D_ARRAY,1,itrc), nciddiags_uv_avg,
     &                               diags_uvVmix2_avg(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxMVmix2+itrc-1))
            write(stdout,1) vname(1,indxMVmix2+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif



      endif
      enddo
!


!---------------------------------------------------------------


  1   format(/1x, 'WRT_DIAGS_UV_AVG while writing variable ''', A,
     &         ''' into diags_uv_avg file.'  /11x, 'Time record:', I6,
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
      ierr=nf_close (nciddiags_uv_avg)
      if (nrpfdiags_uv_avg.gt.0 .and. record.ge.nrpfdiags_uv_avg) nciddiags_uv_avg=-1
#else
      if (nrpfdiags_uv_avg.gt.0 .and. record.ge.nrpfdiags_uv_avg) then
        ierr=nf_close (nciddiags_uv_avg)
        nciddiags_uv_avg=-1
      else
        ierr=nf_sync(nciddiags_uv_avg)
      endif
#endif
      if (ierr .eq. nf_noerr) then
        write(stdout,'(6x,2A,I4,1x,A,I4,1x,A,I3)')  'WRT_DIAGS_UV_AVG -- ',
     &  'wrote fields into record =', record,
     &  '/', nrechis MYID
      else
        write(stdout,'(/1x,2A/)') 'ERROR in wrt_DIAGS_UV_AVG: Can not ',
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
