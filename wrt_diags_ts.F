#include "cppdefs.h"
#if defined DIAGNOSTICS_TS 
                               ! Writes requested model
      subroutine wrt_diags_ts       ! fields at requested levels
                               ! into history netCDF file.
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


      call def_diags_ts (nciddiags_ts, nrecdiags_ts, ierr)
      if (ierr .ne. nf_noerr) goto 99
      lstr=lenstr(diags_tsname)


!      write(*,*) 'is it?', nrecdiags_ts, nrpfdiags_ts


!                                       !!! WARNING: Once time
! Set record within the file.           !!! stepping has been
!                                       !!! started, it is assumed
      nrecdiags_ts=max(nrecdiags_ts,1)            !!! that the global history
      if (nrpfdiags_ts.eq.0) then            !!! record index "nrecdiags_ts"
        record=nrecdiags_ts                  !!! is advanced by main.
      else
        record=1+mod(nrecdiags_ts-1, nrpfdiags_ts)
      endif
!      write(*,*) 'nope', nrecdiags_ts, nrpfdiags_ts, record
 
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
      ibuff(4)=nrecdiags_ts
      start(1)=1
      start(2)=record
      count(1)=4
      count(2)=1
      ierr=nf_put_vara_int (nciddiags_ts, diags_tsTstep, start, count, ibuff)
      if (ierr .ne. nf_noerr) then
        write(stdout,1) 'time_step', record, ierr, nf_strerror(ierr)
     &                   MYID
        goto 99                                           !--> ERROR
      endif
!
! Time
!
      ierr=nf_put_var1_FTYPE (nciddiags_ts, diags_tsTime, record, time)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxTime))
        write(stdout,1) vname(1,indxTime)(1:lvar), record, ierr,
     &                  nf_strerror(ierr) MYID
        goto 99                                           !--> ERROR
      endif
 




!---------------------------------------------------------------
! Tracer diagnostic variables.
!---------------------------------------------------------------
!
!
      do itrc=1,NT
        if (wrtdiags_ts(2)) then


            ivar=r3dvar


!------------------------------
!
! indxTXadv
!

          ierr=nf_fwrite(TXadv(START_2D_ARRAY,1,itrc), nciddiags_ts,
     &                               diags_tsXadv(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxTXadv+itrc-1))
            write(stdout,1) vname(1,indxTXadv+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif   
!
! indxTYadv
!
          ierr=nf_fwrite(TYadv(START_2D_ARRAY,1,itrc), nciddiags_ts,
     &                               diags_tsYadv(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxTYadv+itrc-1))
            write(stdout,1) vname(1,indxTYadv+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif
!
! indxTVadv
!
          ierr=nf_fwrite(TVadv(START_2D_ARRAY,1,itrc), nciddiags_ts,
     &                               diags_tsVadv(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxTVadv+itrc-1))
            write(stdout,1) vname(1,indxTVadv+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif
!
! indxTForc
!
          ierr=nf_fwrite(TForc(START_2D_ARRAY,1,itrc), nciddiags_ts,
     &                               diags_tsForc(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxTForc+itrc-1))
            write(stdout,1) vname(1,indxTForc+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif

!
! indxTHmix
!
          ierr=nf_fwrite(THmix(START_2D_ARRAY,1,itrc), nciddiags_ts,
     &                               diags_tsHmix(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxTHmix+itrc-1))
            write(stdout,1) vname(1,indxTHmix+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif
!
! indxTVmix
!
          ierr=nf_fwrite(TVmix(START_2D_ARRAY,1,itrc), nciddiags_ts,
     &                               diags_tsVmix(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxTVmix+itrc-1))
            write(stdout,1) vname(1,indxTVmix+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif
!
! indxTrate
!
          ierr = nf_fwrite(Trate(START_2D_ARRAY,1,itrc), 
     &        nciddiags_ts, diags_tsrate(itrc), record, ivar)
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



  1   format(/1x, 'WRT_DIAGS_TS ERROR while writing variable ''', A,
     &         ''' into history file.'  /11x, 'Time record:', I6,
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
      ierr=nf_close (nciddiags_ts)
      if (nrpfdiags_ts.gt.0 .and. record.ge.nrpfdiags_ts) nciddiags_ts=-1
#else
      if (nrpfdiags_ts.gt.0 .and. record.ge.nrpfdiags_ts) then
        ierr=nf_close (nciddiags_ts)
        nciddiags_ts=-1
      else
        ierr=nf_sync(nciddiags_ts)
      endif
#endif
      if (ierr .eq. nf_noerr) then
        write(stdout,'(6x,2A,I4,1x,A,I4,1x,A,I3)')  'WRT_DIAGS_TS -- ',
     &               'wrote diag fields into record =', record,
     &               '/', nrecdiags_ts MYID
      else
        write(stdout,'(/1x,2A/)') 'ERROR in wrt_DIAGS_TS: Can not ',
     &                 'synchronize/close history netCDF file.'
        may_day_flag=3
      endif
#if defined MPI && !defined PARALLEL_FILES
      if (mynode .lt. NNODES-1) call MPI_Send (blank, 1, MPI_INTEGER,
     &                            mynode+1, 2, ocean_grid_comm, ierr)
#endif
      return
      end
#endif /* DIAGNOSTICS_TS */
