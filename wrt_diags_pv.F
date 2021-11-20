#include "cppdefs.h"
#if defined DIAGNOSTICS_PV 
                               ! Writes requested model
      subroutine wrt_diags_pv       ! fields at requested levels
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

      call def_diags_pv (nciddiags_pv, nrecdiags_pv, ierr)
      if (ierr .ne. nf_noerr) goto 99
      lstr=lenstr(diags_pvname)



!                                       !!! WARNING: Once time
! Set record within the file.           !!! stepping has been
!                                       !!! started, it is assumed
      nrecdiags_pv=max(nrecdiags_pv,1)            !!! that the global history
      if (nrpfdiags_pv.eq.0) then            !!! record index "nrecdiags_pv"
        record=nrecdiags_pv                  !!! is advanced by main.
      else
        record=1+mod(nrecdiags_pv-1, nrpfdiags_pv)
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
      ibuff(4)=nrecdiags_pv
      start(1)=1
      start(2)=record
      count(1)=4
      count(2)=1
      ierr=nf_put_vara_int (nciddiags_pv, diags_pvTstep, start, count, ibuff)
      if (ierr .ne. nf_noerr) then
        write(stdout,1) 'time_step', record, ierr, nf_strerror(ierr)
     &                   MYID
        goto 99                                           !--> ERROR
      endif
!
! Time
!
      ierr=nf_put_var1_FTYPE (nciddiags_pv, diags_pvTime, record, time)
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

        if (wrtdiags_pv(2)) then


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

          ierr=nf_fwrite(pvForc_sfc(START_2D_ARRAY,itrc), nciddiags_pv,
     &                               diags_pvForc_sfc(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxpvForc_sfc+itrc-1))
            write(stdout,1) vname(1,indxpvForc_sfc+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif   
!
! indxpvForc_bot
!
          ierr=nf_fwrite(pvForc_bot(START_2D_ARRAY,itrc), nciddiags_pv,
     &                               diags_pvForc_bot(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxpvForc_bot+itrc-1))
            write(stdout,1) vname(1,indxpvForc_bot+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif



#ifdef DIAGNOSTICS_PV_SUP


!------------------------------
!
! indxpvForc_0
!

          ierr=nf_fwrite(pvForc_0(START_2D_ARRAY,itrc), nciddiags_pv,
     &                               diags_pvForc_0(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxpvForc_0+itrc-1))
            write(stdout,1) vname(1,indxpvForc_0+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif   


!------------------------------
!
! indxpvForc_1
!

          ierr=nf_fwrite(pvForc_1(START_2D_ARRAY,itrc), nciddiags_pv,
     &                               diags_pvForc_1(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxpvForc_1+itrc-1))
            write(stdout,1) vname(1,indxpvForc_1+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif   


!------------------------------
!
! indxpvForc_2
!

          ierr=nf_fwrite(pvForc_2(START_2D_ARRAY,itrc), nciddiags_pv,
     &                               diags_pvForc_2(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxpvForc_2+itrc-1))
            write(stdout,1) vname(1,indxpvForc_2+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif   


!------------------------------
!
! indxpvForc_3
!

          ierr=nf_fwrite(pvForc_3(START_2D_ARRAY,itrc), nciddiags_pv,
     &                               diags_pvForc_3(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxpvForc_3+itrc-1))
            write(stdout,1) vname(1,indxpvForc_3+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif   


#endif

      enddo




          do itrc=1,NT



            ivar=r2dvar


!
! indxpvDiab_sfc
!
          ierr=nf_fwrite(pvDiab_sfc(START_2D_ARRAY,itrc), nciddiags_pv,
     &                               diags_pvDiab_sfc(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxpvDiab_sfc+itrc-1))
            write(stdout,1) vname(1,indxpvDiab_sfc+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif
!
! indxpvDiab_bot
!
          ierr=nf_fwrite(pvDiab_bot(START_2D_ARRAY,itrc), nciddiags_pv,
     &                               diags_pvDiab_bot(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxpvDiab_bot+itrc-1))
            write(stdout,1) vname(1,indxpvDiab_bot+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif




#ifdef DIAGNOSTICS_PV_SUP


!
! indxpvDiab_0
!
          ierr=nf_fwrite(pvDiab_0(START_2D_ARRAY,itrc), nciddiags_pv,
     &                               diags_pvDiab_0(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxpvDiab_0+itrc-1))
            write(stdout,1) vname(1,indxpvDiab_0+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif


!
! indxpvDiab_1
!
          ierr=nf_fwrite(pvDiab_1(START_2D_ARRAY,itrc), nciddiags_pv,
     &                               diags_pvDiab_1(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxpvDiab_1+itrc-1))
            write(stdout,1) vname(1,indxpvDiab_1+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif


!
! indxpvDiab_2
!
          ierr=nf_fwrite(pvDiab_2(START_2D_ARRAY,itrc), nciddiags_pv,
     &                               diags_pvDiab_2(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxpvDiab_2+itrc-1))
            write(stdout,1) vname(1,indxpvDiab_2+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif

!
! indxpvDiab_3
!
          ierr=nf_fwrite(pvDiab_3(START_2D_ARRAY,itrc), nciddiags_pv,
     &                               diags_pvDiab_3(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxpvDiab_3+itrc-1))
            write(stdout,1) vname(1,indxpvDiab_3+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif

#endif


      enddo
      endif

!

!---------------------------------------------------------------



  1   format(/1x, 'WRT_DIAGS_PV ERROR while writing variable ''', A,
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
      ierr=nf_close (nciddiags_pv)
      if (nrpfdiags_pv.gt.0 .and. record.ge.nrpfdiags_pv) nciddiags_pv=-1
#else
      if (nrpfdiags_pv.gt.0 .and. record.ge.nrpfdiags_pv) then
        ierr=nf_close (nciddiags_pv)
        nciddiags_pv=-1
      else
        ierr=nf_sync(nciddiags_pv)
      endif
#endif
      if (ierr .eq. nf_noerr) then
        write(stdout,'(6x,2A,I4,1x,A,I4,1x,A,I3)')  'WRT_DIAGS_PV -- ',
     &               'wrote diag fields into record =', record,
     &               '/', nrecdiags_pv MYID
      else
        write(stdout,'(/1x,2A/)') 'ERROR in wrt_DIAGS_PV: Can not ',
     &                 'synchronize/close history netCDF file.'
        may_day_flag=3
      endif
#if defined MPI && !defined PARALLEL_FILES
      if (mynode .lt. NNODES-1) call MPI_Send (blank, 1, MPI_INTEGER,
     &                            mynode+1, 2, ocean_grid_comm, ierr)
#endif
      return
      end
#endif /* DIAGNOSTICS_PV */
