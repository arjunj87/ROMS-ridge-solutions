#include "cppdefs.h"
#if defined DIAGNOSTICS_VRT 
                               ! Writes requested model
      subroutine wrt_diags_vrt       ! fields at requested levels
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

      call def_diags_vrt (nciddiags_vrt, nrecdiags_vrt, ierr)
      if (ierr .ne. nf_noerr) goto 99
      lstr=lenstr(diags_vrtname)



!                                       !!! WARNING: Once time
! Set record within the file.           !!! stepping has been
!                                       !!! started, it is assumed
      nrecdiags_vrt=max(nrecdiags_vrt,1)            !!! that the global history
      if (nrpfdiags_vrt.eq.0) then            !!! record index "nrecdiags_vrt"
        record=nrecdiags_vrt                  !!! is advanced by main.
      else
        record=1+mod(nrecdiags_vrt-1, nrpfdiags_vrt)
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
      ibuff(4)=nrecdiags_vrt
      start(1)=1
      start(2)=record
      count(1)=4
      count(2)=1
      ierr=nf_put_vara_int (nciddiags_vrt, diags_vrtTstep, start, count, ibuff)
      if (ierr .ne. nf_noerr) then
        write(stdout,1) 'time_step', record, ierr, nf_strerror(ierr)
     &                   MYID
        goto 99                                           !--> ERROR
      endif
!
! Time
!
      ierr=nf_put_var1_FTYPE (nciddiags_vrt, diags_vrtTime, record, time)
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

        if (wrtdiags_vrt(2)) then

            itrc = 1
            ivar=p2dvar


!------------------------------
!
! indxvrtXadv
!

          ierr=nf_fwrite(vrtXadv, nciddiags_vrt,
     &                               diags_vrtXadv(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxvrtXadv+itrc-1))
            write(stdout,1) vname(1,indxvrtXadv+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif   
!
! indxvrtYadv
!
          ierr=nf_fwrite(vrtYadv(START_2D_ARRAY), nciddiags_vrt,
     &                               diags_vrtYadv(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxvrtYadv+itrc-1))
            write(stdout,1) vname(1,indxvrtYadv+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif

!
! indxvrtCor
!
          ierr=nf_fwrite(vrtCor(START_2D_ARRAY), nciddiags_vrt,
     &                               diags_vrtCor(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxvrtCor+itrc-1))
            write(stdout,1) vname(1,indxvrtCor+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif
!
! indxvrtPrsgrd
!
          ierr=nf_fwrite(vrtPrsgrd(START_2D_ARRAY), nciddiags_vrt,
     &                               diags_vrtPrsgrd(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxvrtPrsgrd+itrc-1))
            write(stdout,1) vname(1,indxvrtPrsgrd+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif
!
! indxvrtHmix
!
          ierr=nf_fwrite(vrtHmix(START_2D_ARRAY), nciddiags_vrt,
     &                               diags_vrtHmix(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxvrtHmix+itrc-1))
            write(stdout,1) vname(1,indxvrtHmix+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif
!
! indxvrtVmix
!
          ierr=nf_fwrite(vrtVmix(START_2D_ARRAY), nciddiags_vrt,
     &                               diags_vrtVmix(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxvrtVmix+itrc-1))
            write(stdout,1) vname(1,indxvrtVmix+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif
!
! indxvrtrate
!
          ierr = nf_fwrite(vrtrate(START_2D_ARRAY), 
     &        nciddiags_vrt, diags_vrtrate(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxvrtrate+itrc-1))
            write(stdout,1) vname(1,indxvrtrate+itrc-1)(1:lvar),
     &           record, ierr, nf_strerror(ierr) MYID 
            goto 99
          endif


!
! indxvrtVmix2
!
          ierr=nf_fwrite(vrtVmix2(START_2D_ARRAY), nciddiags_vrt,
     &                               diags_vrtVmix2(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxvrtVmix2+itrc-1))
            write(stdout,1) vname(1,indxvrtVmix2+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif



!
! indxvrtWind
!
          ierr=nf_fwrite(vrtWind(START_2D_ARRAY), nciddiags_vrt,
     &                               diags_vrtWind(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxvrtWind+itrc-1))
            write(stdout,1) vname(1,indxvrtWind+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif


!
! indxvrtDrag
!
          ierr=nf_fwrite(vrtDrag(START_2D_ARRAY), nciddiags_vrt,
     &                               diags_vrtDrag(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxvrtDrag+itrc-1))
            write(stdout,1) vname(1,indxvrtDrag+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif


      endif

!

!---------------------------------------------------------------



  1   format(/1x, 'WRT_DIAGS_VRT ERROR while writing variable ''', A,
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
      ierr=nf_close (nciddiags_vrt)
      if (nrpfdiags_vrt.gt.0 .and. record.ge.nrpfdiags_vrt) nciddiags_vrt=-1
#else
      if (nrpfdiags_vrt.gt.0 .and. record.ge.nrpfdiags_vrt) then
        ierr=nf_close (nciddiags_vrt)
        nciddiags_vrt=-1
      else
        ierr=nf_sync(nciddiags_vrt)
      endif
#endif
      if (ierr .eq. nf_noerr) then
        write(stdout,'(6x,2A,I4,1x,A,I4,1x,A,I3)')  'WRT_DIAGS_VRT -- ',
     &               'wrote diag fields into record =', record,
     &               '/', nrecdiags_vrt MYID
      else
        write(stdout,'(/1x,2A/)') 'ERROR in wrt_DIAGS_VRT: Can not ',
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
