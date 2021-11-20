#include "cppdefs.h"
#if defined DIAGNOSTICS_UV 
                               ! Writes requested model
      subroutine wrt_diags_uv       ! fields at requested levels
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


      call def_diags_uv (nciddiags_uv, nrecdiags_uv, ierr)
      if (ierr .ne. nf_noerr) goto 99
      lstr=lenstr(diags_uvname)


!      write(*,*) 'is it?', nrecdiags_uv, nrpfdiags_uv


!                                       !!! WARNING: Once time
! Set record within the file.           !!! stepping has been
!                                       !!! started, it is assumed
      nrecdiags_uv=max(nrecdiags_uv,1)            !!! that the global history
      if (nrpfdiags_uv.eq.0) then            !!! record index "nrecdiags_uv"
        record=nrecdiags_uv                  !!! is advanced by main.
      else
        record=1+mod(nrecdiags_uv-1, nrpfdiags_uv)
      endif
!      write(*,*) 'nope', nrecdiags_uv, nrpfdiags_uv, record
 
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
      ibuff(4)=nrecdiags_uv
      start(1)=1
      start(2)=record
      count(1)=4
      count(2)=1
      ierr=nf_put_vara_int (nciddiags_uv, diags_uvTstep, start, count, ibuff)
      if (ierr .ne. nf_noerr) then
        write(stdout,1) 'time_step', record, ierr, nf_strerror(ierr)
     &                   MYID
        goto 99                                           !--> ERROR
      endif
!
! Time
!
      ierr=nf_put_var1_FTYPE (nciddiags_uv, diags_uvTime, record, time)
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
        if (wrtdiags_uv(itrc+1)) then


          if (itrc.eq.1) then
            ivar=u3dvar
          else
            ivar=v3dvar
          endif


!------------------------------
!
! indxMXadv
!

          ierr=nf_fwrite(MXadv(START_2D_ARRAY,1,itrc), nciddiags_uv,
     &                               diags_uvXadv(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxMXadv+itrc-1))
            write(stdout,1) vname(1,indxMXadv+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif   
!
! indxMYadv
!
          ierr=nf_fwrite(MYadv(START_2D_ARRAY,1,itrc), nciddiags_uv,
     &                               diags_uvYadv(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxMYadv+itrc-1))
            write(stdout,1) vname(1,indxMYadv+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif
!
! indxMVadv
!
          ierr=nf_fwrite(MVadv(START_2D_ARRAY,1,itrc), nciddiags_uv,
     &                               diags_uvVadv(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxMVadv+itrc-1))
            write(stdout,1) vname(1,indxMVadv+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif
!
! indxMCor
!
          ierr=nf_fwrite(MCor(START_2D_ARRAY,1,itrc), nciddiags_uv,
     &                               diags_uvCor(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxMCor+itrc-1))
            write(stdout,1) vname(1,indxMCor+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif
!
! indxMPrsgrd
!
          ierr=nf_fwrite(MPrsgrd(START_2D_ARRAY,1,itrc), nciddiags_uv,
     &                               diags_uvPrsgrd(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxMPrsgrd+itrc-1))
            write(stdout,1) vname(1,indxMPrsgrd+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif
!
! indxMHmix
!
          ierr=nf_fwrite(MHmix(START_2D_ARRAY,1,itrc), nciddiags_uv,
     &                               diags_uvHmix(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxMHmix+itrc-1))
            write(stdout,1) vname(1,indxMHmix+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif
!
! indxMVmix
!
          ierr=nf_fwrite(MVmix(START_2D_ARRAY,1,itrc), nciddiags_uv,
     &                               diags_uvVmix(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxMVmix+itrc-1))
            write(stdout,1) vname(1,indxMVmix+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif
!
! indxMrate
!
          ierr = nf_fwrite(Mrate(START_2D_ARRAY,1,itrc), 
     &        nciddiags_uv, diags_uvrate(itrc), record, ivar)
          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxMrate+itrc-1))
            write(stdout,1) vname(1,indxMrate+itrc-1)(1:lvar),
     &           record, ierr, nf_strerror(ierr) MYID 
            goto 99
          endif


!
! indxMVmix2
!
          ierr=nf_fwrite(MVmix2(START_2D_ARRAY,1,itrc), nciddiags_uv,
     &                               diags_uvVmix2(itrc), record, ivar)
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



  1   format(/1x, 'WRT_DIAGS_UV ERROR while writing variable ''', A,
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
      ierr=nf_close (nciddiags_uv)
      if (nrpfdiags_uv.gt.0 .and. record.ge.nrpfdiags_uv) nciddiags_uv=-1
#else
      if (nrpfdiags_uv.gt.0 .and. record.ge.nrpfdiags_uv) then
        ierr=nf_close (nciddiags_uv)
        nciddiags_uv=-1
      else
        ierr=nf_sync(nciddiags_uv)
      endif
#endif
      if (ierr .eq. nf_noerr) then
        write(stdout,'(6x,2A,I4,1x,A,I4,1x,A,I3)')  'WRT_DIAGS_UV -- ',
     &               'wrote diag fields into record =', record,
     &               '/', nrecdiags_uv MYID
      else
        write(stdout,'(/1x,2A/)') 'ERROR in wrt_DIAGS_UV: Can not ',
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
