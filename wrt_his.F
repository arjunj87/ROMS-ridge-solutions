#include "cppdefs.h"
                               ! Writes requested model
      subroutine wrt_his       ! fields at requested levels
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
# ifdef BULK_FLUX
#include "forces.h"
# endif
#ifdef SOLVE3D
      integer tile, itrc,i,j,k
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
      call def_his (ncidhis, nrechis, ierr)
      if (ierr .ne. nf_noerr) goto 99
      lstr=lenstr(hisname)
!                                       !!! WARNING: Once time
! Set record within the file.           !!! stepping has been
!                                       !!! started, it is assumed
      nrechis=max(nrechis,1)            !!! that the global history
      if (nrpfhis.eq.0) then            !!! record index "nrechis"
        record=nrechis                  !!! is advanced by main.
      else
        record=1+mod(nrechis-1, nrpfhis)
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
#ifdef AVERAGES
      ibuff(4)=nrecavg
#else
      ibuff(4)=0
#endif
      start(1)=1
      start(2)=record
      count(1)=4
      count(2)=1
      ierr=nf_put_vara_int (ncidhis, hisTstep, start, count, ibuff)
      if (ierr .ne. nf_noerr) then
        write(stdout,1) 'time_step', record, ierr, nf_strerror(ierr)
     &                   MYID
        goto 99                                           !--> ERROR
      endif
!
! Time
!
      ierr=nf_put_var1_FTYPE (ncidhis, hisTime, record, time)
      if (ierr .ne. nf_noerr) then
        lvar=lenstr(vname(1,indxTime))
        write(stdout,1) vname(1,indxTime)(1:lvar), record, ierr,
     &                  nf_strerror(ierr) MYID
        goto 99                                           !--> ERROR
      endif
 
CR      write(*,*) 'wrt_his:  time  ' MYID
!
! Barotropic mode variables: free-surface and 2D momentum
! components in XI-,ETA-directions.
!
      if (wrthis(indxZ)) then
        ierr=nf_fwrite (zeta(START_2D_ARRAY,knew), ncidhis, hisZ,
     &                                              record, r2dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxZ))
          write(stdout,1) vname(1,indxZ)(1:lvar), record, ierr,
     &                    nf_strerror(ierr) MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his:  wrtZ  ' MYID
      endif
      if (wrthis(indxUb)) then
        ierr=nf_fwrite (ubar(START_2D_ARRAY,knew), ncidhis, hisUb,
     &                                              record, u2dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxUb))
          write(stdout,1) vname(1,indxUb)(1:lvar), record, ierr,
     &                    nf_strerror(ierr) MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his: wrtUBAR' MYID
      endif
      if (wrthis(indxVb)) then
        ierr=nf_fwrite (vbar(START_2D_ARRAY,knew), ncidhis, hisVb,
     &                                              record, v2dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxVb))
          write(stdout,1) vname(1,indxVb)(1:lvar), record, ierr,
     &                    nf_strerror(ierr) MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his: wrtVBAR' MYID
      endif
 
#ifdef SOLVE3D
!
! 3D momentum components in XI- and ETA-directions.
!
      if (wrthis(indxU)) then
        ierr=nf_fwrite (u(START_2D_ARRAY,1,nstp), ncidhis,
     &                             hisU,  record,  u3dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxU))
          write(stdout,1) vname(1,indxU)(1:lvar), record, ierr,
     &                    nf_strerror(ierr) MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his:  wrtU  ' MYID
      endif
      if (wrthis(indxV)) then
        ierr=nf_fwrite (v(START_2D_ARRAY,1,nstp), ncidhis,
     &                             hisV,  record,  v3dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxV))
          write(stdout,1) vname(1,indxV)(1:lvar), record, ierr,
     &                    nf_strerror(ierr) MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his:  wrtV  ' MYID
      endif
!
! Tracer variables.
!
      do itrc=1,NT
        if (wrthis(indxT+itrc-1)) then
          ierr=nf_fwrite (t(START_2D_ARRAY,1,nstp,itrc), ncidhis,
     &                                hisT(itrc), record, r3dvar)


c          write(*,*) 'Replacing T,S with rho1,qp1.'
c          if (itrc.eq.1) then
c            ierr=nf_fwrite (rho1, ncidhis, hisT(itrc), record,r3dvar)
c          elseif (itrc.eq.2) then
c            ierr=nf_fwrite  (qp1, ncidhis, hisT(itrc), record,r3dvar)
c          endif

          if (ierr .ne. nf_noerr) then
            lvar=lenstr(vname(1,indxT+itrc-1))
            write(stdout,1) vname(1,indxT+itrc-1)(1:lvar), record,
     &                      ierr, nf_strerror(ierr) MYID
            goto 99                                       !--> ERROR
          endif
CR      write(*,*) 'wrt_his:  wrtT  ' MYID
        endif
      enddo
!
! Density anomaly.
!        
        do k=1,N
          do j=0,Mm+1
            do i=0,Lm+1
# if defined NONLIN_EOS
              work(i,j,k)=rho1(i,j,k)+qp1(i,j,k)
     &                          *(z_w(i,j,N)-z_r(i,j,k))
# else
              work(i,j,k)=rho(i,j,k)
# endif
            enddo
          enddo
        enddo
      if (wrthis(indxR)) then
        ierr=nf_fwrite (work, ncidhis, hisR, record, r3dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxR))
          write(stdout,1) vname(1,indxR)(1:lvar), record, ierr,
     &                    nf_strerror(ierr) MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his:  wrtRHO' MYID
      endif
!
! S-coordinate omega vertical velocity (m/s).
!
      if (wrthis(indxO)) then
        do k=0,N
          do j=0,Mm+1
            do i=0,Lm+1
              work(i,j,k)=W(i,j,k)*pm(i,j)*pn(i,j)
#ifdef DIAGS_CR
              work(i,j,k)=crzs(i,j,k)
#endif
            enddo
          enddo
        enddo
        ierr=nf_fwrite (work, ncidhis, hisO, record, w3dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxO))
          write(stdout,1) vname(1,indxO)(1:lvar), record, ierr,
     &                    nf_strerror(ierr) MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his:  wrtO  ' MYID
      endif
!
! Write out true vertical velocity (m/s).
!
      if (wrthis(indxW)) then
        do tile=0,NSUB_X*NSUB_E-1
          call Wvlcty (tile, work)
        enddo
#ifdef DIAGS_CR
        do k=0,N
          do j=0,Mm+1
            do i=0,Lm+1
              work(i,j,k)=crzb(i,j,k)
            enddo
          enddo
        enddo
#endif
        ierr=nf_fwrite (work, ncidhis, hisW, record, r3dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxW))
          write(stdout,1) vname(1,indxW)(1:lvar), record, ierr,
     &                    nf_strerror(ierr) MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his:  wrtW  ' MYID
      endif
!
! Write out vertical viscosity coefficient.
!
      if (wrthis(indxAkv)) then
        ierr=nf_fwrite (Akv, ncidhis, hisAkv, record, w3dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxAkv))
          write(stdout,1) vname(1,indxAkv)(1:lvar), record, ierr,
     &                    nf_strerror(ierr) MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his:  wrtAkv' MYID
      endif
!
!  Write out vertical diffusion coefficient for potential temperature.
!
      if (wrthis(indxAkt)) then
        ierr=nf_fwrite (Akt(START_2D_ARRAY,0,itemp), ncidhis,
     &                                hisAkt, record, w3dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxAkt))
          write(stdout,1) vname(1,indxAkt)(1:lvar), record, ierr,
     &                    nf_strerror(ierr) MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his:  wrtAkt' MYID
      endif
# ifdef SALINITY
!
!  Write out vertical diffusion coefficient for salinity.
!
      if (wrthis(indxAks)) then
        ierr=nf_fwrite (Akt(START_2D_ARRAY,0,isalt), ncidhis,
     &                                hisAks, record, w3dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxAks))
          write(stdout,1) vname(1,indxAks)(1:lvar), record, ierr,
     &                    nf_strerror(ierr) MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his:  wrtAks' MYID
      endif
# endif
# ifdef LMD_KPP
!
! Thickness of planetary boundary layer [m].
!
      if (wrthis(indxHbl)) then
        ierr=nf_fwrite (hbls(START_2D_ARRAY,nstp), ncidhis,
     &                              hisHbl, record, r2dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxHbl))
          write(stdout,1) vname(1,indxHbl)(1:lvar), record, ierr,
     &                    nf_strerror(ierr) MYID
          goto 99                                         !--> ERROR
        endif
CR      write(*,*) 'wrt_his:  wrtHBL' MYID
# ifdef LMD_BKPP
        ierr=nf_fwrite (hbbls(START_2D_ARRAY,nstp), ncidhis,
     &                              hisHbbl, record, r2dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxHbbl))
          write(stdout,1) vname(1,indxHbbl)(1:lvar), record, ierr,
     &                    nf_strerror(ierr) MYID
          goto 99                                         !--> ERROR
        endif
# endif
      endif
# endif

#if defined ADV_ISONEUTRAL && defined WRITE_SXSY
        ierr=nf_fwrite (Sx_u(START_2D_ARRAY,1), ncidhis,
     &                                hisSX, record, u3dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxSX))
          write(stdout,1) vname(1,indxSX)(1:lvar), record, ierr,
     &                    nf_strerror(ierr) MYID
          goto 99                                         !--> ERROR
        endif

        ierr=nf_fwrite (Sy_v(START_2D_ARRAY,1), ncidhis,
     &                                hisSY, record, v3dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxSY))
          write(stdout,1) vname(1,indxSY)(1:lvar), record, ierr,
     &                    nf_strerror(ierr) MYID
          goto 99                                         !--> ERROR
        endif
#endif  /* WRITE_SXSY */


# ifdef BULK_FLUX

       if (wrthis(indxSustr_blk)) then
        ierr=nf_fwrite(sustr_blk, ncidhis, hisSustr_blk, record,
# ifdef WND_AT_RHO_POINTS        
     &                                                     r2dvar)
# else
     &                                                     u2dvar)
#endif
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxSustr_blk))
          write(stdout,1) vname(1,indxSustr_blk)(1:lvar), record, ierr
     &                  MYID
          goto 99                                         !--> ERROR
        endif
       endif

       if (wrthis(indxSvstr_blk)) then
        ierr=nf_fwrite(svstr_blk, ncidhis, hisSvstr_blk, record,
# ifdef WND_AT_RHO_POINTS
     &                                                     r2dvar)
# else
     &                                                     v2dvar)
#endif     
         if (ierr .ne. nf_noerr) then
           lvar=lenstr(vname(1,indxSvstr_blk))
           write(stdout,1) vname(1,indxSvstr_blk)(1:lvar), record, ierr
     &                  MYID
           goto 99                                         !--> ERROR
         endif
       endif

!
!  Write out surface heat flux component (degC/s)
!
       if (wrthis(indxShflx_rsw)) then
        ierr=nf_fwrite(shflx_rsw, ncidhis, hisShflx_rsw, record,
     &                                                     r2dvar)
        if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxShflx_rsw))
          write(stdout,1) vname(1,indxShflx_rsw)(1:lvar), record, ierr
     &                  MYID
          goto 99                                         !--> ERROR
        endif
       endif

       if (wrthis(indxShflx_rlw)) then
        ierr=nf_fwrite(shflx_rlw, ncidhis, hisShflx_rlw, record,
     &                                                     r2dvar)
         if (ierr .ne. nf_noerr) then
           lvar=lenstr(vname(1,indxShflx_rlw))
           write(stdout,1) vname(1,indxShflx_rlw)(1:lvar), record, ierr
     &                  MYID
           goto 99                                         !--> ERROR
         endif
       endif

      if (wrthis(indxShflx_lat)) then
         ierr=nf_fwrite(shflx_lat, ncidhis, hisShflx_lat, record,
     &                                                     r2dvar)
         if (ierr .ne. nf_noerr) then
          lvar=lenstr(vname(1,indxShflx_lat))
          write(stdout,1) vname(1,indxShflx_lat)(1:lvar), record, ierr
     &                  MYID
          goto 99                                         !--> ERROR
         endif
      endif

      if (wrthis(indxShflx_sen)) then
        ierr=nf_fwrite(shflx_sen, ncidhis, hisShflx_sen, record,
     &                                                     r2dvar)
          if (ierr .ne. nf_noerr) then
           lvar=lenstr(vname(1,indxShflx_sen))
           write(stdout,1) vname(1,indxShflx_sen)(1:lvar), record, ierr
     &                  MYID
           goto 99                                         !--> ERROR
           endif
      endif
# endif /* BULK_FLUX */



#endif
  1   format(/1x, 'WRT_HIS ERROR while writing variable ''', A,
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
      ierr=nf_close (ncidhis)
      if (nrpfhis.gt.0 .and. record.ge.nrpfhis) ncidhis=-1
#else
      if (nrpfhis.gt.0 .and. record.ge.nrpfhis) then
        ierr=nf_close (ncidhis)
        ncidhis=-1
      else
        ierr=nf_sync(ncidhis)
      endif
#endif
      if (ierr .eq. nf_noerr) then
        write(stdout,'(6x,2A,I4,1x,A,I4,1x,A,I3)')  'WRT_HIS -- ',
     &               'wrote history fields into record =', record,
     &               '/', nrechis MYID
      else
        write(stdout,'(/1x,2A/)') 'ERROR in wrt_his: Can not ',
     &                 'synchronize/close history netCDF file.'
        may_day_flag=3
      endif
#if defined MPI && !defined PARALLEL_FILES
      if (mynode .lt. NNODES-1) call MPI_Send (blank, 1, MPI_INTEGER,
     &                            mynode+1, 2, ocean_grid_comm, ierr)
#endif
      return
      end
 
