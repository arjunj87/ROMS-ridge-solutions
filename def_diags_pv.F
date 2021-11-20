#ifndef AVRH
# include "cppdefs.h"
#endif

#ifdef DIAGNOSTICS_PV

#ifndef AVRH
      subroutine def_diags_pv (ncid, total_rec, ierr)
#else
      subroutine def_diags_pv_avg(ncid, total_rec, ierr)
#endif
!
! Create/open averages/history netCDF file. In the case when a new
! netCDF file is created, define all variables, their dimensions and
! attributes. In the case when a previously existing netCDF file is
! to be opened for addition of new data, verify that all dimensions
! of the file are consistent with the present model configuration
! and all necessary variables exist. Save netCDF IDs for all needed
! variables. Also determine size of the unlimited dimension.
!
! The difference between def_his and def_avg is as follows: they
! have different netCDF file name (hisname/avgname); netCDF file ID
! (passed as argument); time record index (hisindx/avgindx); array
! of switches which variables to write (wrthis/wrtavg); and different
! sets of netCDF variable IDs (hisTime...hisHbl/avgTime...avgHbl);
! and the first attribute of each variable, long_name, has prefix
! 'averaged'. Because most of the code is identical for both
! routines, the second one is generated from the first entirely
! by CPP.
!
      implicit none
      logical create_new_file
      integer ncid, total_rec, ierr, rec, lstr,lvar,lenstr, timedim
     &      , r2dgrd(3),  u2dgrd(3), v2dgrd(3),  auxil(2),  checkdims
#ifdef SOLVE3D
     &      , uv2dgrd(3)
     &      , r3dgrd(4),  u3dgrd(4), v3dgrd(4),  w3dgrd(4), itrc
     &      , p2dgrd(3), p3dgrd(4), pw3dgrd(4)
#endif

      character*60 text
      integer my_nf_def_dim, ifour
      parameter (ifour=4)

#include "param.h"
#include "scalars.h"
#include "ncvars.h"
#include "netcdf.inc"

#include "diags_pv.h"

#ifdef MASK_LAND_DATA
# ifdef HIS_DOUBLE
      real*8 spv_set
      parameter (spv_set=1.D+33)
# else
      real*4 spv_set
      parameter (spv_set=1.E+33)
# endif
#endif


#ifndef AVRH

# define ncname diags_pvname
# define rec_per_file nrpfdiags_pv
# define wrt wrtdiags_pv
# define vidTime diags_pvTime
# define vidTstep diags_pvTstep

#  define viddiags_pvForc_sfc diags_pvForc_sfc
#  define viddiags_pvForc_bot diags_pvForc_bot
#  define viddiags_pvDiab_sfc diags_pvDiab_sfc
#  define viddiags_pvDiab_bot diags_pvDiab_bot

#ifdef DIAGNOSTICS_PV_SUP
#  define viddiags_pvForc_0 diags_pvForc_0
#  define viddiags_pvForc_1 diags_pvForc_1
#  define viddiags_pvForc_2 diags_pvForc_2
#  define viddiags_pvForc_3 diags_pvForc_3
#  define viddiags_pvDiab_0 diags_pvDiab_0
#  define viddiags_pvDiab_1 diags_pvDiab_1
#  define viddiags_pvDiab_2 diags_pvDiab_2
#  define viddiags_pvDiab_3 diags_pvDiab_3
#endif

#else

# define ncname diags_pv_avgname
# define rec_per_file nrpfdiags_pv_avg
# define wrt wrtdiags_pv_avg
# define vidTime diags_pv_avgTime
# define vidTstep diags_pv_avgTstep

#  define viddiags_pvForc_sfc diags_pvForc_sfc_avg
#  define viddiags_pvForc_bot diags_pvForc_bot_avg
#  define viddiags_pvDiab_sfc diags_pvDiab_sfc_avg
#  define viddiags_pvDiab_bot diags_pvDiab_bot_avg
#endif




!
! Put time record index into file name. In  the case when model
! output is to be arranged into sequence of named files, the naming
! convention is as follows: 'his_root.INDEX.[MPI_node.]nc', where
! INDEX is an integer number such that (i) it is divisible by the
! specified number of records per file; and (ii)
!
!      INDEX + record_within_the_file = total_record
!
! where, 1 =< record_within_the_file =< records_per_file, so that
! total_record changes continuously throughout the sequence of files.
!      



      ierr=0
      lstr=lenstr(ncname)
      if (rec_per_file.gt.0) then
        lvar=total_rec-(1+mod(total_rec-1, rec_per_file))
        call insert_time_index (ncname, lstr, lvar, ierr)
        if (ierr .ne. 0) goto 99
      endif



!
! Decide whether to create a new file, or open existing one.
! Overall the whole code below is organized into 3-way switch,
!
! 10  if (create_new_file) then
!        .... create new file, save netCDF ids for all variables;
!     elseif (ncid.eq.-1) then
!        .... try to open existing file and check its dimensions
!       if (cannot be opened or rejected) then
!         create_new_file=.true.
!         goto 10
!       endif   and prepare
!        .... prepare the file for adding new data,
!        .... find and save netCDF ids for all variables
!     else
!        .... just open, no checking, all ids are assumed to be
!        .... already known (MPI single file output only).
!     endif
!
! which is designed to implement flexible opening policy:
! if ldefhis=.true., it forces creation of a new file [if the
! file already exists, it will be overwritten]; on the other hand,
! ldefhis=.false., it is assumed that the file already exists and
! an attempt to open it is made; if the attempt is successful, the
! file is prepared for appending hew data; if it fails, a new file
! is created.
!

       if (lvar.eq.total_rec-1) then
#ifndef AVRH
      create_new_file=ldefdiags_pv
#else
      create_new_file=ldefdiags_pv_avg
#endif
       else
#ifndef AVRH
      create_new_file=.false.
#else
      create_new_file=.false.
#endif
       endif


      if (ncid.ne.-1) create_new_file=.false.
#if defined MPI & !defined PARALLEL_FILES
      if (mynode.gt.0) create_new_file=.false.
#endif

!
! Create new history/averages file:    Put global attributes
! ====== === ======= ======== =====    and define all variables.
!
  10  if (create_new_file) then
c        ierr=nf_create(ncname(1:lstr), nf_clobber, ncid)
        ierr=nf_create(ncname(1:lstr), nf_64bit_offset, ncid)
        if (ierr .ne. nf_noerr) then
          write(stdout,'(/1x,4A/)') 'ERROR in def_his/avg: Can not ',
     &                 'create netCDF file ''', ncname(1:lstr), '''.'
          goto 99
        endif
        if (rec_per_file.eq.0) total_rec=0
!
! Put global attributes.
! --- ------ -----------
!
        call put_global_atts (ncid, ierr)
!
! Define dimensions of staggered fields.
! ------ ---------- -- --------- -------
!
        ierr=my_nf_def_dim (ncid, 'xi_rho',  xi_rho,  r2dgrd(1))
        ierr=my_nf_def_dim (ncid, 'xi_u',    xi_u,    u2dgrd(1))
        ierr=my_nf_def_dim (ncid, 'eta_rho', eta_rho, r2dgrd(2))
        ierr=my_nf_def_dim (ncid, 'eta_v',   eta_v,   v2dgrd(2))
#ifdef SOLVE3D
        ierr=my_nf_def_dim (ncid, 's_rho',   N,       r3dgrd(3))
        ierr=my_nf_def_dim (ncid, 's_w',     N+1,     w3dgrd(3))
#endif
        ierr=my_nf_def_dim (ncid, 'time', nf_unlimited, timedim)
        ierr=my_nf_def_dim (ncid, 'auxil',   ifour,    auxil(1))

        auxil(2)=timedim
 
        r2dgrd(3)=timedim           ! Free surface
 
        u2dgrd(2)=r2dgrd(2)         ! 2D UBAR-type
        u2dgrd(3)=timedim
 
        v2dgrd(1)=r2dgrd(1)         ! 2D VBAR-type
        v2dgrd(3)=timedim

        p2dgrd(1)=u2dgrd(1)         ! 
        p2dgrd(2)=v2dgrd(2)         ! 2D vort-type
        p2dgrd(3)=timedim

#ifdef SOLVE3D

        r3dgrd(1)=r2dgrd(1)         !
        r3dgrd(2)=r2dgrd(2)         ! 3D RHO-type
        r3dgrd(4)=timedim           !
 
        u3dgrd(1)=u2dgrd(1)         !
        u3dgrd(2)=r2dgrd(2)         ! 3D U-type
        u3dgrd(3)=r3dgrd(3)         !
        u3dgrd(4)=timedim
 
        v3dgrd(1)=r2dgrd(1)         !
        v3dgrd(2)=v2dgrd(2)         ! 3D V-type
        v3dgrd(3)=r3dgrd(3)         !
        v3dgrd(4)=timedim
 
        w3dgrd(1)=r2dgrd(1)         !
        w3dgrd(2)=r2dgrd(2)         ! 3D W-type
        w3dgrd(4)=timedim           !

        p3dgrd(1)=u2dgrd(1)         !
        p3dgrd(2)=v2dgrd(2)         ! 3D vort-type
        p3dgrd(3)=r2dgrd(3)         ! 
        p3dgrd(4)=timedim           !

        pw3dgrd(1)=u2dgrd(1)         !
        pw3dgrd(2)=v2dgrd(2)         ! 3D vortW-type
        pw3dgrd(3)=w3dgrd(3)         ! 
        pw3dgrd(4)=timedim         

#endif
#if (defined PUT_GRID_INTO_HISTORY && !defined AVRH)\
 || (defined PUT_GRID_INTO_AVERAGES && defined AVRH)
!
! Define grid variables.
! ------ ---- ----------
!
        if (total_rec.le.1) call def_grid (ncid, r2dgrd)
#endif
 
!
! Define evolving model variables.
! ------ -------- ----- ----------
!
! Time step number and time record indices:
!
        ierr=nf_def_var (ncid, 'time_step', nf_int, 2, auxil,
     &                                                 vidTstep)
        ierr=nf_put_att_text (ncid, vidTstep, 'long_name', 48,
     &       'time step and record numbers from initialization')
!
! Time.
!
        lvar=lenstr(vname(1,indxTime))
        ierr=nf_def_var (ncid, vname(1,indxTime)(1:lvar), NF_FOUT,
     &                                        1, timedim, vidTime)
#ifdef AVRH
        text='averaged '/ /vname(2,indxTime)
#else
        text=vname(2,indxTime)
#endif
        lvar=lenstr(text)
        ierr=nf_put_att_text (ncid, vidTime, 'long_name', lvar,
     &                                             text(1:lvar))
        lvar=lenstr(vname(3,indxTime))
        ierr=nf_put_att_text (ncid, vidTime, 'units',  lvar,
     &                                vname(3,indxTime)(1:lvar))


!------------------------------
! Momentum diagnostics variables.
! -------- ----------- --------



       if (wrt(2)) then


        do itrc=1,2  ! use itrc indice = 1 for u_momentum
!                    !                 = 2 for v_momentum 


          if (itrc.eq.1) then
            uv2dgrd=u2dgrd
          else
            uv2dgrd=v2dgrd
          endif

!------------------------------
!
! indxpvForc_sfc
!          
          lvar=lenstr(vname(1,indxpvForc_sfc+itrc-1))
          ierr=nf_def_var (ncid, vname(1,indxpvForc_sfc+itrc-1)(1:lvar),
     &                     NF_FOUT, 3, uv2dgrd, viddiags_pvForc_sfc(itrc))

!        write(*,*) 'viddiags_pvForc_sfc(itrc)', viddiags_pvForc_sfc(itrc)
# ifdef AVRH
          text='averaged '/ /vname(2,indxpvForc_sfc+itrc-1)
# else
          text=vname(2,indxpvForc_sfc+itrc-1)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, viddiags_pvForc_sfc(itrc), 'long_name', 
     &                          lvar, text(1:lvar))


          lvar=lenstr(vname(3,indxpvForc_sfc+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_pvForc_sfc(itrc), 'units', lvar,
     &                          vname(3,indxpvForc_sfc+itrc-1)(1:lvar))

#ifdef MASK_LAND_DATA 
          ierr=nf_put_att_FOUT (ncid, viddiags_pvForc_sfc(itrc), '_FillValue',
     &                                   NF_FOUT, 1, spv_set)
#endif
!------------------------------
!
! indxpvForc_bot
!          
          lvar=lenstr(vname(1,indxpvForc_bot+itrc-1))
          ierr=nf_def_var (ncid, vname(1,indxpvForc_bot+itrc-1)(1:lvar),
     &                     NF_FOUT, 3, uv2dgrd, viddiags_pvForc_bot(itrc))

!        write(*,*) 'viddiags_pvForc_bot(itrc)', viddiags_pvForc_bot(itrc)


# ifdef AVRH
          text='averaged '/ /vname(2,indxpvForc_bot+itrc-1)
# else
          text=vname(2,indxpvForc_bot+itrc-1)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, viddiags_pvForc_bot(itrc), 'long_name', 
     &                          lvar, text(1:lvar))

          lvar=lenstr(vname(3,indxpvForc_bot+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_pvForc_bot(itrc), 'units', lvar,
     &                          vname(3,indxpvForc_bot+itrc-1)(1:lvar))
#ifdef MASK_LAND_DATA 
          ierr=nf_put_att_FOUT (ncid, viddiags_pvForc_bot(itrc), '_FillValue',
     &                                   NF_FOUT, 1, spv_set)
#endif

#ifdef DIAGNOSTICS_PV_SUP
!------------------------------
#ifndef AVRH


!------------------------------
!
! indxpvForc_0
!          
          lvar=lenstr(vname(1,indxpvForc_0+itrc-1))
          ierr=nf_def_var (ncid, vname(1,indxpvForc_0+itrc-1)(1:lvar),
     &                     NF_FOUT, 3, uv2dgrd, viddiags_pvForc_0(itrc))

!        write(*,*) 'viddiags_pvForc_0(itrc)', viddiags_pvForc_0(itrc)
# ifdef AVRH
          text='averaged '/ /vname(2,indxpvForc_0+itrc-1)
# else
          text=vname(2,indxpvForc_0+itrc-1)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, viddiags_pvForc_0(itrc), 'long_name', 
     &                          lvar, text(1:lvar))


          lvar=lenstr(vname(3,indxpvForc_0+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_pvForc_0(itrc), 'units', lvar,
     &                          vname(3,indxpvForc_0+itrc-1)(1:lvar))

#ifdef MASK_LAND_DATA 
          ierr=nf_put_att_FOUT (ncid, viddiags_pvForc_0(itrc), '_FillValue',
     &                                   NF_FOUT, 1, spv_set)
#endif
!------------------------------
!
! indxpvForc_1
!          
          lvar=lenstr(vname(1,indxpvForc_1+itrc-1))
          ierr=nf_def_var (ncid, vname(1,indxpvForc_1+itrc-1)(1:lvar),
     &                     NF_FOUT, 3, uv2dgrd, viddiags_pvForc_1(itrc))

!        write(*,*) 'viddiags_pvForc_1(itrc)', viddiags_pvForc_1(itrc)


# ifdef AVRH
          text='averaged '/ /vname(2,indxpvForc_1+itrc-1)
# else
          text=vname(2,indxpvForc_1+itrc-1)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, viddiags_pvForc_1(itrc), 'long_name', 
     &                          lvar, text(1:lvar))

          lvar=lenstr(vname(3,indxpvForc_1+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_pvForc_1(itrc), 'units', lvar,
     &                          vname(3,indxpvForc_1+itrc-1)(1:lvar))
#ifdef MASK_LAND_DATA 
          ierr=nf_put_att_FOUT (ncid, viddiags_pvForc_1(itrc), '_FillValue',
     &                                   NF_FOUT, 1, spv_set)
#endif





!------------------------------
!
! indxpvForc_2
!          
          lvar=lenstr(vname(1,indxpvForc_2+itrc-1))
          ierr=nf_def_var (ncid, vname(1,indxpvForc_2+itrc-1)(1:lvar),
     &                     NF_FOUT, 3, uv2dgrd, viddiags_pvForc_2(itrc))

!        write(*,*) 'viddiags_pvForc_2(itrc)', viddiags_pvForc_2(itrc)
# ifdef AVRH
          text='averaged '/ /vname(2,indxpvForc_2+itrc-1)
# else
          text=vname(2,indxpvForc_2+itrc-1)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, viddiags_pvForc_2(itrc), 'long_name', 
     &                          lvar, text(1:lvar))


          lvar=lenstr(vname(3,indxpvForc_2+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_pvForc_2(itrc), 'units', lvar,
     &                          vname(3,indxpvForc_2+itrc-1)(1:lvar))

#ifdef MASK_LAND_DATA 
          ierr=nf_put_att_FOUT (ncid, viddiags_pvForc_2(itrc), '_FillValue',
     &                                   NF_FOUT, 1, spv_set)
#endif
!------------------------------
!
! indxpvForc_3
!          
          lvar=lenstr(vname(1,indxpvForc_3+itrc-1))
          ierr=nf_def_var (ncid, vname(1,indxpvForc_3+itrc-1)(1:lvar),
     &                     NF_FOUT, 3, uv2dgrd, viddiags_pvForc_3(itrc))

!        write(*,*) 'viddiags_pvForc_3(itrc)', viddiags_pvForc_3(itrc)


# ifdef AVRH
          text='averaged '/ /vname(2,indxpvForc_3+itrc-1)
# else
          text=vname(2,indxpvForc_3+itrc-1)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, viddiags_pvForc_3(itrc), 'long_name', 
     &                          lvar, text(1:lvar))

          lvar=lenstr(vname(3,indxpvForc_3+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_pvForc_3(itrc), 'units', lvar,
     &                          vname(3,indxpvForc_3+itrc-1)(1:lvar))
#ifdef MASK_LAND_DATA 
          ierr=nf_put_att_FOUT (ncid, viddiags_pvForc_3(itrc), '_FillValue',
     &                                   NF_FOUT, 1, spv_set)
#endif



#endif
#endif


      enddo


        do itrc=1,NT 


!------------------------------
!
! indxpvDiab_sfc
!          
          lvar=lenstr(vname(1,indxpvDiab_sfc+itrc-1))
          ierr=nf_def_var (ncid, vname(1,indxpvDiab_sfc+itrc-1)(1:lvar),
     &                     NF_FOUT, 3, r2dgrd, viddiags_pvDiab_sfc(itrc))
# ifdef AVRH
          text='averaged '/ /vname(2,indxpvDiab_sfc+itrc-1)
# else
          text=vname(2,indxpvDiab_sfc+itrc-1)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, viddiags_pvDiab_sfc(itrc), 'long_name', 
     &                          lvar, text(1:lvar))

          lvar=lenstr(vname(3,indxpvDiab_sfc+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_pvDiab_sfc(itrc), 'units', lvar,
     &                          vname(3,indxpvDiab_sfc+itrc-1)(1:lvar))
#ifdef MASK_LAND_DATA 
          ierr=nf_put_att_FOUT (ncid, viddiags_pvDiab_sfc(itrc), '_FillValue',
     &                                   NF_FOUT, 1, spv_set)
#endif



!------------------------------
!
! indxpvDiab_bot
!          
          lvar=lenstr(vname(1,indxpvDiab_bot+itrc-1))
          ierr=nf_def_var (ncid, vname(1,indxpvDiab_bot+itrc-1)(1:lvar),
     &                     NF_FOUT, 3, r2dgrd, viddiags_pvDiab_bot(itrc))
# ifdef AVRH
          text='averaged '/ /vname(2,indxpvDiab_bot+itrc-1)
# else
          text=vname(2,indxpvDiab_bot+itrc-1)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, viddiags_pvDiab_bot(itrc), 'long_name', 
     &                          lvar, text(1:lvar))

          lvar=lenstr(vname(3,indxpvDiab_bot+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_pvDiab_bot(itrc), 'units', lvar,
     &                          vname(3,indxpvDiab_bot+itrc-1)(1:lvar))
#ifdef MASK_LAND_DATA 
          ierr=nf_put_att_FOUT (ncid, viddiags_pvDiab_bot(itrc), '_FillValue',
     &                                   NF_FOUT, 1, spv_set)
#endif



#ifdef DIAGNOSTICS_PV_SUP
!------------------------------
#ifndef AVRH
!------------------------------




!------------------------------
!
! indxpvDiab_0
!          
          lvar=lenstr(vname(1,indxpvDiab_0+itrc-1))
          ierr=nf_def_var (ncid, vname(1,indxpvDiab_0+itrc-1)(1:lvar),
     &                     NF_FOUT, 3, r2dgrd, viddiags_pvDiab_0(itrc))
# ifdef AVRH
          text='averaged '/ /vname(2,indxpvDiab_0+itrc-1)
# else
          text=vname(2,indxpvDiab_0+itrc-1)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, viddiags_pvDiab_0(itrc), 'long_name', 
     &                          lvar, text(1:lvar))

          lvar=lenstr(vname(3,indxpvDiab_0+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_pvDiab_0(itrc), 'units', lvar,
     &                          vname(3,indxpvDiab_0+itrc-1)(1:lvar))
#ifdef MASK_LAND_DATA 
          ierr=nf_put_att_FOUT (ncid, viddiags_pvDiab_0(itrc), '_FillValue',
     &                                   NF_FOUT, 1, spv_set)
#endif



!------------------------------
!
! indxpvDiab_1
!          
          lvar=lenstr(vname(1,indxpvDiab_1+itrc-1))
          ierr=nf_def_var (ncid, vname(1,indxpvDiab_1+itrc-1)(1:lvar),
     &                     NF_FOUT, 3, r2dgrd, viddiags_pvDiab_1(itrc))
# ifdef AVRH
          text='averaged '/ /vname(2,indxpvDiab_1+itrc-1)
# else
          text=vname(2,indxpvDiab_1+itrc-1)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, viddiags_pvDiab_1(itrc), 'long_name', 
     &                          lvar, text(1:lvar))

          lvar=lenstr(vname(3,indxpvDiab_1+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_pvDiab_1(itrc), 'units', lvar,
     &                          vname(3,indxpvDiab_1+itrc-1)(1:lvar))
#ifdef MASK_LAND_DATA 
          ierr=nf_put_att_FOUT (ncid, viddiags_pvDiab_1(itrc), '_FillValue',
     &                                   NF_FOUT, 1, spv_set)
#endif





!------------------------------
!
! indxpvDiab_2
!          
          lvar=lenstr(vname(1,indxpvDiab_2+itrc-1))
          ierr=nf_def_var (ncid, vname(1,indxpvDiab_2+itrc-1)(1:lvar),
     &                     NF_FOUT, 3, r2dgrd, viddiags_pvDiab_2(itrc))
# ifdef AVRH
          text='averaged '/ /vname(2,indxpvDiab_2+itrc-1)
# else
          text=vname(2,indxpvDiab_2+itrc-1)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, viddiags_pvDiab_2(itrc), 'long_name', 
     &                          lvar, text(1:lvar))

          lvar=lenstr(vname(3,indxpvDiab_2+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_pvDiab_2(itrc), 'units', lvar,
     &                          vname(3,indxpvDiab_2+itrc-1)(1:lvar))
#ifdef MASK_LAND_DATA 
          ierr=nf_put_att_FOUT (ncid, viddiags_pvDiab_2(itrc), '_FillValue',
     &                                   NF_FOUT, 1, spv_set)
#endif



!------------------------------
!
! indxpvDiab_3
!          
          lvar=lenstr(vname(1,indxpvDiab_3+itrc-1))
          ierr=nf_def_var (ncid, vname(1,indxpvDiab_3+itrc-1)(1:lvar),
     &                     NF_FOUT, 3, r2dgrd, viddiags_pvDiab_3(itrc))
# ifdef AVRH
          text='averaged '/ /vname(2,indxpvDiab_3+itrc-1)
# else
          text=vname(2,indxpvDiab_3+itrc-1)
# endif
          lvar=lenstr(text)
          ierr=nf_put_att_text (ncid, viddiags_pvDiab_3(itrc), 'long_name', 
     &                          lvar, text(1:lvar))

          lvar=lenstr(vname(3,indxpvDiab_3+itrc-1))
          ierr=nf_put_att_text (ncid, viddiags_pvDiab_3(itrc), 'units', lvar,
     &                          vname(3,indxpvDiab_3+itrc-1)(1:lvar))
#ifdef MASK_LAND_DATA 
          ierr=nf_put_att_FOUT (ncid, viddiags_pvDiab_3(itrc), '_FillValue',
     &                                   NF_FOUT, 1, spv_set)
#endif



#endif

#endif

      enddo
        endif









!
! Leave definition mode.
! ----- ---------- -----
!
        ierr=nf_enddef(ncid)
        write(stdout,'(6x,4A,1x,A,i4)') 
#ifdef AVRH
     &        'DEF_DIAGS_PV_AVG - Created ',
#else      
     &        'DEF_DIAGS_PV - Created ',
#endif /*AVRH*/
     &                'new netCDF file ''', 
     &                 ncname(1:lstr), '''.'
     &                 MYID

!
! Open existing file and prepare for appending data:   Inquire about
! ==== ======== ==== === ======= === ========= ====    dimensions and
! variables. Check for consistency with model dimensions. In the case
! when file is rejected (whether it cannot be opened, or something is
! wrong with its dimensions) create a new file.  After that verify
! that all necessary variables are already defined, find and save
! their netCDF IDs.
!
      elseif (ncid.eq.-1) then
        ierr=nf_open (ncname(1:lstr), nf_write, ncid)
        if (ierr. eq. nf_noerr) then
          write(stdout,'(1x,4A,1x,A,I3,A,I4)') 'Opened existing ',
     &                 'netCDF file ''',   ncname(1:lstr),   '''',
     &                 'ncid =', ncid  MYID
          ierr=checkdims (ncid, ncname(1:lstr), lstr, rec)
          if (ierr .eq. nf_noerr) then
            if (rec_per_file.eq.0) then
              ierr=rec+1 - total_rec
            else
              ierr=rec+1 - (1+mod(total_rec-1, rec_per_file))
            endif
            if (ierr.gt.0) then
              MPI_master_only write( stdout,
     &                 '(/1x,A,I5,1x,A/8x,3A,I5,/8x,A,I5,1x,A/)'
     &           ) 'WARNING: def_his/avg: Actual number of records',
     &               rec,  'in netCDF file',  '''',  ncname(1:lstr),
     &             ''' exceeds the record number from restart data',
     &             rec+1-ierr,'/', total_rec,', restart is assumed.'
              rec=rec-ierr
            elseif (rec_per_file.eq.0) then
              total_rec=rec+1           ! <-- set to the next record
#if defined MPI & !defined PARALLEL_FILES
              if (mynode.gt.0) total_rec=total_rec-1
#endif
            endif
            ierr=nf_noerr
          endif
        endif
 
        if (ierr. ne. nf_noerr) then
#if defined MPI & !defined PARALLEL_FILES
          if (mynode.eq.0) then
            create_new_file=.true.
            goto 10
          else
            write(stdout,'(/1x,5A,I4/)') '### ERROR: def_his/avg: ',
     &        'Cannot open file ''', ncname(1:lstr), '''.' MYID
            goto 99                                       !--> ERROR
          endif
#else
          create_new_file=.true.
          goto 10
#endif
        endif
!
! Find netCDF IDs of evolving model variables:
! ---- ------ --- -- -------- ----- ----------
!
! Time step indices:
!
        ierr=nf_inq_varid (ncid, 'time_step', vidTstep)
        if (ierr .ne. nf_noerr) then
          write(stdout,1) 'time_step', ncname(1:lstr)
          goto 99                                         !--> ERROR
        endif
!
! Time.
!
        lvar=lenstr(vname(1,indxTime))
        ierr=nf_inq_varid (ncid,vname(1,indxTime)(1:lvar),vidTime)
        if (ierr .ne. nf_noerr) then
          write(stdout,1) vname(1,indxTime)(1:lvar), ncname(1:lstr)
          goto 99                                         !--> ERROR
        endif




! momentum flux diagnostics variables :
! -------- ---- ----------- ----------


        if (wrt(2)) then
!
        do itrc=1,2  ! use itrc indice = 1 for u_momentum
!                    !                 = 2 for v_momentum 



!  
! indxpvForc_sfc
!
         lvar=lenstr(vname(1,indxpvForc_sfc+itrc-1))
         ierr=nf_inq_varid (ncid, vname(1,indxpvForc_sfc+itrc-1)(1:lvar),
     &                      viddiags_pvForc_sfc(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxpvForc_sfc+itrc-1)(1:lvar),
     &                     ncname(1:lstr) 
           goto 99                                       !--> ERROR
         endif
!
! indxpvForc_bot
!
         lvar=lenstr(vname(1,indxpvForc_bot+itrc-1))
         ierr=nf_inq_varid (ncid, vname(1,indxpvForc_bot+itrc-1)(1:lvar),
     &                      viddiags_pvForc_bot(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxpvForc_bot+itrc-1)(1:lvar),
     &                     ncname(1:lstr) 
           goto 99                                       !--> ERROR
         endif      


#ifdef DIAGNOSTICS_PV_SUP
!------------------------------
#ifndef AVRH
!------------------------------



!  
! indxpvForc_0
!
         lvar=lenstr(vname(1,indxpvForc_0+itrc-1))
         ierr=nf_inq_varid (ncid, vname(1,indxpvForc_0+itrc-1)(1:lvar),
     &                      viddiags_pvForc_0(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxpvForc_0+itrc-1)(1:lvar),
     &                     ncname(1:lstr) 
           goto 99                                       !--> ERROR
         endif
!
! indxpvForc_1
!
         lvar=lenstr(vname(1,indxpvForc_1+itrc-1))
         ierr=nf_inq_varid (ncid, vname(1,indxpvForc_1+itrc-1)(1:lvar),
     &                      viddiags_pvForc_1(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxpvForc_1+itrc-1)(1:lvar),
     &                     ncname(1:lstr) 
           goto 99                                       !--> ERROR
         endif      


!  
! indxpvForc_2
!
         lvar=lenstr(vname(1,indxpvForc_2+itrc-1))
         ierr=nf_inq_varid (ncid, vname(1,indxpvForc_2+itrc-1)(1:lvar),
     &                      viddiags_pvForc_2(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxpvForc_2+itrc-1)(1:lvar),
     &                     ncname(1:lstr) 
           goto 99                                       !--> ERROR
         endif
!
! indxpvForc_3
!
         lvar=lenstr(vname(1,indxpvForc_3+itrc-1))
         ierr=nf_inq_varid (ncid, vname(1,indxpvForc_3+itrc-1)(1:lvar),
     &                      viddiags_pvForc_3(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxpvForc_3+itrc-1)(1:lvar),
     &                     ncname(1:lstr) 
           goto 99                                       !--> ERROR
         endif      







#endif
#endif

        enddo


        do itrc=1,NT 

!
! indxpvDiab_sfc
!
         lvar=lenstr(vname(1,indxpvDiab_sfc+itrc-1))
         ierr=nf_inq_varid (ncid, vname(1,indxpvDiab_sfc+itrc-1)(1:lvar),
     &                      viddiags_pvDiab_sfc(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxpvDiab_sfc+itrc-1)(1:lvar),
     &                     ncname(1:lstr) 
           goto 99                                       !--> ERROR
         endif
!
! indxpvDiab_bot
!
         lvar=lenstr(vname(1,indxpvDiab_bot+itrc-1))
         ierr=nf_inq_varid (ncid, vname(1,indxpvDiab_bot+itrc-1)(1:lvar),
     &                      viddiags_pvDiab_bot(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxpvDiab_bot+itrc-1)(1:lvar),
     &                     ncname(1:lstr)
           goto 99                                       !--> ERROR
         endif


#ifdef DIAGNOSTICS_PV_SUP
!------------------------------
#ifndef AVRH
!------------------------------



!
! indxpvDiab_0
!
         lvar=lenstr(vname(1,indxpvDiab_0+itrc-1))
         ierr=nf_inq_varid (ncid, vname(1,indxpvDiab_0+itrc-1)(1:lvar),
     &                      viddiags_pvDiab_0(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxpvDiab_0+itrc-1)(1:lvar),
     &                     ncname(1:lstr) 
           goto 99                                       !--> ERROR
         endif
!
! indxpvDiab_1
!
         lvar=lenstr(vname(1,indxpvDiab_1+itrc-1))
         ierr=nf_inq_varid (ncid, vname(1,indxpvDiab_1+itrc-1)(1:lvar),
     &                      viddiags_pvDiab_1(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxpvDiab_1+itrc-1)(1:lvar),
     &                     ncname(1:lstr)
           goto 99                                       !--> ERROR
         endif


!
! indxpvDiab_2
!
         lvar=lenstr(vname(1,indxpvDiab_2+itrc-1))
         ierr=nf_inq_varid (ncid, vname(1,indxpvDiab_2+itrc-1)(1:lvar),
     &                      viddiags_pvDiab_2(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxpvDiab_2+itrc-1)(1:lvar),
     &                     ncname(1:lstr) 
           goto 99                                       !--> ERROR
         endif
!
! indxpvDiab_3
!
         lvar=lenstr(vname(1,indxpvDiab_3+itrc-1))
         ierr=nf_inq_varid (ncid, vname(1,indxpvDiab_3+itrc-1)(1:lvar),
     &                      viddiags_pvDiab_3(itrc))
         if (ierr .ne. nf_noerr) then
           write(stdout,1) vname(1,indxpvDiab_3+itrc-1)(1:lvar),
     &                     ncname(1:lstr)
           goto 99                                       !--> ERROR
         endif






#endif
#endif


        enddo
       endif
!




        write(*,'(6x,2A,i4,1x,A,i4)')    'def_his/avg: Opened ',
     &              'existing file from record =', rec MYID 
#if defined MPI & !defined PARALLEL_FILES
      else
        ierr=nf_open (ncname(1:lstr), nf_write, ncid)
        if (ierr .eq. nf_noerr) then
          write(stdout,'(1x,4A,I3)')        'def_his/avg: Opened ',
     &      'existing file ''', ncname(1:lstr), ''' ncid =', ncid
        else
          write(stdout,'(/1x,5A,I4/)')   '### ERROR: def_his/avg: ',
     &      'Cannot open file ''', ncname(1:lstr), '''.' MYID
          goto 99                                         !--> ERROR
        endif
#endif
      endif             !<--  create_new_file
 
      ierr=nf_set_fill (ncid, nf_nofill, lvar)
      if (ierr .ne. nf_noerr) then
        write(*,'(6x,4A,i4)')    '### ERROR: def_his/avg: Cannot ',
     &          'switch to ''nf_nofill'' mode.', nf_strerror(ierr)
     &           MYID
      endif
   1  format(/1x,'### ERROR: def_his/avg: Cannot find variable ''',
     &                       A, ''' in netCDF file ''', A, '''.'/)
#if (defined PUT_GRID_INTO_HISTORY && !defined AVRH)\
 || (defined PUT_GRID_INTO_AVERAGES && defined AVRH)
!
! Write grid variables.
!------ ---- ----------
!
        if (total_rec.le.1) call wrt_grid (ncid, ncname, lstr)
#endif
  99  return
      end
 
#undef ncname
#undef rec_per_file
#undef wrt

#undef vidTime
#undef vidTstep

#undef viddiags_pvForc_sfc
#undef viddiags_pvForc_bot
#undef viddiags_pvDiab_sfc
#undef viddiags_pvDiab_bot

#ifdef DIAGNOSTICS_PV_SUP
#ifndef AVRH
#undef viddiags_pvForc_0
#undef viddiags_pvForc_1
#undef viddiags_pvForc_2
#undef viddiags_pvForc_3
#undef viddiags_pvDiab_0
#undef viddiags_pvDiab_1
#undef viddiags_pvDiab_2
#undef viddiags_pvDiab_3
#endif
#endif

#ifndef AVRH
# define AVRH
# ifdef AVERAGES
#  include "def_diags_pv.F"
# endif
#endif


#else /*DIAGNOSTICS_PV */
      subroutine def_diags_pv_empty()
      return
      end
#endif /* DIAGNOSTICS_PV */


