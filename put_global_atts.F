#include "cppdefs.h"
 
      subroutine put_global_atts (ncid, ierr)
      implicit none
      integer ncid, ierr, nf_ftype, lvar,lenstr
#include "param.h"
#include "scoord.h"
#include "scalars.h"
#include "ncvars.h"
#include "netcdf.inc"
#include "strings.h"
#ifdef FLOATS
# include "floats/init_floats.h"
# include "floats/ncvars_floats.h"
      real tempreal(i_floats)
      integer tempint(i_floats),i
#endif
#ifndef NONLIN_EOS
# include "eos_vars.h"
#endif
 
#if defined MPI && defined PARALLEL_FILES
      integer*4 ibuff(4), npart 
!
! Put global attribute 'partition' which identifies subdomain
! within the processor grid individually for each file.
!
c**    ibuff(1)=inode
c**    ibuff(2)=jnode
c**    ibuff(3)=NP_XI
c**    ibuff(4)=NP_ETA
c**    ibuff(5)=iSW_corn
c**    ibuff(6)=jSW_corn


      npart=4
      ibuff(1)=mynode
      ibuff(2)=nnodes
      if (WESTERN_MPI_EDGE) then
        ibuff(3)=iSW_corn+iwest
      else
        ibuff(3)=iSW_corn+iwest+1
      endif
      if (SOUTHERN_MPI_EDGE) then
        ibuff(4)=jSW_corn+jsouth
      else
        ibuff(4)=jSW_corn+jsouth+1
      endif
      

      ierr=nf_put_att_int (ncid, nf_global, 'partition',
     &                             nf_int, npart, ibuff)
#endif
!
! Decide about output file type and precision for floating point
! variables (restart file always has the same precision as in the
! code, typically double precision, while all others can be made
! single precision.
!
      if (ncid.eq.ncidrst) then
        nf_ftype=NF_FTYPE
      else
        nf_ftype=NF_FOUT
      endif
 
      if (ncid .eq. ncidrst) then
        ierr=nf_put_att_text (ncid, nf_global, 'type',  17,
     &                                 'ROMS restart file')
      elseif (ncid .eq. ncidhis) then
        ierr=nf_put_att_text (ncid, nf_global, 'type',  17,
     &                                 'ROMS history file')
#ifdef AVERAGES
      elseif (ncid .eq. ncidavg) then
        ierr=nf_put_att_text (ncid, nf_global, 'type',  18,
     &                                'ROMS averages file')
#endif
#ifdef STATIONS
      elseif (ncid .eq. ncidsta) then
        ierr=nf_put_att_text(ncid,  nf_global, 'type', 18,
     &                                'ROMS stations file')
#endif
#ifdef FLOATS
      elseif (ncid .eq. ncidflt) then
        ierr=nf_put_att_text(ncid,  nf_global, 'type', 24,
     &                          'ROMS floats history file')
#endif
      endif
 
      lvar=lenstr(title)
      ierr=nf_put_att_text(ncid, nf_global, 'title',   lvar,
     &                                         title(1:lvar))

c      lvar=lenstr(date_str)
c      ierr=nf_put_att_text(ncid, nf_global, 'date',    lvar,
c     &                                      date_str(1:lvar))
c      lvar=lenstr(rstname)
c      ierr=nf_put_att_text(ncid, nf_global, 'rst_file',lvar,
c     &                                       rstname(1:lvar))
c      lvar=lenstr(hisname)
c      ierr=nf_put_att_text(ncid, nf_global, 'his_file',lvar,
c     &                                       hisname(1:lvar))
#ifdef AVERAGES
c      lvar=lenstr(avgname)
c      ierr=nf_put_att_text(ncid, nf_global, 'avg_file',lvar,
c     &                                       avgname(1:lvar))
#endif
#ifdef STATIONS
c      lvar=lenstr(staname)
c      ierr=nf_put_att_text(ncid, nf_global, 'sta_file',lvar,
c     &                                       staname(1:lvar))
c      lvar=lenstr(sposnam)
c      ierr=nf_put_att_text(ncid, nf_global, 'spos_file',lvar,
c     &                                       sposnam(1:lvar))
#endif
#ifdef FLOATS
c      lvar=lenstr(fltname)
c      ierr=nf_put_att_text(ncid, nf_global, 'flt_file',lvar,
c     &                                       fltname(1:lvar))
c      lvar=lenstr(fposnam)
c      ierr=nf_put_att_text(ncid, nf_global, 'fpos_file',lvar,
c     &                                       fposnam(1:lvar))
      if (ncid.eq.ncidflt) then
# ifdef FLOATS_GLOBAL_ATTRIBUTES
        do i=1,i_floats
          tempreal(i)=Ft0(i)
        enddo
        ierr=nf_put_att_FTYPE(ncid, nf_global, 'Ft0',nf_ftype,
     &                                     i_floats, tempreal)
        do i=1,i_floats
          tempreal(i)=Fx0(i)
        enddo
        ierr=nf_put_att_FTYPE(ncid, nf_global, 'Fx0',nf_ftype,
     &                                     i_floats, tempreal)
        do i=1,i_floats
          tempreal(i)=Fy0(i)
        enddo
        ierr=nf_put_att_FTYPE(ncid, nf_global, 'Fy0',nf_ftype,
     &                                     i_floats, tempreal)
        do i=1,i_floats
          tempreal(i)=Fz0(i)
        enddo
        ierr=nf_put_att_FTYPE(ncid, nf_global, 'Fz0',nf_ftype,
     &                                     i_floats, tempreal)
        do i=1,i_floats
          tempint(i)=Fgrd(i)
        enddo
        ierr=nf_put_att_int(ncid, nf_global, 'Fgrd',nf_int,
     &                                   i_floats, tempint)
        do i=1,i_floats
          tempint(i)=Fcoor(i)
        enddo
        ierr=nf_put_att_int(ncid, nf_global, 'Fcoor',nf_int,
     &                                    i_floats, tempint)
        do i=1,i_floats
          tempint(i)=Ftype(i)
        enddo
        ierr=nf_put_att_int(ncid, nf_global, 'Ftype',nf_int,
     &                                    i_floats, tempint)
        do i=1,i_floats
          tempint(i)=Fcount(i)
        enddo
        ierr=nf_put_att_int(ncid, nf_global, 'Fcount',nf_int,
     &                                     i_floats, tempint)
        do i=1,i_floats
          tempreal(i)=Fdt(i)
        enddo
        ierr=nf_put_att_FTYPE(ncid, nf_global, 'Fdt',nf_ftype,
     &                                     i_floats, tempreal)
        do i=1,i_floats
          tempreal(i)=Fdx(i)
        enddo
        ierr=nf_put_att_FTYPE(ncid, nf_global, 'Fdx',nf_ftype,
     &                                     i_floats, tempreal)
        do i=1,i_floats
          tempreal(i)=Fdy(i)
        enddo
        ierr=nf_put_att_FTYPE(ncid, nf_global, 'Fdy',nf_ftype,
     &                                     i_floats, tempreal)
        do i=1,i_floats
          tempreal(i)=Fdz(i)
        enddo
        ierr=nf_put_att_FTYPE(ncid, nf_global, 'Fdz',nf_ftype,
     &                                     i_floats, tempreal)
# endif
      endif
#endif /* FLOATS */

#ifndef ANA_GRID
      lvar=lenstr(grdname)
      ierr=nf_put_att_text(ncid, nf_global, 'grid_file',lvar,
     &                                        grdname(1:lvar))
#endif
#ifndef ANA_INITIAL
      lvar=lenstr(ininame)
      ierr=nf_put_att_text(ncid, nf_global, 'init_file',lvar,
     &                                        ininame(1:lvar))
#endif
#if !defined ANA_SMFLUX || !defined ANA_STFLUX \
 || (defined SG_BBL96 && !defined ANA_BSEDIM)  \
 || (defined SG_BBL96 && !defined ANA_WWAVE)   \
 || (defined SALINITY && !defined ANA_SSFLUX)  \
 || (defined LMD_KPP  && !defined ANA_SRFLUX)
c      lvar=lenstr(frcname)
c      ierr=nf_put_att_text(ncid, nf_global,'frc_file', lvar,
c     &                                       frcname(1:lvar))
#endif

#ifdef SOLVE3D
!
! S-coordinate control parameters "theta_s", "theta_b", and "hc",
! and stretching curves "Cs_w", "Cs_r" at W- and RHO-points.
!

# ifdef  NEW_S_COORD
      ierr=nf_put_att_text (ncid, nf_global, 'VertCoordType',3,'NEW')
# endif
      ierr=nf_put_att_FTYPE (ncid, nf_global, 'theta_s', nf_ftype, 1,
     &                                                       theta_s)
      ierr=nf_put_att_FTYPE (ncid, nf_global, 'theta_b', nf_ftype, 1,
     &                                                       theta_b)
      ierr=nf_put_att_FTYPE(ncid, nf_global, 'hc',nf_ftype, 1, hc)

      ierr=nf_put_att_FTYPE(ncid,nf_global, 'Cs_w',nf_ftype,N+1,Cs_w)

      ierr=nf_put_att_FTYPE(ncid,nf_global, 'Cs_r', nf_ftype, N,Cs_r)
#endif
!
! Time stepping parameters.
!
      ierr=nf_put_att_int(ncid,nf_global,'ntimes',  nf_int,1,ntimes)
      ierr=nf_put_att_int(ncid,nf_global,'ndtfast', nf_int,1,ndtfast)
      ierr=nf_put_att_FTYPE(ncid,nf_global,'dt',    nf_ftype, 1,  dt)
      ierr=nf_put_att_FTYPE(ncid,nf_global,'dtfast',nf_ftype, 1,
     &                                                        dtfast)
      ierr=nf_put_att_int  (ncid,nf_global,'nwrt',  nf_int,  1, nwrt)
#ifdef AVERAGES
      ierr=nf_put_att_int  (ncid,nf_global,'ntsavg',nf_int, 1,ntsavg)


 
      ierr=nf_put_att_int  (ncid,nf_global,'navg',  nf_int, 1,  navg)
#endif
#ifdef STATIONS
      ierr=nf_put_att_int  (ncid,nf_global,'nsta',  nf_int, 1,  nsta)
#endif
!
! Horizontal viscosity and mixing coefficients.
!
#ifdef UV_VIS2
      ierr=nf_put_att_FTYPE(ncid,nf_global,'visc2',nf_ftype,1,visc2)
      ierr=nf_put_att_text (ncid,nf_global,'visc2_units',15,
     &                                             'meter2 second-1')
#endif
#ifdef SOLVE3D
# ifdef TS_DIF2
      ierr=nf_put_att_FTYPE(ncid,nf_global,'tnu2',nf_ftype, 1,tnu2)
      ierr=nf_put_att_text (ncid,nf_global,'tnu2_units',15,
     &                                             'meter2 second-1')
# endif
# ifdef MY25_MIXING
#  ifdef Q2_DIF2
      ierr=nf_put_att_FTYPE(ncid,nf_global,'q2nu2',nf_ftype, 1,q2nu2)
      ierr=nf_put_att_text (ncid,nf_global,'q2nu2_units',15,
     &                                             'meter2 second-1')
#  endif
#  ifdef Q2_DIF4
      ierr=nf_put_att_FTYPE(ncid,nf_global,'q2nu4',nf_ftype, 1,q2nu4)
      ierr=nf_put_att_text (ncid,nf_global,'q2nu4_units',15,
     &                                             'meter4 second-1')
#  endif
# endif
# if !defined LMD_MIXING && !defined BVF_MIXING\
  && !defined MY2_MIXING && !defined MY25_MIXING\
                         && !defined PP_MIXING
!
! Background vertical viscosity and tracer mixing coefficients.
!
      ierr=nf_put_att_FTYPE(ncid,nf_global,'Akv_bak',nf_ftype, 1,
     &                                                      Akv_bak)
      ierr=nf_put_att_text (ncid,nf_global,'Akv_bak',51,
     &         'background vertical mixing coefficient for momentum')
      ierr=nf_put_att_text (ncid,nf_global,'Akv_bak_units',15,
     &                                             'meter2 second-1')
 
      ierr=nf_put_att_FTYPE(ncid,nf_global,'Akt_bak',nf_ftype, NT,
     &                                                       Akt_bak)
      ierr=nf_put_att_text (ncid,nf_global,'Akt_bak_units', 15,
     &                                             'meter2 second-1')
# endif
# ifdef MY25_MIXING
      ierr=nf_put_att_FTYPE(ncid,nf_global,'Akq_bak',nf_ftype,
     &                    1,Akq_bak)
      ierr=nf_put_att_text (ncid,nf_global,'Akq_bak_units',15,
     &                                             'meter2 second-1')
# endif
#endif
!
! Bottom drag coefficients.
!
      ierr=nf_put_att_FTYPE(ncid,nf_global,'rdrg',nf_ftype,1,rdrg)
      ierr=nf_put_att_text (ncid,nf_global,'rdrg_units',14,
     &                                              'meter second-1')
      ierr=nf_put_att_FTYPE(ncid,nf_global,'rdrg2',nf_ftype,1,rdrg2)
      ierr=nf_put_att_text (ncid,nf_global,'rdrg2_units',14,
     &                                              'nondimensional')
#ifdef SOLVE3D
!
! Equation of State parameters.
!
      ierr=nf_put_att_FTYPE(ncid,nf_global, 'rho0',nf_ftype, 1,rho0)
      ierr=nf_put_att_text (ncid,nf_global, 'rho0_units', 16,
     &                                            'kilogram meter-3')
# ifndef NONLIN_EOS
      ierr=nf_put_att_FTYPE(ncid,nf_global, 'R0',  nf_ftype, 1, R0)
      ierr=nf_put_att_text (ncid,nf_global, 'R0_units', 16,
     &                                            'kilogram meter-3')

      ierr=nf_put_att_FTYPE(ncid,nf_global, 'T0',  nf_ftype, 1, T0)
      ierr=nf_put_att_text (ncid,nf_global, 'T0_units', 11,
     &                                                 'deg Celsius')

      ierr=nf_put_att_FTYPE(ncid,nf_global,'Tcoef',nf_ftype, 1,Tcoef)
      ierr=nf_put_att_text (ncid,nf_global,'Tcoef_units', 13,
     &                                               'deg Celsius-1')

#  ifdef SALINITY
      ierr=nf_put_att_FTYPE(ncid,nf_global, 'S0',  nf_ftype, 1, S0)
      ierr=nf_put_att_text (ncid,nf_global, 'S0_units',     3, 'PSU')
 
      ierr=nf_put_att_FTYPE(ncid,nf_global,'Scoef',nf_ftype, 1,Scoef)
      ierr=nf_put_att_text (ncid,nf_global,'Scoef_units', 5, 'PSU-1')
#  endif
# endif
#endif /* SOLVE3D */
!
! Slipperiness parameters.
!
      ierr=nf_put_att_FTYPE(ncid,nf_global,'gamma2',nf_ftype, 1,
     &                                                      gamma2)
!
! List of Source Codes and Activated CPP-switches
!
      lvar=lenstr(srcs)
      ierr=nf_put_att_text (ncid,nf_global, 'SRCS', lvar,
     &                                       srcs(1:lvar))
      lvar=lenstr(cpps)
      ierr=nf_put_att_text (ncid,nf_global, 'CPPS', lvar,
     &                                       cpps(1:lvar))
      return
      end
 
