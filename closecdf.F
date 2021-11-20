#include "cppdefs.h"
                                !
      subroutine closecdf       ! Close output netCDF files.
      implicit none             !
      integer ierr
#include "param.h"
#include "scalars.h"
#include "ncvars.h"
#include "netcdf.inc"
#ifdef FLOATS
# include "floats/ncvars_floats.h"
#endif
 
      if (ncidrst.ne.-1) ierr=nf_close(ncidrst)
      if (ncidhis.ne.-1) ierr=nf_close(ncidhis)
#ifdef AVERAGES
      if (ncidavg.ne.-1) ierr=nf_close(ncidavg)
#endif
#ifdef FLOATS
      if (ncidflt.ne.-1) ierr=nf_close(ncidflt)
#endif

#if defined DIAGNOSTICS_UV
      if (nciddiags_uv.ne.-1) ierr=nf_close(nciddiags_uv)
# ifdef AVERAGES
      if (nciddiags_uv_avg.ne.-1) ierr=nf_close(nciddiags_uv_avg)
# endif
#endif


#if defined DIAGNOSTICS_VRT
      if (nciddiags_vrt.ne.-1) ierr=nf_close(nciddiags_vrt)
# ifdef AVERAGES
      if (nciddiags_vrt_avg.ne.-1) ierr=nf_close(nciddiags_vrt_avg)
# endif
#endif

#if defined DIAGNOSTICS_TS
      if (nciddiags_ts.ne.-1) ierr=nf_close(nciddiags_ts)
# ifdef AVERAGES
      if (nciddiags_ts_avg.ne.-1) ierr=nf_close(nciddiags_ts_avg)
# endif
#endif


#if defined DIAGNOSTICS_PV
      if (nciddiags_pv.ne.-1) ierr=nf_close(nciddiags_pv)
# ifdef AVERAGES
      if (nciddiags_pv_avg.ne.-1) ierr=nf_close(nciddiags_pv_avg)
# endif
#endif

#if defined DIAGNOSTICS_EDDY
      if (nciddiags_eddy.ne.-1) ierr=nf_close(nciddiags_eddy)
# ifdef AVERAGES
      if (nciddiags_eddy_avg.ne.-1) ierr=nf_close(nciddiags_eddy_avg)
# endif
#endif



      write(stdout,'(/1x,2A,I5,4(/8x,A,I5)/)') 'MAIN - '
     &  ,'number of records written into history  file(s):', nrechis
     &  ,'number of records written into restart  file(s):', nrecrst
#ifdef AVERAGES
     &  ,'number of records written into averages file(s):', nrecavg
#endif
#ifdef FLOATS
     &  ,'number of records written into floats file(s):', nrecflt
#endif
#ifdef STATIONS
     &  ,'number of records written into averages file(s):', nrecsta
#endif
 
      if (may_day_flag.eq.0) then
        write(stdout,'(/1x,A/)') 'MAIN: DONE'
      elseif (may_day_flag.eq.1) then
        write(stdout,'(/1x,A/)') 'MAIN: Abnormal termination : BLOWUP'
      elseif (may_day_flag.eq.2) then
        write(stdout,'(/1x,A/)')
     &               '### ERROR: Cause of termination: netCDF input.'
      elseif (may_day_flag.eq.3) then
        write(stdout,'(/1x,A/)')
     &              '### ERROR: Cause of termination: netCDF output.'
      elseif (may_day_flag.eq.4) then
        write(stdout,'(/1x,A/)') '### ERROR: Cannot open netCDF file.'
      elseif (may_day_flag.eq.5) then
        write(stdout,'(/1x,A/)')
     &                     '### ERROR: read_inp: Error in input file.'
      elseif (may_day_flag.eq.6) then
        write(stdout,'(/1x,A/)')
     &          '### ERROR: read_inp: An input file is not available.'
      elseif (may_day_flag.eq.7) then
        write(stdout,'(/1x,A/)') '### ERROR: KPP algoritm failure.'
      elseif (may_day_flag.eq.8) then
        write(stdout,'(/1x,A/)') '### ERROR: Wrong choice of tiling.'
      endif
      return
      end
 
