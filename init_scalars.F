#include "cppdefs.h"
 
       subroutine init_scalars (ierr)
!
! Set initial values for  globally accessable (stored in common
! blocks) scalar variables. Typically this is associated with
! model utility switches, i/o contlol, time indices, global
! summation counters, etc; NONE of the physical parameters of
! the model is set here.
!
      implicit none
      integer ierr, i,j,itrc, lvar, lenstr
      character*20 nametrc, unitt
      character*60 vname1, vname3      
#define ALL_DATA
#include "param.h"
#include "scalars.h"
#include "comm_vars.h"
#include "forces.h"
#include "sources.h"
#include "mixing.h"
#include "ncvars.h"
#include "boundary.h"
#include "diag_vars.h"
# ifdef DIAGNOSTICS_UV
#include "diags_uv.h"
# endif

# ifdef DIAGNOSTICS_VRT
#include "diags_vrt.h"
# endif

# ifdef DIAGNOSTICS_TS
#include "diags_ts.h"
# endif

# ifdef DIAGNOSTICS_PV
#include "diags_pv.h"
# endif

# ifdef DIAGNOSTICS_EDDY
#include "diags_eddy.h"
# endif

c      integer size, kmp_get_stacksize_s

      call system('uname -nmr')

c      size=kmp_get_stacksize_s()
      

cc      write(*,*) 'kmp_stacksize =', size

 
      start_time=0.                  ! Set time and evolving
      nfast=1                        ! time indices to their
      iic=0                          ! default initial values.
      kstp=1
      knew=1
      ntstart=1
#ifdef SOLVE3D
      iif=1
      nstp=1
      nrhs=1
      nnew=1
#endif
#ifdef FLOATS
      nfp1=0     ! set time step indices for
      nf=3       ! 5th order AB4AM4 corrected
      nfm1=2     ! scheme
      nfm2=1
      nfm3=0
#endif
                                     ! Reset intra-thread
      may_day_flag=0                 ! communication signal
      trd_count=0                    ! variables and counters.
      trd_count=0
      do i=1,4
        CPU_time_ALL(i)=0.
      enddo
      do i=1,16
        barr_count(i)=0
      enddo

 
      mcheck=0
      first_time=0
      tile_count=0                   ! Reset global sumation
      bc_count=0                     ! variables and counters
      volume=0.                      ! used for volume/energy
      avke=0.                        ! diagnostics and volume
      avpe=0.                        ! conservation constraint.
      prev_ke=0.
      avke2b=0.
      avke3bc=0.
      avkp=0.
      v2d_max=0.
      v3d_max=0.
      v3bc_max=0.

      Cu_Adv   = 0.
      Cu_W     = 0.
      i_cx_max = 0
      j_cx_max = 0
      k_cx_max = 0
      

 
      hmin=+1.E+20                   ! Set extremal values to
      hmax=-1.E+20                   ! large positive or negative
      grdmin=+1.E+20                 ! unrealistic values.
      grdmax=-1.E+20
      Cu_min=+1.E+20
      Cu_max=-1.E+20
#ifdef SOLVE3D
      rx0=-1.E+20
      rx1=-1.E+20
#endif
 
      bc_crss=QuadZero
#ifdef OBC_VOLCONS
      bc_flux=QuadZero
      ubar_xs=QuadZero
#endif

      ncidrst=-1                     ! Set netCDF file IDs to
      nrecrst=0                      ! unrealistic negative values
      ncidhis=-1                     ! to indicate their initial
      nrechis=0                      ! closed status.
#ifdef AVERAGES
      ncidavg=-1                     ! Also reset record counters for
      nrecavg=0                      ! restart, history, average, and
#endif




#ifdef DIAGNOSTICS_UV
      nciddiags_uv=-1
      nrecdiags_uv=0
# ifdef AVERAGES
      nciddiags_uv_avg=-1
      nrecdiags_uv_avg=0
# endif
# endif

#ifdef DIAGNOSTICS_VRT
      nciddiags_vrt=-1
      nrecdiags_vrt=0
# ifdef AVERAGES
      nciddiags_vrt_avg=-1
      nrecdiags_vrt_avg=0
# endif
# endif

#ifdef DIAGNOSTICS_TS
      nciddiags_ts=-1
      nrecdiags_ts=0
# ifdef AVERAGES
      nciddiags_ts_avg=-1
      nrecdiags_ts_avg=0
# endif
# endif

#ifdef DIAGNOSTICS_PV
      nciddiags_pv=-1
      nrecdiags_pv=0
# ifdef AVERAGES
      nciddiags_pv_avg=-1
      nrecdiags_pv_avg=0
# endif
# endif


#ifdef DIAGNOSTICS_EDDY
      nciddiags_eddy=-1
      nrecdiags_eddy=0
# ifdef AVERAGES
      nciddiags_eddy_avg=-1
      nrecdiags_eddy_avg=0
# endif
# endif


#ifdef STATIONS
      ncidsta=-1                     ! stations netCDF output files.
      nrecsta=0
#endif

#ifdef BULK_FLUX     
      do i=1,max_blk_file
        ncidbulk(i)=-1
      enddo  
#else      
      do i=1,max_frc_file
        ncidfrc(i)=-1
      enddo      
#endif


#ifndef ANA_SMFLUX
      sms_file_id=0
#endif

#ifdef BULK_FLUX
      rad_file_id=0
      tra_file_id=0
      prec_file_id=0
      wnd_file_id=0
#endif

#ifdef SOLVE3D
# if !defined ANA_STFLUX || !defined ANA_SSFLUX
      do i=1,NT
        stf_file_id(i)=0
      enddo
# endif
# ifndef ANA_SRFLUX
      srf_file_id=0
# endif
# if defined QCORRECTION && !defined ANA_SST
      sst_file_id=0
      sss_file_id=0
      dqdt_file_id=0
# endif
#endif
      ncidclm=-1
#ifndef ANA_BRY
       bry_id=-1
#endif
#ifdef PSOURCE
       do j=1,NT                     ! Initialize switch
         do i=1,Msrc                 ! to apply tracer
           Lsrc(i,j)=.false.         ! sources and sinks.
         enddo
       enddo
#endif
 
      call get_date (date_str)
!
! Define names of variables in NetCDF output files.
! Inner dimension is for variable type and outer is for variable
! attributes (name, long-name, units, field).
!
c--#define CR
CR      write(stdout,'(1x,A,I3)') 'indxTime =',indxTime
      vname(1,indxTime)='ocean_time                                '
      vname(2,indxTime)='time since initialization                 '
      vname(3,indxTime)='second                                    '
 
CR      write(stdout,'(1x,A,I3)') 'indxZ =',indxZ
      vname(1,indxZ)='zeta                                      '
      vname(2,indxZ)='free-surface elevation                    '
      vname(3,indxZ)='meter                                     '
 
CR      write(stdout,'(1x,A,I3)') 'indxUb =',indxUb
      vname(1,indxUb)='ubar                                      '
      vname(2,indxUb)='vertically integrated u-momentum component'
      vname(3,indxUb)='meter second-1                            '
 
CR      write(stdout,'(1x,A,I3)') 'indxVb = ',indxVb
      vname(1,indxVb)='vbar                                      '
      vname(2,indxVb)='vertically integrated v-momentum component'
      vname(3,indxVb)='meter second-1                            '
 
#ifdef SOLVE3D
CR      write(stdout,'(1x,A,I3)') 'indxU = ',indxU
      vname(1,indxU)='u                                         '
      vname(2,indxU)='u-momentum component                      '
      vname(3,indxU)='meter second-1                            '
 
CR      write(stdout,'(1x,A,I3)') 'indxV = ',indxV
      vname(1,indxV)='v                                         '
      vname(2,indxV)='v-momentum component                      '
      vname(3,indxV)='meter second-1                            '
 
CR      write(stdout,'(1x,A,I3)') 'indxT = ',indxT
      vname(1,indxT)='temp                                      '
      vname(2,indxT)='potential temperature                     '
      vname(3,indxT)='Celsius                                   '
 
# ifdef SALINITY
CR      write(stdout,'(1x,A,I3)') 'indxS = ',indxS
      vname(1,indxS)='salt                                      '
      vname(2,indxS)='salinity                                  '
      vname(3,indxS)='PSU                                       '
# endif
# ifdef PASSIVE_TRACER
CR      write(stdout,'(1x,A,I3)') 'indxTPAS = ',indxTPAS
      vname(1,indxTPAS)='tpas                                      '
      vname(2,indxTPAS)='passive tracer                            '
      vname(3,indxTPAS)='no unit                                   '
# endif 
CR      write(stdout,'(1x,A,I3)') 'indxO =',indxO
      vname(1,indxO)='omega'
      vname(2,indxO)='S-coordinate vertical velocity'
      vname(3,indxO)='meter second-1'
 
CR      write(stdout,'(1x,A,I3)') 'indxW =',indxW
      vname(1,indxW)='w                                         '
      vname(2,indxW)='vertical momentum component               '
      vname(3,indxW)='meter second-1                            '

#if defined ADV_ISONEUTRAL & WRITE_SXSY
      vname(1,indxSx)='Sx                                        '
      vname(2,indxSx)='Isopycnal Slopes in xi-direction          '
      vname(3,indxSx)='nondimesional                             '

      vname(1,indxSy)='Sy                                        '
      vname(2,indxSy)='Isopycnal Slopes in eta-direction         '
      vname(3,indxSy)='nondimesional                             ' 
# endif
 
CR      write(stdout,'(1x,A,I3)') 'indxR =',indxR
      vname(1,indxR)='rho                                       '
      vname(2,indxR)='density anomaly                           '
      vname(3,indxR)='kilogram meter-3                          '
 
CR      write(stdout,'(1x,A,I3)') 'indxAkv =',indxAkv
      vname(1,indxAkv)='AKv                                       '
      vname(2,indxAkv)='vertical viscosity coefficient            '
      vname(3,indxAkv)='meter2 second-1                           '
 
CR      write(stdout,'(1x,A,I3)') 'indxAkt =',indxAkt
      vname(1,indxAkt)='AKt                                       '
      vname(2,indxAkt)='temperature vertical diffusion coefficient'
      vname(3,indxAkt)='meter2 second-1                           '
# ifdef SALINITY
CR      write(stdout,'(1x,A,I3)') 'indxAkt =',indxAks
      vname(1,indxAks)='AKs                                       '
      vname(2,indxAks)='salinity vertical diffusion coefficient   '
      vname(3,indxAks)='meter2 second-1                           '
# endif
# ifdef LMD_KPP
CR      write(stdout,'(1x,A,I3)') 'indxHbl =',indxHbl
      vname(1,indxHbl)='hbls'
      vname(2,indxHbl)='Thickness of KPP surface boundary layer'
      vname(3,indxHbl)='meter'
# ifdef LMD_BKPP
CR      write(stdout,'(1x,A,I3)') 'indxHbbl =',indxHbbl
      vname(1,indxHbbl)='hbbls'
      vname(2,indxHbbl)='Thickness of KPP bottom boundary layer'
      vname(3,indxHbbl)='meter'
#endif      
#endif

# ifdef BULK_FLUX
CR      write(stdout,'(1x,A,I3)') 'indxSustr_blk =',indxSustr_blk
      vname(1,indxSustr_blk)='sustr_blk                           '
      vname(2,indxSustr_blk)='surface u-momentum stress           '
      vname(3,indxSustr_blk)='Newton meter-2                      '
CR      write(stdout,'(1x,A,I3)') 'indxSvstr_blk =',indxSvstr_blk
      vname(1,indxSvstr_blk)='svstr_blk                           '
      vname(2,indxSvstr_blk)='surface v-momentum stress           '
      vname(3,indxSvstr_blk)='Newton meter-2                      '
CR      write(stdout,'(1x,A,I3)') 'indxShflx_rsw =',indxShflx_rsw
      vname(1,indxShflx_rsw)='shflx_rsw                           '
      vname(2,indxShflx_rsw)='Short-wave surface radiation        '
      vname(3,indxShflx_rsw)='Watts meter-2                       '
CR      write(stdout,'(1x,A,I3)') 'indxShflx_rlw =',indxShflx_rlw
      vname(1,indxShflx_rlw)='shflx_rlw                           '
      vname(2,indxShflx_rlw)='Long-wave surface radiation         '
      vname(3,indxShflx_rlw)='Watts meter-2                       '
CR      write(stdout,'(1x,A,I3)') 'indxShflx_lat =',indxShflx_lat
      vname(1,indxShflx_lat)='shflx_lat                           '
      vname(2,indxShflx_lat)='Latent surface heat flux            '
      vname(3,indxShflx_lat)='Watts meter-2                       '
CR      write(stdout,'(1x,A,I3)') 'indxShflx_sen =',indxShflx_sen
      vname(1,indxShflx_sen)='shflx_sen                           '
      vname(2,indxShflx_sen)='Sensible surface heat flux          '
      vname(3,indxShflx_sen)='Watts meter-2                       '
# endif

#endif
CR      write(stdout,'(1x,A,I3)') 'indxSUSTR =',indxSUSTR
      vname(1,indxSUSTR)='sustr                                     '
      vname(2,indxSUSTR)='surface u-momentum stress                 '
      vname(3,indxSUSTR)='Newton meter-2                            '
 
CR      write(stdout,'(1x,A,I3)') 'indxSVSTR =',indxSVSTR
      vname(1,indxSVSTR)='svstr                                     '
      vname(2,indxSVSTR)='surface v-momentum stress                 '
      vname(3,indxSVSTR)='Newton meter-2                            '
 
#ifdef SOLVE3D
CR      write(stdout,'(1x,A,I3)') 'indxSHFl =',indxSHFl
      vname(1,indxSHFl)='shflux                                    '
      vname(2,indxSHFl)='surface net heat flux                     '
      vname(3,indxSHFl)='Watts meter-2                             '
# ifdef SALINITY
CR      write(stdout,'(1x,A,I3)') 'indxSSFl=',indxSSFl
      vname(1,indxSSFl)='swflux                                    '
      vname(2,indxSSFl)='surface freshwater flux (E-P)             '
      vname(3,indxSSFl)='centimeter day-1                          '
# endif
CR      write(stdout,'(1x,A,I3)') 'indxSWRad =',indxSWRad
      vname(1,indxSWRad)='swrad                                     '
      vname(2,indxSWRad)='solar shortwave radiation                 '
      vname(3,indxSWRad)='Watts meter-2                             '
 
CR      write(stdout,'(1x,A,I3)') 'indxSST =', indxSST
      vname(1,indxSST)='SST                                       '
      vname(2,indxSST)='sea surface temperature                   '
      vname(3,indxSST)='Celsius                                   '
 
CR      write(stdout,'(1x,A,I3)') 'indxdQdSST =', indxdQdSST
      vname(1,indxdQdSST)='dQdSST                                    '
      vname(2,indxdQdSST)='surface net heat flux sensitivity to SST  '
      vname(3,indxdQdSST)='Watts meter-2 Celsius-1                   '

CR      write(stdout,'(1x,A,I3)') 'indxSSS =', indxSSS
      vname(1,indxSSS)='SSS                                       '
      vname(2,indxSSS)='sea surface salinity                      '
      vname(3,indxSSS)='PSU                                       '


#ifdef BULK_FLUX
CR      write(stdout,'(1x,A,I3)') 'indxTAIR =',indxTAIR
      vname(1,indxTAIR)='tair                                     '
      vname(2,indxTAIR)='surface air temperature 10m               '
      vname(3,indxTAIR)='Celsius                                  '
      vname(4,indxTAIR)='surface air temperature, scalar, series  '

CR      write(stdout,'(1x,A,I3)') 'indxQAIR =',indxQAIR
      vname(1,indxQAIR)='qair                                     '
      vname(2,indxQAIR)='surface air specific humidity 10m        '
      vname(3,indxQAIR)='kg/kg                                    '
      vname(4,indxQAIR)='surface relative humidity, scalar, series'

CR      write(stdout,'(1x,A,I3)') 'indxRADLW =',indxRADLW
      vname(1,indxRADLW)='radlw                                   '
      vname(2,indxRADLW)='Surface Downwelling Longwave Flux      '
      vname(3,indxRADLW)='Watts meter-2                           '
      vname(4,indxRADLW)='terrestrial longwave, scalar, series    '

CR      write(stdout,'(1x,A,I3)') 'indxRADSW =',indxRADSW
      vname(1,indxRADSW)='radsw                                   '
      vname(2,indxRADSW)='Surface Downwelling Shortwave Flux      '
      vname(3,indxRADSW)='Watts meter-2                           '
      vname(4,indxRADSW)='solar shortwave, scalar, series         '

CR      write(stdout,'(1x,A,I3)') 'indxPREC =',indxPREC
      vname(1,indxPREC)='rain                                    '
      vname(2,indxPREC)='surface precipitation                   '
      vname(3,indxPREC)='mm-2 second-1                           '
      vname(4,indxPREC)='precipitation rate, scalar, series      '

CR      write(stdout,'(1x,A,I3)') 'indxUWND =',indxUWND
      vname(1,indxUWND)='uwnd                                     '
      vname(2,indxUWND)='surface u-wind speed 10 m                '
      vname(3,indxUWND)='meter second-1                           '
      vname(4,indxUWND)='surface wind speed, scalar, series       '

CR      write(stdout,'(1x,A,I3)') 'indxVWND =',indxVWND
      vname(1,indxVWND)='vwnd                                     '
      vname(2,indxVWND)='surface v-wind speed 10 m                '
      vname(3,indxVWND)='meter second-1                           '
      vname(4,indxVWND)='surface wind speed, scalar, series       '
#endif /* BULK_FLUX */


# ifdef DIAGNOSTICS_UV
      do itrc=1,2
CR       write(stdout,'(1x,A,I3)') 'indxMXadv = ',indxMXadv+itrc-1
       lvar=lenstr(vname(1,indxU+itrc-1))
       nametrc=vname(1,indxU+itrc-1)(1:lvar)
       lvar=lenstr(nametrc)
       unitt=vname(3,indxU+itrc-1)
       write(vname1,*) trim(nametrc),'_xadv                       '
       vname(1,indxMXadv+itrc-1)=vname1
CR       write(stdout,'(20x,A)') vname(1,indxMxadv+itrc-1)
       vname(2,indxMXadv+itrc-1)='Horizontal (xi) advection term  '
       write(vname3,*) trim(unitt),' second-1                     '
       vname(3,indxMXadv+itrc-1)=vname3

CR       write(stdout,'(1x,A,I3)') 'indxMYadv = ',indxMYadv+itrc-1
       nametrc=vname(1,indxU+itrc-1)
       unitt=vname(3,indxU+itrc-1)
       write(vname1,*) trim(nametrc),'_yadv                       '
       vname(1,indxMYadv+itrc-1)=vname1
CR       write(stdout,'(20x,A)') vname(1,indxMYadv+itrc-1)
       vname(2,indxMYadv+itrc-1)='Horizontal (eta) advection term '
       write(vname3,*) trim(unitt),' second-1                     '
       vname(3,indxMYadv+itrc-1)=vname3

CR       write(stdout,'(1x,A,I3)') 'indxMVadv = ',indxMVadv+itrc-1
       nametrc=vname(1,indxU+itrc-1)
       unitt=vname(3,indxU+itrc-1)
       write(vname1,*) trim(nametrc),'_vadv                       ' 
       vname(1,indxMVadv+itrc-1)=vname1
CR       write(stdout,'(20x,A)') vname(1,indxMVadv+itrc-1)
       vname(2,indxMVadv+itrc-1)='Vertical advection term         '
       write(vname3,*) trim(unitt),' second-1                     '
       vname(3,indxMVadv+itrc-1)=vname3

CR       write(stdout,'(1x,A,I3)') 'indxMCor = ',indxMCor+itrc-1
       nametrc=vname(1,indxU+itrc-1)
       unitt=vname(3,indxU+itrc-1)
       write(vname1,*) trim(nametrc),'_cor                        '
       vname(1,indxMCor+itrc-1)=vname1
CR       write(stdout,'(20x,A)') vname(1,indxMcor+itrc-1)
       vname(2,indxMCor+itrc-1)='Coriolis term                    '
       write(vname3,*) trim(unitt),' second-1                     '
       vname(3,indxMCor+itrc-1)=vname3

CR       write(stdout,'(1x,A,I3)') 'indxMPrsgrd = ',indxMPrsgrd+itrc-1
       nametrc=vname(1,indxU+itrc-1)
       unitt=vname(3,indxU+itrc-1)
       write(vname1,*) trim(nametrc),'_Prsgrd                     '
       vname(1,indxMPrsgrd+itrc-1)=vname1
CR       write(stdout,'(20x,A)') vname(1,indxMPrsgrd+itrc-1)
       vname(2,indxMPrsgrd+itrc-1)='Pressure gradient term        '
       write(vname3,*) trim(unitt),' second-1                     '
       vname(3,indxMPrsgrd+itrc-1)=vname3

CR       write(stdout,'(1x,A,I3)') 'indxMHmix = ',indxMHmix+itrc-1
       nametrc=vname(1,indxU+itrc-1)
       unitt=vname(3,indxU+itrc-1)
       write(vname1,*) trim(nametrc),'_hmix                       '
       vname(1,indxMHmix+itrc-1)=vname1
CR       write(stdout,'(20x,A)') vname(1,indxMHmix+itrc-1)
       vname(2,indxMHmix+itrc-1)='Horizontal mixing term          '
       write(vname3,*) trim(unitt),' second-1                     '
       vname(3,indxMHmix+itrc-1)=vname3

CR       write(stdout,'(1x,A,I3)') 'indxMVmix = ',indxMVmix+itrc-1
       nametrc=vname(1,indxU+itrc-1)
       unitt=vname(3,indxU+itrc-1)
       write(vname1,*) trim(nametrc),'_vmix                       ' 
       vname(1,indxMVmix+itrc-1)=vname1
CR       write(stdout,'(20x,A)') vname(1,indxMVmix+itrc-1)
       vname(2,indxMVmix+itrc-1)='Vertical mixing term            '
       write(vname3,*) trim(unitt),' second-1                     '
       vname(3,indxMVmix+itrc-1)=vname3
       
CR       write(stdout,'(1x,A,I3)') 'indxMrate = ',indxMrate+itrc-1
       nametrc=vname(1,indxU+itrc-1)
       unitt=vname(3,indxU+itrc-1)
       write(vname1,*) trim(nametrc),'_rate                       '
       vname(1,indxMrate+itrc-1)=vname1
CR       write(stdout,'(20x,A)') vname(1,indxMrate+itrc-1)
       vname(2,indxMrate+itrc-1)='Time rate of change             '
       write(vname3,*) trim(unitt),' second-1                     '
       vname(3,indxMrate+itrc-1)=vname3



CR       write(stdout,'(1x,A,I3)') 'indxMVmix2 = ',indxMVmix2+itrc-1
       nametrc=vname(1,indxU+itrc-1)
       unitt=vname(3,indxU+itrc-1)
       write(vname1,*) trim(nametrc),'_vmix2                       ' 
       vname(1,indxMVmix2+itrc-1)=vname1
CR       write(stdout,'(20x,A)') vname(1,indxMVmix2+itrc-1)
       vname(2,indxMVmix2+itrc-1)='Vertical mixing term            '
       write(vname3,*) trim(unitt),' second-1                     '
       vname(3,indxMVmix2+itrc-1)=vname3

      enddo


# endif /* DIAGNOSTICS_UV */


# ifdef DIAGNOSTICS_TS
      do itrc=1,NT
CR       write(stdout,'(1x,A,I3)') 'indxTXadv = ',indxTXadv+itrc-1
       lvar=lenstr(vname(1,indxT+itrc-1))
       nametrc=vname(1,indxT+itrc-1)(1:lvar)
       lvar=lenstr(nametrc)
       unitt=vname(3,indxT+itrc-1)
       write(vname1,*) trim(nametrc),'_xadv                       '
       vname(1,indxTXadv+itrc-1)=vname1
CR       write(stdout,'(20x,A)') vname(1,indxTXadv+itrc-1)
       vname(2,indxTXadv+itrc-1)='Horizontal (xi) advection term  '
       write(vname3,*) trim(unitt),' second-1                     '
       vname(3,indxTXadv+itrc-1)=vname3

CR       write(stdout,'(1x,A,I3)') 'indxTYadv = ',indxTYadv+itrc-1
       nametrc=vname(1,indxT+itrc-1)
       unitt=vname(3,indxT+itrc-1)
       write(vname1,*) trim(nametrc),'_yadv                       '
       vname(1,indxTYadv+itrc-1)=vname1
CR       write(stdout,'(20x,A)') vname(1,indxTYadv+itrc-1)
       vname(2,indxTYadv+itrc-1)='Horizontal (eta) advection term '
       write(vname3,*) trim(unitt),' second-1                     '
       vname(3,indxTYadv+itrc-1)=vname3

CR       write(stdout,'(1x,A,I3)') 'indxTVadv = ',indxTVadv+itrc-1
       nametrc=vname(1,indxT+itrc-1)
       unitt=vname(3,indxT+itrc-1)
       write(vname1,*) trim(nametrc),'_vadv                       ' 
       vname(1,indxTVadv+itrc-1)=vname1
CR       write(stdout,'(20x,A)') vname(1,indxTVadv+itrc-1)
       vname(2,indxTVadv+itrc-1)='Vertical advection term         '
       write(vname3,*) trim(unitt),' second-1                     '
       vname(3,indxTVadv+itrc-1)=vname3

CR       write(stdout,'(1x,A,I3)') 'indxTHmix = ',indxTHmix+itrc-1
       nametrc=vname(1,indxT+itrc-1)
       unitt=vname(3,indxT+itrc-1)
       write(vname1,*) trim(nametrc),'_hmix                       '
       vname(1,indxTHmix+itrc-1)=vname1
CR       write(stdout,'(20x,A)') vname(1,indxTHmix+itrc-1)
       vname(2,indxTHmix+itrc-1)='Horizontal mixing term          '
       write(vname3,*) trim(unitt),' second-1                     '
       vname(3,indxTHmix+itrc-1)=vname3

CR       write(stdout,'(1x,A,I3)') 'indxTVmix = ',indxTVmix+itrc-1
       nametrc=vname(1,indxT+itrc-1)
       unitt=vname(3,indxT+itrc-1)
       write(vname1,*) trim(nametrc),'_vmix                       ' 
       vname(1,indxTVmix+itrc-1)=vname1
CR       write(stdout,'(20x,A)') vname(1,indxTVmix+itrc-1)
       vname(2,indxTVmix+itrc-1)='Vertical mixing term            '
       write(vname3,*) trim(unitt),' second-1                     '
       vname(3,indxTVmix+itrc-1)=vname3

CR       write(stdout,'(1x,A,I3)') 'indxTForc = ',indxTForc+itrc-1
       nametrc=vname(1,indxT+itrc-1)
       unitt=vname(3,indxT+itrc-1)
       write(vname1,*) trim(nametrc),'_forc                       '
       vname(1,indxTForc+itrc-1)=vname1
CR       write(stdout,'(20x,A)') vname(1,indxTForc+itrc-1)
       vname(2,indxTForc+itrc-1)='Forcing term (Q & Nudging)      '
       write(vname3,*) trim(unitt),' second-1                     '
       vname(3,indxTForc+itrc-1)=vname3

CR       write(stdout,'(1x,A,I3)') 'indxTrate = ',indxTrate+itrc-1
       nametrc=vname(1,indxT+itrc-1)
       unitt=vname(3,indxT+itrc-1)
       write(vname1,*) trim(nametrc),'_rate                       '
       vname(1,indxTrate+itrc-1)=vname1
CR       write(stdout,'(20x,A)') vname(1,indxTrate+itrc-1)
       vname(2,indxTrate+itrc-1)='Time rate of change             '
       write(vname3,*) trim(unitt),' second-1                     '
       vname(3,indxTrate+itrc-1)=vname3

      enddo


# endif /* DIAGNOSTICS_TS */




# ifdef DIAGNOSTICS_VRT

CR       write(stdout,'(1x,A,I3)') 'indxvrtXadv = ',indxvrtXadv+itrc-1
       nametrc='vrt'
       lvar=lenstr(nametrc)
       unitt=vname(3,indxU)
       write(vname1,*) trim(nametrc),'_xadv                       '
       vname(1,indxvrtXadv)=vname1
CR       write(stdout,'(20x,A)') vname(1,indxvrtxadv)
       vname(2,indxvrtXadv)='Horizontal (xi) advection term  '
       write(vname3,*) trim(unitt),' second-1                     '
       vname(3,indxvrtXadv)=vname3

CR       write(stdout,'(1x,A,I3)') 'indxvrtYadv = ',indxvrtYadv+itrc-1
       nametrc='vrt'
       unitt=vname(3,indxU)
       write(vname1,*) trim(nametrc),'_yadv                       '
       vname(1,indxvrtYadv)=vname1
CR       write(stdout,'(20x,A)') vname(1,indxvrtYadv)
       vname(2,indxvrtYadv)='Horizontal (eta) advection term '
       write(vname3,*) trim(unitt),' second-1                     '
       vname(3,indxvrtYadv)=vname3

CR       write(stdout,'(1x,A,I3)') 'indxvrtCor = ',indxvrtCor+itrc-1
       nametrc='vrt'
       unitt=vname(3,indxU)
       write(vname1,*) trim(nametrc),'_cor                        '
       vname(1,indxvrtCor)=vname1
CR       write(stdout,'(20x,A)') vname(1,indxvrtcor)
       vname(2,indxvrtCor)='Coriolis term                    '
       write(vname3,*) trim(unitt),' second-1                     '
       vname(3,indxvrtCor)=vname3

CR       write(stdout,'(1x,A,I3)') 'indxvrtPrsgrd = ',indxvrtPrsgrd+itrc-1
       nametrc='vrt'
       unitt=vname(3,indxU)
       write(vname1,*) trim(nametrc),'_Prsgrd                     '
       vname(1,indxvrtPrsgrd)=vname1
CR       write(stdout,'(20x,A)') vname(1,indxvrtPrsgrd)
       vname(2,indxvrtPrsgrd)='Pressure gradient term        '
       write(vname3,*) trim(unitt),' second-1                     '
       vname(3,indxvrtPrsgrd)=vname3

CR       write(stdout,'(1x,A,I3)') 'indxvrtHmix = ',indxvrtHmix+itrc-1
       nametrc='vrt'
       unitt=vname(3,indxU)
       write(vname1,*) trim(nametrc),'_hmix                       '
       vname(1,indxvrtHmix)=vname1
CR       write(stdout,'(20x,A)') vname(1,indxvrtHmix)
       vname(2,indxvrtHmix)='Horizontal mixing term          '
       write(vname3,*) trim(unitt),' second-1                     '
       vname(3,indxvrtHmix)=vname3

CR       write(stdout,'(1x,A,I3)') 'indxvrtVmix = ',indxvrtVmix+itrc-1
       nametrc='vrt'
       unitt=vname(3,indxU)
       write(vname1,*) trim(nametrc),'_vmix                       ' 
       vname(1,indxvrtVmix)=vname1
CR       write(stdout,'(20x,A)') vname(1,indxvrtVmix)
       vname(2,indxvrtVmix)='Vertical mixing term            '
       write(vname3,*) trim(unitt),' second-1                     '
       vname(3,indxvrtVmix)=vname3
       
CR       write(stdout,'(1x,A,I3)') 'indxvrtrate = ',indxvrtrate+itrc-1
       nametrc='vrt'
       unitt=vname(3,indxU)
       write(vname1,*) trim(nametrc),'_rate                       '
       vname(1,indxvrtrate)=vname1
CR       write(stdout,'(20x,A)') vname(1,indxvrtrate)
       vname(2,indxvrtrate)='Time rate of change             '
       write(vname3,*) trim(unitt),' second-1                     '
       vname(3,indxvrtrate)=vname3



CR       write(stdout,'(1x,A,I3)') 'indxvrtVmix2 = ',indxvrtVmix2+itrc-1
       nametrc='vrt'
       unitt=vname(3,indxU)
       write(vname1,*) trim(nametrc),'_vmix2                       ' 
       vname(1,indxvrtVmix2)=vname1
CR       write(stdout,'(20x,A)') vname(1,indxvrtVmix2)
       vname(2,indxvrtVmix2)='Vertical mixing term            '
       write(vname3,*) trim(unitt),' second-1                     '
       vname(3,indxvrtVmix2)=vname3

CR       write(stdout,'(1x,A,I3)') 'indxvrtWind = ',indxvrtWind+itrc-1
       nametrc='vrt'
       unitt=vname(3,indxU)
       write(vname1,*) trim(nametrc),'_Wind                       ' 
       vname(1,indxvrtWind)=vname1
CR       write(stdout,'(20x,A)') vname(1,indxvrtWind)
       vname(2,indxvrtWind)='Wind stress curl            '
       vname(3,indxvrtWind)=' meter second-2                     '



CR       write(stdout,'(1x,A,I3)') 'indxvrtDrag = ',indxvrtDrag+itrc-1
       nametrc='vrt'
       unitt=vname(3,indxU)
       write(vname1,*) trim(nametrc),'_Drag                       ' 
       vname(1,indxvrtDrag)=vname1
CR       write(stdout,'(20x,A)') vname(1,indxvrtDrag)
       vname(2,indxvrtDrag)='Bottom drag curl            '
       vname(3,indxvrtDrag)=' meter second-2                     '


# endif /* DIAGNOSTICS_VRT */



# ifdef DIAGNOSTICS_PV

      do itrc=1,2

CR       write(stdout,'(1x,A,I3)') 'indxpvForc_sfc = ',indxpvForc_sfc+itrc-1
       nametrc=vname(1,indxU+itrc-1)
       unitt=vname(3,indxU)
       write(vname1,*) trim(nametrc),'_pvForc_sfc                   ' 
       vname(1,indxpvForc_sfc+itrc-1)=vname1
CR       write(stdout,'(20x,A)') vname(1,indxpvForc_sfc+itrc-1)
       vname(2,indxpvForc_sfc+itrc-1)='Non-cons. Forces (surface)          '
       vname(3,indxpvForc_sfc+itrc-1)=' meter second-2                     '


CR       write(stdout,'(1x,A,I3)') 'indxpvForc_bot = ',indxpvForc_bot+itrc-1
       nametrc=vname(1,indxU+itrc-1)
       unitt=vname(3,indxU)
       write(vname1,*) trim(nametrc),'_pvForc_bot                    ' 
       vname(1,indxpvForc_bot+itrc-1)=vname1
CR       write(stdout,'(20x,A)') vname(1,indxpvForc_bot+itrc-1)
       vname(2,indxpvForc_bot+itrc-1)='Non-cons. Forces (bottom)            '
       vname(3,indxpvForc_bot+itrc-1)=' meter second-2                       ' 




#ifdef DIAGNOSTICS_PV_SUP

CR       write(stdout,'(1x,A,I3)') 'indxpvForc_0 = ',indxpvForc_0+itrc-1
       nametrc=vname(1,indxU+itrc-1)
       unitt=vname(3,indxU)
       write(vname1,*) trim(nametrc),'_pvForc_0                   ' 
       vname(1,indxpvForc_0+itrc-1)=vname1
CR       write(stdout,'(20x,A)') vname(1,indxpvForc_0+itrc-1)
       vname(2,indxpvForc_0+itrc-1)='Non-cons. Forces (surface)          '
       vname(3,indxpvForc_0+itrc-1)=' meter second-2                     '

CR       write(stdout,'(1x,A,I3)') 'indxpvForc_1 = ',indxpvForc_1+itrc-1
       nametrc=vname(1,indxU+itrc-1)
       unitt=vname(3,indxU)
       write(vname1,*) trim(nametrc),'_pvForc_1                   ' 
       vname(1,indxpvForc_1+itrc-1)=vname1
CR       write(stdout,'(20x,A)') vname(1,indxpvForc_1+itrc-1)
       vname(2,indxpvForc_1+itrc-1)='Non-cons. Forces (surface)          '
       vname(3,indxpvForc_1+itrc-1)=' meter second-2                     '


CR       write(stdout,'(1x,A,I3)') 'indxpvForc_2 = ',indxpvForc_2+itrc-1
       nametrc=vname(1,indxU+itrc-1)
       unitt=vname(3,indxU)
       write(vname1,*) trim(nametrc),'_pvForc_2                   ' 
       vname(1,indxpvForc_2+itrc-1)=vname1
CR       write(stdout,'(20x,A)') vname(1,indxpvForc_2+itrc-1)
       vname(2,indxpvForc_2+itrc-1)='Non-cons. Forces (surface)          '
       vname(3,indxpvForc_2+itrc-1)=' meter second-2                     '

CR       write(stdout,'(1x,A,I3)') 'indxpvForc_3 = ',indxpvForc_3+itrc-1
       nametrc=vname(1,indxU+itrc-1)
       unitt=vname(3,indxU)
       write(vname1,*) trim(nametrc),'_pvForc_3                   ' 
       vname(1,indxpvForc_3+itrc-1)=vname1
CR       write(stdout,'(20x,A)') vname(1,indxpvForc_3+itrc-1)
       vname(2,indxpvForc_3+itrc-1)='Non-cons. Forces (surface)          '
       vname(3,indxpvForc_3+itrc-1)=' meter second-2                     '


# endif

      enddo



      do itrc=1,NT


CR       write(stdout,'(1x,A,I3)') 'indxpvDiab_sfc = ',indxpvDiab_sfc+itrc-1
       nametrc=vname(1,indxT+itrc-1)
       unitt=vname(3,indxT+itrc-1)
       write(vname1,*) trim(nametrc),'_pvDiab_sfc     '
       vname(1,indxpvDiab_sfc+itrc-1)=vname1
CR       write(stdout,'(20x,A)') vname(1,indxpvDiab_sfc+itrc-1)
       vname(2,indxpvDiab_sfc+itrc-1)='Diabatic term (surface)    '
       write(vname3,*) trim(unitt),' second-1                     '
       vname(3,indxpvDiab_sfc+itrc-1)=vname3


CR       write(stdout,'(1x,A,I3)') 'indxpvDiab_bot = ',indxpvDiab_bot+itrc-1
       nametrc=vname(1,indxT+itrc-1)
       unitt=vname(3,indxT+itrc-1)
       write(vname1,*) trim(nametrc),'_pvDiab_bot                 '
       vname(1,indxpvDiab_bot+itrc-1)=vname1
CR       write(stdout,'(20x,A)') vname(1,indxpvDiab_bot+itrc-1)
       vname(2,indxpvDiab_bot+itrc-1)='Diabatic term (bottom)     '
       write(vname3,*) trim(unitt),' second-1                     '
       vname(3,indxpvDiab_bot+itrc-1)=vname3

#ifdef DIAGNOSTICS_PV_SUP

CR       write(stdout,'(1x,A,I3)') 'indxpvDiab_0 = ',indxpvDiab_0+itrc-1
       nametrc=vname(1,indxT+itrc-1)
       unitt=vname(3,indxT+itrc-1)
       write(vname1,*) trim(nametrc),'_pvDiab_0     '
       vname(1,indxpvDiab_0+itrc-1)=vname1
CR       write(stdout,'(20x,A)') vname(1,indxpvDiab_0+itrc-1)
       vname(2,indxpvDiab_0+itrc-1)='Diabatic term (surface)    '
       write(vname3,*) trim(unitt),' second-1                     '
       vname(3,indxpvDiab_0+itrc-1)=vname3


CR       write(stdout,'(1x,A,I3)') 'indxpvDiab_1 = ',indxpvDiab_1+itrc-1
       nametrc=vname(1,indxT+itrc-1)
       unitt=vname(3,indxT+itrc-1)
       write(vname1,*) trim(nametrc),'_pvDiab_1     '
       vname(1,indxpvDiab_1+itrc-1)=vname1
CR       write(stdout,'(20x,A)') vname(1,indxpvDiab_1+itrc-1)
       vname(2,indxpvDiab_1+itrc-1)='Diabatic term (surface)    '
       write(vname3,*) trim(unitt),' second-1                     '
       vname(3,indxpvDiab_1+itrc-1)=vname3

CR       write(stdout,'(1x,A,I3)') 'indxpvDiab_2 = ',indxpvDiab_2+itrc-1
       nametrc=vname(1,indxT+itrc-1)
       unitt=vname(3,indxT+itrc-1)
       write(vname1,*) trim(nametrc),'_pvDiab_2     '
       vname(1,indxpvDiab_2+itrc-1)=vname1
CR       write(stdout,'(20x,A)') vname(1,indxpvDiab_2+itrc-1)
       vname(2,indxpvDiab_2+itrc-1)='Diabatic term (surface)    '
       write(vname3,*) trim(unitt),' second-1                     '
       vname(3,indxpvDiab_2+itrc-1)=vname3

CR       write(stdout,'(1x,A,I3)') 'indxpvDiab_3 = ',indxpvDiab_3+itrc-1
       nametrc=vname(1,indxT+itrc-1)
       unitt=vname(3,indxT+itrc-1)
       write(vname1,*) trim(nametrc),'_pvDiab_3     '
       vname(1,indxpvDiab_3+itrc-1)=vname1
CR       write(stdout,'(20x,A)') vname(1,indxpvDiab_3+itrc-1)
       vname(2,indxpvDiab_3+itrc-1)='Diabatic term (surface)    '
       write(vname3,*) trim(unitt),' second-1                     '
       vname(3,indxpvDiab_3+itrc-1)=vname3












# endif


      enddo



# endif /* DIAGNOSTICS_PV */



# ifdef DIAGNOSTICS_EDDY
      itrc=1

CR       write(stdout,'(1x,A,I3)') 'indxeddyuu = ',indxeddyuu+itrc-1
       vname(1,indxeddyuu+itrc-1)='uu'
CR       write(stdout,'(20x,A)') vname(1,indxeddyxadv+itrc-1)
       vname(2,indxeddyuu+itrc-1)='uu'
       vname(3,indxeddyuu+itrc-1)='meter2 second-2'

CR       write(stdout,'(1x,A,I3)') 'indxeddyvv = ',indxeddyvv+itrc-1
       vname(1,indxeddyvv+itrc-1)='vv'
CR       write(stdout,'(20x,A)') vname(1,indxeddyvv+itrc-1)
       vname(2,indxeddyvv+itrc-1)='vv'
       vname(3,indxeddyvv+itrc-1)='meter2 second-2'

CR       write(stdout,'(1x,A,I3)') 'indxeddyuv = ',indxeddyuv+itrc-1
       vname(1,indxeddyuv+itrc-1)='uv'
CR       write(stdout,'(20x,A)') vname(1,indxeddyuv+itrc-1)
       vname(2,indxeddyuv+itrc-1)='uv'
       vname(3,indxeddyuv+itrc-1)='meter2 second-2'

CR       write(stdout,'(1x,A,I3)') 'indxeddyub = ',indxeddyub+itrc-1
       vname(1,indxeddyub+itrc-1)='ub'
CR       write(stdout,'(20x,A)') vname(1,indxeddycor+itrc-1)
       vname(2,indxeddyub+itrc-1)='ub'
       vname(3,indxeddyub+itrc-1)='meter2 second-3'

CR       write(stdout,'(1x,A,I3)') 'indxeddyvb = ',indxeddyvb+itrc-1
       vname(1,indxeddyvb+itrc-1)='vb'
CR       write(stdout,'(20x,A)') vname(1,indxeddyvb+itrc-1)
       vname(2,indxeddyvb+itrc-1)='vb'
       vname(3,indxeddyvb+itrc-1)='meter2 second-3'

CR       write(stdout,'(1x,A,I3)') 'indxeddywb = ',indxeddywb+itrc-1
       vname(1,indxeddywb+itrc-1)='wb'
CR       write(stdout,'(20x,A)') vname(1,indxeddywb+itrc-1)
       vname(2,indxeddywb+itrc-1)='wb'
       vname(3,indxeddywb+itrc-1)='meter2 second-3'

CR       write(stdout,'(1x,A,I3)') 'indxeddynotused1 = ',indxeddynotused1+itrc-1
        vname(1,indxeddynotused1+itrc-1)='uw'
CR       write(stdout,'(20x,A)') vname(1,indxeddynotused1+itrc-1)
       vname(2,indxeddynotused1+itrc-1)='uw'
       vname(3,indxeddynotused1+itrc-1)='meter2 second-2'
       
CR       write(stdout,'(1x,A,I3)') 'indxeddynotused2 = ',indxeddynotused2+itrc-1
       vname(1,indxeddynotused2+itrc-1)='vw'
CR       write(stdout,'(20x,A)') vname(1,indxeddynotused2+itrc-1)
       vname(2,indxeddynotused2+itrc-1)='vw'
       vname(3,indxeddynotused2+itrc-1)='meter2 second-2'

CR       write(stdout,'(1x,A,I3)') 'indxeddynotused3 = ',indxeddynotused3+itrc-1
       vname(1,indxeddynotused3+itrc-1)='notused3'
CR       write(stdout,'(20x,A)') vname(1,indxeddynotused3+itrc-1)
       vname(2,indxeddynotused3+itrc-1)='notused3'
       vname(3,indxeddynotused3+itrc-1)='None'


# endif /* DIAGNOSTICS_UV */



#endif

      ierr=0

      call read_inp (ierr)
      if (ierr.eq.0) then
#ifdef SOLVE3D
        call set_weights
#endif
      endif
      return
      end
 
