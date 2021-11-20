! This is "ncvars.h":   indices for character array "vname" for
!----- -- -----------   names of netCDF variables and attributes
!
! indxTime        time
! indxZ           free-surface
! indxUb,indxVb   vertically integrated 2D U,V-momentum components
!
! indxU,indxV     3D U- and V-momenta.
! indxT,indxS,    tracers (temerature, salinity,
!                 
! indxO,indeW     omega vertical mass flux and true vertical velocity
! indxR           density anomaly
!
! indxAkv,indxAkt,indxAks  vertical viscosity/diffusivity coeffcients
! indxHbl         depth of planetary boundary layer in KPP model
!
! indxSUSTR,indxSVSTR  surface U-, V-momentum stress (wind forcing)
! indxSHFl        net surface heat flux.
! indxSWRad       shortwave radiation flux
! indxSST         sea surface temperature
! indxdQdSST      Q-correction coefficient dQdSST
! indxSSFl        surface fresh water flux
!
! indxAi          fraction of cell covered by ice
! indxUi,indxVi   U,V-components of sea ice velocity
! indxHi,indxHS   depth of ice cover and depth of snow cover
! indxTIsrf       temperature of ice surface
!
! indxBSD,indxBSS bottom sediment grain Density and Size.
! indxWWA,indxWWD,indxWWP   wind induced wave Amplitude,
!                 Direction and Period
!
      integer indxTime, indxZ, indxUb, indxVb
      parameter (indxTime=1, indxZ=2, indxUb=3, indxVb=4)
#ifdef SOLVE3D
      integer indxU, indxV, indxT
      parameter (indxU=5, indxV=6, indxT=7)
# ifdef SALINITY
      integer indxS
      parameter (indxS=indxT+1)
# endif
# if (defined PASSIVE_TRACER && defined SALINITY)
      integer indxTPAS
      parameter (indxTPAS=indxS+1)
# endif
# if (defined PASSIVE_TRACER && !defined SALINITY)
      integer indxTPAS
      parameter (indxTPAS=indxT+1)
# endif
      integer indxO, indxW, indxR, indxAkv, indxAkt
      parameter (indxO=indxT+NT, indxW=indxO+1, indxR=indxO+2,
     &                     indxAkv=indxR+1, indxAkt=indxAkv+1)
# ifdef SALINITY
      integer indxAks
      parameter (indxAks=indxAkt+1)
# endif
# ifdef LMD_KPP
      integer indxHbl
#  ifdef SALINITY
      parameter (indxHbl=indxAks+1)
#  else
      parameter (indxHbl=indxAkt+1)
#  endif
# ifdef LMD_BKPP
      integer indxHbbl
      parameter (indxHbbl=indxHbl+1)      
# endif
# endif
#endif


      integer indxSUSTR, indxSVSTR
#ifdef SOLVE3D
      parameter (indxSUSTR=indxAkt+4, indxSVSTR=indxAkt+5)
#else
      parameter (indxSUSTR=indxVb+1,  indxSVSTR=indxSUSTR+1)
#endif

#ifdef SOLVE3D
      integer indxSHFl, indxSWRad
      parameter (indxSHFl=indxAkt+6)
# ifdef SALINITY
      integer indxSSFl
      parameter (indxSSFl=indxSHFl+1, indxSWRad=indxSHFl+2)
# else
      parameter (indxSWRad=indxSHFl+1)
# endif
      integer indxSST, indxdQdSST, indxSSS
      parameter (indxSST=indxSWRad+1, indxdQdSST=indxSST+1,
     &                                   indxSSS=indxSST+2)
#endif
#if defined ADV_ISONEUTRAL & defined WRITE_SXSY
      integer indxSX,indxSY     
      parameter(indxSX=indxSST+13,indxSY=indxSST+14) 
#endif

# if defined BULK_FLUX
      integer indxTAIR,indxQAIR,indxRADLW,indxRADSW,
     &        indxPREC,indxUWND,indxVWND
      parameter (
     &           indxTAIR = indxSST+15, indxQAIR = indxTAIR+1,
     &           indxRADLW= indxTAIR+2, indxRADSW= indxTAIR+3,
     &           indxPREC = indxTAIR+4,
     &           indxUWND = indxTAIR+5, indxVWND = indxTAIR+6  
     &    )
      integer indxShflx_rsw,indxShflx_rlw,
     & indxShflx_lat,indxShflx_sen,indxSustr_blk,indxSvstr_blk
      parameter (indxSustr_blk=indxTAIR+7, indxSvstr_blk=indxTAIR+8,
     &           indxShflx_rsw=indxTAIR+9, indxShflx_rlw=indxTAIR+10,
     &           indxShflx_lat=indxTAIR+11, indxShflx_sen=indxTAIR+12)
# endif






# ifdef DIAGNOSTICS_UV
      integer indxMXadv,indxMYadv,indxMVadv,indxMCor,
     &        indxMPrsgrd,indxMHmix,indxMVmix,indxMrate,
     &        indxMVmix2
# if defined BULK_FLUX
      parameter (indxMXadv=indxShflx_sen+1,
# else
      parameter (indxMXadv=indxSST+15,
# endif
     &           indxMYadv=indxMXadv+2,
     &           indxMVadv=indxMYadv+2,
     &           indxMCor=indxMVadv+2,
     &           indxMPrsgrd=indxMCor+2,
     &           indxMHmix=indxMPrsgrd+2,
     &           indxMVmix=indxMHmix+2, 
     &           indxMrate=indxMVmix+2,
     &           indxMVmix2=indxMrate+2)
# endif



# ifdef DIAGNOSTICS_VRT
      integer indxvrtXadv,indxvrtYadv,indxvrtCor,
     &        indxvrtPrsgrd,indxvrtHmix,indxvrtVmix,indxvrtrate,
     &        indxvrtVmix2,indxvrtWind,indxvrtDrag
# ifdef DIAGNOSTICS_UV
      parameter (indxvrtXadv=indxMVmix2+2,
# else
# if defined BULK_FLUX
      parameter (indxvrtXadv=indxShflx_sen+1,
# else
      parameter (indxvrtXadv=indxSST+15,
# endif
# endif
     &           indxvrtYadv=indxvrtXadv+1,
     &           indxvrtCor=indxvrtYadv+1,
     &           indxvrtPrsgrd=indxvrtCor+1,
     &           indxvrtHmix=indxvrtPrsgrd+1,
     &           indxvrtVmix=indxvrtHmix+1, 
     &           indxvrtrate=indxvrtVmix+1,
     &           indxvrtVmix2=indxvrtrate+1,
     &           indxvrtWind=indxvrtVmix2+1,
     &           indxvrtDrag=indxvrtWind+1)
# endif

      
      
# ifdef DIAGNOSTICS_TS
      integer indxTXadv,indxTYadv,indxTVadv,indxTForc,
     &        indxTHmix,indxTVmix,indxTrate
# ifdef DIAGNOSTICS_VRT
      parameter (indxTXadv=indxvrtDrag+1,
# else
#ifdef DIAGNOSTICS_UV
      parameter (indxTXadv=indxMVmix2+2,
# else 
# if defined BULK_FLUX
      parameter (indxTXadv=indxShflx_sen+1,
# else
      parameter (indxTXadv=indxSST+15,
# endif
# endif
# endif
     &           indxTYadv=indxTXadv+NT,
     &           indxTVadv=indxTYadv+NT,
     &           indxTForc=indxTVadv+NT,
     &           indxTHmix=indxTForc+NT,
     &           indxTVmix=indxTHmix+NT, 
     &           indxTrate=indxTVmix+NT)
# endif

      
      
      
# ifdef DIAGNOSTICS_PV
      integer indxpvForc_sfc,indxpvForc_bot,
     &        indxpvDiab_sfc,indxpvDiab_bot
# ifdef DIAGNOSTICS_TS     
      parameter (indxpvForc_sfc=indxTrate+NT,
# else                
# ifdef DIAGNOSTICS_VRT
      parameter (indxpvForc_sfc=indxvrtDrag+1,
# else
#ifdef DIAGNOSTICS_UV
      parameter (indxpvForc_sfc=indxMVmix2+2,
# else 
# if defined BULK_FLUX
      parameter (indxpvForc_sfc=indxShflx_sen+1,
# else
      parameter (indxpvForc_sfc=indxSST+15,
# endif
# endif              
# endif
# endif
     &           indxpvForc_bot= indxpvForc_sfc+2,
     &           indxpvDiab_sfc= indxpvForc_bot+NT,
     &           indxpvDiab_bot= indxpvDiab_sfc+NT)
# ifdef DIAGNOSTICS_PV_SUP    
        integer indxpvForc_0,indxpvForc_1,
     &        indxpvForc_2,indxpvForc_3,
     &        indxpvDiab_0,indxpvDiab_1,  
     &        indxpvDiab_2,indxpvDiab_3      
        parameter (indxpvForc_0=indxpvDiab_bot+NT,
     &           indxpvForc_1= indxpvForc_0+2,
     &           indxpvForc_2= indxpvForc_1+2,    
     &           indxpvForc_3= indxpvForc_2+2,  
     &           indxpvDiab_0= indxpvForc_3+2,
     &           indxpvDiab_1= indxpvDiab_0+2,
     &           indxpvDiab_2= indxpvDiab_1+2,
     &           indxpvDiab_3= indxpvDiab_2+2)    
# endif     
# endif


      
# ifdef DIAGNOSTICS_EDDY
      integer indxeddyuu,indxeddyvv,indxeddyuv,indxeddyub,
     &        indxeddyvb,indxeddywb,indxeddynotused1,indxeddynotused2,
     &        indxeddynotused3
# ifdef DIAGNOSTICS_PV_SUP
      parameter (indxeddyuu=indxpvDiab_3+2,
# else  
# ifdef DIAGNOSTICS_PV
      parameter (indxeddyuu=indxpvDiab_bot+NT,
# else  
# ifdef DIAGNOSTICS_TS     
      parameter (indxeddyuu=indxTrate+NT,
# else                
# ifdef DIAGNOSTICS_VRT
      parameter (indxeddyuu=indxvrtDrag+1,
# else
#ifdef DIAGNOSTICS_UV
      parameter (indxeddyuu=indxMVmix2+2,
# else 
# if defined BULK_FLUX
      parameter (indxeddyuu=indxShflx_sen+1,
# else
      parameter (indxeddyuu=indxSST+15,
# endif
# endif
# endif              
# endif
# endif
# endif
     &           indxeddyvv=indxeddyuu+1,
     &           indxeddyuv=indxeddyvv+1,
     &           indxeddyub=indxeddyuv+1,
     &           indxeddyvb=indxeddyub+1,
     &           indxeddywb=indxeddyvb+1,
     &           indxeddynotused1=indxeddywb+1, 
     &           indxeddynotused2=indxeddynotused1+1,
     &           indxeddynotused3=indxeddynotused2+1)
# endif
    
      
      
      
!
! Naming conventions for indices, variable IDs, etc...
!
! prefix ncid_  means netCDF ID for netCDF file
!        nrec_  record number in netCDF file since initialization
!        nrpf_  maximum number of records per file  (output netCDF
!                                                       files only)
! prefix/ending rst_/_rst refers to restart  netCDF file
!               his_/_his           history
!               avg_/_avg           averages
!               sta_/_sta           stations
!                    _frc           forcing
!                    _clm           climatology
!
! endings refer to:  ___Time  time [in seconds]
!                    ___Tstep time step numbers and record numbers
! all objects with   ___Z     free-surface
! these endings are  ___Ub    vertically integrated 2D U-momentum
! either:            ___Vb    vertically integrated 2D V-momentum
!
!  netCDF IDs, if    ___U     3D U-momentum
!  occur with prefix ___V     3D V-momentum
!  rst/his/avg/sta   ___T(NT) tracers
!                    ___R     density anomaly
! or                 ___O     omega vertical velocity
!                    ___W     true vertical velocity
!  parameter indices 
!  if combined with  ___Akv   vertical viscosity coefficient
!  prefix "indx"     ___Akt   vertical T-diffusion coefficient
!  (see above).      ___Aks   vertical S-diffusion coefficient
!                    ___Hbl   depth of mixed layer LMD_KPP.
!
! vname    character array for variable names and attributes;
!
      integer max_frc_file
      parameter (max_frc_file=4)
      integer max_frc, ncidfrc(max_frc_file), nrst, ncidrst, nrecrst,
     &      nrrec, nrpfrst, ncidclm, nwrt, ncidhis, nrechis, nrpfhis
      common /ncvars/       max_frc, ncidfrc, nrst, ncidrst, nrecrst,
     &      nrrec, nrpfrst, ncidclm, nwrt, ncidhis, nrechis, nrpfhis
      








#ifdef AVERAGES
      integer ntsavg, navg
      common /ncvars/ ntsavg, navg
#endif
#ifdef STATIONS
      integer nsta
      common /ncvars/ nsta
#endif
#ifdef FLOATS
      integer nflt
      common /ncvars/ nflt
#endif


# ifdef DIAGNOSTICS_UV
      integer nwrtdiags_uv
      common /ncvars/ nwrtdiags_uv
#ifdef AVERAGES
      integer ntsdiags_uv_avg, nwrtdiags_uv_avg
      common /ncvars/ ntsdiags_uv_avg, nwrtdiags_uv_avg
#endif
#endif


# ifdef DIAGNOSTICS_VRT
      integer nwrtdiags_vrt
      common /ncvars/ nwrtdiags_vrt
#ifdef AVERAGES
      integer ntsdiags_vrt_avg, nwrtdiags_vrt_avg
      common /ncvars/ ntsdiags_vrt_avg, nwrtdiags_vrt_avg
#endif
#endif



# ifdef DIAGNOSTICS_TS
      integer nwrtdiags_ts
      common /ncvars/ nwrtdiags_ts
#ifdef AVERAGES
      integer ntsdiags_ts_avg, nwrtdiags_ts_avg
      common /ncvars/ ntsdiags_ts_avg, nwrtdiags_ts_avg
#endif
#endif


# ifdef DIAGNOSTICS_PV
      integer nwrtdiags_pv
      common /ncvars/ nwrtdiags_pv
#ifdef AVERAGES
      integer ntsdiags_pv_avg, nwrtdiags_pv_avg
      common /ncvars/ ntsdiags_pv_avg, nwrtdiags_pv_avg
#endif
#endif

# ifdef DIAGNOSTICS_EDDY
      integer nwrtdiags_eddy
      common /ncvars/ nwrtdiags_eddy
#ifdef AVERAGES
      integer ntsdiags_eddy_avg, nwrtdiags_eddy_avg
      common /ncvars/ ntsdiags_eddy_avg, nwrtdiags_eddy_avg
#endif
#endif

      integer rstTime, rstTstep,      rstZ,   rstUb,  rstVb,
     &        hisTime, hisTstep,      hisZ,   hisUb,  hisVb
      common /ncvars/
     &        rstTime, rstTstep,      rstZ,   rstUb,  rstVb,
     &        hisTime, hisTstep,      hisZ,   hisUb,  hisVb
#ifdef SOLVE3D
# ifdef EXACT_RESTART
      integer rst_DU_avg2, rst_DV_avg2
      common /ncvars/ rst_DU_avg2, rst_DV_avg2
# endif
      integer rstU, rstV, rstT(NT+1), hisO,   hisW,   hisR,
     &        hisU, hisV, hisT(NT+1), hisAkv, hisAkt, hisAks
      common /ncvars/
     &        rstU, rstV, rstT,       hisO,   hisW,   hisR,
     &        hisU, hisV, hisT,       hisAkv, hisAkt, hisAks

# ifdef LMD_KPP
      integer rstHbl, hisHbl
      common /ncvars/ rstHbl, hisHbl
# endif
# ifdef LMD_BKPP
      integer rstHbbl, hisHbbl
      common /ncvars/ rstHbbl, hisHbbl
# endif
#if defined ADV_ISONEUTRAL & defined WRITE_SXSY
      integer hisSX,hisSY
      common /ncvars/ hisSX, hisSY
# endif
# ifdef BULK_FLUX
      integer hisSustr_blk,   hisSvstr_blk
     &      , hisShflx_rsw,   hisShflx_rlw
     &      , hisShflx_lat,   hisShflx_sen
      common /ncvars/ hisSustr_blk,   hisSvstr_blk
     &      , hisShflx_rsw,   hisShflx_rlw
     &      , hisShflx_lat,   hisShflx_sen
# endif
# ifdef DIAGNOSTICS_UV
      real timediags_uv
      common /ncvars/
     &        timediags_uv


      integer nciddiags_uv, nrecdiags_uv, nrpfdiags_uv 
     &      , diags_uvTime, diags_uvTstep
     &      , diags_uvXadv(2), diags_uvYadv(2), diags_uvVadv(2)
     &      , diags_uvCor(2), diags_uvPrsgrd(2), diags_uvHmix(2)
     &      , diags_uvVmix(2), diags_uvrate(2)
     &      , diags_uvVmix2(2)
      common /ncvars/   nciddiags_uv, nrecdiags_uv, nrpfdiags_uv
     &      , diags_uvTime, diags_uvTstep
     &      , diags_uvXadv, diags_uvYadv, diags_uvVadv, diags_uvCor
     &      , diags_uvPrsgrd, diags_uvHmix, diags_uvVmix, diags_uvrate
     &      , diags_uvVmix2
# endif


# ifdef DIAGNOSTICS_VRT
      real timediags_vrt
      common /ncvars/
     &        timediags_vrt


      integer nciddiags_vrt, nrecdiags_vrt, nrpfdiags_vrt 
     &      , diags_vrtTime, diags_vrtTstep
     &      , diags_vrtXadv(2), diags_vrtYadv(2)
     &      , diags_vrtCor(2), diags_vrtPrsgrd(2), diags_vrtHmix(2)
     &      , diags_vrtVmix(2), diags_vrtrate(2)
     &      , diags_vrtVmix2(2), diags_vrtWind(2), diags_vrtDrag(2)
      common /ncvars/   nciddiags_vrt, nrecdiags_vrt, nrpfdiags_vrt
     &      , diags_vrtTime, diags_vrtTstep
     &      , diags_vrtXadv, diags_vrtYadv, diags_vrtCor
     &      , diags_vrtPrsgrd, diags_vrtHmix, diags_vrtVmix, diags_vrtrate
     &      , diags_vrtVmix2, diags_vrtWind, diags_vrtDrag
# endif

# ifdef DIAGNOSTICS_TS
      real timediags_ts
      common /ncvars/
     &        timediags_ts


      integer nciddiags_ts, nrecdiags_ts, nrpfdiags_ts 
     &      , diags_tsTime, diags_tsTstep
     &      , diags_tsXadv(2), diags_tsYadv(2), diags_tsVadv(2)
     &      , diags_tsForc(2),  diags_tsHmix(2)
     &      , diags_tsVmix(2), diags_tsrate(2)
      common /ncvars/   nciddiags_ts, nrecdiags_ts, nrpfdiags_ts
     &      , diags_tsTime, diags_tsTstep
     &      , diags_tsXadv, diags_tsYadv, diags_tsVadv, diags_tsForc
     &      , diags_tsHmix, diags_tsVmix, diags_tsrate
# endif


# ifdef DIAGNOSTICS_PV
      real timediags_pv
      common /ncvars/
     &        timediags_pv


      integer nciddiags_pv, nrecdiags_pv, nrpfdiags_pv 
     &      , diags_pvTime, diags_pvTstep
     &      , diags_pvForc_sfc(2), diags_pvForc_bot(2)
     &      , diags_pvDiab_sfc(NT), diags_pvDiab_bot(NT)
      common /ncvars/   nciddiags_pv, nrecdiags_pv, nrpfdiags_pv
     &      , diags_pvTime, diags_pvTstep
     &      , diags_pvForc_sfc, diags_pvForc_bot
     &      , diags_pvDiab_sfc, diags_pvDiab_bot
# ifdef DIAGNOSTICS_PV_SUP
         integer diags_pvForc_0(2), diags_pvForc_1(2)
     &      , diags_pvForc_2(2), diags_pvForc_3(2)  
     &      , diags_pvDiab_0(NT), diags_pvDiab_1(NT)
     &      , diags_pvDiab_2(NT), diags_pvDiab_3(NT)    
      common /ncvars/  diags_pvForc_0, diags_pvForc_1
     &      , diags_pvForc_2, diags_pvForc_3
     &      , diags_pvDiab_0, diags_pvDiab_1      
     &      , diags_pvDiab_2, diags_pvDiab_3 

# endif    
# endif

# ifdef DIAGNOSTICS_EDDY
      real timediags_eddy
      common /ncvars/
     &        timediags_eddy


      integer nciddiags_eddy, nrecdiags_eddy, nrpfdiags_eddy 
     &      , diags_eddyTime, diags_eddyTstep
     &      , diags_eddyuu(2), diags_eddyvv(2), diags_eddyuv(2)
     &      , diags_eddyub(2), diags_eddyvb(2), diags_eddywb(2)
     &      , diags_eddynotused1(2), diags_eddynotused2(2)
     &      , diags_eddynotused3(2)
      common /ncvars/   nciddiags_eddy, nrecdiags_eddy, nrpfdiags_eddy
     &      , diags_eddyTime, diags_eddyTstep
     &      , diags_eddyuu, diags_eddyvv, diags_eddyuv, diags_eddyub
     &      , diags_eddyvb, diags_eddywb, diags_eddynotused1, diags_eddynotused2
     &      , diags_eddynotused3
# endif

# endif

#ifdef AVERAGES
      integer ncidavg, nrecavg,  nrpfavg,
     &        avgTime, avgTstep, avgZ,    avgUb, avgVb
      common /ncvars/  ncidavg,  nrecavg, nrpfavg,
     &        avgTime, avgTstep, avgZ,    avgUb, avgVb
# ifdef SOLVE3D
      integer avgU,  avgV,  avgT(NT+1), avgR,
     &        avgO,  avgW,  avgAkv,     avgAkt,  avgAks
      common /ncvars/ avgU, avgV,       avgT,    avgR, 
     &        avgO,  avgW,  avgAkv,     avgAkt,  avgAks
#  ifdef LMD_KPP
      integer avgHbl
      common /ncvars/ avgHbl
#  endif
#  ifdef LMD_BKPP
      integer avgHbbl
      common /ncvars/ avgHbbl
#  endif
#if defined ADV_ISONEUTRAL & defined WRITE_SXSY
      integer avgSX,avgSY
      common /ncvars/ avgSX, avgSY
# endif

#  ifdef BULK_FLUX
      integer avgSustr_blk,   avgSvstr_blk
     &      , avgShflx_rsw,   avgShflx_rlw
     &      , avgShflx_lat,   avgShflx_sen
      common /ncvars/ avgSustr_blk,   avgSvstr_blk
     &      , avgShflx_rsw,   avgShflx_rlw
     &      , avgShflx_lat,   avgShflx_sen 
#  endif


#  ifdef DIAGNOSTICS_UV
      real timediags_uv_avg
      common /ncvars/
     &        timediags_uv_avg

       integer nciddiags_uv_avg, nrecdiags_uv_avg, nrpfdiags_uv_avg 
     &      , diags_uv_avgTime, diags_uv_avgTstep
     &      , diags_uvXadv_avg(2), diags_uvYadv_avg(2), diags_uvVadv_avg(2)
     &      , diags_uvCor_avg(2), diags_uvPrsgrd_avg(2), diags_uvHmix_avg(2)
     &      , diags_uvVmix_avg(2), diags_uvrate_avg(2)
     &      , diags_uvVmix2_avg(2)
      common /ncvars/   nciddiags_uv_avg, nrecdiags_uv_avg, nrpfdiags_uv_avg
     &      , diags_uv_avgTime, diags_uv_avgTstep
     &      , diags_uvXadv_avg, diags_uvYadv_avg, diags_uvVadv_avg
     &      , diags_uvCor_avg, diags_uvPrsgrd_avg, diags_uvHmix_avg
     &      , diags_uvVmix_avg, diags_uvrate_avg
     &      , diags_uvVmix2_avg
#  endif

#  ifdef DIAGNOSTICS_VRT
      real timediags_vrt_avg
      common /ncvars/
     &        timediags_vrt_avg

       integer nciddiags_vrt_avg, nrecdiags_vrt_avg, nrpfdiags_vrt_avg 
     &      , diags_vrt_avgTime, diags_vrt_avgTstep
     &      , diags_vrtXadv_avg(2), diags_vrtYadv_avg(2)
     &      , diags_vrtCor_avg(2), diags_vrtPrsgrd_avg(2), diags_vrtHmix_avg(2)
     &      , diags_vrtVmix_avg(2), diags_vrtrate_avg(2)
     &      , diags_vrtVmix2_avg(2), diags_vrtWind_avg(2), diags_vrtDrag_avg(2)
      common /ncvars/   nciddiags_vrt_avg, nrecdiags_vrt_avg, nrpfdiags_vrt_avg
     &      , diags_vrt_avgTime, diags_vrt_avgTstep
     &      , diags_vrtXadv_avg, diags_vrtYadv_avg
     &      , diags_vrtCor_avg, diags_vrtPrsgrd_avg, diags_vrtHmix_avg
     &      , diags_vrtVmix_avg, diags_vrtrate_avg
     &      , diags_vrtVmix2_avg, diags_vrtWind_avg, diags_vrtDrag_avg
#  endif


# ifdef DIAGNOSTICS_TS
      real timediags_ts_avg
      common /ncvars/
     &        timediags_ts_avg


      integer nciddiags_ts_avg, nrecdiags_ts_avg, nrpfdiags_ts _avg
     &      , diags_ts_avgTime, diags_ts_avgTstep
     &      , diags_tsXadv_avg(2), diags_tsYadv_avg(2), diags_tsVadv_avg(2)
     &      , diags_tsForc_avg(2),  diags_tsHmix_avg(2)
     &      , diags_tsVmix_avg(2), diags_tsrate_avg(2)
      common /ncvars/   nciddiags_ts_avg, nrecdiags_ts_avg, nrpfdiags_ts_avg
     &      , diags_ts_avgTime, diags_ts_avgTstep
     &      , diags_tsXadv_avg, diags_tsYadv_avg, diags_tsVadv_avg, diags_tsForc_avg
     &      , diags_tsHmix_avg, diags_tsVmix_avg, diags_tsrate_avg
# endif


# ifdef DIAGNOSTICS_PV
      real timediags_pv_avg
      common /ncvars/
     &        timediags_pv_avg


      integer nciddiags_pv_avg, nrecdiags_pv_avg, nrpfdiags_pv _avg
     &      , diags_pv_avgTime, diags_pv_avgTstep
     &      , diags_pvForc_sfc_avg(2), diags_pvForc_bot_avg(2)
     &      , diags_pvDiab_sfc_avg(NT), diags_pvDiab_bot_avg(NT)
      common /ncvars/   nciddiags_pv_avg, nrecdiags_pv_avg, nrpfdiags_pv_avg
     &      , diags_pv_avgTime, diags_pv_avgTstep
     &      , diags_pvForc_sfc_avg, diags_pvForc_bot_avg
     &      , diags_pvDiab_sfc_avg, diags_pvDiab_bot_avg
# endif


#  ifdef DIAGNOSTICS_EDDY
      real timediags_eddy_avg
      common /ncvars/
     &        timediags_eddy_avg

       integer nciddiags_eddy_avg, nrecdiags_eddy_avg, nrpfdiags_eddy_avg 
     &      , diags_eddy_avgTime, diags_eddy_avgTstep
     &      , diags_eddyuu_avg(2), diags_eddyvv_avg(2), diags_eddyuv_avg(2)
     &      , diags_eddyub_avg(2), diags_eddyvb_avg(2), diags_eddywb_avg(2)
     &      , diags_eddynotused1_avg(2), diags_eddynotused2_avg(2)
     &      , diags_eddynotused3_avg(2)
      common /ncvars/   nciddiags_eddy_avg, nrecdiags_eddy_avg, nrpfdiags_eddy_avg
     &      , diags_eddy_avgTime, diags_eddy_avgTstep
     &      , diags_eddyuu_avg, diags_eddyvv_avg, diags_eddyuv_avg
     &      , diags_eddyub_avg, diags_eddyvb_avg, diags_eddywb_avg
     &      , diags_eddynotused1_avg, diags_eddynotused2_avg
     &      , diags_eddynotused3_avg
#  endif


# endif
#endif


#ifdef SOLVE3D
# define NWRTHIS 150+NT
#else
# define NWRTHIS 14      
#endif
      logical ldefhis, wrthis(NWRTHIS)
      common /ncvars/ ldefhis, wrthis
#ifdef DIAGNOSTICS_UV
      logical ldefdiags_uv, wrtdiags_uv(3)
      common /ncvars/ ldefdiags_uv, wrtdiags_uv
#endif

#ifdef DIAGNOSTICS_VRT
      logical ldefdiags_vrt, wrtdiags_vrt(3)
      common /ncvars/ ldefdiags_vrt, wrtdiags_vrt
#endif

#ifdef DIAGNOSTICS_TS
      logical ldefdiags_ts, wrtdiags_ts(3)
      common /ncvars/ ldefdiags_ts, wrtdiags_ts
#endif

#ifdef DIAGNOSTICS_PV
      logical ldefdiags_pv, wrtdiags_pv(3)
      common /ncvars/ ldefdiags_pv, wrtdiags_pv
#endif

#ifdef DIAGNOSTICS_EDDY
      logical ldefdiags_eddy, wrtdiags_eddy(3)
      common /ncvars/ ldefdiags_eddy, wrtdiags_eddy
#endif


#ifdef AVERAGES
      logical wrtavg(NWRTHIS)
      common /ncvars/ wrtavg
      
#ifdef DIAGNOSTICS_UV
      logical ldefdiags_uv_avg, wrtdiags_uv_avg(3)
      common /ncvars/ ldefdiags_uv_avg, wrtdiags_uv_avg
#endif

#ifdef DIAGNOSTICS_VRT
      logical ldefdiags_vrt_avg, wrtdiags_vrt_avg(3)
      common /ncvars/ ldefdiags_vrt_avg, wrtdiags_vrt_avg
#endif

#ifdef DIAGNOSTICS_TS
      logical ldefdiags_ts_avg, wrtdiags_ts_avg(3)
      common /ncvars/ ldefdiags_ts_avg, wrtdiags_ts_avg
#endif

#ifdef DIAGNOSTICS_PV
      logical ldefdiags_pv_avg, wrtdiags_pv_avg(3)
      common /ncvars/ ldefdiags_pv_avg, wrtdiags_pv_avg
#endif

#ifdef DIAGNOSTICS_EDDY
      logical ldefdiags_eddy_avg, wrtdiags_eddy_avg(3)
      common /ncvars/ ldefdiags_eddy_avg, wrtdiags_eddy_avg
#endif

#endif


#ifdef FLOATS
      logical ldefflt
      common /ncvars/ ldefflt
#endif
#ifdef STATIONS
      logical wsta(NWRTHIS)
      common /ncvars/ wsta
#endif
!
! Grid Type Codes:  r2dvar....w3hvar are codes for array types.
! ==== ==== ======  The codes are set according to the rule:
!                     horiz_grid_type+4*vert_grid_type
!    where horiz_grid_type=0,1,2,3 for RHO-,U-,V-,PSI-points
!    respectively and vert_grid_type=0 for 2D fields; 1,2 for
!    3D-RHO- and W-vertical points.
!
      integer r2dvar, u2dvar, v2dvar, p2dvar, r3dvar,
     &                u3dvar, v3dvar, p3dvar, w3dvar,
     &                pw3dvar
      parameter (r2dvar=0, u2dvar=1, v2dvar=2, p2dvar=3,
     & r3dvar=4, u3dvar=5, v3dvar=6, p3dvar=7, w3dvar=8,
     & pw3dvar=9)
!
!            Horizontal array dimensions in netCDF files.
! xi_rho     NOTE: In MPI mode using PARALLEL_FILES these
! xi_u       depend on corresonding sizes of individual MPI
! eta_rho    subdomains rather than the whole physical grid, 
! eta_v      and therefore become live variables placed in 
!            common block here and set in mpi_setup.
!
      integer xi_rho,xi_u, eta_rho,eta_v
#if defined MPI && defined PARALLEL_FILES
      common /ncvars/ xi_rho,xi_u, eta_rho,eta_v
#else
      parameter (xi_rho=LLm+2,   eta_rho=MMm+2,
     &           xi_u=xi_rho-1,  eta_v=eta_rho-1)
#endif

      integer max_name_size
      parameter (max_name_size=64)
      character date_str*44, title*80
      character*(max_name_size) ininame, grdname,
     &                 hisname, rstname, 
#if defined BULK_FLUX     
     &    blkfile(max_blk_file) 
#else
     &    frcfile(max_frc_file)     
# endif      
      common /cncvars/ date_str, title,  ininame,
     &        grdname, hisname, rstname
#if defined BULK_FLUX     
     & , blkfile 
#else
     & , frcfile      
#endif      
#ifdef AVERAGES
      character*(max_name_size) avgname
      common /cncvars/ avgname
#endif

#ifdef DIAGNOSTICS_UV
      character*(max_name_size) diags_uvname
      common /cncvars/ diags_uvname
#ifdef AVERAGES
      character*(max_name_size) diags_uv_avgname
      common /cncvars/ diags_uv_avgname
#endif
#endif

#ifdef DIAGNOSTICS_VRT
      character*(max_name_size) diags_vrtname
      common /cncvars/ diags_vrtname
#ifdef AVERAGES
      character*(max_name_size) diags_vrt_avgname
      common /cncvars/ diags_vrt_avgname
#endif
#endif

#ifdef DIAGNOSTICS_TS
      character*(max_name_size) diags_tsname
      common /cncvars/ diags_tsname
#ifdef AVERAGES
      character*(max_name_size) diags_ts_avgname
      common /cncvars/ diags_ts_avgname
#endif
#endif

#ifdef DIAGNOSTICS_PV
      character*(max_name_size) diags_pvname
      common /cncvars/ diags_pvname
#ifdef AVERAGES
      character*(max_name_size) diags_pv_avgname
      common /cncvars/ diags_pv_avgname
#endif
#endif

#ifdef DIAGNOSTICS_EDDY
      character*(max_name_size) diags_eddyname
      common /cncvars/ diags_eddyname
#ifdef AVERAGES
      character*(max_name_size) diags_eddy_avgname
      common /cncvars/ diags_eddy_avgname
#endif
#endif

#if (defined TCLIMATOLOGY && !defined ANA_TCLIMA) || !defined ANA_SSH
      character*(max_name_size) clm_file
      common /cncvars/ clm_file
#endif
#if defined T_FRC_BRY  || defined M2_FRC_BRY || \
    defined M3_FRC_BRY || defined Z_FRC_BRY
      character*(max_name_size) bry_file 
      common /cncvars/ bry_file
#endif

      character*42  vname(3,
     &                       199)
      common /cncvars/ vname
 
