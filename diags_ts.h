! This is include file "diags_vrt.h"
!  ==== == ======= ==== ==========
!

#ifdef DIAGNOSTICS_TS
      real TXadv(GLOBAL_2D_ARRAY,N,NT)
CSDISTRIBUTE_RESHAPE TXadv(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real TYadv(GLOBAL_2D_ARRAY,N,NT)
CSDISTRIBUTE_RESHAPE TYadv(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real TVadv(GLOBAL_2D_ARRAY,N,NT)
CSDISTRIBUTE_RESHAPE TVadv(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real TForc(GLOBAL_2D_ARRAY,N,NT)
CSDISTRIBUTE_RESHAPE TForc(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real THmix(GLOBAL_2D_ARRAY,N,NT)
CSDISTRIBUTE_RESHAPE THmix(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real TVmix(GLOBAL_2D_ARRAY,N,NT)
CSDISTRIBUTE_RESHAPE TVmix(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real Trate(GLOBAL_2D_ARRAY,N,NT)
CSDISTRIBUTE_RESHAPE Trate(BLOCK_PATTERN,*) BLOCK_CLAUSE


# ifdef AVERAGES
      real TXadv_avg(GLOBAL_2D_ARRAY,N,NT)
CSDISTRIBUTE_RESHAPE TXadv_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real TYadv_avg(GLOBAL_2D_ARRAY,N,NT)
CSDISTRIBUTE_RESHAPE TYadv_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real TVadv_avg(GLOBAL_2D_ARRAY,N,NT)
CSDISTRIBUTE_RESHAPE TVadv_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real TForc_avg(GLOBAL_2D_ARRAY,N,NT)
CSDISTRIBUTE_RESHAPE TForc_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real THmix_avg(GLOBAL_2D_ARRAY,N,NT)
CSDISTRIBUTE_RESHAPE THmix_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real TVmix_avg(GLOBAL_2D_ARRAY,N,NT)
CSDISTRIBUTE_RESHAPE TVmix_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real Trate_avg(GLOBAL_2D_ARRAY,N,NT)
CSDISTRIBUTE_RESHAPE Trate_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
# endif



      common /diag_TXadv/TXadv   
     &       /diag_TYadv/TYadv
     &       /diag_TVadv/TVadv  
     &       /diag_TForc/TForc
     &       /diag_THmix/THmix
     &       /diag_TVmix/TVmix
     &       /diag_Trate/Trate
# ifdef AVERAGES
      common /diag_TXadv_avg/TXadv_avg
     &       /diag_TYadv_avg/TYadv_avg
     &       /diag_TVadv_avg/TVadv_avg
     &       /diag_TForc_avg/TForc_avg
     &       /diag_THmix_avg/THmix_avg
     &       /diag_TVmix_avg/TVmix_avg
     &       /diag_Trate_avg/Trate_avg
# endif      
 

#endif /* DIAGNOSTICS_TS */



 
