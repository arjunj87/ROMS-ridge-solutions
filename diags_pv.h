! This is include file "diags_pv.h"
!  ==== == ======= ==== ==========
!

#ifdef DIAGNOSTICS_PV

      real pvForc_sfc(GLOBAL_2D_ARRAY,2)
CSDISTRIBUTE_RESHAPE pvForc_sfc(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real pvForc_bot(GLOBAL_2D_ARRAY,2)
CSDISTRIBUTE_RESHAPE pvForc_bot(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real pvDiab_sfc(GLOBAL_2D_ARRAY,NT)
CSDISTRIBUTE_RESHAPE pvDiab_sfc(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real pvDiab_bot(GLOBAL_2D_ARRAY,NT)
CSDISTRIBUTE_RESHAPE pvDiab_bot(BLOCK_PATTERN,*) BLOCK_CLAUSE


# ifdef AVERAGES
      real pvForc_sfc_avg(GLOBAL_2D_ARRAY,2)
CSDISTRIBUTE_RESHAPE pvForc_sfc_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real pvForc_bot_avg(GLOBAL_2D_ARRAY,2)
CSDISTRIBUTE_RESHAPE pvForc_bot_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real pvDiab_sfc_avg(GLOBAL_2D_ARRAY,NT)
CSDISTRIBUTE_RESHAPE pvDiab_sfc_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real pvDiab_bot_avg(GLOBAL_2D_ARRAY,NT)
CSDISTRIBUTE_RESHAPE pvDiab_bot_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
# endif



      common /diag_pvForc_sfc/pvForc_sfc 
     &       /diag_pvForc_bot/pvForc_bot
     &       /diag_pvDiab_sfc/pvDiab_sfc
     &       /diag_pvDiab_bot/pvDiab_bot

# ifdef AVERAGES
      common /diag_pvForc_sfc_avg /pvForc_sfc_avg 
     &       /diag_pvForc_bot_avg/pvForc_bot_avg
     &       /diag_pvDiab_sfc_avg/pvDiab_sfc_avg
     &       /diag_pvDiab_bot_avg/pvDiab_bot_avg
# endif      


#ifdef DIAGNOSTICS_PV_SUP

      real pvForc_0(GLOBAL_2D_ARRAY,2)
CSDISTRIBUTE_RESHAPE pvForc_0(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real pvForc_1(GLOBAL_2D_ARRAY,2)
CSDISTRIBUTE_RESHAPE pvForc_1(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real pvForc_2(GLOBAL_2D_ARRAY,2)
CSDISTRIBUTE_RESHAPE pvForc_2(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real pvForc_3(GLOBAL_2D_ARRAY,2)
CSDISTRIBUTE_RESHAPE pvForc_3(BLOCK_PATTERN,*) BLOCK_CLAUSE

      real pvDiab_0(GLOBAL_2D_ARRAY,NT)
CSDISTRIBUTE_RESHAPE pvDiab_0(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real pvDiab_1(GLOBAL_2D_ARRAY,NT)
CSDISTRIBUTE_RESHAPE pvDiab_1(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real pvDiab_2(GLOBAL_2D_ARRAY,NT)
CSDISTRIBUTE_RESHAPE pvDiab_2(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real pvDiab_3(GLOBAL_2D_ARRAY,NT)
CSDISTRIBUTE_RESHAPE pvDiab_3(BLOCK_PATTERN,*) BLOCK_CLAUSE

      common /diag_pvForc_0/pvForc_0
     &       /diag_pvForc_1/pvForc_1
     &       /diag_pvForc_2/pvForc_2   
     &       /diag_pvForc_3/pvForc_3   
     &       /diag_pvDiab_0/pvDiab_0
     &       /diag_pvDiab_1/pvDiab_1
     &       /diag_pvDiab_2/pvDiab_2
     &       /diag_pvDiab_3/pvDiab_3
     
#endif /* DIAGNOSTICS_VRT */



#endif /* DIAGNOSTICS_VRT */



 
