! This is include file "diags_eddy.h"
!  ==== == ======= ==== ==========
!

#ifdef DIAGNOSTICS_EDDY
# ifdef AVERAGES
      real eddyuu_avg(GLOBAL_2D_ARRAY,N)
CSDISTRIBUTE_RESHAPE eddyuu_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real eddyvv_avg(GLOBAL_2D_ARRAY,N)
CSDISTRIBUTE_RESHAPE eddyvv_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real eddyuv_avg(GLOBAL_2D_ARRAY,N)
CSDISTRIBUTE_RESHAPE eddyuv_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real eddyub_avg(GLOBAL_2D_ARRAY,N)
CSDISTRIBUTE_RESHAPE eddyub_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real eddyvb_avg(GLOBAL_2D_ARRAY,N)
CSDISTRIBUTE_RESHAPE eddyvb_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real eddywb_avg(GLOBAL_2D_ARRAY,N)
CSDISTRIBUTE_RESHAPE eddywb_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real eddynotused1_avg(GLOBAL_2D_ARRAY,N)
CSDISTRIBUTE_RESHAPE eddynotused1_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real eddynotused2_avg(GLOBAL_2D_ARRAY,N)
CSDISTRIBUTE_RESHAPE eddynotused2_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real eddynotused3_avg(GLOBAL_2D_ARRAY,N)
CSDISTRIBUTE_RESHAPE eddynotused3_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
# endif

# ifdef AVERAGES
      common /diag_eddyuu_avg/eddyuu_avg
     &       /diag_eddyvv_avg/eddyvv_avg
     &       /diag_eddyuv_avg/eddyuv_avg
     &       /diag_eddyub_avg/eddyub_avg
     &       /diag_eddyvb_avg/eddyvb_avg
     &       /diag_eddywb_avg/eddywb_avg
     &       /diag_eddynotused1_avg/eddynotused1_avg
     &       /diag_eddynotused2_avg/eddynotused2_avg
     &       /diag_eddynotused3_avg/eddynotused3_avg
# endif      
#endif /* DIAGNOSTICS_EDDY*/



 
