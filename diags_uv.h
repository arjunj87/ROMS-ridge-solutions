! This is include file "diags_vrt.h"
!  ==== == ======= ==== ==========
!

#ifdef DIAGNOSTICS_UV

      real MXadv(GLOBAL_2D_ARRAY,N,2)
CSDISTRIBUTE_RESHAPE MXadv(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real MYadv(GLOBAL_2D_ARRAY,N,2)
CSDISTRIBUTE_RESHAPE MYadv(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real MVadv(GLOBAL_2D_ARRAY,N,2)
CSDISTRIBUTE_RESHAPE MVadv(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real MCor(GLOBAL_2D_ARRAY,N,2)
CSDISTRIBUTE_RESHAPE MCor(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real MPrsgrd(GLOBAL_2D_ARRAY,N,2)
CSDISTRIBUTE_RESHAPE MPrsgrd(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real MHmix(GLOBAL_2D_ARRAY,N,2)
CSDISTRIBUTE_RESHAPE MHmix(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real MVmix(GLOBAL_2D_ARRAY,N,2)
CSDISTRIBUTE_RESHAPE MVmix(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real Mrate(GLOBAL_2D_ARRAY,N,2)
CSDISTRIBUTE_RESHAPE Mrate(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real MVmix2(GLOBAL_2D_ARRAY,N,2)
CSDISTRIBUTE_RESHAPE MVmix2(BLOCK_PATTERN,*) BLOCK_CLAUSE

# ifdef AVERAGES
      real MXadv_avg(GLOBAL_2D_ARRAY,N,2)
CSDISTRIBUTE_RESHAPE MXadv_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real MYadv_avg(GLOBAL_2D_ARRAY,N,2)
CSDISTRIBUTE_RESHAPE MYadv_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real MVadv_avg(GLOBAL_2D_ARRAY,N,2)
CSDISTRIBUTE_RESHAPE MVadv_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real MCor_avg(GLOBAL_2D_ARRAY,N,2)
CSDISTRIBUTE_RESHAPE MCor_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real MPrsgrd_avg(GLOBAL_2D_ARRAY,N,2)
CSDISTRIBUTE_RESHAPE MPrsgrd_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real MHmix_avg(GLOBAL_2D_ARRAY,N,2)
CSDISTRIBUTE_RESHAPE MHmix_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real MVmix_avg(GLOBAL_2D_ARRAY,N,2)
CSDISTRIBUTE_RESHAPE MVmix_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real Mrate_avg(GLOBAL_2D_ARRAY,N,2)
CSDISTRIBUTE_RESHAPE Mrate_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real MVmix2_avg(GLOBAL_2D_ARRAY,N,2)
CSDISTRIBUTE_RESHAPE MVmix2_avg(BLOCK_PATTERN,*) BLOCK_CLAUSE
# endif



      common /diag_MXadv/MXadv   
     &       /diag_MYadv/MYadv
     &       /diag_MVadv/MVadv  
     &       /diag_MCor/MCor
     &       /diag_MPrsgrd/MPrsgrd
     &       /diag_MHmix/MHmix
     &       /diag_MVmix/MVmix
     &       /diag_Mrate/Mrate
     &       /diag_MVmix2/MVmix2
# ifdef AVERAGES
      common /diag_MXadv_avg/MXadv_avg
     &       /diag_MYadv_avg/MYadv_avg
     &       /diag_MVadv_avg/MVadv_avg
     &       /diag_MCor_avg/MCor_avg
     &       /diag_MPrsgrd_avg/MPrsgrd_avg
     &       /diag_MHmix_avg/MHmix_avg
     &       /diag_MVmix_avg/MVmix_avg
     &       /diag_Mrate_avg/Mrate_avg
     &       /diag_MVmix2_avg/MVmix2_avg
# endif      
 

#endif /* DIAGNOSTICS_UV */



 
