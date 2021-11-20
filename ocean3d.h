! This is include file "ocean3d.h".
!----- -- ------- ---- ------------
!
#ifdef SOLVE3D
      real u(GLOBAL_2D_ARRAY,N,3)
CSDISTRIBUTE_RESHAPE u(BLOCK_PATTERN,*,*) BLOCK_CLAUSE
      real v(GLOBAL_2D_ARRAY,N,3)
CSDISTRIBUTE_RESHAPE v(BLOCK_PATTERN,*,*) BLOCK_CLAUSE
      real t(GLOBAL_2D_ARRAY,N,3,NT)
CSDISTRIBUTE_RESHAPE t(BLOCK_PATTERN,*,*,*) BLOCK_CLAUSE
      common /ocean_u/u /ocean_v/v /ocean_t/t
      
#ifdef ADV_ISONEUTRAL            
      real diff3u(GLOBAL_2D_ARRAY,N)
CSDISTRIBUTE_RESHAPE diff3u(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real diff3v(GLOBAL_2D_ARRAY,N)
CSDISTRIBUTE_RESHAPE diff3v(BLOCK_PATTERN,*) BLOCK_CLAUSE
      common /ocean_diffu/diff3u /ocean_diffv/diff3v    
      real Sx_u(GLOBAL_2D_ARRAY,1:N)
CSDISTRIBUTE_RESHAPE Sx_u(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real Sy_v(GLOBAL_2D_ARRAY,1:N)
CSDISTRIBUTE_RESHAPE Sy_v(BLOCK_PATTERN,*) BLOCK_CLAUSE
      common /ocean_slopex/Sx_u /ocean_slopey/Sy_v
      real DrDx(GLOBAL_2D_ARRAY,1:N)
CSDISTRIBUTE_RESHAPE dRdX(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real DrDe(GLOBAL_2D_ARRAY,1:N)
CSDISTRIBUTE_RESHAPE dRde(BLOCK_PATTERN,*) BLOCK_CLAUSE
      common /ocean_drhoX/dRdX /ocean_dRhoY/dRdE      
      real DzDr(GLOBAL_2D_ARRAY,0:N)
CSDISTRIBUTE_RESHAPE dZdR(BLOCK_PATTERN,*) BLOCK_CLAUSE
      common /ocean_drhoZ/dZdR      
#endif /* ADV_ISONEUTRAL */

      real Hz(GLOBAL_2D_ARRAY,N)
CSDISTRIBUTE_RESHAPE Hz(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real Hz_bak(GLOBAL_2D_ARRAY,N)
CSDISTRIBUTE_RESHAPE Hz_bak(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real z_r(GLOBAL_2D_ARRAY,N)
CSDISTRIBUTE_RESHAPE z_r(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real z_w(GLOBAL_2D_ARRAY,0:N)
CSDISTRIBUTE_RESHAPE z_w(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real FlxU(GLOBAL_2D_ARRAY,N)
CSDISTRIBUTE_RESHAPE  FlxU(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real FlxV(GLOBAL_2D_ARRAY,N)
CSDISTRIBUTE_RESHAPE  FlxV(BLOCK_PATTERN,*) BLOCK_CLAUSE
      real W(GLOBAL_2D_ARRAY,0:N)
CSDISTRIBUTE_RESHAPE W(BLOCK_PATTERN,*) BLOCK_CLAUSE
      common /grid_Hz/Hz    /grid_zr/z_r  /grid_W/W
     &  /grid_Hz_bak/Hz_bak /grid_zw/z_w  /grid_FlxU/FlxU
     &                                    /grid_FlxV/FlxV
#endif  /* SOLVE3D */
 
