! Dimensions of Physical Grid and array dimensions:
!----------- -- -------- ---- --- ----- -----------
! LLm   Number of the internal points of the PHYSICAL grid in XI-
! MMm   and ETA-directions, excluding physical side boundary points,
!       peroodic ghost points, and MPI-margins (if any). 
!
! Lm    Number of the internal points [see above] of array covering
! Mm    a single MPI-subdomain.  These are identical to LLm,MMm there 
!       is no MPI-partitioning. 
!
      integer  LLm, MMm, N
# if defined GRID_LEVEL && GRID_LEVEL < MAX_GRID_LEVEL
     &    , imin_child,imax_child, jmin_child,jmax_child
# endif
      parameter (
#if defined BASIN
     &                LLm=180, MMm=140, N=12
#elif defined ATLBIG
     &                LLm=1500, MMm=2000, N=50 
#elif defined NWAT
     &                LLm=1600, MMm=1000, N=50
#elif defined CHABU
     &                LLm=1200, MMm=1400, N=50
#elif defined SEAMOUNT
c     &               LLm=720,  MMm=864,  N=70
c     &               LLm=360,  MMm=432,  N=70
c     &               LLm=1002,  MMm=1206,  N=70
c     &               LLm=448,  MMm=576,  N=70
c     &              LLm=768,  MMm=768,  N=35
c     &              LLm=264,  MMm=264,  N=70    
#elif defined CONICAL_ISLAND
c     &               LLm=500,  MMm=600,  N=70
c     &               LLm=1440,  MMm=1728,  N=70
c     &               LLm=120,  MMm=270,  N=110
     &               LLm=300,  MMm=800,  N=110
c     &               LLm=500,  MMm=768,  N=200
c     &               LLm=576,  MMm=4,  N=200
c     &               LLm=250,  MMm=300,  N=40
#else
     &                LLm=??, MMm=??, N=??
#endif
     &                                      )
!
! Domain subdivision parameters:
!------- ----------- -----------
! NNODES          total number of MPI processes (nodes);
! NP_XI,  NP_ETA  number of MPI subdomains in XI-, ETA-directions;
! NSUB_X, NSUB_E  number of shared memory subdomains (tiles) in
!                                      XI- and ETA-directions;
      integer NSUB_X, NSUB_E
#ifdef MPI
      integer NP_XI, NP_ETA, NNODES
      parameter (NP_XI=24, NP_ETA= 8, NSUB_X=1, NSUB_E=1)
      parameter (NNODES=NP_XI*NP_ETA)
#else
        parameter (NSUB_X=6, NSUB_E=24)
#endif
!
! Array dimensions and bounds of the used portions of sub-arrays
!
      integer Lm,Mm, irc
#ifdef MPI
      parameter (Lm=(LLm+NP_XI-1)/NP_XI, Mm=(MMm+NP_ETA-1)/NP_ETA)
      logical west_inter, east_inter, south_inter, north_inter
      integer ocean_grid_comm, mynode, iwest, ieast, jsouth, jnorth,
     &                                 iSW_corn,     jSW_corn
      common /mpi_comm_vars/  ocean_grid_comm, mynode, iwest, ieast,
     &        jsouth, jnorth, iSW_corn, jSW_corn, irc,
     &        west_inter, east_inter, south_inter, north_inter
#else
      parameter (Lm=LLm, Mm=MMm, irc=1)
#endif


!
! Number of tracers and tracer identification indices:
!------- -- ------- --- ------ -------------- --------
!
#ifdef SOLVE3D
      integer NT, itemp
# ifdef SALINITY
     &          , isalt
# endif
# ifdef PASSIVE_TRACER
     &          , itpas
# endif


      parameter (itemp=1,
# ifdef SALINITY
     &           isalt=2,
# ifdef PASSIVE_TRACER     
     &           itpas=3,
     &           NT=3     
# else
     &           NT=2     
# endif   /* PASSIVE_TRACER */
# else
# ifdef PASSIVE_TRACER
     &           itpas=2,
     &           NT=2
# else
     &           NT=1
# endif     /* PASSIVE_TRACER */

# endif     /* SALINITY */
     &           )



#endif /*SOLVE3D */
#ifdef PSOURCE
      integer Msrc         ! Number of point sources, if any
      parameter (Msrc=10)
#endif
#ifdef STATIONS
      integer NS           ! Number of output stations (if any).
      parameter (NS=5)
#endif
#ifdef FLOATS
       integer Mfloats          ! Maximum number of floats
       parameter (Mfloats=32000)! ====== == ====== ========
#endif
!
! Derived dimension parameters.
!
      integer stdout, padd_X,padd_E
      parameter (stdout=6, padd_X=(Lm+2)/2-(Lm+1)/2,
     &                     padd_E=(Mm+2)/2-(Mm+1)/2)
 
