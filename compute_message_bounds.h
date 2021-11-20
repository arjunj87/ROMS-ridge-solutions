      integer imin,imax,jmin,jmax,  ishft,jshft, itg,jtg

! Determine east-west extents for messages exchanged in the northern
! and southern directions: 

# ifdef EW_PERIODIC
      if (NP_XI.eq.1) then            ! periodic margins must be
        imin=istr-2                   ! included into the message
        imax=iend+2                   ! if there is no partitioning
      else                            ! in XI-direction; otherwise
        imin=istr                     ! strip them out. 
        imax=iend
      endif
# else
      if (WESTERN_EDGE) then          ! include an extra point on
        imin=istr-1                   ! either side for ghost points
      else                            ! associated with PHYSICAL
        imin=istr                     ! boundaries.
      endif
      if (EASTERN_EDGE) then
        imax=iend+1
      else
        imax=iend
      endif
# endif

! Determine north-south extents for messages exchanged in the eastern
! and western directions: this leads to two variants of code, depending
! which exchange method is used.  In the case of 8-message exchange the
! corners are sent/received directly by special small messages, so
! extents are computed exactly the same way as east-west extents above.
! Conversely, the 4-message exchane is a two-stage procedure. At first,
! the messages sent/received to/from northern and southern neighbours
! are unpacked first and the data is placed into the respective ghost
! rows.  Once this is complete, messages going into eastern and western
! directions are packed in such a way that the corner points withing
! the computational margins of already filled are included into the
! eastern and western bound exchanges.  This eliminates the need for
! special corner messages and makes the north-south extents slightly
! larger than in the case of 8-message approach. 
 
# ifdef NS_PERIODIC
      if (NP_ETA.eq.1) then           ! if there is no partition in
        jmin=jstr-2                   ! ETA-direction include periodic
        jmax=jend+2                   ! margins into the message;
      else
#  ifdef FOUR_MESSAGES
        if (SOUTH_INTER) then         ! similarly, and unlike in
          jmin=jstr-2                 ! XI-direction above, include
        else                          ! computational margins
          jmin=jstr                   ! (4-message exchange only)
        endif
        if (NORTH_INTER) then
          jmax=jend+2
        else
          jmax=jend
        endif
#  else
        jmin=jstr                     ! strip out periodic margins 
        jmax=jend                     ! 8-message exchange, NP_ETA>1
#  endif
      endif
# else
      if (SOUTHERN_EDGE) then         ! add one extra row of points
        jmin=jstr-1                   ! for PHYSICAL boundaries or
#  ifdef FOUR_MESSAGES
      elseif (SOUTH_INTER) then       ! two for computational margins
        jmin=jstr-2                   ! (4-message exchange only)
#  endif
      else
        jmin=jstr
      endif
      if (NORTHERN_EDGE) then
        jmax=jend+1
#  ifdef FOUR_MESSAGES
      elseif (NORTH_INTER) then
        jmax=jend+2
#  endif
      else
        jmax=jend
      endif
# endif

      ishft=imax-imin+1               ! one-dimensional sizes 
      jshft=jmax-jmin+1               ! of messages to be send 

! Create unique labels for MPI tags to ensure that messages sent
! or received by different tiles of the same MPI subdomain can be
! distinguished.  Essentially translate tile bounds back into tile
! indices and and multiply it by 8.   The actual tags in MPI cals
! are then obtained by adding directional bias "idir", e.g.,
! tag=itg+idir or tag=jtg+idir where idir=1,2,...,8


      itg=(istr+iend-2*iwest)/(2*(iend-istr+1))   ! two-dimensional
      jtg=(jstr+jend-2*jsouth)/(2*(jend-jstr+1))  ! indices of tile

      itg=8*itg
      jtg=8*jtg

      if (nmax.gt.N) then
        itg=itg +10000
        jtg=jtg +10000
      elseif (nmax.gt.1) then
        itg=itg +1000
        jtg=jtg +1000
      endif
 
