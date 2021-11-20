#include "cppdefs.h"
#ifdef MPI
 
      subroutine mpi_exchange4_tile (istr,iend,jstr,jend, A, nmax)
      implicit none
# include "mpif.h"
# include "param.h"
# include "mess_buffers.h"
# include "hidden_mpi_vars.h"
      integer istr,iend,jstr,jend, nmax
      real A(GLOBAL_2D_ARRAY,nmax)
CSDISTRIBUTE_RESHAPE A(BLOCK_PATTERN) BLOCK_CLAUSE

      integer kshft, isize,jsize, md_XI,md_ETA, i,j,k, iter,ierr,
     &        stts1(MPI_STATUS_SIZE),  stts2(MPI_STATUS_SIZE),
     &        stts3(MPI_STATUS_SIZE),  stts4(MPI_STATUS_SIZE)

! Everything is the same as in "mpi_exchange4.F", except that now
! MPI_SendRecv is used everywhere instead of paired Send-Recv and 
! Recv-Send sequences. 


# define FOUR_MESSAGES
# include "compute_message_bounds.h"

      isize=2*ishft*nmax
      jsize=2*jshft*nmax

c*    write(*,*) mynode, proc(2), ' itg,jtg =', itg,jtg
c*    write(*,'(2(6x,A,I2,2x,A,I3,2x,A,I3))')
c*   &        'inode=',inode, 'imin=',imin, 'imax=',imax,
c*   &        'jnode=',jnode, 'jmin=',jmin, 'jmax=',jmax

#define write !

      do iter=0,1
        md_ETA=mod(jnode+iter,2)
        if (md_ETA.eq.0) then
          if (SOUTH_INTER) then
            write(*,*) 'mpi_exchange4SR: 1.3', mynode
            do k=1,nmax
              kshft=2*(k-1)*ishft
              do i=imin,imax
                sendS(i-imin       +kshft)=A(i,jsouth  ,k)
                sendS(i-imin+ishft +kshft)=A(i,jsouth+1,k)
              enddo
            enddo


            call MPI_SendRecv (sendS, isize, MPI_DOUBLE_PRECISION,
     &                         p_S, itg+1,
     &                         recvS, isize, MPI_DOUBLE_PRECISION,
     &                         p_S, itg+2,
     &                              ocean_grid_comm, stts2, ierr)

            do k=1,nmax
              kshft=2*(k-1)*ishft
              do i=imin,imax
                A(i,jsouth-2,k)=recvS(i-imin       +kshft)
                A(i,jsouth-1,k)=recvS(i-imin+ishft +kshft)
              enddo
            enddo
          endif
        else
          if (NORTH_INTER) then
            write(*,*) 'mpi_exchange4SR: 1.4', mynode
            do k=1,nmax
              kshft=2*(k-1)*ishft
              do i=imin,imax
                sendN(i-imin       +kshft)=A(i,jnorth-1,k)
                sendN(i-imin+ishft +kshft)=A(i,jnorth  ,k)
              enddo
            enddo

            call MPI_SendRecv (sendN, isize, MPI_DOUBLE_PRECISION,
     &                         p_N, itg+2,
     &                         recvN, isize, MPI_DOUBLE_PRECISION,
     &                         p_N, itg+1,
     &                              ocean_grid_comm, stts1, ierr)
            do k=1,nmax
              kshft=2*(k-1)*ishft
              do i=imin,imax
                A(i,jnorth+1,k)=recvN(i-imin       +kshft)
                A(i,jnorth+2,k)=recvN(i-imin+ishft +kshft)
              enddo
            enddo
          endif
        endif
      enddo  !<-- iter


      do iter=0,1
        md_XI=mod(inode+iter,2)
        if (md_XI.eq.0) then
          if (WEST_INTER) then
            write(*,*) 'mpi_exchange4SR: 1.1', mynode
            do k=1,nmax
              kshft=2*(k-1)*jshft
              do j=jmin,jmax
                sendW(j-jmin       +kshft)=A(iwest  ,j,k)
                sendW(j-jmin+jshft +kshft)=A(iwest+1,j,k)
              enddo
            enddo

            call MPI_SendRecv (sendW, jsize, MPI_DOUBLE_PRECISION,
     &                         p_W, jtg+3,
     &                         recvW, jsize, MPI_DOUBLE_PRECISION,
     &                         p_W, jtg+4,
     &                              ocean_grid_comm, stts4, ierr)

            do k=1,nmax
              kshft=2*(k-1)*jshft
              do j=jmin,jmax
                A(iwest-2,j,k)=recvW(j-jmin       +kshft)
                A(iwest-1,j,k)=recvW(j-jmin+jshft +kshft)
             enddo
            enddo
          endif
        else
          if (EAST_INTER) then
            write(*,*) 'mpi_exchange4SR: 1.2', mynode
            do k=1,nmax
              kshft=2*(k-1)*jshft
              do j=jmin,jmax
                sendE(j-jmin       +kshft)=A(ieast-1,j,k)
                sendE(j-jmin+jshft +kshft)=A(ieast  ,j,k)
              enddo
            enddo

            call MPI_SendRecv (sendE, jsize, MPI_DOUBLE_PRECISION,
     &                         p_E, jtg+4,
     &                         recvE, jsize, MPI_DOUBLE_PRECISION,
     &                         p_E, jtg+3,
     &                              ocean_grid_comm, stts3, ierr)
            do k=1,nmax
              kshft=2*(k-1)*jshft
              do j=jmin,jmax
                A(ieast+1,j,k)=recvE(j-jmin       +kshft)
                A(ieast+2,j,k)=recvE(j-jmin+jshft +kshft)
              enddo
            enddo
          endif
        endif
      enddo

      write(*,*) 'mpi_exchange4SR out', mynode
      return
      end
#else
      subroutine mpi_exchange_empty
      end
#endif
