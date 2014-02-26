!Crown Copyright 2012 AWE.
!
! This file is part of CloverLeaf.
!
! CloverLeaf is free software: you can redistribute it and/or modify it under 
! the terms of the GNU General Public License as published by the 
! Free Software Foundation, either version 3 of the License, or (at your option) 
! any later version.
!
! CloverLeaf is distributed in the hope that it will be useful, but 
! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
! FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
! details.
!
! You should have received a copy of the GNU General Public License along with 
! CloverLeaf. If not, see http://www.gnu.org/licenses/.

!>  @brief Communication Utilities
!>  @author Wayne Gaudin
!>  @details Contains all utilities required to run CloverLeaf in a distributed
!>  environment, including initialisation, mesh decompostion, reductions and
!>  halo exchange using explicit buffers.
!>
!>  Note the halo exchange is currently coded as simply as possible and no 
!>  optimisations have been implemented, such as post receives before sends or packing
!>  buffers with multiple data fields. This is intentional so the effect of these
!>  optimisations can be measured on large systems, as and when they are added.
!>
!>  Even without these modifications CloverLeaf weak scales well on moderately sized
!>  systems of the order of 10K cores.

MODULE clover_module

  USE data_module
  USE definitions_module
  !USE iso_c_binding
  USE pack_kernel_module

  IMPLICIT NONE

  INCLUDE 'mpp/shmem.fh'

  REAL(KIND=8) :: sum_total, sum_value
  REAL(KIND=8) :: min_value, min_final
  REAL(KIND=8) :: max_value, max_final
  INTEGER      :: error_value, error_final

  REAL(KIND=8) :: pWrk_sum(MAX(1/2+1, SHMEM_REDUCE_MIN_WRKDATA_SIZE))
  INTEGER :: pSync_sum(SHMEM_REDUCE_SYNC_SIZE)

  REAL(KIND=8) :: pWrk_min(MAX(1/2+1, SHMEM_REDUCE_MIN_WRKDATA_SIZE))
  INTEGER :: pSync_min(SHMEM_REDUCE_SYNC_SIZE)

  REAL(KIND=8) :: pWrk_max(MAX(1/2+1, SHMEM_REDUCE_MIN_WRKDATA_SIZE))
  INTEGER :: pSync_max(SHMEM_REDUCE_SYNC_SIZE)

  INTEGER :: pWrk_error(MAX(1/2+1, SHMEM_REDUCE_MIN_WRKDATA_SIZE))
  INTEGER :: pSync_error(SHMEM_REDUCE_SYNC_SIZE)

  INTEGER :: pSync_collect(SHMEM_COLLECT_SYNC_SIZE)

  INTEGER(KIND=4) :: left_rcv_flag, right_rcv_flag, left_write_flag, right_write_flag
  INTEGER(KIND=4) :: top_rcv_flag, bottom_rcv_flag, top_write_flag, bottom_write_flag
  INTEGER(KIND=4) :: left_top_rcv_flag, right_top_rcv_flag, right_bottom_rcv_flag, left_bottom_rcv_flag
  INTEGER(KIND=4) :: left_top_write_flag, right_top_write_flag, right_bottom_write_flag, left_bottom_write_flag

  COMMON/FLAG/left_rcv_flag, right_rcv_flag, left_write_flag, right_write_flag, top_rcv_flag, bottom_rcv_flag, & 
              top_write_flag, bottom_write_flag, left_top_rcv_flag, right_top_rcv_flag, right_bottom_rcv_flag, & 
              left_bottom_rcv_flag, left_top_write_flag, right_top_write_flag, right_bottom_write_flag, left_bottom_write_flag
CONTAINS

SUBROUTINE clover_barrier

  CALL SHMEM_BARRIER_ALL()

END SUBROUTINE clover_barrier

SUBROUTINE clover_abort

  CALL SHMEM_FINALIZE

END SUBROUTINE clover_abort

SUBROUTINE clover_finalize

  CLOSE(g_out)
  CALL FLUSH(0)
  CALL FLUSH(6)
  CALL FLUSH(g_out)

  CALL SHMEM_FINALIZE

END SUBROUTINE clover_finalize

SUBROUTINE clover_init_comms

  IMPLICIT NONE

  INTEGER :: err,rank,size

  left_rcv_flag = 0
  right_rcv_flag = 0
  left_write_flag = 1 
  right_write_flag = 1

  bottom_rcv_flag = 0
  top_rcv_flag = 0
  bottom_write_flag = 1 
  top_write_flag = 1

  left_top_rcv_flag = 0
  right_top_rcv_flag = 0
  right_bottom_rcv_flag = 0
  left_bottom_rcv_flag = 0

  left_top_write_flag = 1
  right_top_write_flag = 1
  right_bottom_write_flag = 1
  left_bottom_write_flag = 1

  rank=0
  size=1

  CALL START_PES(0)

  rank=SHMEM_MY_PE()
  size=SHMEM_N_PES()

  parallel%parallel=.TRUE.
  parallel%task=rank

  IF(rank.EQ.0) THEN
    parallel%boss=.TRUE.
  ENDIF

  parallel%boss_task=0
  parallel%max_task=size


END SUBROUTINE clover_init_comms

SUBROUTINE clover_get_num_chunks(count)

  IMPLICIT NONE

  INTEGER :: count

! Should be changed so there can be more than one chunk per mpi task

  count=parallel%max_task

END SUBROUTINE clover_get_num_chunks

SUBROUTINE clover_decompose(x_cells,y_cells,left,right,bottom,top)

  ! This decomposes the mesh into a number of chunks.
  ! The number of chunks may be a multiple of the number of mpi tasks
  ! Doesn't always return the best split if there are few factors
  ! All factors need to be stored and the best picked. But its ok for now

  IMPLICIT NONE

  INTEGER :: x_cells,y_cells,left(:),right(:),top(:),bottom(:)
  INTEGER :: c,delta_x,delta_y

  REAL(KIND=8) :: mesh_ratio,factor_x,factor_y
  INTEGER  :: chunk_x,chunk_y,mod_x,mod_y,split_found

  INTEGER  :: cx,cy,chunk,add_x,add_y,add_x_prev,add_y_prev

  ! 2D Decomposition of the mesh

  mesh_ratio=real(x_cells)/real(y_cells)

  chunk_x=number_of_chunks
  chunk_y=1

  split_found=0 ! Used to detect 1D decomposition
  DO c=1,number_of_chunks
    IF (MOD(number_of_chunks,c).EQ.0) THEN
      factor_x=number_of_chunks/real(c)
      factor_y=c
      !Compare the factor ratio with the mesh ratio
      IF(factor_x/factor_y.LE.mesh_ratio) THEN
        chunk_y=c
        chunk_x=number_of_chunks/c
        split_found=1
        EXIT
      ENDIF
    ENDIF
  ENDDO

  IF(split_found.EQ.0.OR.chunk_y.EQ.number_of_chunks) THEN ! Prime number or 1D decomp detected
    IF(mesh_ratio.GE.1.0) THEN
      chunk_x=number_of_chunks
      chunk_y=1
    ELSE
      chunk_x=1
      chunk_y=number_of_chunks
    ENDIF
  ENDIF

  delta_x=x_cells/chunk_x
  delta_y=y_cells/chunk_y
  mod_x=MOD(x_cells,chunk_x)
  mod_y=MOD(y_cells,chunk_y)

  ! Set up chunk mesh ranges and chunk connectivity

  add_x_prev=0
  add_y_prev=0
  chunk=1
  DO cy=1,chunk_y
    DO cx=1,chunk_x
      add_x=0
      add_y=0
      IF(cx.LE.mod_x)add_x=1
      IF(cy.LE.mod_y)add_y=1
      left(chunk)=(cx-1)*delta_x+1+add_x_prev
      right(chunk)=left(chunk)+delta_x-1+add_x
      bottom(chunk)=(cy-1)*delta_y+1+add_y_prev
      top(chunk)=bottom(chunk)+delta_y-1+add_y

      chunks(chunk)%chunk_neighbours(chunk_left)=chunk_x*(cy-1)+cx-1
      chunks(chunk)%chunk_neighbours(chunk_right)=chunk_x*(cy-1)+cx+1
      chunks(chunk)%chunk_neighbours(chunk_bottom)=chunk_x*(cy-2)+cx
      chunks(chunk)%chunk_neighbours(chunk_top)=chunk_x*(cy)+cx

      chunks(chunk)%chunk_neighbours(CHUNK_LEFT_TOP) = chunk_x*cy+cx-1
      chunks(chunk)%chunk_neighbours(CHUNK_LEFT_BOTTOM) = chunk_x*(cy-2)+cx-1
      chunks(chunk)%chunk_neighbours(CHUNK_RIGHT_TOP) = chunk_x*cy+cx+1
      chunks(chunk)%chunk_neighbours(CHUNK_RIGHT_BOTTOM) = chunk_x*(cy-2)+cx+1

      IF(cx.EQ.1) THEN
        chunks(chunk)%chunk_neighbours(chunk_left)=external_face
        chunks(chunk)%chunk_neighbours(CHUNK_LEFT_TOP)=external_face
        chunks(chunk)%chunk_neighbours(CHUNK_LEFT_BOTTOM)=external_face
      ENDIF
      IF(cx.EQ.chunk_x) THEN
        chunks(chunk)%chunk_neighbours(chunk_right)=external_face
        chunks(chunk)%chunk_neighbours(CHUNK_RIGHT_TOP)=external_face
        chunks(chunk)%chunk_neighbours(CHUNK_RIGHT_BOTTOM)=external_face
      ENDIF
      IF(cy.EQ.1) THEN 
        chunks(chunk)%chunk_neighbours(chunk_bottom)=external_face
        chunks(chunk)%chunk_neighbours(CHUNK_LEFT_BOTTOM)=external_face
        chunks(chunk)%chunk_neighbours(CHUNK_RIGHT_BOTTOM)=external_face
      ENDIF
      IF(cy.EQ.chunk_y) THEN 
        chunks(chunk)%chunk_neighbours(chunk_top)=external_face
        chunks(chunk)%chunk_neighbours(CHUNK_LEFT_TOP)=external_face
        chunks(chunk)%chunk_neighbours(CHUNK_RIGHT_TOP)=external_face
      ENDIF

      IF(cx.LE.mod_x)add_x_prev=add_x_prev+1
      chunk=chunk+1

    ENDDO
    add_x_prev=0
    IF(cy.LE.mod_y)add_y_prev=add_y_prev+1
  ENDDO

  num_chunks_x = chunk_x
  num_chunks_y = chunk_y

  IF(parallel%boss)THEN
    WRITE(g_out,*)
    WRITE(g_out,*)"Mesh ratio of ",mesh_ratio
    WRITE(g_out,*)"Decomposing the mesh into ",chunk_x," by ",chunk_y," chunks"
    WRITE(g_out,*)
  ENDIF

END SUBROUTINE clover_decompose

SUBROUTINE clover_allocate_buffers(chunk)

  IMPLICIT NONE

  INTEGER      :: chunk, err, idefault
  INTEGER      :: buffer_size_x, buffer_size_y, element_size, diag_buff_size
  REAL(KIND=8) :: r8default 
  
  ! Unallocated buffers for external boundaries caused issues on some systems so they are now
  !  all allocated
  IF(parallel%task.EQ.chunks(chunk)%task)THEN

    element_size = kind(r8default)/kind(idefault)
    buffer_size_x = element_size*2*(chunks(chunk)%field%x_max+5)
    buffer_size_y = element_size*2*(chunks(chunk)%field%y_max+5)
    diag_buff_size = element_size*6

    CALL SHPALLOC(density0_pls,  buffer_size_y, err, 1)
    CALL SHPALLOC(density0_plr,  buffer_size_y, err, 1)
    CALL SHPALLOC(density0_prs,  buffer_size_y, err, 1)
    CALL SHPALLOC(density0_prr,  buffer_size_y, err, 1)
    CALL SHPALLOC(density0_pbs,  buffer_size_x, err, 1)
    CALL SHPALLOC(density0_pbr,  buffer_size_x, err, 1)
    CALL SHPALLOC(density0_pts,  buffer_size_x, err, 1)
    CALL SHPALLOC(density0_ptr,  buffer_size_x, err, 1)
    CALL SHPALLOC(density0_pltr, diag_buff_size, err, 1)
    CALL SHPALLOC(density0_prtr, diag_buff_size, err, 1)
    CALL SHPALLOC(density0_plbr, diag_buff_size, err, 1)
    CALL SHPALLOC(density0_prbr, diag_buff_size, err, 1)
    CALL SHPALLOC(density0_plts, diag_buff_size, err, 1)
    CALL SHPALLOC(density0_prts, diag_buff_size, err, 1)
    CALL SHPALLOC(density0_plbs, diag_buff_size, err, 1)
    CALL SHPALLOC(density0_prbs, diag_buff_size, err, 1)

    CALL SHPALLOC(density1_pls,  buffer_size_y, err, 1)
    CALL SHPALLOC(density1_plr,  buffer_size_y, err, 1)
    CALL SHPALLOC(density1_prs,  buffer_size_y, err, 1)
    CALL SHPALLOC(density1_prr,  buffer_size_y, err, 1)
    CALL SHPALLOC(density1_pbs,  buffer_size_x, err, 1)
    CALL SHPALLOC(density1_pbr,  buffer_size_x, err, 1)
    CALL SHPALLOC(density1_pts,  buffer_size_x, err, 1)
    CALL SHPALLOC(density1_ptr,  buffer_size_x, err, 1)
    CALL SHPALLOC(density1_pltr, diag_buff_size, err, 1)
    CALL SHPALLOC(density1_prtr, diag_buff_size, err, 1)
    CALL SHPALLOC(density1_plbr, diag_buff_size, err, 1)
    CALL SHPALLOC(density1_prbr, diag_buff_size, err, 1)
    CALL SHPALLOC(density1_plts, diag_buff_size, err, 1)
    CALL SHPALLOC(density1_prts, diag_buff_size, err, 1)
    CALL SHPALLOC(density1_plbs, diag_buff_size, err, 1)
    CALL SHPALLOC(density1_prbs, diag_buff_size, err, 1)

    CALL SHPALLOC(energy0_pls,  buffer_size_y, err, 1)
    CALL SHPALLOC(energy0_plr,  buffer_size_y, err, 1)
    CALL SHPALLOC(energy0_prs,  buffer_size_y, err, 1)
    CALL SHPALLOC(energy0_prr,  buffer_size_y, err, 1)
    CALL SHPALLOC(energy0_pbs,  buffer_size_x, err, 1)
    CALL SHPALLOC(energy0_pbr,  buffer_size_x, err, 1)
    CALL SHPALLOC(energy0_pts,  buffer_size_x, err, 1)
    CALL SHPALLOC(energy0_ptr,  buffer_size_x, err, 1)
    CALL SHPALLOC(energy0_pltr, diag_buff_size, err, 1)
    CALL SHPALLOC(energy0_prtr, diag_buff_size, err, 1)
    CALL SHPALLOC(energy0_plbr, diag_buff_size, err, 1)
    CALL SHPALLOC(energy0_prbr, diag_buff_size, err, 1)
    CALL SHPALLOC(energy0_plts, diag_buff_size, err, 1)
    CALL SHPALLOC(energy0_prts, diag_buff_size, err, 1)
    CALL SHPALLOC(energy0_plbs, diag_buff_size, err, 1)
    CALL SHPALLOC(energy0_prbs, diag_buff_size, err, 1)

    CALL SHPALLOC(energy1_pls,  buffer_size_y, err, 1)
    CALL SHPALLOC(energy1_plr,  buffer_size_y, err, 1)
    CALL SHPALLOC(energy1_prs,  buffer_size_y, err, 1)
    CALL SHPALLOC(energy1_prr,  buffer_size_y, err, 1)
    CALL SHPALLOC(energy1_pbs,  buffer_size_x, err, 1)
    CALL SHPALLOC(energy1_pbr,  buffer_size_x, err, 1)
    CALL SHPALLOC(energy1_pts,  buffer_size_x, err, 1)
    CALL SHPALLOC(energy1_ptr,  buffer_size_x, err, 1)
    CALL SHPALLOC(energy1_pltr, diag_buff_size, err, 1)
    CALL SHPALLOC(energy1_prtr, diag_buff_size, err, 1)
    CALL SHPALLOC(energy1_plbr, diag_buff_size, err, 1)
    CALL SHPALLOC(energy1_prbr, diag_buff_size, err, 1)
    CALL SHPALLOC(energy1_plts, diag_buff_size, err, 1)
    CALL SHPALLOC(energy1_prts, diag_buff_size, err, 1)
    CALL SHPALLOC(energy1_plbs, diag_buff_size, err, 1)
    CALL SHPALLOC(energy1_prbs, diag_buff_size, err, 1)

    CALL SHPALLOC(pressure_pls,  buffer_size_y, err, 1)
    CALL SHPALLOC(pressure_plr,  buffer_size_y, err, 1)
    CALL SHPALLOC(pressure_prs,  buffer_size_y, err, 1)
    CALL SHPALLOC(pressure_prr,  buffer_size_y, err, 1)
    CALL SHPALLOC(pressure_pbs,  buffer_size_x, err, 1)
    CALL SHPALLOC(pressure_pbr,  buffer_size_x, err, 1)
    CALL SHPALLOC(pressure_pts,  buffer_size_x, err, 1)
    CALL SHPALLOC(pressure_ptr,  buffer_size_x, err, 1)
    CALL SHPALLOC(pressure_pltr, diag_buff_size, err, 1)
    CALL SHPALLOC(pressure_prtr, diag_buff_size, err, 1)
    CALL SHPALLOC(pressure_plbr, diag_buff_size, err, 1)
    CALL SHPALLOC(pressure_prbr, diag_buff_size, err, 1)
    CALL SHPALLOC(pressure_plts, diag_buff_size, err, 1)
    CALL SHPALLOC(pressure_prts, diag_buff_size, err, 1)
    CALL SHPALLOC(pressure_plbs, diag_buff_size, err, 1)
    CALL SHPALLOC(pressure_prbs, diag_buff_size, err, 1)

    CALL SHPALLOC(viscosity_pls,  buffer_size_y, err, 1)
    CALL SHPALLOC(viscosity_plr,  buffer_size_y, err, 1)
    CALL SHPALLOC(viscosity_prs,  buffer_size_y, err, 1)
    CALL SHPALLOC(viscosity_prr,  buffer_size_y, err, 1)
    CALL SHPALLOC(viscosity_pbs,  buffer_size_x, err, 1)
    CALL SHPALLOC(viscosity_pbr,  buffer_size_x, err, 1)
    CALL SHPALLOC(viscosity_pts,  buffer_size_x, err, 1)
    CALL SHPALLOC(viscosity_ptr,  buffer_size_x, err, 1)
    CALL SHPALLOC(viscosity_pltr, diag_buff_size, err, 1)
    CALL SHPALLOC(viscosity_prtr, diag_buff_size, err, 1)
    CALL SHPALLOC(viscosity_plbr, diag_buff_size, err, 1)
    CALL SHPALLOC(viscosity_prbr, diag_buff_size, err, 1)
    CALL SHPALLOC(viscosity_plts, diag_buff_size, err, 1)
    CALL SHPALLOC(viscosity_prts, diag_buff_size, err, 1)
    CALL SHPALLOC(viscosity_plbs, diag_buff_size, err, 1)
    CALL SHPALLOC(viscosity_prbs, diag_buff_size, err, 1)

    CALL SHPALLOC(soundspeed_pls,  buffer_size_y, err, 1)
    CALL SHPALLOC(soundspeed_plr,  buffer_size_y, err, 1)
    CALL SHPALLOC(soundspeed_prs,  buffer_size_y, err, 1)
    CALL SHPALLOC(soundspeed_prr,  buffer_size_y, err, 1)
    CALL SHPALLOC(soundspeed_pbs,  buffer_size_x, err, 1)
    CALL SHPALLOC(soundspeed_pbr,  buffer_size_x, err, 1)
    CALL SHPALLOC(soundspeed_pts,  buffer_size_x, err, 1)
    CALL SHPALLOC(soundspeed_ptr,  buffer_size_x, err, 1)
    CALL SHPALLOC(soundspeed_pltr, diag_buff_size, err, 1)
    CALL SHPALLOC(soundspeed_prtr, diag_buff_size, err, 1)
    CALL SHPALLOC(soundspeed_plbr, diag_buff_size, err, 1)
    CALL SHPALLOC(soundspeed_prbr, diag_buff_size, err, 1)
    CALL SHPALLOC(soundspeed_plts, diag_buff_size, err, 1)
    CALL SHPALLOC(soundspeed_prts, diag_buff_size, err, 1)
    CALL SHPALLOC(soundspeed_plbs, diag_buff_size, err, 1)
    CALL SHPALLOC(soundspeed_prbs, diag_buff_size, err, 1)

    CALL SHPALLOC(xvel0_pls,  buffer_size_y, err, 1)
    CALL SHPALLOC(xvel0_plr,  buffer_size_y, err, 1)
    CALL SHPALLOC(xvel0_prs,  buffer_size_y, err, 1)
    CALL SHPALLOC(xvel0_prr,  buffer_size_y, err, 1)
    CALL SHPALLOC(xvel0_pbs,  buffer_size_x, err, 1)
    CALL SHPALLOC(xvel0_pbr,  buffer_size_x, err, 1)
    CALL SHPALLOC(xvel0_pts,  buffer_size_x, err, 1)
    CALL SHPALLOC(xvel0_ptr,  buffer_size_x, err, 1)
    CALL SHPALLOC(xvel0_pltr, diag_buff_size, err, 1)
    CALL SHPALLOC(xvel0_prtr, diag_buff_size, err, 1)
    CALL SHPALLOC(xvel0_plbr, diag_buff_size, err, 1)
    CALL SHPALLOC(xvel0_prbr, diag_buff_size, err, 1)
    CALL SHPALLOC(xvel0_plts, diag_buff_size, err, 1)
    CALL SHPALLOC(xvel0_prts, diag_buff_size, err, 1)
    CALL SHPALLOC(xvel0_plbs, diag_buff_size, err, 1)
    CALL SHPALLOC(xvel0_prbs, diag_buff_size, err, 1)

    CALL SHPALLOC(xvel1_pls,  buffer_size_y, err, 1)
    CALL SHPALLOC(xvel1_plr,  buffer_size_y, err, 1)
    CALL SHPALLOC(xvel1_prs,  buffer_size_y, err, 1)
    CALL SHPALLOC(xvel1_prr,  buffer_size_y, err, 1)
    CALL SHPALLOC(xvel1_pbs,  buffer_size_x, err, 1)
    CALL SHPALLOC(xvel1_pbr,  buffer_size_x, err, 1)
    CALL SHPALLOC(xvel1_pts,  buffer_size_x, err, 1)
    CALL SHPALLOC(xvel1_ptr,  buffer_size_x, err, 1)
    CALL SHPALLOC(xvel1_pltr, diag_buff_size, err, 1)
    CALL SHPALLOC(xvel1_prtr, diag_buff_size, err, 1)
    CALL SHPALLOC(xvel1_plbr, diag_buff_size, err, 1)
    CALL SHPALLOC(xvel1_prbr, diag_buff_size, err, 1)
    CALL SHPALLOC(xvel1_plts, diag_buff_size, err, 1)
    CALL SHPALLOC(xvel1_prts, diag_buff_size, err, 1)
    CALL SHPALLOC(xvel1_plbs, diag_buff_size, err, 1)
    CALL SHPALLOC(xvel1_prbs, diag_buff_size, err, 1)

    CALL SHPALLOC(yvel0_pls,  buffer_size_y, err, 1)
    CALL SHPALLOC(yvel0_plr,  buffer_size_y, err, 1)
    CALL SHPALLOC(yvel0_prs,  buffer_size_y, err, 1)
    CALL SHPALLOC(yvel0_prr,  buffer_size_y, err, 1)
    CALL SHPALLOC(yvel0_pbs,  buffer_size_x, err, 1)
    CALL SHPALLOC(yvel0_pbr,  buffer_size_x, err, 1)
    CALL SHPALLOC(yvel0_pts,  buffer_size_x, err, 1)
    CALL SHPALLOC(yvel0_ptr,  buffer_size_x, err, 1)
    CALL SHPALLOC(yvel0_pltr, diag_buff_size, err, 1)
    CALL SHPALLOC(yvel0_prtr, diag_buff_size, err, 1)
    CALL SHPALLOC(yvel0_plbr, diag_buff_size, err, 1)
    CALL SHPALLOC(yvel0_prbr, diag_buff_size, err, 1)
    CALL SHPALLOC(yvel0_plts, diag_buff_size, err, 1)
    CALL SHPALLOC(yvel0_prts, diag_buff_size, err, 1)
    CALL SHPALLOC(yvel0_plbs, diag_buff_size, err, 1)
    CALL SHPALLOC(yvel0_prbs, diag_buff_size, err, 1)

    CALL SHPALLOC(yvel1_pls,  buffer_size_y, err, 1)
    CALL SHPALLOC(yvel1_plr,  buffer_size_y, err, 1)
    CALL SHPALLOC(yvel1_prs,  buffer_size_y, err, 1)
    CALL SHPALLOC(yvel1_prr,  buffer_size_y, err, 1)
    CALL SHPALLOC(yvel1_pbs,  buffer_size_x, err, 1)
    CALL SHPALLOC(yvel1_pbr,  buffer_size_x, err, 1)
    CALL SHPALLOC(yvel1_pts,  buffer_size_x, err, 1)
    CALL SHPALLOC(yvel1_ptr,  buffer_size_x, err, 1)
    CALL SHPALLOC(yvel1_pltr, diag_buff_size, err, 1)
    CALL SHPALLOC(yvel1_prtr, diag_buff_size, err, 1)
    CALL SHPALLOC(yvel1_plbr, diag_buff_size, err, 1)
    CALL SHPALLOC(yvel1_prbr, diag_buff_size, err, 1)
    CALL SHPALLOC(yvel1_plts, diag_buff_size, err, 1)
    CALL SHPALLOC(yvel1_prts, diag_buff_size, err, 1)
    CALL SHPALLOC(yvel1_plbs, diag_buff_size, err, 1)
    CALL SHPALLOC(yvel1_prbs, diag_buff_size, err, 1)

    CALL SHPALLOC(volflux_x_pls,  buffer_size_y, err, 1)
    CALL SHPALLOC(volflux_x_plr,  buffer_size_y, err, 1)
    CALL SHPALLOC(volflux_x_prs,  buffer_size_y, err, 1)
    CALL SHPALLOC(volflux_x_prr,  buffer_size_y, err, 1)
    CALL SHPALLOC(volflux_x_pbs,  buffer_size_x, err, 1)
    CALL SHPALLOC(volflux_x_pbr,  buffer_size_x, err, 1)
    CALL SHPALLOC(volflux_x_pts,  buffer_size_x, err, 1)
    CALL SHPALLOC(volflux_x_ptr,  buffer_size_x, err, 1)
    CALL SHPALLOC(volflux_x_pltr, diag_buff_size, err, 1)
    CALL SHPALLOC(volflux_x_prtr, diag_buff_size, err, 1)
    CALL SHPALLOC(volflux_x_plbr, diag_buff_size, err, 1)
    CALL SHPALLOC(volflux_x_prbr, diag_buff_size, err, 1)
    CALL SHPALLOC(volflux_x_plts, diag_buff_size, err, 1)
    CALL SHPALLOC(volflux_x_prts, diag_buff_size, err, 1)
    CALL SHPALLOC(volflux_x_plbs, diag_buff_size, err, 1)
    CALL SHPALLOC(volflux_x_prbs, diag_buff_size, err, 1)

    CALL SHPALLOC(volflux_y_pls,  buffer_size_y, err, 1)
    CALL SHPALLOC(volflux_y_plr,  buffer_size_y, err, 1)
    CALL SHPALLOC(volflux_y_prs,  buffer_size_y, err, 1)
    CALL SHPALLOC(volflux_y_prr,  buffer_size_y, err, 1)
    CALL SHPALLOC(volflux_y_pbs,  buffer_size_x, err, 1)
    CALL SHPALLOC(volflux_y_pbr,  buffer_size_x, err, 1)
    CALL SHPALLOC(volflux_y_pts,  buffer_size_x, err, 1)
    CALL SHPALLOC(volflux_y_ptr,  buffer_size_x, err, 1)
    CALL SHPALLOC(volflux_y_pltr, diag_buff_size, err, 1)
    CALL SHPALLOC(volflux_y_prtr, diag_buff_size, err, 1)
    CALL SHPALLOC(volflux_y_plbr, diag_buff_size, err, 1)
    CALL SHPALLOC(volflux_y_prbr, diag_buff_size, err, 1)
    CALL SHPALLOC(volflux_y_plts, diag_buff_size, err, 1)
    CALL SHPALLOC(volflux_y_prts, diag_buff_size, err, 1)
    CALL SHPALLOC(volflux_y_plbs, diag_buff_size, err, 1)
    CALL SHPALLOC(volflux_y_prbs, diag_buff_size, err, 1)

    CALL SHPALLOC(massflux_x_pls,  buffer_size_y, err, 1)
    CALL SHPALLOC(massflux_x_plr,  buffer_size_y, err, 1)
    CALL SHPALLOC(massflux_x_prs,  buffer_size_y, err, 1)
    CALL SHPALLOC(massflux_x_prr,  buffer_size_y, err, 1)
    CALL SHPALLOC(massflux_x_pbs,  buffer_size_x, err, 1)
    CALL SHPALLOC(massflux_x_pbr,  buffer_size_x, err, 1)
    CALL SHPALLOC(massflux_x_pts,  buffer_size_x, err, 1)
    CALL SHPALLOC(massflux_x_ptr,  buffer_size_x, err, 1)
    CALL SHPALLOC(massflux_x_pltr, diag_buff_size, err, 1)
    CALL SHPALLOC(massflux_x_prtr, diag_buff_size, err, 1)
    CALL SHPALLOC(massflux_x_plbr, diag_buff_size, err, 1)
    CALL SHPALLOC(massflux_x_prbr, diag_buff_size, err, 1)
    CALL SHPALLOC(massflux_x_plts, diag_buff_size, err, 1)
    CALL SHPALLOC(massflux_x_prts, diag_buff_size, err, 1)
    CALL SHPALLOC(massflux_x_plbs, diag_buff_size, err, 1)
    CALL SHPALLOC(massflux_x_prbs, diag_buff_size, err, 1)

    CALL SHPALLOC(massflux_y_pls,  buffer_size_y, err, 1)
    CALL SHPALLOC(massflux_y_plr,  buffer_size_y, err, 1)
    CALL SHPALLOC(massflux_y_prs,  buffer_size_y, err, 1)
    CALL SHPALLOC(massflux_y_prr,  buffer_size_y, err, 1)
    CALL SHPALLOC(massflux_y_pbs,  buffer_size_x, err, 1)
    CALL SHPALLOC(massflux_y_pbr,  buffer_size_x, err, 1)
    CALL SHPALLOC(massflux_y_pts,  buffer_size_x, err, 1)
    CALL SHPALLOC(massflux_y_ptr,  buffer_size_x, err, 1)
    CALL SHPALLOC(massflux_y_pltr, diag_buff_size, err, 1)
    CALL SHPALLOC(massflux_y_prtr, diag_buff_size, err, 1)
    CALL SHPALLOC(massflux_y_plbr, diag_buff_size, err, 1)
    CALL SHPALLOC(massflux_y_prbr, diag_buff_size, err, 1)
    CALL SHPALLOC(massflux_y_plts, diag_buff_size, err, 1)
    CALL SHPALLOC(massflux_y_prts, diag_buff_size, err, 1)
    CALL SHPALLOC(massflux_y_plbs, diag_buff_size, err, 1)
    CALL SHPALLOC(massflux_y_prbs, diag_buff_size, err, 1)



    WRITE(*,*) "Successfully shpalloc'd the comms buffers"

  ENDIF

END SUBROUTINE clover_allocate_buffers


SUBROUTINE clover_exchange(fields,depth)

    IMPLICIT NONE

    INTEGER      :: fields(:),depth

    CALL clover_exchange_send_async(parallel%task+1, depth, fields)

    CALL clover_exchange_receive_async(parallel%task+1, depth, fields)

END SUBROUTINE clover_exchange


SUBROUTINE clover_exchange_send_async(chunk, depth, fields)

    IMPLICIT NONE

    INTEGER :: chunk, depth, fields(NUM_FIELDS), receiver 

    IF(chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) THEN
        !wait until can write
        IF (left_write_flag .EQ. 0) THEN
          CALL SHMEM_INT4_WAIT_UNTIL(left_write_flag, SHMEM_CMP_EQ, 1)
        ENDIF
        
        !call method to write all buffers in this direction 
        CALL clover_exchange_write_all_buffers_left(chunk, depth, fields)

        left_write_flag = 0 
    ENDIF
    IF(chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) THEN
        !wait until can write
        IF (right_write_flag .EQ. 0) THEN
          CALL SHMEM_INT4_WAIT_UNTIL(right_write_flag, SHMEM_CMP_EQ, 1)
        ENDIF

        !call method to write all buffers in this direction 
        CALL clover_exchange_write_all_buffers_right(chunk, depth, fields)

        right_write_flag = 0 
    ENDIF
    IF(chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) THEN
        !wait until can write
        IF (bottom_write_flag .EQ. 0) THEN
          CALL SHMEM_INT4_WAIT_UNTIL(bottom_write_flag, SHMEM_CMP_EQ, 1)
        ENDIF

        !call method to write all buffers in this direction 
        CALL clover_exchange_write_all_buffers_bottom(chunk, depth, fields)

        bottom_write_flag = 0
    ENDIF
    IF(chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) THEN
        !wait until can write
        IF (top_write_flag .EQ. 0) THEN
            CALL SHMEM_INT4_WAIT_UNTIL(top_write_flag, SHMEM_CMP_EQ, 1)
        ENDIF

        !call method to write all buffers in this direction 
        CALL clover_exchange_write_all_buffers_top(chunk, depth, fields)

        top_write_flag = 0 
    ENDIF
    IF ( (chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) ) THEN
        !wait until can write
        IF (left_top_write_flag .EQ. 0) THEN
            CALL SHMEM_INT4_WAIT_UNTIL(left_top_write_flag, SHMEM_CMP_EQ, 1)
        ENDIF

        !call method to write all buffers in this direction 
        CALL clover_exchange_write_all_buffers_left_top(chunk, depth, fields)

        left_top_write_flag = 0
    ENDIF
    IF ( (chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) ) THEN
        !wait until can write
        IF (right_top_write_flag .EQ. 0) THEN
            CALL SHMEM_INT4_WAIT_UNTIL(right_top_write_flag, SHMEM_CMP_EQ, 1)
        ENDIF

        !call method to write all buffers in this direction 
        CALL clover_exchange_write_all_buffers_right_top(chunk, depth, fields)

        right_top_write_flag = 0
    ENDIF
    IF ( (chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) ) THEN
        !wait until can write
        IF (right_bottom_write_flag .EQ. 0) THEN
            CALL SHMEM_INT4_WAIT_UNTIL(right_bottom_write_flag, SHMEM_CMP_EQ, 1)
        ENDIF

        !call method to write all buffers in this direction 
        CALL clover_exchange_write_all_buffers_right_bottom(chunk, depth, fields)

        right_bottom_write_flag = 0
    ENDIF
    IF ( (chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) ) THEN
        !wait until can write
        IF (left_bottom_write_flag .EQ. 0) THEN
            CALL SHMEM_INT4_WAIT_UNTIL(left_bottom_write_flag, SHMEM_CMP_EQ, 1)
        ENDIF

        !call method to write all buffers in this direction 
        CALL clover_exchange_write_all_buffers_left_bottom(chunk, depth, fields)

        left_bottom_write_flag = 0
    ENDIF


#ifdef FENCE_NOT_QUIET
    CALL SHMEM_FENCE
#else
    CALL SHMEM_QUIET
#endif

    IF(chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) THEN
        receiver=chunks(chunks(chunk)%chunk_neighbours(chunk_left))%task
        CALL SHMEM_PUT4_NB(right_rcv_flag, 1, 1, receiver)
    ENDIF
    IF(chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) THEN
        receiver=chunks(chunks(chunk)%chunk_neighbours(chunk_right))%task
        CALL SHMEM_PUT4_NB(left_rcv_flag, 1, 1, receiver)
    ENDIF
    IF(chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) THEN
        receiver=chunks(chunks(chunk)%chunk_neighbours(chunk_bottom))%task
        CALL SHMEM_PUT4_NB(top_rcv_flag, 1, 1, receiver)
    ENDIF
    IF(chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) THEN
        receiver=chunks(chunks(chunk)%chunk_neighbours(chunk_top))%task
        CALL SHMEM_PUT4_NB(bottom_rcv_flag, 1, 1, receiver)
    ENDIF
    IF ( (chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) ) THEN
        receiver = chunks(chunks(chunk)%chunk_neighbours(chunk_left_top))%task
        CALL SHMEM_PUT4_NB(right_bottom_rcv_flag, 1, 1, receiver)
    ENDIF
    IF ( (chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) ) THEN
        receiver = chunks(chunks(chunk)%chunk_neighbours(chunk_right_top))%task
        CALL SHMEM_PUT4_NB(left_bottom_rcv_flag, 1, 1, receiver)
    ENDIF
    IF ( (chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) ) THEN
        receiver = chunks(chunks(chunk)%chunk_neighbours(chunk_right_bottom))%task
        CALL SHMEM_PUT4_NB(left_top_rcv_flag, 1, 1, receiver)
    ENDIF
    IF ( (chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) ) THEN
        receiver = chunks(chunks(chunk)%chunk_neighbours(chunk_left_bottom))%task
        CALL SHMEM_PUT4_NB(right_top_rcv_flag, 1, 1, receiver)
    ENDIF

END SUBROUTINE clover_exchange_send_async


SUBROUTINE clover_exchange_receive_async(chunk, depth, fields)

    IMPLICIT NONE

    INTEGER :: chunk, depth, fields(NUM_FIELDS), receiver 

    IF(chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) THEN
        !wait for flag to be set
        IF (left_rcv_flag .EQ. 0) THEN
            CALL SHMEM_INT4_WAIT_UNTIL(left_rcv_flag, SHMEM_CMP_EQ, 1)
        ENDIF
        !unpack the buffer
        CALL clover_exchange_unpack_all_buffers_left(chunk, depth, fields)
        left_rcv_flag = 0
    ENDIF
    IF(chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) THEN
        !wait for flag to be set
        IF (right_rcv_flag .EQ. 0) THEN
            CALL SHMEM_INT4_WAIT_UNTIL(right_rcv_flag, SHMEM_CMP_EQ, 1)
        ENDIF
        !unpack the buffer
        CALL clover_exchange_unpack_all_buffers_right(chunk, depth, fields)
        right_rcv_flag = 0
    ENDIF
    IF(chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) THEN
        !wait for flag to be set
        IF (bottom_rcv_flag .EQ. 0) THEN
            CALL SHMEM_INT4_WAIT_UNTIL(bottom_rcv_flag, SHMEM_CMP_EQ, 1)
        ENDIF
        !unpack the buffer
        CALL clover_exchange_unpack_all_buffers_bottom(chunk, depth, fields)
        bottom_rcv_flag = 0
    ENDIF
    IF(chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) THEN
        !wait for flag to be set
        IF (top_rcv_flag .EQ. 0) THEN
            CALL SHMEM_INT4_WAIT_UNTIL(top_rcv_flag, SHMEM_CMP_EQ, 1)
        ENDIF
        !unpack the buffer
        CALL clover_exchange_unpack_all_buffers_top(chunk, depth, fields)
        top_rcv_flag = 0
    ENDIF
    IF ( (chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) ) THEN
        !wait for flag to be set
        IF (left_top_rcv_flag .EQ. 0) THEN
            CALL SHMEM_INT4_WAIT_UNTIL(left_top_rcv_flag, SHMEM_CMP_EQ, 1)
        ENDIF
        !unpack the buffer
        CALL clover_exchange_unpack_all_buffers_left_top(chunk, depth, fields)
        left_top_rcv_flag = 0
    ENDIF
    IF ( (chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) ) THEN
        !wait for flag to be set
        IF (right_top_rcv_flag .EQ. 0) THEN
            CALL SHMEM_INT4_WAIT_UNTIL(right_top_rcv_flag, SHMEM_CMP_EQ, 1)
        ENDIF
        !unpack the buffer
        CALL clover_exchange_unpack_all_buffers_right_top(chunk, depth, fields)
        right_top_rcv_flag = 0
    ENDIF
    IF ( (chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) ) THEN
        !wait for flag to be set
        IF (right_bottom_rcv_flag .EQ. 0) THEN
            CALL SHMEM_INT4_WAIT_UNTIL(right_bottom_rcv_flag, SHMEM_CMP_EQ, 1)
        ENDIF
        !unpack the buffer
        CALL clover_exchange_unpack_all_buffers_right_bottom(chunk, depth, fields)
        right_bottom_rcv_flag = 0
    ENDIF
    IF ( (chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) ) THEN
        !wait for flag to be set
        IF (left_bottom_rcv_flag .EQ. 0) THEN
            CALL SHMEM_INT4_WAIT_UNTIL(left_bottom_rcv_flag, SHMEM_CMP_EQ, 1)
        ENDIF
        !unpack the buffer
        CALL clover_exchange_unpack_all_buffers_left_bottom(chunk, depth, fields)
        left_bottom_rcv_flag = 0
    ENDIF


    IF(chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) THEN
        receiver=chunks(chunks(chunk)%chunk_neighbours(chunk_bottom))%task
        CALL SHMEM_PUT4_NB(top_write_flag, 1, 1, receiver)
    ENDIF
    IF(chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) THEN
        receiver=chunks(chunks(chunk)%chunk_neighbours(chunk_top))%task
        CALL SHMEM_PUT4_NB(bottom_write_flag, 1, 1, receiver)
    ENDIF
    IF(chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) THEN
        receiver=chunks(chunks(chunk)%chunk_neighbours(chunk_left))%task
        CALL SHMEM_PUT4_NB(right_write_flag, 1, 1, receiver)
    ENDIF
    IF(chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) THEN
        receiver=chunks(chunks(chunk)%chunk_neighbours(chunk_right))%task
        CALL SHMEM_PUT4_NB(left_write_flag, 1, 1, receiver)
    ENDIF
    IF ( (chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) ) THEN
        receiver = chunks(chunks(chunk)%chunk_neighbours(chunk_left_top))%task
        CALL SHMEM_PUT4_NB(right_bottom_write_flag, 1, 1, receiver)
    ENDIF
    IF ( (chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) ) THEN
        receiver = chunks(chunks(chunk)%chunk_neighbours(chunk_right_top))%task
        CALL SHMEM_PUT4_NB(left_bottom_write_flag, 1, 1, receiver)
    ENDIF
    IF ( (chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) ) THEN
        receiver = chunks(chunks(chunk)%chunk_neighbours(chunk_right_bottom))%task
        CALL SHMEM_PUT4_NB(left_top_write_flag, 1, 1, receiver)
    ENDIF
    IF ( (chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) ) THEN
        receiver = chunks(chunks(chunk)%chunk_neighbours(chunk_left_bottom))%task
        CALL SHMEM_PUT4_NB(right_top_write_flag, 1, 1, receiver)
    ENDIF

END SUBROUTINE clover_exchange_receive_async



SUBROUTINE clover_exchange_write_all_buffers_left(chunk, depth, fields)

    IMPLICIT NONE 

    INTEGER :: chunk, depth, receiver, topedge, bottomedge, fields(NUM_FIELDS)
    
    receiver=chunks(chunks(chunk)%chunk_neighbours(chunk_left))%task

    topedge = 0
    bottomedge = 0
    IF (chunks(chunk)%chunk_neighbours(chunk_top).EQ.external_face) THEN
        topedge = depth
    ENDIF
    IF (chunks(chunk)%chunk_neighbours(chunk_bottom).EQ.external_face) THEN
        bottomedge = depth
    ENDIF

    IF(fields(FIELD_DENSITY0).EQ.1) THEN
        CALL clover_exchange_write_message_left(chunk, depth, receiver, topedge, bottomedge, CELL_DATA, & 
                                                density0_left_snd_buffer, density0_right_rcv_buffer, chunks(chunk)%field%density0)
    ENDIF
    IF(fields(FIELD_DENSITY1).EQ.1) THEN
        CALL clover_exchange_write_message_left(chunk, depth, receiver, topedge, bottomedge, CELL_DATA, & 
                                                density1_left_snd_buffer, density1_right_rcv_buffer, chunks(chunk)%field%density1)
    ENDIF
    IF(fields(FIELD_ENERGY0).EQ.1) THEN
        CALL clover_exchange_write_message_left(chunk, depth, receiver, topedge, bottomedge, CELL_DATA, & 
                                                energy0_left_snd_buffer, energy0_right_rcv_buffer, chunks(chunk)%field%energy0)
    ENDIF
    IF(fields(FIELD_ENERGY1).EQ.1) THEN
        CALL clover_exchange_write_message_left(chunk, depth, receiver, topedge, bottomedge, CELL_DATA, & 
                                                energy1_left_snd_buffer, energy1_right_rcv_buffer, chunks(chunk)%field%energy1)
    ENDIF
    IF(fields(FIELD_PRESSURE).EQ.1) THEN
        CALL clover_exchange_write_message_left(chunk, depth, receiver, topedge, bottomedge, CELL_DATA, & 
                                                pressure_left_snd_buffer, pressure_right_rcv_buffer, chunks(chunk)%field%pressure)
    ENDIF
    IF(fields(FIELD_VISCOSITY).EQ.1) THEN
        CALL clover_exchange_write_message_left(chunk, depth, receiver, topedge, bottomedge, CELL_DATA, & 
                                                viscosity_left_snd_buffer, viscosity_right_rcv_buffer, chunks(chunk)%field%viscosity)
    ENDIF
    IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
        CALL clover_exchange_write_message_left(chunk, depth, receiver, topedge, bottomedge, CELL_DATA, & 
                                                soundspeed_left_snd_buffer, soundspeed_right_rcv_buffer, chunks(chunk)%field%soundspeed)
    ENDIF
    IF(fields(FIELD_XVEL0).EQ.1) THEN
        CALL clover_exchange_write_message_left(chunk, depth, receiver, topedge, bottomedge, VERTEX_DATA, & 
                                                xvel0_left_snd_buffer, xvel0_right_rcv_buffer, chunks(chunk)%field%xvel0)
    ENDIF
    IF(fields(FIELD_XVEL1).EQ.1) THEN
        CALL clover_exchange_write_message_left(chunk, depth, receiver, topedge, bottomedge, VERTEX_DATA, & 
                                                xvel1_left_snd_buffer, xvel1_right_rcv_buffer, chunks(chunk)%field%xvel1)
    ENDIF
    IF(fields(FIELD_YVEL0).EQ.1) THEN
        CALL clover_exchange_write_message_left(chunk, depth, receiver, topedge, bottomedge, VERTEX_DATA, & 
                                                yvel0_left_snd_buffer, yvel0_right_rcv_buffer, chunks(chunk)%field%yvel0)
    ENDIF
    IF(fields(FIELD_YVEL1).EQ.1) THEN
        CALL clover_exchange_write_message_left(chunk, depth, receiver, topedge, bottomedge, VERTEX_DATA, & 
                                                yvel1_left_snd_buffer, yvel1_right_rcv_buffer, chunks(chunk)%field%yvel1)
    ENDIF
    IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
        CALL clover_exchange_write_message_left(chunk, depth, receiver, topedge, bottomedge, X_FACE_DATA, & 
                                                volflux_x_left_snd_buffer, volflux_x_right_rcv_buffer, chunks(chunk)%field%vol_flux_x)
    ENDIF
    IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
        CALL clover_exchange_write_message_left(chunk, depth, receiver, topedge, bottomedge, Y_FACE_DATA, & 
                                                volflux_y_left_snd_buffer, volflux_y_right_rcv_buffer, chunks(chunk)%field%vol_flux_y)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
        CALL clover_exchange_write_message_left(chunk, depth, receiver, topedge, bottomedge, X_FACE_DATA, & 
                                                massflux_x_left_snd_buffer, massflux_x_right_rcv_buffer, chunks(chunk)%field%mass_flux_x)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
        CALL clover_exchange_write_message_left(chunk, depth, receiver, topedge, bottomedge, Y_FACE_DATA, & 
                                                massflux_y_left_snd_buffer, massflux_y_right_rcv_buffer, chunks(chunk)%field%mass_flux_y)
    ENDIF

END SUBROUTINE clover_exchange_write_all_buffers_left

SUBROUTINE clover_exchange_write_all_buffers_right(chunk, depth, fields)

    IMPLICIT NONE 

    INTEGER :: chunk, depth, receiver, topedge, bottomedge, fields(NUM_FIELDS)

    receiver=chunks(chunks(chunk)%chunk_neighbours(chunk_right))%task

    topedge = 0
    bottomedge = 0
    IF (chunks(chunk)%chunk_neighbours(chunk_top).EQ.external_face) THEN
        topedge = depth
    ENDIF
    IF (chunks(chunk)%chunk_neighbours(chunk_bottom).EQ.external_face) THEN
        bottomedge = depth
    ENDIF

    IF(fields(FIELD_DENSITY0).EQ.1) THEN
        CALL clover_exchange_write_message_right(chunk, depth, receiver, topedge, bottomedge, CELL_DATA, & 
                                                 density0_right_snd_buffer, density0_left_rcv_buffer, chunks(chunk)%field%density0)
    ENDIF
    IF(fields(FIELD_DENSITY1).EQ.1) THEN
        CALL clover_exchange_write_message_right(chunk, depth, receiver, topedge, bottomedge, CELL_DATA, & 
                                                 density1_right_snd_buffer, density1_left_rcv_buffer, chunks(chunk)%field%density1)
    ENDIF
    IF(fields(FIELD_ENERGY0).EQ.1) THEN
        CALL clover_exchange_write_message_right(chunk, depth, receiver, topedge, bottomedge, CELL_DATA, & 
                                                 energy0_right_snd_buffer, energy0_left_rcv_buffer, chunks(chunk)%field%energy0)
    ENDIF
    IF(fields(FIELD_ENERGY1).EQ.1) THEN
        CALL clover_exchange_write_message_right(chunk, depth, receiver, topedge, bottomedge, CELL_DATA, & 
                                                 energy1_right_snd_buffer, energy1_left_rcv_buffer, chunks(chunk)%field%energy1)
    ENDIF
    IF(fields(FIELD_PRESSURE).EQ.1) THEN
        CALL clover_exchange_write_message_right(chunk, depth, receiver, topedge, bottomedge, CELL_DATA, & 
                                                 pressure_right_snd_buffer, pressure_left_rcv_buffer, chunks(chunk)%field%pressure)
    ENDIF
    IF(fields(FIELD_VISCOSITY).EQ.1) THEN
        CALL clover_exchange_write_message_right(chunk, depth, receiver, topedge, bottomedge, CELL_DATA, & 
                                                 viscosity_right_snd_buffer, viscosity_left_rcv_buffer, chunks(chunk)%field%viscosity)
    ENDIF
    IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
        CALL clover_exchange_write_message_right(chunk, depth, receiver, topedge, bottomedge, CELL_DATA, & 
                                                 soundspeed_right_snd_buffer, soundspeed_left_rcv_buffer, chunks(chunk)%field%soundspeed)
    ENDIF
    IF(fields(FIELD_XVEL0).EQ.1) THEN
        CALL clover_exchange_write_message_right(chunk, depth, receiver, topedge, bottomedge, VERTEX_DATA, & 
                                                 xvel0_right_snd_buffer, xvel0_left_rcv_buffer, chunks(chunk)%field%xvel0)
    ENDIF
    IF(fields(FIELD_XVEL1).EQ.1) THEN
        CALL clover_exchange_write_message_right(chunk, depth, receiver, topedge, bottomedge, VERTEX_DATA, & 
                                                 xvel1_right_snd_buffer, xvel1_left_rcv_buffer, chunks(chunk)%field%xvel1)
    ENDIF
    IF(fields(FIELD_YVEL0).EQ.1) THEN
        CALL clover_exchange_write_message_right(chunk, depth, receiver, topedge, bottomedge, VERTEX_DATA, & 
                                                 yvel0_right_snd_buffer, yvel0_left_rcv_buffer, chunks(chunk)%field%yvel0)
    ENDIF
    IF(fields(FIELD_YVEL1).EQ.1) THEN
        CALL clover_exchange_write_message_right(chunk, depth, receiver, topedge, bottomedge, VERTEX_DATA, & 
                                                 yvel1_right_snd_buffer, yvel1_left_rcv_buffer, chunks(chunk)%field%yvel1)
    ENDIF
    IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
        CALL clover_exchange_write_message_right(chunk, depth, receiver, topedge, bottomedge, X_FACE_DATA, & 
                                                 volflux_x_right_snd_buffer, volflux_x_left_rcv_buffer, chunks(chunk)%field%vol_flux_x)
    ENDIF
    IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
        CALL clover_exchange_write_message_right(chunk, depth, receiver, topedge, bottomedge, Y_FACE_DATA, & 
                                                 volflux_y_right_snd_buffer, volflux_y_left_rcv_buffer, chunks(chunk)%field%vol_flux_y)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
        CALL clover_exchange_write_message_right(chunk, depth, receiver, topedge, bottomedge, X_FACE_DATA, & 
                                                 massflux_x_right_snd_buffer, massflux_x_left_rcv_buffer, chunks(chunk)%field%mass_flux_x)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
        CALL clover_exchange_write_message_right(chunk, depth, receiver, topedge, bottomedge, Y_FACE_DATA, & 
                                                 massflux_y_right_snd_buffer, massflux_y_left_rcv_buffer, chunks(chunk)%field%mass_flux_y)
    ENDIF

END SUBROUTINE clover_exchange_write_all_buffers_right

SUBROUTINE clover_exchange_write_all_buffers_bottom(chunk, depth, fields)

    IMPLICIT NONE 

    INTEGER :: chunk, depth, receiver, topedge, bottomedge, fields(NUM_FIELDS), leftedge, rightedge

    receiver=chunks(chunks(chunk)%chunk_neighbours(chunk_bottom))%task

    leftedge= 0
    rightedge= 0
    IF (chunks(chunk)%chunk_neighbours(chunk_left).EQ.external_face) THEN
        leftedge = depth
    ENDIF
    IF (chunks(chunk)%chunk_neighbours(chunk_right).EQ.external_face) THEN
        rightedge = depth
    ENDIF

    IF(fields(FIELD_DENSITY0).EQ.1) THEN
        CALL clover_exchange_write_message_bottom(chunk, depth, receiver, leftedge, rightedge, CELL_DATA, & 
                                                  density0_bottom_snd_buffer, density0_top_rcv_buffer, chunks(chunk)%field%density0)
    ENDIF
    IF(fields(FIELD_DENSITY1).EQ.1) THEN
        CALL clover_exchange_write_message_bottom(chunk, depth, receiver, leftedge, rightedge, CELL_DATA, & 
                                                  density1_bottom_snd_buffer, density1_top_rcv_buffer, chunks(chunk)%field%density1)
    ENDIF
    IF(fields(FIELD_ENERGY0).EQ.1) THEN
        CALL clover_exchange_write_message_bottom(chunk, depth, receiver, leftedge, rightedge, CELL_DATA, & 
                                                  energy0_bottom_snd_buffer, energy0_top_rcv_buffer, chunks(chunk)%field%energy0)
    ENDIF
    IF(fields(FIELD_ENERGY1).EQ.1) THEN
        CALL clover_exchange_write_message_bottom(chunk, depth, receiver, leftedge, rightedge, CELL_DATA, & 
                                                  energy1_bottom_snd_buffer, energy1_top_rcv_buffer, chunks(chunk)%field%energy1)
    ENDIF
    IF(fields(FIELD_PRESSURE).EQ.1) THEN
        CALL clover_exchange_write_message_bottom(chunk, depth, receiver, leftedge, rightedge, CELL_DATA, & 
                                                  pressure_bottom_snd_buffer, pressure_top_rcv_buffer, chunks(chunk)%field%pressure)
    ENDIF
    IF(fields(FIELD_VISCOSITY).EQ.1) THEN
        CALL clover_exchange_write_message_bottom(chunk, depth, receiver, leftedge, rightedge, CELL_DATA, & 
                                                  viscosity_bottom_snd_buffer, viscosity_top_rcv_buffer, chunks(chunk)%field%viscosity)
    ENDIF
    IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
        CALL clover_exchange_write_message_bottom(chunk, depth, receiver, leftedge, rightedge, CELL_DATA, & 
                                                  soundspeed_bottom_snd_buffer, soundspeed_top_rcv_buffer, chunks(chunk)%field%soundspeed)
    ENDIF
    IF(fields(FIELD_XVEL0).EQ.1) THEN
        CALL clover_exchange_write_message_bottom(chunk, depth, receiver, leftedge, rightedge, VERTEX_DATA, & 
                                                  xvel0_bottom_snd_buffer, xvel0_top_rcv_buffer, chunks(chunk)%field%xvel0)
    ENDIF
    IF(fields(FIELD_XVEL1).EQ.1) THEN
        CALL clover_exchange_write_message_bottom(chunk, depth, receiver, leftedge, rightedge, VERTEX_DATA, & 
                                                  xvel1_bottom_snd_buffer, xvel1_top_rcv_buffer, chunks(chunk)%field%xvel1)
    ENDIF
    IF(fields(FIELD_YVEL0).EQ.1) THEN
        CALL clover_exchange_write_message_bottom(chunk, depth, receiver, leftedge, rightedge, VERTEX_DATA, & 
                                                  yvel0_bottom_snd_buffer, yvel0_top_rcv_buffer, chunks(chunk)%field%yvel0)
    ENDIF
    IF(fields(FIELD_YVEL1).EQ.1) THEN
        CALL clover_exchange_write_message_bottom(chunk, depth, receiver, leftedge, rightedge, VERTEX_DATA, & 
                                                  yvel1_bottom_snd_buffer, yvel1_top_rcv_buffer, chunks(chunk)%field%yvel1)
    ENDIF
    IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
        CALL clover_exchange_write_message_bottom(chunk, depth, receiver, leftedge, rightedge, X_FACE_DATA, & 
                                                  volflux_x_bottom_snd_buffer, volflux_x_top_rcv_buffer, chunks(chunk)%field%vol_flux_x)
    ENDIF
    IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
        CALL clover_exchange_write_message_bottom(chunk, depth, receiver, leftedge, rightedge, Y_FACE_DATA, & 
                                                  volflux_y_bottom_snd_buffer, volflux_y_top_rcv_buffer, chunks(chunk)%field%vol_flux_y)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
        CALL clover_exchange_write_message_bottom(chunk, depth, receiver, leftedge, rightedge, X_FACE_DATA, & 
                                                  massflux_x_bottom_snd_buffer, massflux_x_top_rcv_buffer, chunks(chunk)%field%mass_flux_x)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
        CALL clover_exchange_write_message_bottom(chunk, depth, receiver, leftedge, rightedge, Y_FACE_DATA, & 
                                                  massflux_y_bottom_snd_buffer, massflux_y_top_rcv_buffer, chunks(chunk)%field%mass_flux_y)
    ENDIF

END SUBROUTINE clover_exchange_write_all_buffers_bottom


SUBROUTINE clover_exchange_write_all_buffers_top(chunk, depth, fields)

    IMPLICIT NONE 

    INTEGER :: chunk, depth, receiver, topedge, bottomedge, fields(NUM_FIELDS), leftedge, rightedge

    receiver=chunks(chunks(chunk)%chunk_neighbours(chunk_top))%task

    leftedge= 0
    rightedge= 0
    IF (chunks(chunk)%chunk_neighbours(chunk_left).EQ.external_face) THEN
        leftedge = depth
    ENDIF
    IF (chunks(chunk)%chunk_neighbours(chunk_right).EQ.external_face) THEN
        rightedge = depth
    ENDIF

    IF(fields(FIELD_DENSITY0).EQ.1) THEN
        CALL clover_exchange_write_message_top(chunk, depth, receiver, leftedge, rightedge, CELL_DATA, & 
                                               density0_top_snd_buffer, density0_bottom_rcv_buffer, chunks(chunk)%field%density0)
    ENDIF
    IF(fields(FIELD_DENSITY1).EQ.1) THEN
        CALL clover_exchange_write_message_top(chunk, depth, receiver, leftedge, rightedge, CELL_DATA, & 
                                               density1_top_snd_buffer, density1_bottom_rcv_buffer, chunks(chunk)%field%density1)
    ENDIF
    IF(fields(FIELD_ENERGY0).EQ.1) THEN
        CALL clover_exchange_write_message_top(chunk, depth, receiver, leftedge, rightedge, CELL_DATA, & 
                                               energy0_top_snd_buffer, energy0_bottom_rcv_buffer, chunks(chunk)%field%energy0)
    ENDIF
    IF(fields(FIELD_ENERGY1).EQ.1) THEN
        CALL clover_exchange_write_message_top(chunk, depth, receiver, leftedge, rightedge, CELL_DATA, & 
                                               energy1_top_snd_buffer, energy1_bottom_rcv_buffer, chunks(chunk)%field%energy1)
    ENDIF
    IF(fields(FIELD_PRESSURE).EQ.1) THEN
        CALL clover_exchange_write_message_top(chunk, depth, receiver, leftedge, rightedge, CELL_DATA, & 
                                               pressure_top_snd_buffer, pressure_bottom_rcv_buffer, chunks(chunk)%field%pressure)
    ENDIF
    IF(fields(FIELD_VISCOSITY).EQ.1) THEN
        CALL clover_exchange_write_message_top(chunk, depth, receiver, leftedge, rightedge, CELL_DATA, & 
                                               viscosity_top_snd_buffer, viscosity_bottom_rcv_buffer, chunks(chunk)%field%viscosity)
    ENDIF
    IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
        CALL clover_exchange_write_message_top(chunk, depth, receiver, leftedge, rightedge, CELL_DATA, & 
                                               soundspeed_top_snd_buffer, soundspeed_bottom_rcv_buffer, chunks(chunk)%field%soundspeed)
    ENDIF
    IF(fields(FIELD_XVEL0).EQ.1) THEN
        CALL clover_exchange_write_message_top(chunk, depth, receiver, leftedge, rightedge, VERTEX_DATA, & 
                                               xvel0_top_snd_buffer, xvel0_bottom_rcv_buffer, chunks(chunk)%field%xvel0)
    ENDIF
    IF(fields(FIELD_XVEL1).EQ.1) THEN
        CALL clover_exchange_write_message_top(chunk, depth, receiver, leftedge, rightedge, VERTEX_DATA, & 
                                               xvel1_top_snd_buffer, xvel1_bottom_rcv_buffer, chunks(chunk)%field%xvel1)
    ENDIF
    IF(fields(FIELD_YVEL0).EQ.1) THEN
        CALL clover_exchange_write_message_top(chunk, depth, receiver, leftedge, rightedge, VERTEX_DATA, & 
                                               yvel0_top_snd_buffer, yvel0_bottom_rcv_buffer, chunks(chunk)%field%yvel0)
    ENDIF
    IF(fields(FIELD_YVEL1).EQ.1) THEN
        CALL clover_exchange_write_message_top(chunk, depth, receiver, leftedge, rightedge, VERTEX_DATA, & 
                                               yvel1_top_snd_buffer, yvel1_bottom_rcv_buffer, chunks(chunk)%field%yvel1)
    ENDIF
    IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
        CALL clover_exchange_write_message_top(chunk, depth, receiver, leftedge, rightedge, X_FACE_DATA, & 
                                               volflux_x_top_snd_buffer, volflux_x_bottom_rcv_buffer, chunks(chunk)%field%vol_flux_x)
    ENDIF
    IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
        CALL clover_exchange_write_message_top(chunk, depth, receiver, leftedge, rightedge, Y_FACE_DATA, & 
                                               volflux_y_top_snd_buffer, volflux_y_bottom_rcv_buffer, chunks(chunk)%field%vol_flux_y)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
        CALL clover_exchange_write_message_top(chunk, depth, receiver, leftedge, rightedge, X_FACE_DATA, & 
                                               massflux_x_top_snd_buffer, massflux_x_bottom_rcv_buffer, chunks(chunk)%field%mass_flux_x)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
        CALL clover_exchange_write_message_top(chunk, depth, receiver, leftedge, rightedge, Y_FACE_DATA, & 
                                               massflux_y_top_snd_buffer, massflux_y_bottom_rcv_buffer, chunks(chunk)%field%mass_flux_y)
    ENDIF

END SUBROUTINE clover_exchange_write_all_buffers_top


SUBROUTINE clover_exchange_write_all_buffers_left_top(chunk, depth, fields)

    IMPLICIT NONE 

    INTEGER :: size, chunk, depth, receiver, topedge, bottomedge, fields(NUM_FIELDS)

    receiver = chunks(chunks(chunk)%chunk_neighbours(chunk_left_top))%task

    size = depth*depth

    IF(fields(FIELD_DENSITY0).EQ.1) THEN
        CALL clover_exchange_write_message_left_top(chunk, depth, receiver, size, CELL_DATA, & 
                                                    density0_left_top_snd_buffer, density0_right_bottom_rcv_buffer, chunks(chunk)%field%density0)
    ENDIF
    IF(fields(FIELD_DENSITY1).EQ.1) THEN
        CALL clover_exchange_write_message_left_top(chunk, depth, receiver, size, CELL_DATA, & 
                                                    density1_left_top_snd_buffer, density1_right_bottom_rcv_buffer, chunks(chunk)%field%density1)
    ENDIF
    IF(fields(FIELD_ENERGY0).EQ.1) THEN
        CALL clover_exchange_write_message_left_top(chunk, depth, receiver, size, CELL_DATA, & 
                                                    energy0_left_top_snd_buffer, energy0_right_bottom_rcv_buffer, chunks(chunk)%field%energy0)
    ENDIF
    IF(fields(FIELD_ENERGY1).EQ.1) THEN
        CALL clover_exchange_write_message_left_top(chunk, depth, receiver, size, CELL_DATA, & 
                                                    energy1_left_top_snd_buffer, energy1_right_bottom_rcv_buffer, chunks(chunk)%field%energy1)
    ENDIF
    IF(fields(FIELD_PRESSURE).EQ.1) THEN
        CALL clover_exchange_write_message_left_top(chunk, depth, receiver, size, CELL_DATA, & 
                                                    pressure_left_top_snd_buffer, pressure_right_bottom_rcv_buffer, chunks(chunk)%field%pressure)
    ENDIF
    IF(fields(FIELD_VISCOSITY).EQ.1) THEN
        CALL clover_exchange_write_message_left_top(chunk, depth, receiver, size, CELL_DATA, & 
                                                    viscosity_left_top_snd_buffer, viscosity_right_bottom_rcv_buffer, chunks(chunk)%field%viscosity)
    ENDIF
    IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
        CALL clover_exchange_write_message_left_top(chunk, depth, receiver, size, CELL_DATA, & 
                                                    soundspeed_left_top_snd_buffer, soundspeed_right_bottom_rcv_buffer, chunks(chunk)%field%soundspeed)
    ENDIF
    IF(fields(FIELD_XVEL0).EQ.1) THEN
        CALL clover_exchange_write_message_left_top(chunk, depth, receiver, size, VERTEX_DATA, & 
                                                    xvel0_left_top_snd_buffer, xvel0_right_bottom_rcv_buffer, chunks(chunk)%field%xvel0)
    ENDIF
    IF(fields(FIELD_XVEL1).EQ.1) THEN
        CALL clover_exchange_write_message_left_top(chunk, depth, receiver, size, VERTEX_DATA, & 
                                                    xvel1_left_top_snd_buffer, xvel1_right_bottom_rcv_buffer, chunks(chunk)%field%xvel1)
    ENDIF
    IF(fields(FIELD_YVEL0).EQ.1) THEN
        CALL clover_exchange_write_message_left_top(chunk, depth, receiver, size, VERTEX_DATA, & 
                                                    yvel0_left_top_snd_buffer, yvel0_right_bottom_rcv_buffer, chunks(chunk)%field%yvel0)
    ENDIF
    IF(fields(FIELD_YVEL1).EQ.1) THEN
        CALL clover_exchange_write_message_left_top(chunk, depth, receiver, size, VERTEX_DATA, & 
                                                    yvel1_left_top_snd_buffer, yvel1_right_bottom_rcv_buffer, chunks(chunk)%field%yvel1)
    ENDIF
    IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
        CALL clover_exchange_write_message_left_top(chunk, depth, receiver, size, X_FACE_DATA, & 
                                                    volflux_x_left_top_snd_buffer, volflux_x_right_bottom_rcv_buffer, chunks(chunk)%field%vol_flux_x)
    ENDIF
    IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
        CALL clover_exchange_write_message_left_top(chunk, depth, receiver, size, Y_FACE_DATA, & 
                                                    volflux_y_left_top_snd_buffer, volflux_y_right_bottom_rcv_buffer, chunks(chunk)%field%vol_flux_y)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
        CALL clover_exchange_write_message_left_top(chunk, depth, receiver, size, X_FACE_DATA, & 
                                                    massflux_x_left_top_snd_buffer, massflux_x_right_bottom_rcv_buffer, chunks(chunk)%field%mass_flux_x)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
        CALL clover_exchange_write_message_left_top(chunk, depth, receiver, size, Y_FACE_DATA, & 
                                                    massflux_y_left_top_snd_buffer, massflux_y_right_bottom_rcv_buffer, chunks(chunk)%field%mass_flux_y)
    ENDIF

END SUBROUTINE clover_exchange_write_all_buffers_left_top

SUBROUTINE clover_exchange_write_all_buffers_right_top(chunk, depth, fields)

    IMPLICIT NONE 

    INTEGER :: size, chunk, depth, receiver, topedge, bottomedge, fields(NUM_FIELDS)

    receiver = chunks(chunks(chunk)%chunk_neighbours(chunk_right_top))%task
    size = depth*depth

    IF(fields(FIELD_DENSITY0).EQ.1) THEN
        CALL clover_exchange_write_message_right_top(chunk, depth, receiver, size, CELL_DATA, & 
                                                     density0_right_top_snd_buffer, density0_left_bottom_rcv_buffer, chunks(chunk)%field%density0)
    ENDIF
    IF(fields(FIELD_DENSITY1).EQ.1) THEN
        CALL clover_exchange_write_message_right_top(chunk, depth, receiver, size, CELL_DATA, & 
                                                     density1_right_top_snd_buffer, density1_left_bottom_rcv_buffer, chunks(chunk)%field%density1)
    ENDIF
    IF(fields(FIELD_ENERGY0).EQ.1) THEN
        CALL clover_exchange_write_message_right_top(chunk, depth, receiver, size, CELL_DATA, & 
                                                     energy0_right_top_snd_buffer, energy0_left_bottom_rcv_buffer, chunks(chunk)%field%energy0)
    ENDIF
    IF(fields(FIELD_ENERGY1).EQ.1) THEN
        CALL clover_exchange_write_message_right_top(chunk, depth, receiver, size, CELL_DATA, & 
                                                     energy1_right_top_snd_buffer, energy1_left_bottom_rcv_buffer, chunks(chunk)%field%energy1)
    ENDIF
    IF(fields(FIELD_PRESSURE).EQ.1) THEN
        CALL clover_exchange_write_message_right_top(chunk, depth, receiver, size, CELL_DATA, & 
                                                     pressure_right_top_snd_buffer, pressure_left_bottom_rcv_buffer, chunks(chunk)%field%pressure)
    ENDIF
    IF(fields(FIELD_VISCOSITY).EQ.1) THEN
        CALL clover_exchange_write_message_right_top(chunk, depth, receiver, size, CELL_DATA, & 
                                                     viscosity_right_top_snd_buffer, viscosity_left_bottom_rcv_buffer, chunks(chunk)%field%viscosity)
    ENDIF
    IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
        CALL clover_exchange_write_message_right_top(chunk, depth, receiver, size, CELL_DATA, & 
                                                     soundspeed_right_top_snd_buffer, soundspeed_left_bottom_rcv_buffer, chunks(chunk)%field%soundspeed)
    ENDIF
    IF(fields(FIELD_XVEL0).EQ.1) THEN
        CALL clover_exchange_write_message_right_top(chunk, depth, receiver, size, VERTEX_DATA, & 
                                                     xvel0_right_top_snd_buffer, xvel0_left_bottom_rcv_buffer, chunks(chunk)%field%xvel0)
    ENDIF
    IF(fields(FIELD_XVEL1).EQ.1) THEN
        CALL clover_exchange_write_message_right_top(chunk, depth, receiver, size, VERTEX_DATA, & 
                                                     xvel1_right_top_snd_buffer, xvel1_left_bottom_rcv_buffer, chunks(chunk)%field%xvel1)
    ENDIF
    IF(fields(FIELD_YVEL0).EQ.1) THEN
        CALL clover_exchange_write_message_right_top(chunk, depth, receiver, size, VERTEX_DATA, & 
                                                     yvel0_right_top_snd_buffer, yvel0_left_bottom_rcv_buffer, chunks(chunk)%field%yvel0)
    ENDIF
    IF(fields(FIELD_YVEL1).EQ.1) THEN
        CALL clover_exchange_write_message_right_top(chunk, depth, receiver, size, VERTEX_DATA, & 
                                                     yvel1_right_top_snd_buffer, yvel1_left_bottom_rcv_buffer, chunks(chunk)%field%yvel1)
    ENDIF
    IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
        CALL clover_exchange_write_message_right_top(chunk, depth, receiver, size, X_FACE_DATA, & 
                                                     volflux_x_right_top_snd_buffer, volflux_x_left_bottom_rcv_buffer, chunks(chunk)%field%vol_flux_x)
    ENDIF
    IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
        CALL clover_exchange_write_message_right_top(chunk, depth, receiver, size, Y_FACE_DATA, & 
                                                     volflux_y_right_top_snd_buffer, volflux_y_left_bottom_rcv_buffer, chunks(chunk)%field%vol_flux_y)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
        CALL clover_exchange_write_message_right_top(chunk, depth, receiver, size, X_FACE_DATA, & 
                                                     massflux_x_right_top_snd_buffer, massflux_x_left_bottom_rcv_buffer, chunks(chunk)%field%mass_flux_x)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
        CALL clover_exchange_write_message_right_top(chunk, depth, receiver, size, Y_FACE_DATA, & 
                                                     massflux_y_right_top_snd_buffer, massflux_y_left_bottom_rcv_buffer, chunks(chunk)%field%mass_flux_y)
    ENDIF

END SUBROUTINE clover_exchange_write_all_buffers_right_top

SUBROUTINE clover_exchange_write_all_buffers_right_bottom(chunk, depth, fields)

    IMPLICIT NONE 

    INTEGER :: size, chunk, depth, receiver, topedge, bottomedge, fields(NUM_FIELDS)

    receiver = chunks(chunks(chunk)%chunk_neighbours(chunk_right_bottom))%task
    size = depth*depth

    IF(fields(FIELD_DENSITY0).EQ.1) THEN
        CALL clover_exchange_write_message_right_bottom(chunk, depth, receiver, size, CELL_DATA, & 
                                                        density0_right_bottom_snd_buffer, density0_left_top_rcv_buffer, chunks(chunk)%field%density0)
    ENDIF
    IF(fields(FIELD_DENSITY1).EQ.1) THEN
        CALL clover_exchange_write_message_right_bottom(chunk, depth, receiver, size, CELL_DATA, & 
                                                        density1_right_bottom_snd_buffer, density1_left_top_rcv_buffer, chunks(chunk)%field%density1)
    ENDIF
    IF(fields(FIELD_ENERGY0).EQ.1) THEN
        CALL clover_exchange_write_message_right_bottom(chunk, depth, receiver, size, CELL_DATA, & 
                                                        energy0_right_bottom_snd_buffer, energy0_left_top_rcv_buffer, chunks(chunk)%field%energy0)
    ENDIF
    IF(fields(FIELD_ENERGY1).EQ.1) THEN
        CALL clover_exchange_write_message_right_bottom(chunk, depth, receiver, size, CELL_DATA, & 
                                                        energy1_right_bottom_snd_buffer, energy1_left_top_rcv_buffer, chunks(chunk)%field%energy1)
    ENDIF
    IF(fields(FIELD_PRESSURE).EQ.1) THEN
        CALL clover_exchange_write_message_right_bottom(chunk, depth, receiver, size, CELL_DATA, & 
                                                        pressure_right_bottom_snd_buffer, pressure_left_top_rcv_buffer, chunks(chunk)%field%pressure)
    ENDIF
    IF(fields(FIELD_VISCOSITY).EQ.1) THEN
        CALL clover_exchange_write_message_right_bottom(chunk, depth, receiver, size, CELL_DATA, & 
                                                        viscosity_right_bottom_snd_buffer, viscosity_left_top_rcv_buffer, chunks(chunk)%field%viscosity)
    ENDIF
    IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
        CALL clover_exchange_write_message_right_bottom(chunk, depth, receiver, size, CELL_DATA, & 
                                                        soundspeed_right_bottom_snd_buffer, soundspeed_left_top_rcv_buffer, chunks(chunk)%field%soundspeed)
    ENDIF
    IF(fields(FIELD_XVEL0).EQ.1) THEN
        CALL clover_exchange_write_message_right_bottom(chunk, depth, receiver, size, VERTEX_DATA, & 
                                                        xvel0_right_bottom_snd_buffer, xvel0_left_top_rcv_buffer, chunks(chunk)%field%xvel0)
    ENDIF
    IF(fields(FIELD_XVEL1).EQ.1) THEN
        CALL clover_exchange_write_message_right_bottom(chunk, depth, receiver, size, VERTEX_DATA, & 
                                                        xvel1_right_bottom_snd_buffer, xvel1_left_top_rcv_buffer, chunks(chunk)%field%xvel1)
    ENDIF
    IF(fields(FIELD_YVEL0).EQ.1) THEN
        CALL clover_exchange_write_message_right_bottom(chunk, depth, receiver, size, VERTEX_DATA, & 
                                                        yvel0_right_bottom_snd_buffer, yvel0_left_top_rcv_buffer, chunks(chunk)%field%yvel0)
    ENDIF
    IF(fields(FIELD_YVEL1).EQ.1) THEN
        CALL clover_exchange_write_message_right_bottom(chunk, depth, receiver, size, VERTEX_DATA, & 
                                                        yvel1_right_bottom_snd_buffer, yvel1_left_top_rcv_buffer, chunks(chunk)%field%yvel1)
    ENDIF
    IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
        CALL clover_exchange_write_message_right_bottom(chunk, depth, receiver, size, X_FACE_DATA, & 
                                                        volflux_x_right_bottom_snd_buffer, volflux_x_left_top_rcv_buffer, chunks(chunk)%field%vol_flux_x)
    ENDIF
    IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
        CALL clover_exchange_write_message_right_bottom(chunk, depth, receiver, size, Y_FACE_DATA, & 
                                                        volflux_y_right_bottom_snd_buffer, volflux_y_left_top_rcv_buffer, chunks(chunk)%field%vol_flux_y)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
        CALL clover_exchange_write_message_right_bottom(chunk, depth, receiver, size, X_FACE_DATA, & 
                                                        massflux_x_right_bottom_snd_buffer, massflux_x_left_top_rcv_buffer, chunks(chunk)%field%mass_flux_x)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
        CALL clover_exchange_write_message_right_bottom(chunk, depth, receiver, size, Y_FACE_DATA, & 
                                                        massflux_y_right_bottom_snd_buffer, massflux_y_left_top_rcv_buffer, chunks(chunk)%field%mass_flux_y)
    ENDIF

END SUBROUTINE clover_exchange_write_all_buffers_right_bottom

SUBROUTINE clover_exchange_write_all_buffers_left_bottom(chunk, depth, fields)

    IMPLICIT NONE 

    INTEGER :: size, chunk, depth, receiver, topedge, bottomedge, fields(NUM_FIELDS)

    receiver = chunks(chunks(chunk)%chunk_neighbours(chunk_left_bottom))%task
    size = depth*depth

    IF(fields(FIELD_DENSITY0).EQ.1) THEN
        CALL clover_exchange_write_message_left_bottom(chunk, depth, receiver, size, CELL_DATA, & 
                                                       density0_left_bottom_snd_buffer, density0_right_top_rcv_buffer, chunks(chunk)%field%density0)
    ENDIF
    IF(fields(FIELD_DENSITY1).EQ.1) THEN
        CALL clover_exchange_write_message_left_bottom(chunk, depth, receiver, size, CELL_DATA, & 
                                                       density1_left_bottom_snd_buffer, density1_right_top_rcv_buffer, chunks(chunk)%field%density1)
    ENDIF
    IF(fields(FIELD_ENERGY0).EQ.1) THEN
        CALL clover_exchange_write_message_left_bottom(chunk, depth, receiver, size, CELL_DATA, & 
                                                       energy0_left_bottom_snd_buffer, energy0_right_top_rcv_buffer, chunks(chunk)%field%energy0)
    ENDIF
    IF(fields(FIELD_ENERGY1).EQ.1) THEN
        CALL clover_exchange_write_message_left_bottom(chunk, depth, receiver, size, CELL_DATA, & 
                                                       energy1_left_bottom_snd_buffer, energy1_right_top_rcv_buffer, chunks(chunk)%field%energy1)
    ENDIF
    IF(fields(FIELD_PRESSURE).EQ.1) THEN
        CALL clover_exchange_write_message_left_bottom(chunk, depth, receiver, size, CELL_DATA, & 
                                                       pressure_left_bottom_snd_buffer, pressure_right_top_rcv_buffer, chunks(chunk)%field%pressure)
    ENDIF
    IF(fields(FIELD_VISCOSITY).EQ.1) THEN
        CALL clover_exchange_write_message_left_bottom(chunk, depth, receiver, size, CELL_DATA, & 
                                                       viscosity_left_bottom_snd_buffer, viscosity_right_top_rcv_buffer, chunks(chunk)%field%viscosity)
    ENDIF
    IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
        CALL clover_exchange_write_message_left_bottom(chunk, depth, receiver, size, CELL_DATA, & 
                                                       soundspeed_left_bottom_snd_buffer, soundspeed_right_top_rcv_buffer, chunks(chunk)%field%soundspeed)
    ENDIF
    IF(fields(FIELD_XVEL0).EQ.1) THEN
        CALL clover_exchange_write_message_left_bottom(chunk, depth, receiver, size, VERTEX_DATA, & 
                                                       xvel0_left_bottom_snd_buffer, xvel0_right_top_rcv_buffer, chunks(chunk)%field%xvel0)
    ENDIF
    IF(fields(FIELD_XVEL1).EQ.1) THEN
        CALL clover_exchange_write_message_left_bottom(chunk, depth, receiver, size, VERTEX_DATA, & 
                                                       xvel1_left_bottom_snd_buffer, xvel1_right_top_rcv_buffer, chunks(chunk)%field%xvel1)
    ENDIF
    IF(fields(FIELD_YVEL0).EQ.1) THEN
        CALL clover_exchange_write_message_left_bottom(chunk, depth, receiver, size, VERTEX_DATA, & 
                                                       yvel0_left_bottom_snd_buffer, yvel0_right_top_rcv_buffer, chunks(chunk)%field%yvel0)
    ENDIF
    IF(fields(FIELD_YVEL1).EQ.1) THEN
        CALL clover_exchange_write_message_left_bottom(chunk, depth, receiver, size, VERTEX_DATA, & 
                                                       yvel1_left_bottom_snd_buffer, yvel1_right_top_rcv_buffer, chunks(chunk)%field%yvel1)
    ENDIF
    IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
        CALL clover_exchange_write_message_left_bottom(chunk, depth, receiver, size, X_FACE_DATA, & 
                                                       volflux_x_left_bottom_snd_buffer, volflux_x_right_top_rcv_buffer, chunks(chunk)%field%vol_flux_x)
    ENDIF
    IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
        CALL clover_exchange_write_message_left_bottom(chunk, depth, receiver, size, Y_FACE_DATA, & 
                                                       volflux_y_left_bottom_snd_buffer, volflux_y_right_top_rcv_buffer, chunks(chunk)%field%vol_flux_y)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
        CALL clover_exchange_write_message_left_bottom(chunk, depth, receiver, size, X_FACE_DATA, & 
                                                       massflux_x_left_bottom_snd_buffer, massflux_x_right_top_rcv_buffer, chunks(chunk)%field%mass_flux_x)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
        CALL clover_exchange_write_message_left_bottom(chunk, depth, receiver, size, Y_FACE_DATA, & 
                                                       massflux_y_left_bottom_snd_buffer, massflux_y_right_top_rcv_buffer, chunks(chunk)%field%mass_flux_y)
    ENDIF

END SUBROUTINE clover_exchange_write_all_buffers_left_bottom





!SUBROUTINE clover_exchange(fields,depth)
!
!  IMPLICIT NONE
!
!  INTEGER      :: fields(:),depth
!
!  ! Assuming 1 patch per task, this will be changed
!  ! Also, not packing all fields for each communication, doing one at a time
!
!  IF(fields(FIELD_DENSITY0).EQ.1) THEN
!    CALL clover_exchange_message(parallel%task+1,chunks(parallel%task+1)%field%density0,      &
!                                 left_snd_buffer,                                             &
!                                 left_rcv_buffer,                                             &
!                                 right_snd_buffer,                                            &
!                                 right_rcv_buffer,                                            &
!                                 bottom_snd_buffer,                                           &
!                                 bottom_rcv_buffer,                                           &
!                                 top_snd_buffer,                                              &
!                                 top_rcv_buffer,                                              &
!                                 left_top_snd_buffer,                                         &
!                                 left_top_rcv_buffer,                                         &
!                                 right_top_snd_buffer,                                        &
!                                 right_top_rcv_buffer,                                        &
!                                 right_bottom_snd_buffer,                                     &            
!                                 right_bottom_rcv_buffer,                                     &
!                                 left_bottom_snd_buffer,                                      &
!                                 left_bottom_rcv_buffer,                                      &
!                                 depth,CELL_DATA)
!  ENDIF
!
!  IF(fields(FIELD_DENSITY1).EQ.1) THEN
!    CALL clover_exchange_message(parallel%task+1,chunks(parallel%task+1)%field%density1,      &
!                                 left_snd_buffer,                                             &
!                                 left_rcv_buffer,                                             &
!                                 right_snd_buffer,                                            &
!                                 right_rcv_buffer,                                            &
!                                 bottom_snd_buffer,                                           &
!                                 bottom_rcv_buffer,                                           &
!                                 top_snd_buffer,                                              &
!                                 top_rcv_buffer,                                              &
!                                 left_top_snd_buffer,                                         &
!                                 left_top_rcv_buffer,                                         &
!                                 right_top_snd_buffer,                                        &
!                                 right_top_rcv_buffer,                                        &
!                                 right_bottom_snd_buffer,                                     &            
!                                 right_bottom_rcv_buffer,                                     &
!                                 left_bottom_snd_buffer,                                      &
!                                 left_bottom_rcv_buffer,                                      &
!                                 depth,CELL_DATA)
!  ENDIF
!
!  IF(fields(FIELD_ENERGY0).EQ.1) THEN
!    CALL clover_exchange_message(parallel%task+1,chunks(parallel%task+1)%field%energy0,       &
!                                 left_snd_buffer,                                             &
!                                 left_rcv_buffer,                                             &
!                                 right_snd_buffer,                                            &
!                                 right_rcv_buffer,                                            &
!                                 bottom_snd_buffer,                                           &
!                                 bottom_rcv_buffer,                                           &
!                                 top_snd_buffer,                                              &
!                                 top_rcv_buffer,                                              &
!                                 left_top_snd_buffer,                                         &
!                                 left_top_rcv_buffer,                                         &
!                                 right_top_snd_buffer,                                        &
!                                 right_top_rcv_buffer,                                        &
!                                 right_bottom_snd_buffer,                                     &            
!                                 right_bottom_rcv_buffer,                                     &
!                                 left_bottom_snd_buffer,                                      &
!                                 left_bottom_rcv_buffer,                                      &
!                                 depth,CELL_DATA)
!  ENDIF
!
!  IF(fields(FIELD_ENERGY1).EQ.1) THEN
!    CALL clover_exchange_message(parallel%task+1,chunks(parallel%task+1)%field%energy1,       &
!                                 left_snd_buffer,                                             &
!                                 left_rcv_buffer,                                             &
!                                 right_snd_buffer,                                            &
!                                 right_rcv_buffer,                                            &
!                                 bottom_snd_buffer,                                           &
!                                 bottom_rcv_buffer,                                           &
!                                 top_snd_buffer,                                              &
!                                 top_rcv_buffer,                                              &
!                                 left_top_snd_buffer,                                         &
!                                 left_top_rcv_buffer,                                         &
!                                 right_top_snd_buffer,                                        &
!                                 right_top_rcv_buffer,                                        &
!                                 right_bottom_snd_buffer,                                     &            
!                                 right_bottom_rcv_buffer,                                     &
!                                 left_bottom_snd_buffer,                                      &
!                                 left_bottom_rcv_buffer,                                      &
!                                 depth,CELL_DATA)
!  ENDIF
!
!  IF(fields(FIELD_PRESSURE).EQ.1) THEN
!    CALL clover_exchange_message(parallel%task+1,chunks(parallel%task+1)%field%pressure,      &
!                                 left_snd_buffer,                                             &
!                                 left_rcv_buffer,                                             &
!                                 right_snd_buffer,                                            &
!                                 right_rcv_buffer,                                            &
!                                 bottom_snd_buffer,                                           &
!                                 bottom_rcv_buffer,                                           &
!                                 top_snd_buffer,                                              &
!                                 top_rcv_buffer,                                              &
!                                 left_top_snd_buffer,                                         &
!                                 left_top_rcv_buffer,                                         &
!                                 right_top_snd_buffer,                                        &
!                                 right_top_rcv_buffer,                                        &
!                                 right_bottom_snd_buffer,                                     &            
!                                 right_bottom_rcv_buffer,                                     &
!                                 left_bottom_snd_buffer,                                      &
!                                 left_bottom_rcv_buffer,                                      &
!                                 depth,CELL_DATA)
!  ENDIF
!
!  IF(fields(FIELD_VISCOSITY).EQ.1) THEN
!    CALL clover_exchange_message(parallel%task+1,chunks(parallel%task+1)%field%viscosity,     &
!                                 left_snd_buffer,                                             &
!                                 left_rcv_buffer,                                             &
!                                 right_snd_buffer,                                            &
!                                 right_rcv_buffer,                                            &
!                                 bottom_snd_buffer,                                           &
!                                 bottom_rcv_buffer,                                           &
!                                 top_snd_buffer,                                              &
!                                 top_rcv_buffer,                                              &
!                                 left_top_snd_buffer,                                         &
!                                 left_top_rcv_buffer,                                         &
!                                 right_top_snd_buffer,                                        &
!                                 right_top_rcv_buffer,                                        &
!                                 right_bottom_snd_buffer,                                     &            
!                                 right_bottom_rcv_buffer,                                     &
!                                 left_bottom_snd_buffer,                                      &
!                                 left_bottom_rcv_buffer,                                      &
!                                 depth,CELL_DATA)
!  ENDIF
!
!  IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
!    CALL clover_exchange_message(parallel%task+1,chunks(parallel%task+1)%field%soundspeed,    &
!                                 left_snd_buffer,                                             &
!                                 left_rcv_buffer,                                             &
!                                 right_snd_buffer,                                            &
!                                 right_rcv_buffer,                                            &
!                                 bottom_snd_buffer,                                           &
!                                 bottom_rcv_buffer,                                           &
!                                 top_snd_buffer,                                              &
!                                 top_rcv_buffer,                                              &
!                                 left_top_snd_buffer,                                         &
!                                 left_top_rcv_buffer,                                         &
!                                 right_top_snd_buffer,                                        &
!                                 right_top_rcv_buffer,                                        &
!                                 right_bottom_snd_buffer,                                     &            
!                                 right_bottom_rcv_buffer,                                     &
!                                 left_bottom_snd_buffer,                                      &
!                                 left_bottom_rcv_buffer,                                      &
!                                 depth,CELL_DATA)
!  ENDIF
!
!  IF(fields(FIELD_XVEL0).EQ.1) THEN
!    CALL clover_exchange_message(parallel%task+1,chunks(parallel%task+1)%field%xvel0,         &
!                                 left_snd_buffer,                                             &
!                                 left_rcv_buffer,                                             &
!                                 right_snd_buffer,                                            &
!                                 right_rcv_buffer,                                            &
!                                 bottom_snd_buffer,                                           &
!                                 bottom_rcv_buffer,                                           &
!                                 top_snd_buffer,                                              &
!                                 top_rcv_buffer,                                              &
!                                 left_top_snd_buffer,                                         &
!                                 left_top_rcv_buffer,                                         &
!                                 right_top_snd_buffer,                                        &
!                                 right_top_rcv_buffer,                                        &
!                                 right_bottom_snd_buffer,                                     &            
!                                 right_bottom_rcv_buffer,                                     &
!                                 left_bottom_snd_buffer,                                      &
!                                 left_bottom_rcv_buffer,                                      &
!                                 depth,VERTEX_DATA)
!  ENDIF
!
!  IF(fields(FIELD_XVEL1).EQ.1) THEN
!    CALL clover_exchange_message(parallel%task+1,chunks(parallel%task+1)%field%xvel1,         &
!                                 left_snd_buffer,                                             &
!                                 left_rcv_buffer,                                             &
!                                 right_snd_buffer,                                            &
!                                 right_rcv_buffer,                                            &
!                                 bottom_snd_buffer,                                           &
!                                 bottom_rcv_buffer,                                           &
!                                 top_snd_buffer,                                              &
!                                 top_rcv_buffer,                                              &
!                                 left_top_snd_buffer,                                         &
!                                 left_top_rcv_buffer,                                         &
!                                 right_top_snd_buffer,                                        &
!                                 right_top_rcv_buffer,                                        &
!                                 right_bottom_snd_buffer,                                     &            
!                                 right_bottom_rcv_buffer,                                     &
!                                 left_bottom_snd_buffer,                                      &
!                                 left_bottom_rcv_buffer,                                      &
!                                 depth,VERTEX_DATA)
!  ENDIF
!
!  IF(fields(FIELD_YVEL0).EQ.1) THEN
!    CALL clover_exchange_message(parallel%task+1,chunks(parallel%task+1)%field%yvel0,         &
!                                 left_snd_buffer,                                             &
!                                 left_rcv_buffer,                                             &
!                                 right_snd_buffer,                                            &
!                                 right_rcv_buffer,                                            &
!                                 bottom_snd_buffer,                                           &
!                                 bottom_rcv_buffer,                                           &
!                                 top_snd_buffer,                                              &
!                                 top_rcv_buffer,                                              &
!                                 left_top_snd_buffer,                                         &
!                                 left_top_rcv_buffer,                                         &
!                                 right_top_snd_buffer,                                        &
!                                 right_top_rcv_buffer,                                        &
!                                 right_bottom_snd_buffer,                                     &            
!                                 right_bottom_rcv_buffer,                                     &
!                                 left_bottom_snd_buffer,                                      &
!                                 left_bottom_rcv_buffer,                                      &
!                                 depth,VERTEX_DATA)
!  ENDIF
!
!  IF(fields(FIELD_YVEL1).EQ.1) THEN
!    CALL clover_exchange_message(parallel%task+1,chunks(parallel%task+1)%field%yvel1,         &
!                                 left_snd_buffer,                                             &
!                                 left_rcv_buffer,                                             &
!                                 right_snd_buffer,                                            &
!                                 right_rcv_buffer,                                            &
!                                 bottom_snd_buffer,                                           &
!                                 bottom_rcv_buffer,                                           &
!                                 top_snd_buffer,                                              &
!                                 top_rcv_buffer,                                              &
!                                 left_top_snd_buffer,                                         &
!                                 left_top_rcv_buffer,                                         &
!                                 right_top_snd_buffer,                                        &
!                                 right_top_rcv_buffer,                                        &
!                                 right_bottom_snd_buffer,                                     &            
!                                 right_bottom_rcv_buffer,                                     &
!                                 left_bottom_snd_buffer,                                      &
!                                 left_bottom_rcv_buffer,                                      &
!                                 depth,VERTEX_DATA)
!  ENDIF
!
!  IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
!    CALL clover_exchange_message(parallel%task+1,chunks(parallel%task+1)%field%vol_flux_x,    &
!                                 left_snd_buffer,                                             &
!                                 left_rcv_buffer,                                             &
!                                 right_snd_buffer,                                            &
!                                 right_rcv_buffer,                                            &
!                                 bottom_snd_buffer,                                           &
!                                 bottom_rcv_buffer,                                           &
!                                 top_snd_buffer,                                              &
!                                 top_rcv_buffer,                                              &
!                                 left_top_snd_buffer,                                         &
!                                 left_top_rcv_buffer,                                         &
!                                 right_top_snd_buffer,                                        &
!                                 right_top_rcv_buffer,                                        &
!                                 right_bottom_snd_buffer,                                     &            
!                                 right_bottom_rcv_buffer,                                     &
!                                 left_bottom_snd_buffer,                                      &
!                                 left_bottom_rcv_buffer,                                      &
!                                 depth,X_FACE_DATA)
!  ENDIF
!
!  IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
!    CALL clover_exchange_message(parallel%task+1,chunks(parallel%task+1)%field%vol_flux_y,    &
!                                 left_snd_buffer,                                             &
!                                 left_rcv_buffer,                                             &
!                                 right_snd_buffer,                                            &
!                                 right_rcv_buffer,                                            &
!                                 bottom_snd_buffer,                                           &
!                                 bottom_rcv_buffer,                                           &
!                                 top_snd_buffer,                                              &
!                                 top_rcv_buffer,                                              &
!                                 left_top_snd_buffer,                                         &
!                                 left_top_rcv_buffer,                                         &
!                                 right_top_snd_buffer,                                        &
!                                 right_top_rcv_buffer,                                        &
!                                 right_bottom_snd_buffer,                                     &            
!                                 right_bottom_rcv_buffer,                                     &
!                                 left_bottom_snd_buffer,                                      &
!                                 left_bottom_rcv_buffer,                                      &
!                                 depth,Y_FACE_DATA)
!  ENDIF
!
!  IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
!    CALL clover_exchange_message(parallel%task+1,chunks(parallel%task+1)%field%mass_flux_x,   &
!                                 left_snd_buffer,                                             &
!                                 left_rcv_buffer,                                             &
!                                 right_snd_buffer,                                            &
!                                 right_rcv_buffer,                                            &
!                                 bottom_snd_buffer,                                           &
!                                 bottom_rcv_buffer,                                           &
!                                 top_snd_buffer,                                              &
!                                 top_rcv_buffer,                                              &
!                                 left_top_snd_buffer,                                         &
!                                 left_top_rcv_buffer,                                         &
!                                 right_top_snd_buffer,                                        &
!                                 right_top_rcv_buffer,                                        &
!                                 right_bottom_snd_buffer,                                     &            
!                                 right_bottom_rcv_buffer,                                     &
!                                 left_bottom_snd_buffer,                                      &
!                                 left_bottom_rcv_buffer,                                      &
!                                 depth,X_FACE_DATA)
!  ENDIF
!
!  IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
!    CALL clover_exchange_message(parallel%task+1,chunks(parallel%task+1)%field%mass_flux_y,   &
!                                 left_snd_buffer,                                             &
!                                 left_rcv_buffer,                                             &
!                                 right_snd_buffer,                                            &
!                                 right_rcv_buffer,                                            &
!                                 bottom_snd_buffer,                                           &
!                                 bottom_rcv_buffer,                                           &
!                                 top_snd_buffer,                                              &
!                                 top_rcv_buffer,                                              &
!                                 left_top_snd_buffer,                                         &
!                                 left_top_rcv_buffer,                                         &
!                                 right_top_snd_buffer,                                        &
!                                 right_top_rcv_buffer,                                        &
!                                 right_bottom_snd_buffer,                                     &            
!                                 right_bottom_rcv_buffer,                                     &
!                                 left_bottom_snd_buffer,                                      &
!                                 left_bottom_rcv_buffer,                                      &
!                                 depth,Y_FACE_DATA)
!  ENDIF
!
!  CALL SHMEM_BARRIER_ALL()
!
!
!END SUBROUTINE clover_exchange


SUBROUTINE clover_exchange_write_message_left(chunk, depth, receiver, topedge, bottomedge, field_type, left_snd_buffer, right_rcv_buffer, field)

    IMPLICIT NONE 

    INTEGER :: chunk, depth, receiver, size, field_type, x_inc, y_inc, topedge, bottomedge
    REAL(KIND=8) :: left_snd_buffer(:), right_rcv_buffer(:) 
    REAL(KIND=8) :: field(-1:,-1:)

    IF(field_type.EQ.CELL_DATA) THEN
        x_inc=0
        y_inc=0
    ENDIF
    IF(field_type.EQ.VERTEX_DATA) THEN
        x_inc=1
        y_inc=1
    ENDIF
    IF(field_type.EQ.X_FACE_DATA) THEN
        x_inc=1
        y_inc=0
    ENDIF
    IF(field_type.EQ.Y_FACE_DATA) THEN
        x_inc=0
        y_inc=1
    ENDIF


    size=(1+(chunks(chunk)%field%y_max+y_inc+topedge)-(chunks(chunk)%field%y_min-bottomedge))*depth

    ! pack buffer 
    CALL pack_left_buffer_seq(chunk, depth, x_inc, y_inc, left_snd_buffer, field)

    ! send buffer 
    CALL SHMEM_PUT64_NB(right_rcv_buffer, left_snd_buffer, size, receiver)

END SUBROUTINE clover_exchange_write_message_left

SUBROUTINE clover_exchange_write_message_right(chunk, depth, receiver, topedge, bottomedge, field_type, right_snd_buffer, left_rcv_buffer, field)

    IMPLICIT NONE 

    INTEGER :: chunk, depth, receiver, size, field_type, x_inc, y_inc, topedge, bottomedge
    REAL(KIND=8) :: right_snd_buffer(:), left_rcv_buffer(:) 
    REAL(KIND=8) :: field(-1:,-1:)

    IF(field_type.EQ.CELL_DATA) THEN
        x_inc=0
        y_inc=0
    ENDIF
    IF(field_type.EQ.VERTEX_DATA) THEN
        x_inc=1
        y_inc=1
    ENDIF
    IF(field_type.EQ.X_FACE_DATA) THEN
        x_inc=1
        y_inc=0
    ENDIF
    IF(field_type.EQ.Y_FACE_DATA) THEN
        x_inc=0
        y_inc=1
    ENDIF

    size=(1+(chunks(chunk)%field%y_max+y_inc+topedge)-(chunks(chunk)%field%y_min-bottomedge))*depth

    ! pack buffer 
    CALL pack_right_buffer_seq(chunk, depth, x_inc, y_inc, right_snd_buffer, field)

    ! send buffer 
    CALL SHMEM_PUT64_NB(left_rcv_buffer, right_snd_buffer, size, receiver)

END SUBROUTINE clover_exchange_write_message_right

SUBROUTINE clover_exchange_write_message_bottom(chunk, depth, receiver, leftedge, rightedge, field_type, bottom_snd_buffer, top_rcv_buffer, field)

    IMPLICIT NONE 

    INTEGER :: chunk, depth, receiver, size, field_type, x_inc, y_inc, leftedge, rightedge
    REAL(KIND=8) :: bottom_snd_buffer(:), top_rcv_buffer(:) 
    REAL(KIND=8) :: field(-1:,-1:)

    IF(field_type.EQ.CELL_DATA) THEN
        x_inc=0
        y_inc=0
    ENDIF
    IF(field_type.EQ.VERTEX_DATA) THEN
        x_inc=1
        y_inc=1
    ENDIF
    IF(field_type.EQ.X_FACE_DATA) THEN
        x_inc=1
        y_inc=0
    ENDIF
    IF(field_type.EQ.Y_FACE_DATA) THEN
        x_inc=0
        y_inc=1
    ENDIF

    size=(1+(chunks(chunk)%field%x_max+x_inc+rightedge)-(chunks(chunk)%field%x_min-leftedge))*depth

    ! pack buffer 
    CALL pack_bottom_buffer_seq(chunk, depth, x_inc, y_inc, bottom_snd_buffer, field)

    ! send buffer 
    CALL SHMEM_PUT64_NB(top_rcv_buffer, bottom_snd_buffer, size, receiver)

END SUBROUTINE clover_exchange_write_message_bottom

SUBROUTINE clover_exchange_write_message_top(chunk, depth, receiver, leftedge, rightedge, field_type, top_snd_buffer, bottom_rcv_buffer, field)

    IMPLICIT NONE 

    INTEGER :: chunk, depth, receiver, size, field_type, x_inc, y_inc, leftedge, rightedge
    REAL(KIND=8) :: top_snd_buffer(:), bottom_rcv_buffer(:) 
    REAL(KIND=8) :: field(-1:,-1:)

    IF(field_type.EQ.CELL_DATA) THEN
        x_inc=0
        y_inc=0
    ENDIF
    IF(field_type.EQ.VERTEX_DATA) THEN
        x_inc=1
        y_inc=1
    ENDIF
    IF(field_type.EQ.X_FACE_DATA) THEN
        x_inc=1
        y_inc=0
    ENDIF
    IF(field_type.EQ.Y_FACE_DATA) THEN
        x_inc=0
        y_inc=1
    ENDIF

    size=(1+(chunks(chunk)%field%x_max+x_inc+rightedge)-(chunks(chunk)%field%x_min-leftedge))*depth

    ! pack buffer 
    CALL pack_top_buffer_seq(chunk, depth, x_inc, y_inc, top_snd_buffer, field)

    ! send buffer 
    CALL SHMEM_PUT64_NB(bottom_rcv_buffer, top_snd_buffer, size, receiver)

END SUBROUTINE clover_exchange_write_message_top


SUBROUTINE clover_exchange_write_message_left_top(chunk, depth, receiver, size, field_type, left_top_snd_buffer, right_bottom_rcv_buffer, field)

    IMPLICIT NONE 

    INTEGER :: chunk, depth, receiver, size, field_type, x_inc, y_inc
    REAL(KIND=8) :: left_top_snd_buffer(:), right_bottom_rcv_buffer(:) 
    REAL(KIND=8) :: field(-1:,-1:)

    IF(field_type.EQ.CELL_DATA) THEN
        x_inc=0
        y_inc=0
    ENDIF
    IF(field_type.EQ.VERTEX_DATA) THEN
        x_inc=1
        y_inc=1
    ENDIF
    IF(field_type.EQ.X_FACE_DATA) THEN
        x_inc=1
        y_inc=0
    ENDIF
    IF(field_type.EQ.Y_FACE_DATA) THEN
        x_inc=0
        y_inc=1
    ENDIF

    ! pack buffer 
    CALL pack_left_top_buffer_seq(chunk, depth, x_inc, y_inc, left_top_snd_buffer, field)

    ! send buffer 
    CALL SHMEM_PUT64_NB(right_bottom_rcv_buffer, left_top_snd_buffer, size, receiver)

END SUBROUTINE clover_exchange_write_message_left_top

SUBROUTINE clover_exchange_write_message_right_top(chunk, depth, receiver, size, field_type, right_top_snd_buffer, left_bottom_rcv_buffer, field)

    IMPLICIT NONE 

    INTEGER :: chunk, depth, receiver, size, field_type, x_inc, y_inc
    REAL(KIND=8) :: right_top_snd_buffer(:), left_bottom_rcv_buffer(:) 
    REAL(KIND=8) :: field(-1:,-1:)

    IF(field_type.EQ.CELL_DATA) THEN
        x_inc=0
        y_inc=0
    ENDIF
    IF(field_type.EQ.VERTEX_DATA) THEN
        x_inc=1
        y_inc=1
    ENDIF
    IF(field_type.EQ.X_FACE_DATA) THEN
        x_inc=1
        y_inc=0
    ENDIF
    IF(field_type.EQ.Y_FACE_DATA) THEN
        x_inc=0
        y_inc=1
    ENDIF

    ! pack buffer 
    CALL pack_right_top_buffer_seq(chunk, depth, x_inc, y_inc, right_top_snd_buffer, field)

    ! send buffer 
    CALL SHMEM_PUT64_NB(left_bottom_rcv_buffer, right_top_snd_buffer, size, receiver)

END SUBROUTINE clover_exchange_write_message_right_top

SUBROUTINE clover_exchange_write_message_right_bottom(chunk, depth, receiver, size, field_type, right_bottom_snd_buffer, left_top_rcv_buffer, field)

    IMPLICIT NONE 

    INTEGER :: chunk, depth, receiver, size, field_type, x_inc, y_inc
    REAL(KIND=8) :: right_bottom_snd_buffer(:), left_top_rcv_buffer(:) 
    REAL(KIND=8) :: field(-1:,-1:)

    IF(field_type.EQ.CELL_DATA) THEN
        x_inc=0
        y_inc=0
    ENDIF
    IF(field_type.EQ.VERTEX_DATA) THEN
        x_inc=1
        y_inc=1
    ENDIF
    IF(field_type.EQ.X_FACE_DATA) THEN
        x_inc=1
        y_inc=0
    ENDIF
    IF(field_type.EQ.Y_FACE_DATA) THEN
        x_inc=0
        y_inc=1
    ENDIF

    ! pack buffer 
    CALL pack_right_bottom_buffer_seq(chunk, depth, x_inc, y_inc, right_bottom_snd_buffer, field)

    ! send buffer 
    CALL SHMEM_PUT64_NB(left_top_rcv_buffer, right_bottom_snd_buffer, size, receiver)

END SUBROUTINE clover_exchange_write_message_right_bottom

SUBROUTINE clover_exchange_write_message_left_bottom(chunk, depth, receiver, size, field_type, left_bottom_snd_buffer, right_top_rcv_buffer, field)

    IMPLICIT NONE 

    INTEGER :: chunk, depth, receiver, size, field_type, x_inc, y_inc
    REAL(KIND=8) :: left_bottom_snd_buffer(:), right_top_rcv_buffer(:) 
    REAL(KIND=8) :: field(-1:,-1:)

    IF(field_type.EQ.CELL_DATA) THEN
        x_inc=0
        y_inc=0
    ENDIF
    IF(field_type.EQ.VERTEX_DATA) THEN
        x_inc=1
        y_inc=1
    ENDIF
    IF(field_type.EQ.X_FACE_DATA) THEN
        x_inc=1
        y_inc=0
    ENDIF
    IF(field_type.EQ.Y_FACE_DATA) THEN
        x_inc=0
        y_inc=1
    ENDIF

    ! pack buffer 
    CALL pack_left_bottom_buffer_seq(chunk, depth, x_inc, y_inc, left_bottom_snd_buffer, field)

    ! send buffer 
    CALL SHMEM_PUT64_NB(right_top_rcv_buffer, left_bottom_snd_buffer, size, receiver)

END SUBROUTINE clover_exchange_write_message_left_bottom





SUBROUTINE clover_exchange_unpack_all_buffers_left(chunk, depth, fields)

    IMPLICIT NONE 

    INTEGER :: chunk, depth, fields(NUM_FIELDS)

    IF(fields(FIELD_DENSITY0).EQ.1) THEN
        CALL unpack_left_buffer_seq(chunk, depth, CELL_DATA, density0_left_rcv_buffer, chunks(chunk)%field%density0)
    ENDIF
    IF(fields(FIELD_DENSITY1).EQ.1) THEN
        CALL unpack_left_buffer_seq(chunk, depth, CELL_DATA, density1_left_rcv_buffer, chunks(chunk)%field%density1)
    ENDIF
    IF(fields(FIELD_ENERGY0).EQ.1) THEN
        CALL unpack_left_buffer_seq(chunk, depth, CELL_DATA, energy0_left_rcv_buffer, chunks(chunk)%field%energy0)
    ENDIF
    IF(fields(FIELD_ENERGY1).EQ.1) THEN
        CALL unpack_left_buffer_seq(chunk, depth, CELL_DATA, energy1_left_rcv_buffer, chunks(chunk)%field%energy1)
    ENDIF
    IF(fields(FIELD_PRESSURE).EQ.1) THEN
        CALL unpack_left_buffer_seq(chunk, depth, CELL_DATA, pressure_left_rcv_buffer, chunks(chunk)%field%pressure)
    ENDIF
    IF(fields(FIELD_VISCOSITY).EQ.1) THEN
        CALL unpack_left_buffer_seq(chunk, depth, CELL_DATA, viscosity_left_rcv_buffer, chunks(chunk)%field%viscosity)
    ENDIF
    IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
        CALL unpack_left_buffer_seq(chunk, depth, CELL_DATA, soundspeed_left_rcv_buffer, chunks(chunk)%field%soundspeed)
    ENDIF
    IF(fields(FIELD_XVEL0).EQ.1) THEN
        CALL unpack_left_buffer_seq(chunk, depth, VERTEX_DATA, xvel0_left_rcv_buffer, chunks(chunk)%field%xvel0)
    ENDIF
    IF(fields(FIELD_XVEL1).EQ.1) THEN
        CALL unpack_left_buffer_seq(chunk, depth, VERTEX_DATA, xvel1_left_rcv_buffer, chunks(chunk)%field%xvel1)
    ENDIF
    IF(fields(FIELD_YVEL0).EQ.1) THEN
        CALL unpack_left_buffer_seq(chunk, depth, VERTEX_DATA, yvel0_left_rcv_buffer, chunks(chunk)%field%yvel0)
    ENDIF
    IF(fields(FIELD_YVEL1).EQ.1) THEN
        CALL unpack_left_buffer_seq(chunk, depth, VERTEX_DATA, yvel1_left_rcv_buffer, chunks(chunk)%field%yvel1)
    ENDIF
    IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
        CALL unpack_left_buffer_seq(chunk, depth, X_FACE_DATA, volflux_x_left_rcv_buffer, chunks(chunk)%field%vol_flux_x)
    ENDIF
    IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
        CALL unpack_left_buffer_seq(chunk, depth, Y_FACE_DATA, volflux_y_left_rcv_buffer, chunks(chunk)%field%vol_flux_y)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
        CALL unpack_left_buffer_seq(chunk, depth, X_FACE_DATA, massflux_x_left_rcv_buffer, chunks(chunk)%field%mass_flux_x)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
        CALL unpack_left_buffer_seq(chunk, depth, Y_FACE_DATA, massflux_y_left_rcv_buffer, chunks(chunk)%field%mass_flux_y)
    ENDIF

END SUBROUTINE clover_exchange_unpack_all_buffers_left

SUBROUTINE clover_exchange_unpack_all_buffers_right(chunk, depth, fields)

    IMPLICIT NONE 

    INTEGER :: chunk, depth, fields(NUM_FIELDS)

    IF(fields(FIELD_DENSITY0).EQ.1) THEN
        CALL unpack_right_buffer_seq(chunk, depth, CELL_DATA, density0_right_rcv_buffer, chunks(chunk)%field%density0)
    ENDIF
    IF(fields(FIELD_DENSITY1).EQ.1) THEN
        CALL unpack_right_buffer_seq(chunk, depth, CELL_DATA, density1_right_rcv_buffer, chunks(chunk)%field%density1)
    ENDIF
    IF(fields(FIELD_ENERGY0).EQ.1) THEN
        CALL unpack_right_buffer_seq(chunk, depth, CELL_DATA, energy0_right_rcv_buffer, chunks(chunk)%field%energy0)
    ENDIF
    IF(fields(FIELD_ENERGY1).EQ.1) THEN
        CALL unpack_right_buffer_seq(chunk, depth, CELL_DATA, energy1_right_rcv_buffer, chunks(chunk)%field%energy1)
    ENDIF
    IF(fields(FIELD_PRESSURE).EQ.1) THEN
        CALL unpack_right_buffer_seq(chunk, depth, CELL_DATA, pressure_right_rcv_buffer, chunks(chunk)%field%pressure)
    ENDIF
    IF(fields(FIELD_VISCOSITY).EQ.1) THEN
        CALL unpack_right_buffer_seq(chunk, depth, CELL_DATA, viscosity_right_rcv_buffer, chunks(chunk)%field%viscosity)
    ENDIF
    IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
        CALL unpack_right_buffer_seq(chunk, depth, CELL_DATA, soundspeed_right_rcv_buffer, chunks(chunk)%field%soundspeed)
    ENDIF
    IF(fields(FIELD_XVEL0).EQ.1) THEN
        CALL unpack_right_buffer_seq(chunk, depth, VERTEX_DATA, xvel0_right_rcv_buffer, chunks(chunk)%field%xvel0)
    ENDIF
    IF(fields(FIELD_XVEL1).EQ.1) THEN
        CALL unpack_right_buffer_seq(chunk, depth, VERTEX_DATA, xvel1_right_rcv_buffer, chunks(chunk)%field%xvel1)
    ENDIF
    IF(fields(FIELD_YVEL0).EQ.1) THEN
        CALL unpack_right_buffer_seq(chunk, depth, VERTEX_DATA, yvel0_right_rcv_buffer, chunks(chunk)%field%yvel0)
    ENDIF
    IF(fields(FIELD_YVEL1).EQ.1) THEN
        CALL unpack_right_buffer_seq(chunk, depth, VERTEX_DATA, yvel1_right_rcv_buffer, chunks(chunk)%field%yvel1)
    ENDIF
    IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
        CALL unpack_right_buffer_seq(chunk, depth, X_FACE_DATA, volflux_x_right_rcv_buffer, chunks(chunk)%field%vol_flux_x)
    ENDIF
    IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
        CALL unpack_right_buffer_seq(chunk, depth, Y_FACE_DATA, volflux_y_right_rcv_buffer, chunks(chunk)%field%vol_flux_y)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
        CALL unpack_right_buffer_seq(chunk, depth, X_FACE_DATA, massflux_x_right_rcv_buffer, chunks(chunk)%field%mass_flux_x)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
        CALL unpack_right_buffer_seq(chunk, depth, Y_FACE_DATA, massflux_y_right_rcv_buffer, chunks(chunk)%field%mass_flux_y)
    ENDIF

END SUBROUTINE clover_exchange_unpack_all_buffers_right


SUBROUTINE clover_exchange_unpack_all_buffers_bottom(chunk, depth, fields)

    IMPLICIT NONE 

    INTEGER :: chunk, depth, fields(NUM_FIELDS)

    IF(fields(FIELD_DENSITY0).EQ.1) THEN
        CALL unpack_bottom_buffer_seq(chunk, depth, CELL_DATA, density0_bottom_rcv_buffer, chunks(chunk)%field%density0)
    ENDIF
    IF(fields(FIELD_DENSITY1).EQ.1) THEN
        CALL unpack_bottom_buffer_seq(chunk, depth, CELL_DATA, density1_bottom_rcv_buffer, chunks(chunk)%field%density1)
    ENDIF
    IF(fields(FIELD_ENERGY0).EQ.1) THEN
        CALL unpack_bottom_buffer_seq(chunk, depth, CELL_DATA, energy0_bottom_rcv_buffer, chunks(chunk)%field%energy0)
    ENDIF
    IF(fields(FIELD_ENERGY1).EQ.1) THEN
        CALL unpack_bottom_buffer_seq(chunk, depth, CELL_DATA, energy1_bottom_rcv_buffer, chunks(chunk)%field%energy1)
    ENDIF
    IF(fields(FIELD_PRESSURE).EQ.1) THEN
        CALL unpack_bottom_buffer_seq(chunk, depth, CELL_DATA, pressure_bottom_rcv_buffer, chunks(chunk)%field%pressure)
    ENDIF
    IF(fields(FIELD_VISCOSITY).EQ.1) THEN
        CALL unpack_bottom_buffer_seq(chunk, depth, CELL_DATA, viscosity_bottom_rcv_buffer, chunks(chunk)%field%viscosity)
    ENDIF
    IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
        CALL unpack_bottom_buffer_seq(chunk, depth, CELL_DATA, soundspeed_bottom_rcv_buffer, chunks(chunk)%field%soundspeed)
    ENDIF
    IF(fields(FIELD_XVEL0).EQ.1) THEN
        CALL unpack_bottom_buffer_seq(chunk, depth, VERTEX_DATA, xvel0_bottom_rcv_buffer, chunks(chunk)%field%xvel0)
    ENDIF
    IF(fields(FIELD_XVEL1).EQ.1) THEN
        CALL unpack_bottom_buffer_seq(chunk, depth, VERTEX_DATA, xvel1_bottom_rcv_buffer, chunks(chunk)%field%xvel1)
    ENDIF
    IF(fields(FIELD_YVEL0).EQ.1) THEN
        CALL unpack_bottom_buffer_seq(chunk, depth, VERTEX_DATA, yvel0_bottom_rcv_buffer, chunks(chunk)%field%yvel0)
    ENDIF
    IF(fields(FIELD_YVEL1).EQ.1) THEN
        CALL unpack_bottom_buffer_seq(chunk, depth, VERTEX_DATA, yvel1_bottom_rcv_buffer, chunks(chunk)%field%yvel1)
    ENDIF
    IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
        CALL unpack_bottom_buffer_seq(chunk, depth, X_FACE_DATA, volflux_x_bottom_rcv_buffer, chunks(chunk)%field%vol_flux_x)
    ENDIF
    IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
        CALL unpack_bottom_buffer_seq(chunk, depth, Y_FACE_DATA, volflux_y_bottom_rcv_buffer, chunks(chunk)%field%vol_flux_y)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
        CALL unpack_bottom_buffer_seq(chunk, depth, X_FACE_DATA, massflux_x_bottom_rcv_buffer, chunks(chunk)%field%mass_flux_x)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
        CALL unpack_bottom_buffer_seq(chunk, depth, Y_FACE_DATA, massflux_y_bottom_rcv_buffer, chunks(chunk)%field%mass_flux_y)
    ENDIF

END SUBROUTINE clover_exchange_unpack_all_buffers_bottom

SUBROUTINE clover_exchange_unpack_all_buffers_top(chunk, depth, fields)

    IMPLICIT NONE 

    INTEGER :: chunk, depth, fields(NUM_FIELDS)

    IF(fields(FIELD_DENSITY0).EQ.1) THEN
        CALL unpack_top_buffer_seq(chunk, depth, CELL_DATA, density0_top_rcv_buffer, chunks(chunk)%field%density0)
    ENDIF
    IF(fields(FIELD_DENSITY1).EQ.1) THEN
        CALL unpack_top_buffer_seq(chunk, depth, CELL_DATA, density1_top_rcv_buffer, chunks(chunk)%field%density1)
    ENDIF
    IF(fields(FIELD_ENERGY0).EQ.1) THEN
        CALL unpack_top_buffer_seq(chunk, depth, CELL_DATA, energy0_top_rcv_buffer, chunks(chunk)%field%energy0)
    ENDIF
    IF(fields(FIELD_ENERGY1).EQ.1) THEN
        CALL unpack_top_buffer_seq(chunk, depth, CELL_DATA, energy1_top_rcv_buffer, chunks(chunk)%field%energy1)
    ENDIF
    IF(fields(FIELD_PRESSURE).EQ.1) THEN
        CALL unpack_top_buffer_seq(chunk, depth, CELL_DATA, pressure_top_rcv_buffer, chunks(chunk)%field%pressure)
    ENDIF
    IF(fields(FIELD_VISCOSITY).EQ.1) THEN
        CALL unpack_top_buffer_seq(chunk, depth, CELL_DATA, viscosity_top_rcv_buffer, chunks(chunk)%field%viscosity)
    ENDIF
    IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
        CALL unpack_top_buffer_seq(chunk, depth, CELL_DATA, soundspeed_top_rcv_buffer, chunks(chunk)%field%soundspeed)
    ENDIF
    IF(fields(FIELD_XVEL0).EQ.1) THEN
        CALL unpack_top_buffer_seq(chunk, depth, VERTEX_DATA, xvel0_top_rcv_buffer, chunks(chunk)%field%xvel0)
    ENDIF
    IF(fields(FIELD_XVEL1).EQ.1) THEN
        CALL unpack_top_buffer_seq(chunk, depth, VERTEX_DATA, xvel1_top_rcv_buffer, chunks(chunk)%field%xvel1)
    ENDIF
    IF(fields(FIELD_YVEL0).EQ.1) THEN
        CALL unpack_top_buffer_seq(chunk, depth, VERTEX_DATA, yvel0_top_rcv_buffer, chunks(chunk)%field%yvel0)
    ENDIF
    IF(fields(FIELD_YVEL1).EQ.1) THEN
        CALL unpack_top_buffer_seq(chunk, depth, VERTEX_DATA, yvel1_top_rcv_buffer, chunks(chunk)%field%yvel1)
    ENDIF
    IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
        CALL unpack_top_buffer_seq(chunk, depth, X_FACE_DATA, volflux_x_top_rcv_buffer, chunks(chunk)%field%vol_flux_x)
    ENDIF
    IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
        CALL unpack_top_buffer_seq(chunk, depth, Y_FACE_DATA, volflux_y_top_rcv_buffer, chunks(chunk)%field%vol_flux_y)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
        CALL unpack_top_buffer_seq(chunk, depth, X_FACE_DATA, massflux_x_top_rcv_buffer, chunks(chunk)%field%mass_flux_x)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
        CALL unpack_top_buffer_seq(chunk, depth, Y_FACE_DATA, massflux_y_top_rcv_buffer, chunks(chunk)%field%mass_flux_y)
    ENDIF

END SUBROUTINE clover_exchange_unpack_all_buffers_top

SUBROUTINE clover_exchange_unpack_all_buffers_left_top(chunk, depth, fields)

    IMPLICIT NONE 

    INTEGER :: chunk, depth, fields(NUM_FIELDS)

    IF(fields(FIELD_DENSITY0).EQ.1) THEN
        CALL unpack_left_top_buffer_seq(chunk, depth, CELL_DATA, density0_left_top_rcv_buffer, chunks(chunk)%field%density0)
    ENDIF
    IF(fields(FIELD_DENSITY1).EQ.1) THEN
        CALL unpack_left_top_buffer_seq(chunk, depth, CELL_DATA, density1_left_top_rcv_buffer, chunks(chunk)%field%density1)
    ENDIF
    IF(fields(FIELD_ENERGY0).EQ.1) THEN
        CALL unpack_left_top_buffer_seq(chunk, depth, CELL_DATA, energy0_left_top_rcv_buffer, chunks(chunk)%field%energy0)
    ENDIF
    IF(fields(FIELD_ENERGY1).EQ.1) THEN
        CALL unpack_left_top_buffer_seq(chunk, depth, CELL_DATA, energy1_left_top_rcv_buffer, chunks(chunk)%field%energy1)
    ENDIF
    IF(fields(FIELD_PRESSURE).EQ.1) THEN
        CALL unpack_left_top_buffer_seq(chunk, depth, CELL_DATA, pressure_left_top_rcv_buffer, chunks(chunk)%field%pressure)
    ENDIF
    IF(fields(FIELD_VISCOSITY).EQ.1) THEN
        CALL unpack_left_top_buffer_seq(chunk, depth, CELL_DATA, viscosity_left_top_rcv_buffer, chunks(chunk)%field%viscosity)
    ENDIF
    IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
        CALL unpack_left_top_buffer_seq(chunk, depth, CELL_DATA, soundspeed_left_top_rcv_buffer, chunks(chunk)%field%soundspeed)
    ENDIF
    IF(fields(FIELD_XVEL0).EQ.1) THEN
        CALL unpack_left_top_buffer_seq(chunk, depth, VERTEX_DATA, xvel0_left_top_rcv_buffer, chunks(chunk)%field%xvel0)
    ENDIF
    IF(fields(FIELD_XVEL1).EQ.1) THEN
        CALL unpack_left_top_buffer_seq(chunk, depth, VERTEX_DATA, xvel1_left_top_rcv_buffer, chunks(chunk)%field%xvel1)
    ENDIF
    IF(fields(FIELD_YVEL0).EQ.1) THEN
        CALL unpack_left_top_buffer_seq(chunk, depth, VERTEX_DATA, yvel0_left_top_rcv_buffer, chunks(chunk)%field%yvel0)
    ENDIF
    IF(fields(FIELD_YVEL1).EQ.1) THEN
        CALL unpack_left_top_buffer_seq(chunk, depth, VERTEX_DATA, yvel1_left_top_rcv_buffer, chunks(chunk)%field%yvel1)
    ENDIF
    IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
        CALL unpack_left_top_buffer_seq(chunk, depth, X_FACE_DATA, volflux_x_left_top_rcv_buffer, chunks(chunk)%field%vol_flux_x)
    ENDIF
    IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
        CALL unpack_left_top_buffer_seq(chunk, depth, Y_FACE_DATA, volflux_y_left_top_rcv_buffer, chunks(chunk)%field%vol_flux_y)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
        CALL unpack_left_top_buffer_seq(chunk, depth, X_FACE_DATA, massflux_x_left_top_rcv_buffer, chunks(chunk)%field%mass_flux_x)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
        CALL unpack_left_top_buffer_seq(chunk, depth, Y_FACE_DATA, massflux_y_left_top_rcv_buffer, chunks(chunk)%field%mass_flux_y)
    ENDIF

END SUBROUTINE clover_exchange_unpack_all_buffers_left_top

SUBROUTINE clover_exchange_unpack_all_buffers_right_top(chunk, depth, fields)

    IMPLICIT NONE 

    INTEGER :: chunk, depth, fields(NUM_FIELDS)

    IF(fields(FIELD_DENSITY0).EQ.1) THEN
        CALL unpack_right_top_buffer_seq(chunk, depth, CELL_DATA, density0_right_top_rcv_buffer, chunks(chunk)%field%density0)
    ENDIF
    IF(fields(FIELD_DENSITY1).EQ.1) THEN
        CALL unpack_right_top_buffer_seq(chunk, depth, CELL_DATA, density1_right_top_rcv_buffer, chunks(chunk)%field%density1)
    ENDIF
    IF(fields(FIELD_ENERGY0).EQ.1) THEN
        CALL unpack_right_top_buffer_seq(chunk, depth, CELL_DATA, energy0_right_top_rcv_buffer, chunks(chunk)%field%energy0)
    ENDIF
    IF(fields(FIELD_ENERGY1).EQ.1) THEN
        CALL unpack_right_top_buffer_seq(chunk, depth, CELL_DATA, energy1_right_top_rcv_buffer, chunks(chunk)%field%energy1)
    ENDIF
    IF(fields(FIELD_PRESSURE).EQ.1) THEN
        CALL unpack_right_top_buffer_seq(chunk, depth, CELL_DATA, pressure_right_top_rcv_buffer, chunks(chunk)%field%pressure)
    ENDIF
    IF(fields(FIELD_VISCOSITY).EQ.1) THEN
        CALL unpack_right_top_buffer_seq(chunk, depth, CELL_DATA, viscosity_right_top_rcv_buffer, chunks(chunk)%field%viscosity)
    ENDIF
    IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
        CALL unpack_right_top_buffer_seq(chunk, depth, CELL_DATA, soundspeed_right_top_rcv_buffer, chunks(chunk)%field%soundspeed)
    ENDIF
    IF(fields(FIELD_XVEL0).EQ.1) THEN
        CALL unpack_right_top_buffer_seq(chunk, depth, VERTEX_DATA, xvel0_right_top_rcv_buffer, chunks(chunk)%field%xvel0)
    ENDIF
    IF(fields(FIELD_XVEL1).EQ.1) THEN
        CALL unpack_right_top_buffer_seq(chunk, depth, VERTEX_DATA, xvel1_right_top_rcv_buffer, chunks(chunk)%field%xvel1)
    ENDIF
    IF(fields(FIELD_YVEL0).EQ.1) THEN
        CALL unpack_right_top_buffer_seq(chunk, depth, VERTEX_DATA, yvel0_right_top_rcv_buffer, chunks(chunk)%field%yvel0)
    ENDIF
    IF(fields(FIELD_YVEL1).EQ.1) THEN
        CALL unpack_right_top_buffer_seq(chunk, depth, VERTEX_DATA, yvel1_right_top_rcv_buffer, chunks(chunk)%field%yvel1)
    ENDIF
    IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
        CALL unpack_right_top_buffer_seq(chunk, depth, X_FACE_DATA, volflux_x_right_top_rcv_buffer, chunks(chunk)%field%vol_flux_x)
    ENDIF
    IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
        CALL unpack_right_top_buffer_seq(chunk, depth, Y_FACE_DATA, volflux_y_right_top_rcv_buffer, chunks(chunk)%field%vol_flux_y)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
        CALL unpack_right_top_buffer_seq(chunk, depth, X_FACE_DATA, massflux_x_right_top_rcv_buffer, chunks(chunk)%field%mass_flux_x)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
        CALL unpack_right_top_buffer_seq(chunk, depth, Y_FACE_DATA, massflux_y_right_top_rcv_buffer, chunks(chunk)%field%mass_flux_y)
    ENDIF

END SUBROUTINE clover_exchange_unpack_all_buffers_right_top

SUBROUTINE clover_exchange_unpack_all_buffers_right_bottom(chunk, depth, fields)

    IMPLICIT NONE 

    INTEGER :: chunk, depth, fields(NUM_FIELDS)

    IF(fields(FIELD_DENSITY0).EQ.1) THEN
        CALL unpack_right_bottom_buffer_seq(chunk, depth, CELL_DATA, density0_right_bottom_rcv_buffer, chunks(chunk)%field%density0)
    ENDIF
    IF(fields(FIELD_DENSITY1).EQ.1) THEN
        CALL unpack_right_bottom_buffer_seq(chunk, depth, CELL_DATA, density1_right_bottom_rcv_buffer, chunks(chunk)%field%density1)
    ENDIF
    IF(fields(FIELD_ENERGY0).EQ.1) THEN
        CALL unpack_right_bottom_buffer_seq(chunk, depth, CELL_DATA, energy0_right_bottom_rcv_buffer, chunks(chunk)%field%energy0)
    ENDIF
    IF(fields(FIELD_ENERGY1).EQ.1) THEN
        CALL unpack_right_bottom_buffer_seq(chunk, depth, CELL_DATA, energy1_right_bottom_rcv_buffer, chunks(chunk)%field%energy1)
    ENDIF
    IF(fields(FIELD_PRESSURE).EQ.1) THEN
        CALL unpack_right_bottom_buffer_seq(chunk, depth, CELL_DATA, pressure_right_bottom_rcv_buffer, chunks(chunk)%field%pressure)
    ENDIF
    IF(fields(FIELD_VISCOSITY).EQ.1) THEN
        CALL unpack_right_bottom_buffer_seq(chunk, depth, CELL_DATA, viscosity_right_bottom_rcv_buffer, chunks(chunk)%field%viscosity)
    ENDIF
    IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
        CALL unpack_right_bottom_buffer_seq(chunk, depth, CELL_DATA, soundspeed_right_bottom_rcv_buffer, chunks(chunk)%field%soundspeed)
    ENDIF
    IF(fields(FIELD_XVEL0).EQ.1) THEN
        CALL unpack_right_bottom_buffer_seq(chunk, depth, VERTEX_DATA, xvel0_right_bottom_rcv_buffer, chunks(chunk)%field%xvel0)
    ENDIF
    IF(fields(FIELD_XVEL1).EQ.1) THEN
        CALL unpack_right_bottom_buffer_seq(chunk, depth, VERTEX_DATA, xvel1_right_bottom_rcv_buffer, chunks(chunk)%field%xvel1)
    ENDIF
    IF(fields(FIELD_YVEL0).EQ.1) THEN
        CALL unpack_right_bottom_buffer_seq(chunk, depth, VERTEX_DATA, yvel0_right_bottom_rcv_buffer, chunks(chunk)%field%yvel0)
    ENDIF
    IF(fields(FIELD_YVEL1).EQ.1) THEN
        CALL unpack_right_bottom_buffer_seq(chunk, depth, VERTEX_DATA, yvel1_right_bottom_rcv_buffer, chunks(chunk)%field%yvel1)
    ENDIF
    IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
        CALL unpack_right_bottom_buffer_seq(chunk, depth, X_FACE_DATA, volflux_x_right_bottom_rcv_buffer, chunks(chunk)%field%vol_flux_x)
    ENDIF
    IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
        CALL unpack_right_bottom_buffer_seq(chunk, depth, Y_FACE_DATA, volflux_y_right_bottom_rcv_buffer, chunks(chunk)%field%vol_flux_y)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
        CALL unpack_right_bottom_buffer_seq(chunk, depth, X_FACE_DATA, massflux_x_right_bottom_rcv_buffer, chunks(chunk)%field%mass_flux_x)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
        CALL unpack_right_bottom_buffer_seq(chunk, depth, Y_FACE_DATA, massflux_y_right_bottom_rcv_buffer, chunks(chunk)%field%mass_flux_y)
    ENDIF

END SUBROUTINE clover_exchange_unpack_all_buffers_right_bottom

SUBROUTINE clover_exchange_unpack_all_buffers_left_bottom(chunk, depth, fields)

    IMPLICIT NONE 

    INTEGER :: chunk, depth, fields(NUM_FIELDS)

    IF(fields(FIELD_DENSITY0).EQ.1) THEN
        CALL unpack_left_bottom_buffer_seq(chunk, depth, CELL_DATA, density0_left_bottom_rcv_buffer, chunks(chunk)%field%density0)
    ENDIF
    IF(fields(FIELD_DENSITY1).EQ.1) THEN
        CALL unpack_left_bottom_buffer_seq(chunk, depth, CELL_DATA, density1_left_bottom_rcv_buffer, chunks(chunk)%field%density1)
    ENDIF
    IF(fields(FIELD_ENERGY0).EQ.1) THEN
        CALL unpack_left_bottom_buffer_seq(chunk, depth, CELL_DATA, energy0_left_bottom_rcv_buffer, chunks(chunk)%field%energy0)
    ENDIF
    IF(fields(FIELD_ENERGY1).EQ.1) THEN
        CALL unpack_left_bottom_buffer_seq(chunk, depth, CELL_DATA, energy1_left_bottom_rcv_buffer, chunks(chunk)%field%energy1)
    ENDIF
    IF(fields(FIELD_PRESSURE).EQ.1) THEN
        CALL unpack_left_bottom_buffer_seq(chunk, depth, CELL_DATA, pressure_left_bottom_rcv_buffer, chunks(chunk)%field%pressure)
    ENDIF
    IF(fields(FIELD_VISCOSITY).EQ.1) THEN
        CALL unpack_left_bottom_buffer_seq(chunk, depth, CELL_DATA, viscosity_left_bottom_rcv_buffer, chunks(chunk)%field%viscosity)
    ENDIF
    IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
        CALL unpack_left_bottom_buffer_seq(chunk, depth, CELL_DATA, soundspeed_left_bottom_rcv_buffer, chunks(chunk)%field%soundspeed)
    ENDIF
    IF(fields(FIELD_XVEL0).EQ.1) THEN
        CALL unpack_left_bottom_buffer_seq(chunk, depth, VERTEX_DATA, xvel0_left_bottom_rcv_buffer, chunks(chunk)%field%xvel0)
    ENDIF
    IF(fields(FIELD_XVEL1).EQ.1) THEN
        CALL unpack_left_bottom_buffer_seq(chunk, depth, VERTEX_DATA, xvel1_left_bottom_rcv_buffer, chunks(chunk)%field%xvel1)
    ENDIF
    IF(fields(FIELD_YVEL0).EQ.1) THEN
        CALL unpack_left_bottom_buffer_seq(chunk, depth, VERTEX_DATA, yvel0_left_bottom_rcv_buffer, chunks(chunk)%field%yvel0)
    ENDIF
    IF(fields(FIELD_YVEL1).EQ.1) THEN
        CALL unpack_left_bottom_buffer_seq(chunk, depth, VERTEX_DATA, yvel1_left_bottom_rcv_buffer, chunks(chunk)%field%yvel1)
    ENDIF
    IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
        CALL unpack_left_bottom_buffer_seq(chunk, depth, X_FACE_DATA, volflux_x_left_bottom_rcv_buffer, chunks(chunk)%field%vol_flux_x)
    ENDIF
    IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
        CALL unpack_left_bottom_buffer_seq(chunk, depth, Y_FACE_DATA, volflux_y_left_bottom_rcv_buffer, chunks(chunk)%field%vol_flux_y)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
        CALL unpack_left_bottom_buffer_seq(chunk, depth, X_FACE_DATA, massflux_x_left_bottom_rcv_buffer, chunks(chunk)%field%mass_flux_x)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
        CALL unpack_left_bottom_buffer_seq(chunk, depth, Y_FACE_DATA, massflux_y_left_bottom_rcv_buffer, chunks(chunk)%field%mass_flux_y)
    ENDIF

END SUBROUTINE clover_exchange_unpack_all_buffers_left_bottom




!SUBROUTINE clover_exchange_message(chunk,field,                            &
!                                   left_snd_buffer,                        &
!                                   left_rcv_buffer,                        &
!                                   right_snd_buffer,                       &
!                                   right_rcv_buffer,                       &
!                                   bottom_snd_buffer,                      &
!                                   bottom_rcv_buffer,                      &
!                                   top_snd_buffer,                         &
!                                   top_rcv_buffer,                         &
!                                   left_top_snd_buffer,                                         &
!                                   left_top_rcv_buffer,                                         &
!                                   right_top_snd_buffer,                                        &
!                                   right_top_rcv_buffer,                                        &
!                                   right_bottom_snd_buffer,                                     &            
!                                   right_bottom_rcv_buffer,                                     &
!                                   left_bottom_snd_buffer,                                      &
!                                   left_bottom_rcv_buffer,                                      &
!                                   depth,field_type)
!
!  !USE pack_kernel_module
!
!  IMPLICIT NONE
!
!  REAL(KIND=8) :: field(-1:,-1:) ! This seems to work for any type of mesh data
!  REAL(KIND=8) :: left_snd_buffer(:),left_rcv_buffer(:),right_snd_buffer(:),right_rcv_buffer(:)
!  REAL(KIND=8) :: bottom_snd_buffer(:),bottom_rcv_buffer(:),top_snd_buffer(:),top_rcv_buffer(:)
!
!  REAL(KIND=8) :: left_top_snd_buffer(:),left_top_rcv_buffer(:),right_top_snd_buffer(:),right_top_rcv_buffer(:)
!  REAL(KIND=8) :: right_bottom_snd_buffer(:),right_bottom_rcv_buffer(:),left_bottom_snd_buffer(:),left_bottom_rcv_buffer(:)
!
!  INTEGER      :: chunk,depth,field_type
!
!  INTEGER      :: size,x_inc,y_inc
!  INTEGER      :: receiver,sender
!  INTEGER      :: topedge, bottomedge, leftedge, rightedge
!
!  ! Field type will either be cell, vertex, x_face or y_face to get the message limits correct
!
!  ! I am packing my own buffers. I am sure this could be improved with MPI data types
!  !  but this will do for now
!
!  ! I am also sending buffers to chunks with the same task id for now.
!  ! This can be improved in the future but at the moment there is just 1 chunk per task anyway
!
!  ! The tag will be a function of the sending chunk and the face it is coming from
!  !  like chunk 6 sending the left face
!
!  ! No open mp in here either. May be beneficial will packing and unpacking in the future, though I am not sure.
!
!  ! Change this so it will allow more than 1 chunk per task
!
!  ! Pack and send
!
!  ! These array modifications still need to be added on, plus the donor data location changes as in update_halo
!  IF(field_type.EQ.CELL_DATA) THEN
!    x_inc=0
!    y_inc=0
!  ENDIF
!  IF(field_type.EQ.VERTEX_DATA) THEN
!    x_inc=1
!    y_inc=1
!  ENDIF
!  IF(field_type.EQ.X_FACE_DATA) THEN
!    x_inc=1
!    y_inc=0
!  ENDIF
!  IF(field_type.EQ.Y_FACE_DATA) THEN
!    x_inc=0
!    y_inc=1
!  ENDIF
!
!  ! Pack real data into buffers
!  IF(parallel%task.EQ.chunks(chunk)%task) THEN
!
!    !pack all the data buffers required 
!    IF(chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) THEN
!        CALL pack_left_buffer_seq(chunk, depth, x_inc, y_inc, left_snd_buffer, field)
!    ENDIF
!    IF(chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) THEN
!        CALL pack_right_buffer_seq(chunk, depth, x_inc, y_inc, right_snd_buffer, field)
!    ENDIF
!    IF(chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) THEN
!        CALL pack_bottom_buffer_seq(chunk, depth, x_inc, y_inc, bottom_snd_buffer, field)
!    ENDIF
!    IF(chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) THEN
!        CALL pack_top_buffer_seq(chunk, depth, x_inc, y_inc, top_snd_buffer, field)
!    ENDIF
!    IF ( (chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) ) THEN
!        CALL pack_left_top_buffer_seq(chunk, depth, x_inc, y_inc, left_top_snd_buffer, field)
!    ENDIF
!    IF ( (chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) ) THEN
!        CALL pack_right_top_buffer_seq(chunk, depth, x_inc, y_inc, right_top_snd_buffer, field)
!    ENDIF
!    IF ( (chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) ) THEN
!        CALL pack_right_bottom_buffer_seq(chunk, depth, x_inc, y_inc, right_bottom_snd_buffer, field)
!    ENDIF
!    IF ( (chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) ) THEN
!        CALL pack_left_bottom_buffer_seq(chunk, depth, x_inc, y_inc, left_bottom_snd_buffer, field)
!    ENDIF
!
!
!    topedge = 0
!    bottomedge = 0
!    IF (chunks(chunk)%chunk_neighbours(chunk_top).EQ.external_face) THEN
!      topedge = depth
!    ENDIF
!    IF (chunks(chunk)%chunk_neighbours(chunk_bottom).EQ.external_face) THEN
!      bottomedge = depth
!    ENDIF
!    size=(1+(chunks(chunk)%field%y_max+y_inc+topedge)-(chunks(chunk)%field%y_min-bottomedge))*depth
!
!
!    ! Send/receive the data
!    IF(chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) THEN
!      receiver=chunks(chunks(chunk)%chunk_neighbours(chunk_left))%task
!
!      IF (left_write_flag .EQ. 0) THEN
!        CALL SHMEM_INT4_WAIT_UNTIL(left_write_flag, SHMEM_CMP_EQ, 1)
!      ENDIF
!      CALL SHMEM_PUT64_NB(right_rcv_buffer, left_snd_buffer, size, receiver)
!      
!      left_write_flag = 0 
!
!    ENDIF
!
!    IF(chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) THEN
!      receiver=chunks(chunks(chunk)%chunk_neighbours(chunk_right))%task
!
!      IF (right_write_flag .EQ. 0) THEN
!        CALL SHMEM_INT4_WAIT_UNTIL(right_write_flag, SHMEM_CMP_EQ, 1)
!      ENDIF
!      CALL SHMEM_PUT64_NB(left_rcv_buffer, right_snd_buffer, size, receiver)
!
!      right_write_flag = 0 
!
!    ENDIF
!
!    leftedge= 0
!    rightedge= 0
!    IF (chunks(chunk)%chunk_neighbours(chunk_left).EQ.external_face) THEN
!      leftedge = depth
!    ENDIF
!    IF (chunks(chunk)%chunk_neighbours(chunk_right).EQ.external_face) THEN
!      rightedge = depth
!    ENDIF
!    size=(1+(chunks(chunk)%field%x_max+x_inc+rightedge)-(chunks(chunk)%field%x_min-leftedge))*depth
!
!    ! Send/receive the data
!    IF(chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) THEN
!      receiver=chunks(chunks(chunk)%chunk_neighbours(chunk_bottom))%task
!
!      IF (bottom_write_flag .EQ. 0) THEN
!        CALL SHMEM_INT4_WAIT_UNTIL(bottom_write_flag, SHMEM_CMP_EQ, 1)
!      ENDIF
!
!      CALL SHMEM_PUT64_NB(top_rcv_buffer, bottom_snd_buffer, size, receiver)
!      bottom_write_flag = 0
!    ENDIF
!
!    IF(chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) THEN
!        receiver=chunks(chunks(chunk)%chunk_neighbours(chunk_top))%task
!
!        IF (top_write_flag .EQ. 0) THEN
!            CALL SHMEM_INT4_WAIT_UNTIL(top_write_flag, SHMEM_CMP_EQ, 1)
!        ENDIF
!
!        CALL SHMEM_PUT64_NB(bottom_rcv_buffer, top_snd_buffer, size, receiver)
!        top_write_flag = 0 
!    ENDIF
!
!
!    size = depth*depth
!
!    IF ( (chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) ) THEN
!        receiver = chunks(chunks(chunk)%chunk_neighbours(chunk_left_top))%task
!
!        IF (left_top_write_flag .EQ. 0) THEN
!            CALL SHMEM_INT4_WAIT_UNTIL(left_top_write_flag, SHMEM_CMP_EQ, 1)
!        ENDIF
!
!        CALL SHMEM_PUT64_NB(right_bottom_rcv_buffer, left_top_snd_buffer, size, receiver)
!        left_top_write_flag = 0
!    ENDIF
!    
!    IF ( (chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) ) THEN
!        receiver = chunks(chunks(chunk)%chunk_neighbours(chunk_right_top))%task
!
!        IF (right_top_write_flag .EQ. 0) THEN
!            CALL SHMEM_INT4_WAIT_UNTIL(right_top_write_flag, SHMEM_CMP_EQ, 1)
!        ENDIF
!
!        CALL SHMEM_PUT64_NB(left_bottom_rcv_buffer, right_top_snd_buffer, size, receiver)
!        right_top_write_flag = 0
!    ENDIF
!    
!    IF ( (chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) ) THEN
!        receiver = chunks(chunks(chunk)%chunk_neighbours(chunk_right_bottom))%task
!
!        IF (right_bottom_write_flag .EQ. 0) THEN
!            CALL SHMEM_INT4_WAIT_UNTIL(right_bottom_write_flag, SHMEM_CMP_EQ, 1)
!        ENDIF
!
!        CALL SHMEM_PUT64_NB(left_top_rcv_buffer, right_bottom_snd_buffer, size, receiver)
!        right_bottom_write_flag = 0
!    ENDIF
!    
!    IF ( (chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) ) THEN
!       receiver = chunks(chunks(chunk)%chunk_neighbours(chunk_left_bottom))%task
!
!       IF (left_bottom_write_flag .EQ. 0) THEN
!           CALL SHMEM_INT4_WAIT_UNTIL(left_bottom_write_flag, SHMEM_CMP_EQ, 1)
!       ENDIF
!
!       CALL SHMEM_PUT64_NB(right_top_rcv_buffer, left_bottom_snd_buffer, size, receiver)
!       left_bottom_write_flag = 0
!    ENDIF
!
!  ENDIF
!
!  ! Wait for the messages
!#ifdef FENCE_NOT_QUIET
!        CALL SHMEM_FENCE
!#else
!        CALL SHMEM_QUIET
!#endif
!
!
!
!
!  IF(parallel%task.EQ.chunks(chunk)%task) THEN
!
!    ! Send/receive the data
!    IF(chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) THEN
!      receiver=chunks(chunks(chunk)%chunk_neighbours(chunk_left))%task
!      CALL SHMEM_PUT4_NB(right_rcv_flag, 1, 1, receiver)
!    ENDIF
!
!    IF(chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) THEN
!      receiver=chunks(chunks(chunk)%chunk_neighbours(chunk_right))%task
!      CALL SHMEM_PUT4_NB(left_rcv_flag, 1, 1, receiver)
!    ENDIF
!
!    IF(chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) THEN
!      receiver=chunks(chunks(chunk)%chunk_neighbours(chunk_bottom))%task
!      CALL SHMEM_PUT4_NB(top_rcv_flag, 1, 1, receiver)
!    ENDIF
!
!    IF(chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) THEN
!      receiver=chunks(chunks(chunk)%chunk_neighbours(chunk_top))%task
!      CALL SHMEM_PUT4_NB(bottom_rcv_flag, 1, 1, receiver)
!    ENDIF
!    IF ( (chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) ) THEN
!        receiver = chunks(chunks(chunk)%chunk_neighbours(chunk_left_top))%task
!        CALL SHMEM_PUT4_NB(right_bottom_rcv_flag, 1, 1, receiver)
!    ENDIF
!    IF ( (chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) ) THEN
!        receiver = chunks(chunks(chunk)%chunk_neighbours(chunk_right_top))%task
!        CALL SHMEM_PUT4_NB(left_bottom_rcv_flag, 1, 1, receiver)
!    ENDIF
!    IF ( (chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) ) THEN
!        receiver = chunks(chunks(chunk)%chunk_neighbours(chunk_right_bottom))%task
!        CALL SHMEM_PUT4_NB(left_top_rcv_flag, 1, 1, receiver)
!    ENDIF
!    IF ( (chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) ) THEN
!        receiver = chunks(chunks(chunk)%chunk_neighbours(chunk_left_bottom))%task
!        CALL SHMEM_PUT4_NB(right_top_rcv_flag, 1, 1, receiver)
!    ENDIF
!  ENDIF
!
!  IF(parallel%task.EQ.chunks(chunk)%task) THEN
!
!    ! Send/receive the data
!    IF(chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) THEN
!      IF (left_rcv_flag .EQ. 0) THEN
!        CALL SHMEM_INT4_WAIT_UNTIL(left_rcv_flag, SHMEM_CMP_EQ, 1)
!      ENDIF
!      left_rcv_flag = 0
!    ENDIF
!
!    IF(chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) THEN
!      IF (right_rcv_flag .EQ. 0) THEN
!        CALL SHMEM_INT4_WAIT_UNTIL(right_rcv_flag, SHMEM_CMP_EQ, 1)
!      ENDIF
!      right_rcv_flag = 0
!    ENDIF
!
!    IF(chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) THEN
!      IF (bottom_rcv_flag .EQ. 0) THEN
!        CALL SHMEM_INT4_WAIT_UNTIL(bottom_rcv_flag, SHMEM_CMP_EQ, 1)
!      ENDIF
!      bottom_rcv_flag = 0
!    ENDIF
!
!    IF(chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) THEN
!      IF (top_rcv_flag .EQ. 0) THEN
!        CALL SHMEM_INT4_WAIT_UNTIL(top_rcv_flag, SHMEM_CMP_EQ, 1)
!      ENDIF
!      top_rcv_flag = 0
!    ENDIF
!    IF ( (chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) ) THEN
!        IF (left_top_rcv_flag .EQ. 0) THEN
!            CALL SHMEM_INT4_WAIT_UNTIL(left_top_rcv_flag, SHMEM_CMP_EQ, 1)
!        ENDIF
!        left_top_rcv_flag = 0
!    ENDIF
!    IF ( (chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) ) THEN
!        IF (right_top_rcv_flag .EQ. 0) THEN
!            CALL SHMEM_INT4_WAIT_UNTIL(right_top_rcv_flag, SHMEM_CMP_EQ, 1)
!        ENDIF
!        right_top_rcv_flag = 0
!    ENDIF
!    IF ( (chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) ) THEN
!        IF (right_bottom_rcv_flag .EQ. 0) THEN
!            CALL SHMEM_INT4_WAIT_UNTIL(right_bottom_rcv_flag, SHMEM_CMP_EQ, 1)
!        ENDIF
!        right_bottom_rcv_flag = 0
!    ENDIF
!    IF ( (chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) ) THEN
!        IF (left_bottom_rcv_flag .EQ. 0) THEN
!            CALL SHMEM_INT4_WAIT_UNTIL(left_bottom_rcv_flag, SHMEM_CMP_EQ, 1)
!        ENDIF
!        left_bottom_rcv_flag = 0
!    ENDIF
!  ENDIF
!
!
!
!
!
!  ! Unpack buffers in halo cells
!  IF(parallel%task.EQ.chunks(chunk)%task) THEN
!
!    IF(chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) THEN
!        CALL unpack_left_buffer_seq(chunk, depth, x_inc, y_inc, left_rcv_buffer, field)
!    ENDIF
!    IF(chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) THEN
!        CALL unpack_right_buffer_seq(chunk, depth, x_inc, y_inc, right_rcv_buffer, field)
!    ENDIF
!    IF(chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) THEN
!        CALL unpack_bottom_buffer_seq(chunk, depth, x_inc, y_inc, bottom_rcv_buffer, field)
!    ENDIF
!    IF(chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) THEN
!        CALL unpack_top_buffer_seq(chunk, depth, x_inc, y_inc, top_rcv_buffer, field)
!    ENDIF
!    IF ( (chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) ) THEN
!        CALL unpack_left_top_buffer_seq(chunk, depth, x_inc, y_inc, left_top_rcv_buffer, field)
!    ENDIF
!    IF ( (chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) ) THEN
!        CALL unpack_right_top_buffer_seq(chunk, depth, x_inc, y_inc, right_top_rcv_buffer, field)
!    ENDIF
!    IF ( (chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) ) THEN
!        CALL unpack_right_bottom_buffer_seq(chunk, depth, x_inc, y_inc, right_bottom_rcv_buffer, field)
!    ENDIF
!    IF ( (chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) ) THEN
!        CALL unpack_left_bottom_buffer_seq(chunk, depth, x_inc, y_inc, left_bottom_rcv_buffer, field)
!    ENDIF
!
!  ENDIF
!
!  IF(parallel%task.EQ.chunks(chunk)%task) THEN
!
!    IF(chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) THEN
!      receiver=chunks(chunks(chunk)%chunk_neighbours(chunk_bottom))%task
!      CALL SHMEM_PUT4_NB(top_write_flag, 1, 1, receiver)
!    ENDIF
!    IF(chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) THEN
!      receiver=chunks(chunks(chunk)%chunk_neighbours(chunk_top))%task
!      CALL SHMEM_PUT4_NB(bottom_write_flag, 1, 1, receiver)
!    ENDIF
!    IF(chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) THEN
!      receiver=chunks(chunks(chunk)%chunk_neighbours(chunk_left))%task
!      CALL SHMEM_PUT4_NB(right_write_flag, 1, 1, receiver)
!    ENDIF
!    IF(chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) THEN
!      receiver=chunks(chunks(chunk)%chunk_neighbours(chunk_right))%task
!      CALL SHMEM_PUT4_NB(left_write_flag, 1, 1, receiver)
!    ENDIF
!    IF ( (chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) ) THEN
!        receiver = chunks(chunks(chunk)%chunk_neighbours(chunk_left_top))%task
!        CALL SHMEM_PUT4_NB(right_bottom_write_flag, 1, 1, receiver)
!    ENDIF
!    IF ( (chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) ) THEN
!        receiver = chunks(chunks(chunk)%chunk_neighbours(chunk_right_top))%task
!        CALL SHMEM_PUT4_NB(left_bottom_write_flag, 1, 1, receiver)
!    ENDIF
!    IF ( (chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) ) THEN
!        receiver = chunks(chunks(chunk)%chunk_neighbours(chunk_right_bottom))%task
!        CALL SHMEM_PUT4_NB(left_top_write_flag, 1, 1, receiver)
!    ENDIF
!    IF ( (chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) ) THEN
!        receiver = chunks(chunks(chunk)%chunk_neighbours(chunk_left_bottom))%task
!        CALL SHMEM_PUT4_NB(right_top_write_flag, 1, 1, receiver)
!    ENDIF
!
!  ENDIF
!
!END SUBROUTINE clover_exchange_message






SUBROUTINE clover_sum(value)

  ! Only sums to the master

  IMPLICIT NONE

  REAL(KIND=8) :: value

  REAL(KIND=8) :: total

  pSync_sum = SHMEM_SYNC_VALUE
  sum_value=value
  sum_total=0

  CALL SHMEM_REAL8_SUM_TO_ALL(sum_total,sum_value,1,0,0,parallel%max_task,pWrk_sum,pSync_sum)

  value=sum_total

END SUBROUTINE clover_sum

SUBROUTINE clover_min(value)

  IMPLICIT NONE

  REAL(KIND=8) :: value

  pSync_min = SHMEM_SYNC_VALUE
  min_value = value

  CALL SHMEM_REAL8_MIN_TO_ALL(min_final, min_value, 1, 0, 0, parallel%max_task, pWrk_min, pSync_min)

  value = min_final

END SUBROUTINE clover_min

SUBROUTINE clover_max(value)

  IMPLICIT NONE

  REAL(KIND=8) :: value

  pSync_max = SHMEM_SYNC_VALUE
  max_value = value

  CALL SHMEM_REAL8_MAX_TO_ALL(max_final, max_value, 1, 0, 0, parallel%max_task, pWrk_max, pSync_max)

  value = max_final

END SUBROUTINE clover_max

SUBROUTINE clover_allgather(value,values)

  IMPLICIT NONE

  REAL(KIND=8) :: value, values(parallel%max_task)

  values(1)=value ! Just to ensure it will work in serial

  CALL SHMEM_FCOLLECT8(values, value, 1, 0, 0, parallel%max_task, pSync_collect)

END SUBROUTINE clover_allgather

SUBROUTINE clover_check_error(error)

  IMPLICIT NONE

  INTEGER :: error

  pSync_error = SHMEM_SYNC_VALUE
  error_value = error

  CALL SHMEM_INT4_MAX_TO_ALL(error_final, error_value, 1, 0, 0, parallel%max_task, pWrk_error, pSync_error)

  error = error_final

END SUBROUTINE clover_check_error


END MODULE clover_module
