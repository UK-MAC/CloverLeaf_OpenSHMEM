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

    INTEGER, PARAMETER :: NR = 1

    REAL(KIND=8) :: vol, mass, press, ie, ke, dt

    REAL(KIND=8) :: pWrk_pri(MAX(NR/2+1, SHMEM_REDUCE_MIN_WRKDATA_SIZE))
    INTEGER :: pSync_pri(SHMEM_REDUCE_SYNC_SIZE)
    DATA pSync_pri /SHMEM_REDUCE_SYNC_SIZE*SHMEM_SYNC_VALUE/

    REAL(KIND=8) :: pWrk_sec(MAX(NR/2+1, SHMEM_REDUCE_MIN_WRKDATA_SIZE))
    INTEGER :: pSync_sec(SHMEM_REDUCE_SYNC_SIZE)
    DATA pSync_sec /SHMEM_REDUCE_SYNC_SIZE*SHMEM_SYNC_VALUE/

    LOGICAL :: use_primary 

    INTEGER :: pSync_collect(SHMEM_COLLECT_SYNC_SIZE)

    COMMON /COLL/ vol, mass, press, ie, ke, dt

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

#ifdef ENABLE_SHMEM_FINALIZE
  CALL SHMEM_FINALIZE
#endif

END SUBROUTINE clover_abort

SUBROUTINE clover_finalize

  CLOSE(g_out)
  CALL FLUSH(0)
  CALL FLUSH(6)
  CALL FLUSH(g_out)

#ifdef ENABLE_SHMEM_FINALIZE
  CALL SHMEM_FINALIZE
#endif

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

    use_primary = .TRUE.

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

            IF (chunk .EQ. parallel%task+1) THEN

                left(1)   = (cx-1)*delta_x+1+add_x_prev
                right(1)  = left(1)+delta_x-1+add_x
                bottom(1) = (cy-1)*delta_y+1+add_y_prev
                top(1)    = bottom(1)+delta_y-1+add_y

                chunks(1)%chunk_neighbours(chunk_left)=chunk_x*(cy-1)+cx-1
                chunks(1)%chunk_neighbours(chunk_right)=chunk_x*(cy-1)+cx+1
                chunks(1)%chunk_neighbours(chunk_bottom)=chunk_x*(cy-2)+cx
                chunks(1)%chunk_neighbours(chunk_top)=chunk_x*(cy)+cx

                chunks(1)%chunk_neighbours(CHUNK_LEFT_TOP) = chunk_x*cy+cx-1
                chunks(1)%chunk_neighbours(CHUNK_LEFT_BOTTOM) = chunk_x*(cy-2)+cx-1
                chunks(1)%chunk_neighbours(CHUNK_RIGHT_TOP) = chunk_x*cy+cx+1
                chunks(1)%chunk_neighbours(CHUNK_RIGHT_BOTTOM) = chunk_x*(cy-2)+cx+1

                IF(cx.EQ.1) THEN
                  chunks(1)%chunk_neighbours(chunk_left)=external_face
                  chunks(1)%chunk_neighbours(CHUNK_LEFT_TOP)=external_face
                  chunks(1)%chunk_neighbours(CHUNK_LEFT_BOTTOM)=external_face
                ENDIF
                IF(cx.EQ.chunk_x) THEN
                  chunks(1)%chunk_neighbours(chunk_right)=external_face
                  chunks(1)%chunk_neighbours(CHUNK_RIGHT_TOP)=external_face
                  chunks(1)%chunk_neighbours(CHUNK_RIGHT_BOTTOM)=external_face
                ENDIF
                IF(cy.EQ.1) THEN 
                  chunks(1)%chunk_neighbours(chunk_bottom)=external_face
                  chunks(1)%chunk_neighbours(CHUNK_LEFT_BOTTOM)=external_face
                  chunks(1)%chunk_neighbours(CHUNK_RIGHT_BOTTOM)=external_face
                ENDIF
                IF(cy.EQ.chunk_y) THEN 
                  chunks(1)%chunk_neighbours(chunk_top)=external_face
                  chunks(1)%chunk_neighbours(CHUNK_LEFT_TOP)=external_face
                  chunks(1)%chunk_neighbours(CHUNK_RIGHT_TOP)=external_face
                ENDIF
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

    CALL SHPALLOC( pls, buffer_size_y, err, 1)
    CALL SHPALLOC( plr, buffer_size_y, err, 1)
    CALL SHPALLOC( prs, buffer_size_y, err, 1)
    CALL SHPALLOC( prr, buffer_size_y, err, 1)
    CALL SHPALLOC( pbs, buffer_size_x, err, 1)
    CALL SHPALLOC( pbr, buffer_size_x, err, 1)
    CALL SHPALLOC( pts, buffer_size_x, err, 1)
    CALL SHPALLOC( ptr, buffer_size_x, err, 1)

    CALL SHPALLOC(pltr, diag_buff_size, err, 1)
    CALL SHPALLOC(prtr, diag_buff_size, err, 1)
    CALL SHPALLOC(plbr, diag_buff_size, err, 1)
    CALL SHPALLOC(prbr, diag_buff_size, err, 1)

    CALL SHPALLOC(plts, diag_buff_size, err, 1)
    CALL SHPALLOC(prts, diag_buff_size, err, 1)
    CALL SHPALLOC(plbs, diag_buff_size, err, 1)
    CALL SHPALLOC(prbs, diag_buff_size, err, 1)

    WRITE(*,*) "Successfully shpalloc'd the comms buffers"


  ENDIF

END SUBROUTINE clover_allocate_buffers

SUBROUTINE clover_exchange(fields,depth)

  IMPLICIT NONE

  INTEGER      :: fields(:),depth,location_of_tasks_chunk

  ! Assuming 1 patch per task, this will be changed
  ! Also, not packing all fields for each communication, doing one at a time

    location_of_tasks_chunk = 1

  IF(fields(FIELD_DENSITY0).EQ.1) THEN
    CALL clover_exchange_message(location_of_tasks_chunk,chunks(location_of_tasks_chunk)%field%density0,      &
                                 left_snd_buffer,                                             &
                                 left_rcv_buffer,                                             &
                                 right_snd_buffer,                                            &
                                 right_rcv_buffer,                                            &
                                 bottom_snd_buffer,                                           &
                                 bottom_rcv_buffer,                                           &
                                 top_snd_buffer,                                              &
                                 top_rcv_buffer,                                              &
                                 left_top_snd_buffer,                                         &
                                 left_top_rcv_buffer,                                         &
                                 right_top_snd_buffer,                                        &
                                 right_top_rcv_buffer,                                        &
                                 right_bottom_snd_buffer,                                     &            
                                 right_bottom_rcv_buffer,                                     &
                                 left_bottom_snd_buffer,                                      &
                                 left_bottom_rcv_buffer,                                      &
                                 depth,CELL_DATA)
  ENDIF

  IF(fields(FIELD_DENSITY1).EQ.1) THEN
    CALL clover_exchange_message(location_of_tasks_chunk,chunks(location_of_tasks_chunk)%field%density1,      &
                                 left_snd_buffer,                                             &
                                 left_rcv_buffer,                                             &
                                 right_snd_buffer,                                            &
                                 right_rcv_buffer,                                            &
                                 bottom_snd_buffer,                                           &
                                 bottom_rcv_buffer,                                           &
                                 top_snd_buffer,                                              &
                                 top_rcv_buffer,                                              &
                                 left_top_snd_buffer,                                         &
                                 left_top_rcv_buffer,                                         &
                                 right_top_snd_buffer,                                        &
                                 right_top_rcv_buffer,                                        &
                                 right_bottom_snd_buffer,                                     &            
                                 right_bottom_rcv_buffer,                                     &
                                 left_bottom_snd_buffer,                                      &
                                 left_bottom_rcv_buffer,                                      &
                                 depth,CELL_DATA)
  ENDIF

  IF(fields(FIELD_ENERGY0).EQ.1) THEN
    CALL clover_exchange_message(location_of_tasks_chunk,chunks(location_of_tasks_chunk)%field%energy0,       &
                                 left_snd_buffer,                                             &
                                 left_rcv_buffer,                                             &
                                 right_snd_buffer,                                            &
                                 right_rcv_buffer,                                            &
                                 bottom_snd_buffer,                                           &
                                 bottom_rcv_buffer,                                           &
                                 top_snd_buffer,                                              &
                                 top_rcv_buffer,                                              &
                                 left_top_snd_buffer,                                         &
                                 left_top_rcv_buffer,                                         &
                                 right_top_snd_buffer,                                        &
                                 right_top_rcv_buffer,                                        &
                                 right_bottom_snd_buffer,                                     &            
                                 right_bottom_rcv_buffer,                                     &
                                 left_bottom_snd_buffer,                                      &
                                 left_bottom_rcv_buffer,                                      &
                                 depth,CELL_DATA)
  ENDIF

  IF(fields(FIELD_ENERGY1).EQ.1) THEN
    CALL clover_exchange_message(location_of_tasks_chunk,chunks(location_of_tasks_chunk)%field%energy1,       &
                                 left_snd_buffer,                                             &
                                 left_rcv_buffer,                                             &
                                 right_snd_buffer,                                            &
                                 right_rcv_buffer,                                            &
                                 bottom_snd_buffer,                                           &
                                 bottom_rcv_buffer,                                           &
                                 top_snd_buffer,                                              &
                                 top_rcv_buffer,                                              &
                                 left_top_snd_buffer,                                         &
                                 left_top_rcv_buffer,                                         &
                                 right_top_snd_buffer,                                        &
                                 right_top_rcv_buffer,                                        &
                                 right_bottom_snd_buffer,                                     &            
                                 right_bottom_rcv_buffer,                                     &
                                 left_bottom_snd_buffer,                                      &
                                 left_bottom_rcv_buffer,                                      &
                                 depth,CELL_DATA)
  ENDIF

  IF(fields(FIELD_PRESSURE).EQ.1) THEN
    CALL clover_exchange_message(location_of_tasks_chunk,chunks(location_of_tasks_chunk)%field%pressure,      &
                                 left_snd_buffer,                                             &
                                 left_rcv_buffer,                                             &
                                 right_snd_buffer,                                            &
                                 right_rcv_buffer,                                            &
                                 bottom_snd_buffer,                                           &
                                 bottom_rcv_buffer,                                           &
                                 top_snd_buffer,                                              &
                                 top_rcv_buffer,                                              &
                                 left_top_snd_buffer,                                         &
                                 left_top_rcv_buffer,                                         &
                                 right_top_snd_buffer,                                        &
                                 right_top_rcv_buffer,                                        &
                                 right_bottom_snd_buffer,                                     &            
                                 right_bottom_rcv_buffer,                                     &
                                 left_bottom_snd_buffer,                                      &
                                 left_bottom_rcv_buffer,                                      &
                                 depth,CELL_DATA)
  ENDIF

  IF(fields(FIELD_VISCOSITY).EQ.1) THEN
    CALL clover_exchange_message(location_of_tasks_chunk,chunks(location_of_tasks_chunk)%field%viscosity,     &
                                 left_snd_buffer,                                             &
                                 left_rcv_buffer,                                             &
                                 right_snd_buffer,                                            &
                                 right_rcv_buffer,                                            &
                                 bottom_snd_buffer,                                           &
                                 bottom_rcv_buffer,                                           &
                                 top_snd_buffer,                                              &
                                 top_rcv_buffer,                                              &
                                 left_top_snd_buffer,                                         &
                                 left_top_rcv_buffer,                                         &
                                 right_top_snd_buffer,                                        &
                                 right_top_rcv_buffer,                                        &
                                 right_bottom_snd_buffer,                                     &            
                                 right_bottom_rcv_buffer,                                     &
                                 left_bottom_snd_buffer,                                      &
                                 left_bottom_rcv_buffer,                                      &
                                 depth,CELL_DATA)
  ENDIF

  IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
    CALL clover_exchange_message(location_of_tasks_chunk,chunks(location_of_tasks_chunk)%field%soundspeed,    &
                                 left_snd_buffer,                                             &
                                 left_rcv_buffer,                                             &
                                 right_snd_buffer,                                            &
                                 right_rcv_buffer,                                            &
                                 bottom_snd_buffer,                                           &
                                 bottom_rcv_buffer,                                           &
                                 top_snd_buffer,                                              &
                                 top_rcv_buffer,                                              &
                                 left_top_snd_buffer,                                         &
                                 left_top_rcv_buffer,                                         &
                                 right_top_snd_buffer,                                        &
                                 right_top_rcv_buffer,                                        &
                                 right_bottom_snd_buffer,                                     &            
                                 right_bottom_rcv_buffer,                                     &
                                 left_bottom_snd_buffer,                                      &
                                 left_bottom_rcv_buffer,                                      &
                                 depth,CELL_DATA)
  ENDIF

  IF(fields(FIELD_XVEL0).EQ.1) THEN
    CALL clover_exchange_message(location_of_tasks_chunk,chunks(location_of_tasks_chunk)%field%xvel0,         &
                                 left_snd_buffer,                                             &
                                 left_rcv_buffer,                                             &
                                 right_snd_buffer,                                            &
                                 right_rcv_buffer,                                            &
                                 bottom_snd_buffer,                                           &
                                 bottom_rcv_buffer,                                           &
                                 top_snd_buffer,                                              &
                                 top_rcv_buffer,                                              &
                                 left_top_snd_buffer,                                         &
                                 left_top_rcv_buffer,                                         &
                                 right_top_snd_buffer,                                        &
                                 right_top_rcv_buffer,                                        &
                                 right_bottom_snd_buffer,                                     &            
                                 right_bottom_rcv_buffer,                                     &
                                 left_bottom_snd_buffer,                                      &
                                 left_bottom_rcv_buffer,                                      &
                                 depth,VERTEX_DATA)
  ENDIF

  IF(fields(FIELD_XVEL1).EQ.1) THEN
    CALL clover_exchange_message(location_of_tasks_chunk,chunks(location_of_tasks_chunk)%field%xvel1,         &
                                 left_snd_buffer,                                             &
                                 left_rcv_buffer,                                             &
                                 right_snd_buffer,                                            &
                                 right_rcv_buffer,                                            &
                                 bottom_snd_buffer,                                           &
                                 bottom_rcv_buffer,                                           &
                                 top_snd_buffer,                                              &
                                 top_rcv_buffer,                                              &
                                 left_top_snd_buffer,                                         &
                                 left_top_rcv_buffer,                                         &
                                 right_top_snd_buffer,                                        &
                                 right_top_rcv_buffer,                                        &
                                 right_bottom_snd_buffer,                                     &            
                                 right_bottom_rcv_buffer,                                     &
                                 left_bottom_snd_buffer,                                      &
                                 left_bottom_rcv_buffer,                                      &
                                 depth,VERTEX_DATA)
  ENDIF

  IF(fields(FIELD_YVEL0).EQ.1) THEN
    CALL clover_exchange_message(location_of_tasks_chunk,chunks(location_of_tasks_chunk)%field%yvel0,         &
                                 left_snd_buffer,                                             &
                                 left_rcv_buffer,                                             &
                                 right_snd_buffer,                                            &
                                 right_rcv_buffer,                                            &
                                 bottom_snd_buffer,                                           &
                                 bottom_rcv_buffer,                                           &
                                 top_snd_buffer,                                              &
                                 top_rcv_buffer,                                              &
                                 left_top_snd_buffer,                                         &
                                 left_top_rcv_buffer,                                         &
                                 right_top_snd_buffer,                                        &
                                 right_top_rcv_buffer,                                        &
                                 right_bottom_snd_buffer,                                     &            
                                 right_bottom_rcv_buffer,                                     &
                                 left_bottom_snd_buffer,                                      &
                                 left_bottom_rcv_buffer,                                      &
                                 depth,VERTEX_DATA)
  ENDIF

  IF(fields(FIELD_YVEL1).EQ.1) THEN
    CALL clover_exchange_message(location_of_tasks_chunk,chunks(location_of_tasks_chunk)%field%yvel1,         &
                                 left_snd_buffer,                                             &
                                 left_rcv_buffer,                                             &
                                 right_snd_buffer,                                            &
                                 right_rcv_buffer,                                            &
                                 bottom_snd_buffer,                                           &
                                 bottom_rcv_buffer,                                           &
                                 top_snd_buffer,                                              &
                                 top_rcv_buffer,                                              &
                                 left_top_snd_buffer,                                         &
                                 left_top_rcv_buffer,                                         &
                                 right_top_snd_buffer,                                        &
                                 right_top_rcv_buffer,                                        &
                                 right_bottom_snd_buffer,                                     &            
                                 right_bottom_rcv_buffer,                                     &
                                 left_bottom_snd_buffer,                                      &
                                 left_bottom_rcv_buffer,                                      &
                                 depth,VERTEX_DATA)
  ENDIF

  IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
    CALL clover_exchange_message(location_of_tasks_chunk,chunks(location_of_tasks_chunk)%field%vol_flux_x,    &
                                 left_snd_buffer,                                             &
                                 left_rcv_buffer,                                             &
                                 right_snd_buffer,                                            &
                                 right_rcv_buffer,                                            &
                                 bottom_snd_buffer,                                           &
                                 bottom_rcv_buffer,                                           &
                                 top_snd_buffer,                                              &
                                 top_rcv_buffer,                                              &
                                 left_top_snd_buffer,                                         &
                                 left_top_rcv_buffer,                                         &
                                 right_top_snd_buffer,                                        &
                                 right_top_rcv_buffer,                                        &
                                 right_bottom_snd_buffer,                                     &            
                                 right_bottom_rcv_buffer,                                     &
                                 left_bottom_snd_buffer,                                      &
                                 left_bottom_rcv_buffer,                                      &
                                 depth,X_FACE_DATA)
  ENDIF

  IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
    CALL clover_exchange_message(location_of_tasks_chunk,chunks(location_of_tasks_chunk)%field%vol_flux_y,    &
                                 left_snd_buffer,                                             &
                                 left_rcv_buffer,                                             &
                                 right_snd_buffer,                                            &
                                 right_rcv_buffer,                                            &
                                 bottom_snd_buffer,                                           &
                                 bottom_rcv_buffer,                                           &
                                 top_snd_buffer,                                              &
                                 top_rcv_buffer,                                              &
                                 left_top_snd_buffer,                                         &
                                 left_top_rcv_buffer,                                         &
                                 right_top_snd_buffer,                                        &
                                 right_top_rcv_buffer,                                        &
                                 right_bottom_snd_buffer,                                     &            
                                 right_bottom_rcv_buffer,                                     &
                                 left_bottom_snd_buffer,                                      &
                                 left_bottom_rcv_buffer,                                      &
                                 depth,Y_FACE_DATA)
  ENDIF

  IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
    CALL clover_exchange_message(location_of_tasks_chunk,chunks(location_of_tasks_chunk)%field%mass_flux_x,   &
                                 left_snd_buffer,                                             &
                                 left_rcv_buffer,                                             &
                                 right_snd_buffer,                                            &
                                 right_rcv_buffer,                                            &
                                 bottom_snd_buffer,                                           &
                                 bottom_rcv_buffer,                                           &
                                 top_snd_buffer,                                              &
                                 top_rcv_buffer,                                              &
                                 left_top_snd_buffer,                                         &
                                 left_top_rcv_buffer,                                         &
                                 right_top_snd_buffer,                                        &
                                 right_top_rcv_buffer,                                        &
                                 right_bottom_snd_buffer,                                     &            
                                 right_bottom_rcv_buffer,                                     &
                                 left_bottom_snd_buffer,                                      &
                                 left_bottom_rcv_buffer,                                      &
                                 depth,X_FACE_DATA)
  ENDIF

  IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
    CALL clover_exchange_message(location_of_tasks_chunk,chunks(location_of_tasks_chunk)%field%mass_flux_y,   &
                                 left_snd_buffer,                                             &
                                 left_rcv_buffer,                                             &
                                 right_snd_buffer,                                            &
                                 right_rcv_buffer,                                            &
                                 bottom_snd_buffer,                                           &
                                 bottom_rcv_buffer,                                           &
                                 top_snd_buffer,                                              &
                                 top_rcv_buffer,                                              &
                                 left_top_snd_buffer,                                         &
                                 left_top_rcv_buffer,                                         &
                                 right_top_snd_buffer,                                        &
                                 right_top_rcv_buffer,                                        &
                                 right_bottom_snd_buffer,                                     &            
                                 right_bottom_rcv_buffer,                                     &
                                 left_bottom_snd_buffer,                                      &
                                 left_bottom_rcv_buffer,                                      &
                                 depth,Y_FACE_DATA)
  ENDIF

END SUBROUTINE clover_exchange

SUBROUTINE clover_exchange_message(chunk,field,                            &
                                   left_snd_buffer,                        &
                                   left_rcv_buffer,                        &
                                   right_snd_buffer,                       &
                                   right_rcv_buffer,                       &
                                   bottom_snd_buffer,                      &
                                   bottom_rcv_buffer,                      &
                                   top_snd_buffer,                         &
                                   top_rcv_buffer,                         &
                                   left_top_snd_buffer,                                         &
                                   left_top_rcv_buffer,                                         &
                                   right_top_snd_buffer,                                        &
                                   right_top_rcv_buffer,                                        &
                                   right_bottom_snd_buffer,                                     &            
                                   right_bottom_rcv_buffer,                                     &
                                   left_bottom_snd_buffer,                                      &
                                   left_bottom_rcv_buffer,                                      &
                                   depth,field_type)

  !USE pack_kernel_module

  IMPLICIT NONE

  REAL(KIND=8) :: field(-1:,-1:) ! This seems to work for any type of mesh data
  REAL(KIND=8) :: left_snd_buffer(:),left_rcv_buffer(:),right_snd_buffer(:),right_rcv_buffer(:)
  REAL(KIND=8) :: bottom_snd_buffer(:),bottom_rcv_buffer(:),top_snd_buffer(:),top_rcv_buffer(:)

  REAL(KIND=8) :: left_top_snd_buffer(:),left_top_rcv_buffer(:),right_top_snd_buffer(:),right_top_rcv_buffer(:)
  REAL(KIND=8) :: right_bottom_snd_buffer(:),right_bottom_rcv_buffer(:),left_bottom_snd_buffer(:),left_bottom_rcv_buffer(:)

  INTEGER      :: chunk,depth,field_type

  INTEGER      :: size,x_inc,y_inc
  INTEGER      :: receiver,sender
  INTEGER      :: topedge, bottomedge, leftedge, rightedge

  ! Field type will either be cell, vertex, x_face or y_face to get the message limits correct

  ! I am packing my own buffers. I am sure this could be improved with MPI data types
  !  but this will do for now

  ! I am also sending buffers to chunks with the same task id for now.
  ! This can be improved in the future but at the moment there is just 1 chunk per task anyway

  ! The tag will be a function of the sending chunk and the face it is coming from
  !  like chunk 6 sending the left face

  ! No open mp in here either. May be beneficial will packing and unpacking in the future, though I am not sure.

  ! Change this so it will allow more than 1 chunk per task

  ! Pack and send

  ! These array modifications still need to be added on, plus the donor data location changes as in update_halo
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

  ! Pack real data into buffers
  IF(parallel%task.EQ.chunks(chunk)%task) THEN

    !pack all the data buffers required 
    IF(chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) THEN
        CALL pack_left_buffer_seq(chunk, depth, x_inc, y_inc, left_snd_buffer, field)
    ENDIF
    IF(chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) THEN
        CALL pack_right_buffer_seq(chunk, depth, x_inc, y_inc, right_snd_buffer, field)
    ENDIF
    IF(chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) THEN
        CALL pack_bottom_buffer_seq(chunk, depth, x_inc, y_inc, bottom_snd_buffer, field)
    ENDIF
    IF(chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) THEN
        CALL pack_top_buffer_seq(chunk, depth, x_inc, y_inc, top_snd_buffer, field)
    ENDIF
    IF ( (chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) ) THEN
        CALL pack_left_top_buffer_seq(chunk, depth, x_inc, y_inc, left_top_snd_buffer, field)
    ENDIF
    IF ( (chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) ) THEN
        CALL pack_right_top_buffer_seq(chunk, depth, x_inc, y_inc, right_top_snd_buffer, field)
    ENDIF
    IF ( (chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) ) THEN
        CALL pack_right_bottom_buffer_seq(chunk, depth, x_inc, y_inc, right_bottom_snd_buffer, field)
    ENDIF
    IF ( (chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) ) THEN
        CALL pack_left_bottom_buffer_seq(chunk, depth, x_inc, y_inc, left_bottom_snd_buffer, field)
    ENDIF


    topedge = 0
    bottomedge = 0
    IF (chunks(chunk)%chunk_neighbours(chunk_top).EQ.external_face) THEN
      topedge = depth
    ENDIF
    IF (chunks(chunk)%chunk_neighbours(chunk_bottom).EQ.external_face) THEN
      bottomedge = depth
    ENDIF
    size=(1+(chunks(chunk)%field%y_max+y_inc+topedge)-(chunks(chunk)%field%y_min-bottomedge))*depth


    ! Send/receive the data
    IF(chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) THEN
      receiver=chunks(chunk)%chunk_neighbours(chunk_left)-1

      IF (left_write_flag .EQ. 0) THEN
        CALL SHMEM_INT4_WAIT_UNTIL(left_write_flag, SHMEM_CMP_EQ, 1)
      ENDIF
      CALL SHMEM_PUT64(right_rcv_buffer, left_snd_buffer, size, receiver)
      
      left_write_flag = 0 

    ENDIF

    IF(chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) THEN
      receiver=chunks(chunk)%chunk_neighbours(chunk_right)-1

      IF (right_write_flag .EQ. 0) THEN
        CALL SHMEM_INT4_WAIT_UNTIL(right_write_flag, SHMEM_CMP_EQ, 1)
      ENDIF
      CALL SHMEM_PUT64(left_rcv_buffer, right_snd_buffer, size, receiver)

      right_write_flag = 0 

    ENDIF

    leftedge= 0
    rightedge= 0
    IF (chunks(chunk)%chunk_neighbours(chunk_left).EQ.external_face) THEN
      leftedge = depth
    ENDIF
    IF (chunks(chunk)%chunk_neighbours(chunk_right).EQ.external_face) THEN
      rightedge = depth
    ENDIF
    size=(1+(chunks(chunk)%field%x_max+x_inc+rightedge)-(chunks(chunk)%field%x_min-leftedge))*depth

    ! Send/receive the data
    IF(chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) THEN
      receiver=chunks(chunk)%chunk_neighbours(chunk_bottom) - 1

      IF (bottom_write_flag .EQ. 0) THEN
        CALL SHMEM_INT4_WAIT_UNTIL(bottom_write_flag, SHMEM_CMP_EQ, 1)
      ENDIF

      CALL SHMEM_PUT64(top_rcv_buffer, bottom_snd_buffer, size, receiver)
      bottom_write_flag = 0
    ENDIF

    IF(chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) THEN
        receiver=chunks(chunk)%chunk_neighbours(chunk_top) - 1

        IF (top_write_flag .EQ. 0) THEN
            CALL SHMEM_INT4_WAIT_UNTIL(top_write_flag, SHMEM_CMP_EQ, 1)
        ENDIF

        CALL SHMEM_PUT64(bottom_rcv_buffer, top_snd_buffer, size, receiver)
        top_write_flag = 0 
    ENDIF


    size = depth*depth

    IF ( (chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) ) THEN
        receiver = chunks(chunk)%chunk_neighbours(chunk_left_top) - 1

        IF (left_top_write_flag .EQ. 0) THEN
            CALL SHMEM_INT4_WAIT_UNTIL(left_top_write_flag, SHMEM_CMP_EQ, 1)
        ENDIF

        CALL SHMEM_PUT64(right_bottom_rcv_buffer, left_top_snd_buffer, size, receiver)
        left_top_write_flag = 0
    ENDIF
    
    IF ( (chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) ) THEN
        receiver = chunks(chunk)%chunk_neighbours(chunk_right_top) - 1

        IF (right_top_write_flag .EQ. 0) THEN
            CALL SHMEM_INT4_WAIT_UNTIL(right_top_write_flag, SHMEM_CMP_EQ, 1)
        ENDIF

        CALL SHMEM_PUT64(left_bottom_rcv_buffer, right_top_snd_buffer, size, receiver)
        right_top_write_flag = 0
    ENDIF
    
    IF ( (chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) ) THEN
        receiver = chunks(chunk)%chunk_neighbours(chunk_right_bottom) - 1

        IF (right_bottom_write_flag .EQ. 0) THEN
            CALL SHMEM_INT4_WAIT_UNTIL(right_bottom_write_flag, SHMEM_CMP_EQ, 1)
        ENDIF

        CALL SHMEM_PUT64(left_top_rcv_buffer, right_bottom_snd_buffer, size, receiver)
        right_bottom_write_flag = 0
    ENDIF
    
    IF ( (chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) ) THEN
       receiver = chunks(chunk)%chunk_neighbours(chunk_left_bottom) - 1

       IF (left_bottom_write_flag .EQ. 0) THEN
           CALL SHMEM_INT4_WAIT_UNTIL(left_bottom_write_flag, SHMEM_CMP_EQ, 1)
       ENDIF

       CALL SHMEM_PUT64(right_top_rcv_buffer, left_bottom_snd_buffer, size, receiver)
       left_bottom_write_flag = 0
    ENDIF

  ENDIF

  ! Wait for the messages
#ifdef FENCE_NOT_QUIET
        CALL SHMEM_FENCE
#else
        CALL SHMEM_QUIET
#endif




  IF(parallel%task.EQ.chunks(chunk)%task) THEN

    ! Send/receive the data
    IF(chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) THEN
      receiver=chunks(chunk)%chunk_neighbours(chunk_left) - 1
      CALL SHMEM_PUT4(right_rcv_flag, 1, 1, receiver)
    ENDIF

    IF(chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) THEN
      receiver=chunks(chunk)%chunk_neighbours(chunk_right) - 1
      CALL SHMEM_PUT4(left_rcv_flag, 1, 1, receiver)
    ENDIF

    IF(chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) THEN
      receiver=chunks(chunk)%chunk_neighbours(chunk_bottom) - 1
      CALL SHMEM_PUT4(top_rcv_flag, 1, 1, receiver)
    ENDIF

    IF(chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) THEN
      receiver=chunks(chunk)%chunk_neighbours(chunk_top) - 1
      CALL SHMEM_PUT4(bottom_rcv_flag, 1, 1, receiver)
    ENDIF
    IF ( (chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) ) THEN
        receiver = chunks(chunk)%chunk_neighbours(chunk_left_top) - 1
        CALL SHMEM_PUT4(right_bottom_rcv_flag, 1, 1, receiver)
    ENDIF
    IF ( (chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) ) THEN
        receiver = chunks(chunk)%chunk_neighbours(chunk_right_top) - 1
        CALL SHMEM_PUT4(left_bottom_rcv_flag, 1, 1, receiver)
    ENDIF
    IF ( (chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) ) THEN
        receiver = chunks(chunk)%chunk_neighbours(chunk_right_bottom) - 1
        CALL SHMEM_PUT4(left_top_rcv_flag, 1, 1, receiver)
    ENDIF
    IF ( (chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) ) THEN
        receiver = chunks(chunk)%chunk_neighbours(chunk_left_bottom) - 1
        CALL SHMEM_PUT4(right_top_rcv_flag, 1, 1, receiver)
    ENDIF
  ENDIF

  IF(parallel%task.EQ.chunks(chunk)%task) THEN

    ! Send/receive the data
    IF(chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) THEN
      IF (left_rcv_flag .EQ. 0) THEN
        CALL SHMEM_INT4_WAIT_UNTIL(left_rcv_flag, SHMEM_CMP_EQ, 1)
      ENDIF
      left_rcv_flag = 0
    ENDIF

    IF(chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) THEN
      IF (right_rcv_flag .EQ. 0) THEN
        CALL SHMEM_INT4_WAIT_UNTIL(right_rcv_flag, SHMEM_CMP_EQ, 1)
      ENDIF
      right_rcv_flag = 0
    ENDIF

    IF(chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) THEN
      IF (bottom_rcv_flag .EQ. 0) THEN
        CALL SHMEM_INT4_WAIT_UNTIL(bottom_rcv_flag, SHMEM_CMP_EQ, 1)
      ENDIF
      bottom_rcv_flag = 0
    ENDIF

    IF(chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) THEN
      IF (top_rcv_flag .EQ. 0) THEN
        CALL SHMEM_INT4_WAIT_UNTIL(top_rcv_flag, SHMEM_CMP_EQ, 1)
      ENDIF
      top_rcv_flag = 0
    ENDIF
    IF ( (chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) ) THEN
        IF (left_top_rcv_flag .EQ. 0) THEN
            CALL SHMEM_INT4_WAIT_UNTIL(left_top_rcv_flag, SHMEM_CMP_EQ, 1)
        ENDIF
        left_top_rcv_flag = 0
    ENDIF
    IF ( (chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) ) THEN
        IF (right_top_rcv_flag .EQ. 0) THEN
            CALL SHMEM_INT4_WAIT_UNTIL(right_top_rcv_flag, SHMEM_CMP_EQ, 1)
        ENDIF
        right_top_rcv_flag = 0
    ENDIF
    IF ( (chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) ) THEN
        IF (right_bottom_rcv_flag .EQ. 0) THEN
            CALL SHMEM_INT4_WAIT_UNTIL(right_bottom_rcv_flag, SHMEM_CMP_EQ, 1)
        ENDIF
        right_bottom_rcv_flag = 0
    ENDIF
    IF ( (chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) ) THEN
        IF (left_bottom_rcv_flag .EQ. 0) THEN
            CALL SHMEM_INT4_WAIT_UNTIL(left_bottom_rcv_flag, SHMEM_CMP_EQ, 1)
        ENDIF
        left_bottom_rcv_flag = 0
    ENDIF
  ENDIF





  ! Unpack buffers in halo cells
  IF(parallel%task.EQ.chunks(chunk)%task) THEN

    IF(chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) THEN
        CALL unpack_left_buffer_seq(chunk, depth, x_inc, y_inc, left_rcv_buffer, field)
    ENDIF
    IF(chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) THEN
        CALL unpack_right_buffer_seq(chunk, depth, x_inc, y_inc, right_rcv_buffer, field)
    ENDIF
    IF(chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) THEN
        CALL unpack_bottom_buffer_seq(chunk, depth, x_inc, y_inc, bottom_rcv_buffer, field)
    ENDIF
    IF(chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) THEN
        CALL unpack_top_buffer_seq(chunk, depth, x_inc, y_inc, top_rcv_buffer, field)
    ENDIF
    IF ( (chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) ) THEN
        CALL unpack_left_top_buffer_seq(chunk, depth, x_inc, y_inc, left_top_rcv_buffer, field)
    ENDIF
    IF ( (chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) ) THEN
        CALL unpack_right_top_buffer_seq(chunk, depth, x_inc, y_inc, right_top_rcv_buffer, field)
    ENDIF
    IF ( (chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) ) THEN
        CALL unpack_right_bottom_buffer_seq(chunk, depth, x_inc, y_inc, right_bottom_rcv_buffer, field)
    ENDIF
    IF ( (chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) ) THEN
        CALL unpack_left_bottom_buffer_seq(chunk, depth, x_inc, y_inc, left_bottom_rcv_buffer, field)
    ENDIF

  ENDIF

  IF(parallel%task.EQ.chunks(chunk)%task) THEN

    IF(chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) THEN
      receiver=chunks(chunk)%chunk_neighbours(chunk_bottom) - 1
      CALL SHMEM_PUT4(top_write_flag, 1, 1, receiver)
    ENDIF
    IF(chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) THEN
      receiver=chunks(chunk)%chunk_neighbours(chunk_top) - 1
      CALL SHMEM_PUT4(bottom_write_flag, 1, 1, receiver)
    ENDIF
    IF(chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) THEN
      receiver=chunks(chunk)%chunk_neighbours(chunk_left) - 1
      CALL SHMEM_PUT4(right_write_flag, 1, 1, receiver)
    ENDIF
    IF(chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) THEN
      receiver=chunks(chunk)%chunk_neighbours(chunk_right) - 1
      CALL SHMEM_PUT4(left_write_flag, 1, 1, receiver)
    ENDIF
    IF ( (chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) ) THEN
        receiver = chunks(chunk)%chunk_neighbours(chunk_left_top) - 1
        CALL SHMEM_PUT4(right_bottom_write_flag, 1, 1, receiver)
    ENDIF
    IF ( (chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) ) THEN
        receiver = chunks(chunk)%chunk_neighbours(chunk_right_top) - 1
        CALL SHMEM_PUT4(left_bottom_write_flag, 1, 1, receiver)
    ENDIF
    IF ( (chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) ) THEN
        receiver = chunks(chunk)%chunk_neighbours(chunk_right_bottom) - 1
        CALL SHMEM_PUT4(left_top_write_flag, 1, 1, receiver)
    ENDIF
    IF ( (chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) ) THEN
        receiver = chunks(chunk)%chunk_neighbours(chunk_left_bottom) - 1
        CALL SHMEM_PUT4(right_top_write_flag, 1, 1, receiver)
    ENDIF

  ENDIF

END SUBROUTINE clover_exchange_message

SUBROUTINE clover_sum(value)

    IMPLICIT NONE

    REAL(KIND=8) :: value

    !value will be different each time as called with a different variable within field summary 

    IF (use_primary) THEN
        CALL SHMEM_REAL8_SUM_TO_ALL(value,value,1,0,0,parallel%max_task,pWrk_pri,pSync_pri)
    ELSE
        CALL SHMEM_REAL8_SUM_TO_ALL(value,value,1,0,0,parallel%max_task,pWrk_sec,pSync_sec)
    ENDIF
    
    use_primary = .NOT. use_primary

END SUBROUTINE clover_sum

SUBROUTINE clover_min(value)

    IMPLICIT NONE

    REAL(KIND=8) :: value
    ! there will always ne another collective between calls to min so ok to use same min value 

    IF (use_primary) THEN
        CALL SHMEM_REAL8_MIN_TO_ALL(value, value, 1, 0, 0, parallel%max_task, pWrk_pri, pSync_pri)
    ELSE
        CALL SHMEM_REAL8_MIN_TO_ALL(value, value, 1, 0, 0, parallel%max_task, pWrk_sec, pSync_sec)
    ENDIF

    use_primary = .NOT. use_primary

END SUBROUTINE clover_min

SUBROUTINE clover_max(value)

    IMPLICIT NONE

    REAL(KIND=8) :: value

    IF (use_primary) THEN
        CALL SHMEM_REAL8_MAX_TO_ALL(value, value, 1, 0, 0, parallel%max_task, pWrk_pri, pSync_pri)
    ELSE
        CALL SHMEM_REAL8_MAX_TO_ALL(value, value, 1, 0, 0, parallel%max_task, pWrk_sec, pSync_sec)
    ENDIF

    use_primary = .NOT. use_primary

END SUBROUTINE clover_max

SUBROUTINE clover_allgather(value,values)

    IMPLICIT NONE

    REAL(KIND=8) :: value, values(parallel%max_task)

    values(1)=value ! Just to ensure it will work in serial

    CALL SHMEM_FCOLLECT8(values, value, 1, 0, 0, parallel%max_task, pSync_collect)

END SUBROUTINE clover_allgather

SUBROUTINE clover_check_error(error)

    IMPLICIT NONE

    REAL(KIND=8) :: error

    IF (use_primary) THEN
        CALL SHMEM_REAL8_MAX_TO_ALL(error, error, 1, 0, 0, parallel%max_task, pWrk_pri, pSync_pri)
    ELSE
        CALL SHMEM_REAL8_MAX_TO_ALL(error, error, 1, 0, 0, parallel%max_task, pWrk_sec, pSync_sec)
    ENDIF

    use_primary = .NOT. use_primary
    
END SUBROUTINE clover_check_error


END MODULE clover_module
