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
  USE iso_c_binding

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

  COMMON/FLAG/left_rcv_flag, right_rcv_flag, left_write_flag, right_write_flag, top_rcv_flag, bottom_rcv_flag, top_write_flag, bottom_write_flag

CONTAINS

SUBROUTINE clover_barrier

  CALL SHMEM_BARRIER_ALL()

END SUBROUTINE clover_barrier

SUBROUTINE clover_abort

  CALL SHMEM_FINALIZE

END SUBROUTINE clover_abort

SUBROUTINE clover_finalize

  INTEGER :: err

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
      IF(cx.EQ.1)chunks(chunk)%chunk_neighbours(chunk_left)=external_face
      IF(cx.EQ.chunk_x)chunks(chunk)%chunk_neighbours(chunk_right)=external_face
      IF(cy.EQ.1)chunks(chunk)%chunk_neighbours(chunk_bottom)=external_face
      IF(cy.EQ.chunk_y)chunks(chunk)%chunk_neighbours(chunk_top)=external_face
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
  INTEGER      :: buffer_size_x, buffer_size_y
  REAL(KIND=8) :: r8default 
  
  ! Unallocated buffers for external boundaries caused issues on some systems so they are now
  !  all allocated
  IF(parallel%task.EQ.chunks(chunk)%task)THEN

    buffer_size_x = kind(r8default)/kind(idefault)*2*(chunks(chunk)%field%x_max+5)
    buffer_size_y = kind(r8default)/kind(idefault)*2*(chunks(chunk)%field%y_max+5)

    CALL SHPALLOC( pls, buffer_size_y, err, 0)
    CALL SHPALLOC( plr, buffer_size_y, err, 0)
    CALL SHPALLOC( prs, buffer_size_y, err, 0)
    CALL SHPALLOC( prr, buffer_size_y, err, 0)
    CALL SHPALLOC( pbs, buffer_size_x, err, 0)
    CALL SHPALLOC( pbr, buffer_size_x, err, 0)
    CALL SHPALLOC( pts, buffer_size_x, err, 0)
    CALL SHPALLOC( ptr, buffer_size_x, err, 0)

    !!IF(chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) THEN
    !  ALLOCATE(chunks(chunk)%left_snd_buffer(2*(chunks(chunk)%field%y_max+5)))
    !  ALLOCATE(chunks(chunk)%left_rcv_buffer(2*(chunks(chunk)%field%y_max+5)))
    !!ENDIF
    !!IF(chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) THEN
    !  ALLOCATE(chunks(chunk)%right_snd_buffer(2*(chunks(chunk)%field%y_max+5)))
    !  ALLOCATE(chunks(chunk)%right_rcv_buffer(2*(chunks(chunk)%field%y_max+5)))
    !!ENDIF
    !!IF(chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) THEN
    !  ALLOCATE(chunks(chunk)%bottom_snd_buffer(2*(chunks(chunk)%field%x_max+5)))
    !  ALLOCATE(chunks(chunk)%bottom_rcv_buffer(2*(chunks(chunk)%field%x_max+5)))
    !!ENDIF
    !!IF(chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) THEN
    !  ALLOCATE(chunks(chunk)%top_snd_buffer(2*(chunks(chunk)%field%x_max+5)))
    !  ALLOCATE(chunks(chunk)%top_rcv_buffer(2*(chunks(chunk)%field%x_max+5)))
    !!ENDIF
  ENDIF

END SUBROUTINE clover_allocate_buffers

SUBROUTINE clover_exchange(fields,depth)

  IMPLICIT NONE

  INTEGER      :: fields(:),depth
    
  ! Assuming 1 patch per task, this will be changed
  ! Also, not packing all fields for each communication, doing one at a time

  IF(fields(FIELD_DENSITY0).EQ.1) THEN
    CALL clover_exchange_message(parallel%task+1,density0,      &
                                 depth,CELL_DATA)
  ENDIF

  IF(fields(FIELD_DENSITY1).EQ.1) THEN
    CALL clover_exchange_message(parallel%task+1,density1,      &
                                 depth,CELL_DATA)
  ENDIF

  IF(fields(FIELD_ENERGY0).EQ.1) THEN
    CALL clover_exchange_message(parallel%task+1,energy0,       &
                                 depth,CELL_DATA)
  ENDIF

  IF(fields(FIELD_ENERGY1).EQ.1) THEN
    CALL clover_exchange_message(parallel%task+1,energy1,       &
                                 depth,CELL_DATA)
  ENDIF

  IF(fields(FIELD_PRESSURE).EQ.1) THEN
    CALL clover_exchange_message(parallel%task+1,pressure,      &
                                 depth,CELL_DATA)
  ENDIF

  IF(fields(FIELD_VISCOSITY).EQ.1) THEN
    CALL clover_exchange_message(parallel%task+1,viscosity,     &
                                 depth,CELL_DATA)
  ENDIF

  IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
    CALL clover_exchange_message(parallel%task+1,soundspeed,    &
                                 depth,CELL_DATA)
  ENDIF

  IF(fields(FIELD_XVEL0).EQ.1) THEN
    CALL clover_exchange_message(parallel%task+1,xvel0,         &
                                 depth,VERTEX_DATA)
  ENDIF

  IF(fields(FIELD_XVEL1).EQ.1) THEN
    CALL clover_exchange_message(parallel%task+1,xvel1,         &
                                 depth,VERTEX_DATA)
  ENDIF

  IF(fields(FIELD_YVEL0).EQ.1) THEN
    CALL clover_exchange_message(parallel%task+1,yvel0,         &
                                 depth,VERTEX_DATA)
  ENDIF

  IF(fields(FIELD_YVEL1).EQ.1) THEN
    CALL clover_exchange_message(parallel%task+1,yvel1,         &
                                 depth,VERTEX_DATA)
  ENDIF

  IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
    CALL clover_exchange_message(parallel%task+1,vol_flux_x,    &
                                 depth,X_FACE_DATA)
  ENDIF

  IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
    CALL clover_exchange_message(parallel%task+1,vol_flux_y,    &
                                 depth,Y_FACE_DATA)
  ENDIF

  IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
    CALL clover_exchange_message(parallel%task+1,mass_flux_x,   &
                                 depth,X_FACE_DATA)
  ENDIF

  IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
    CALL clover_exchange_message(parallel%task+1,mass_flux_y,   &
                                 depth,Y_FACE_DATA)
  ENDIF


END SUBROUTINE clover_exchange

SUBROUTINE clover_exchange_message(chunk,field,                            &
                                   depth,field_type)

  !USE pack_kernel_module

  IMPLICIT NONE

  REAL(KIND=8) :: field(-1:,-1:) ! This seems to work for any type of mesh data
  !REAL(KIND=8) :: left_snd_buffer(:),left_rcv_buffer(:),right_snd_buffer(:),right_rcv_buffer(:)
  !REAL(KIND=8) :: bottom_snd_buffer(:),bottom_rcv_buffer(:),top_snd_buffer(:),top_rcv_buffer(:)

  INTEGER      :: chunk,depth,field_type

  INTEGER      :: size,x_inc,y_inc
  INTEGER      :: receiver,sender, num_elements, array_width

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

    CALL SHMEM_BARRIER_ALL

  ! Pack real data into buffers
  IF(parallel%task.EQ.chunks(chunk)%task) THEN

        num_elements = 1+(chunks(chunk)%field%y_max+y_inc+depth)-(chunks(chunk)%field%y_min-depth)
        !size=(1+(chunks(chunk)%field%y_max+y_inc+depth)-(chunks(chunk)%field%y_min-depth))*depth
        WRITE(*,*) "Task: ", parallel%task, " num elements comms: ", num_elements

        array_width = chunks(chunk)%field%x_max+4+x_inc
        WRITE(*,*) "Task: ", parallel%task, " array width: ", array_width

    ! Send/receive the data
        IF(chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) THEN
            receiver=chunks(chunks(chunk)%chunk_neighbours(chunk_left))%task

            WRITE(*,*) "Task: ", parallel%task, "sending left to: ", receiver, " with depth: ", depth

            IF (depth .EQ. 1) THEN

                WRITE(*,*) "Task : ", parallel%task, " send to left to: ", receiver, " stride: ", array_width, " depth: ", depth, " num elements: ", num_elements, &
                                                     " source x: ", chunks(chunk)%field%x_min+x_inc, &
                                                     " source y: ", chunks(chunk)%field%y_min-depth, &
                                                     " target x: ", chunks(chunk)%field%x_max+x_inc+1, &
                                                     " target y: ", chunks(chunk)%field%y_min-depth 


                CALL SHMEM_IPUT64(field(chunks(chunk)%field%x_max+x_inc+1,chunks(chunk)%field%y_min-depth), &
                                  field(chunks(chunk)%field%x_min+x_inc,chunks(chunk)%field%y_min-depth), &
                                  array_width, array_width,                                  &
                                  num_elements, receiver)

            ELSE

                CALL SHMEM_IPUT64(field(chunks(chunk)%field%x_max+x_inc+1,chunks(chunk)%field%y_min-depth), &
                                  field(chunks(chunk)%field%x_min+x_inc,chunks(chunk)%field%y_min-depth), &
                                  array_width, array_width,                                  &
                                  num_elements, receiver)
                WRITE(*,*) "Task: ", parallel%task, " first send to left to: ", receiver,  " stride: ", array_width, " depth: ", depth,  &
                                                                       " source x: ", chunks(chunk)%field%x_min+x_inc,  &
                                                                       " source y: ", chunks(chunk)%field%y_min-depth,  &
                                                                       " target x: ", chunks(chunk)%field%x_max+x_inc+1,& 
                                                                       " target y: ", chunks(chunk)%field%y_min-depth  

                CALL SHMEM_IPUT64(field(chunks(chunk)%field%x_max+x_inc+2,chunks(chunk)%field%y_min-depth), &
                                  field(chunks(chunk)%field%x_min+x_inc+1,chunks(chunk)%field%y_min-depth), &
                                  array_width, array_width,                                  &
                                  num_elements, receiver)
                WRITE(*,*) "Task: ", parallel%task, " second sec to left to: ", receiver, " stride: ", array_width, " depth: ", depth, &
                                                                       " source x: ", chunks(chunk)%field%x_min+x_inc+1, &
                                                                       " source y: ", chunks(chunk)%field%y_min-depth, &
                                                                       " target x: ", chunks(chunk)%field%x_max+x_inc+2, &
                                                                       " target y: ", chunks(chunk)%field%y_min-depth
            ENDIF


          !IF (left_write_flag .EQ. 0) THEN
          !  CALL SHMEM_INT4_WAIT_UNTIL(left_write_flag, SHMEM_CMP_EQ, 1)
          !ENDIF
          !CALL SHMEM_PUT64_NB(right_rcv_buffer, left_snd_buffer, size, receiver)
          !
          !left_write_flag = 0 

        ENDIF

        IF(chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) THEN
            receiver=chunks(chunks(chunk)%chunk_neighbours(chunk_right))%task
            WRITE(*,*) "Task: ", parallel%task, " sending right to: ", receiver, " with depth: ", depth

            IF (depth .EQ. 1) THEN

                WRITE(*,*) "Task: ", parallel%task, " send right to: ", receiver, " depth: ", depth, "stride: ", array_width, " num elements: ", num_elements, &
                                                    " source x: ", chunks(chunk)%field%x_max, &
                                                    " source y: ", chunks(chunk)%field%y_min-depth, &
                                                    " target x: ", chunks(chunk)%field%x_min-depth, &    
                                                    " target y: ", chunks(chunk)%field%y_min-depth 

                CALL SHMEM_IPUT64(field(chunks(chunk)%field%x_min-depth,chunks(chunk)%field%y_min-depth), &
                                  field(chunks(chunk)%field%x_max,chunks(chunk)%field%y_min-depth), &
                                  array_width, array_width,                                  &
                                  num_elements, receiver)

            ELSE
                CALL SHMEM_IPUT64(field(chunks(chunk)%field%x_min-1,chunks(chunk)%field%y_min-depth), &
                                     field(chunks(chunk)%field%x_max,chunks(chunk)%field%y_min-depth), &
                                     array_width, array_width,                                  &
                                     num_elements, receiver)
                WRITE(*,*) "Task: ", parallel%task, " first send to right to: ", receiver, " stride: ", array_width, " depth: ", depth, &
                                                                " source x: ", chunks(chunk)%field%x_max, &
                                                                " source y: ", chunks(chunk)%field%y_min-depth, &
                                                                " target x: ", chunks(chunk)%field%x_min-1, &
                                                                " target y: ", chunks(chunk)%field%y_min-depth 

                CALL SHMEM_IPUT64(field(chunks(chunk)%field%x_min-2,chunks(chunk)%field%y_min-depth), &
                                  field(chunks(chunk)%field%x_max-1,chunks(chunk)%field%y_min-depth), &
                                  array_width, array_width,                                  &
                                  num_elements, receiver)
                WRITE(*,*) "Task: ", parallel%task, " second send to right to: ", receiver, " stride: ", array_width, " depth: ", depth, &
                                                                " source x: ", chunks(chunk)%field%x_max-1, &
                                                                " source y: ", chunks(chunk)%field%y_min-depth, &
                                                                " target x: ", chunks(chunk)%field%x_min-2, &
                                                                " target y: ", chunks(chunk)%field%y_min-depth 
            ENDIF



          !IF (right_write_flag .EQ. 0) THEN
          !  CALL SHMEM_INT4_WAIT_UNTIL(right_write_flag, SHMEM_CMP_EQ, 1)
          !ENDIF
          !CALL SHMEM_PUT64_NB(left_rcv_buffer, right_snd_buffer, size, receiver)

          !right_write_flag = 0 

        ENDIF
  ENDIF

  !CALL SHMEM_QUIET
    CALL SHMEM_BARRIER_ALL

  !IF(parallel%task.EQ.chunks(chunk)%task) THEN

  !  ! Send/receive the data
  !  IF(chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) THEN
  !    receiver=chunks(chunks(chunk)%chunk_neighbours(chunk_left))%task

  !    CALL SHMEM_PUT4_NB(right_rcv_flag, 1, 1, receiver)

  !  ENDIF

  !  IF(chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) THEN
  !    receiver=chunks(chunks(chunk)%chunk_neighbours(chunk_right))%task

  !    CALL SHMEM_PUT4_NB(left_rcv_flag, 1, 1, receiver)

  !  ENDIF
  !ENDIF

  !IF(parallel%task.EQ.chunks(chunk)%task) THEN

  !  ! Send/receive the data
  !  IF(chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) THEN

  !    IF (left_rcv_flag .EQ. 0) THEN
  !      CALL SHMEM_INT4_WAIT_UNTIL(left_rcv_flag, SHMEM_CMP_EQ, 1)
  !    ENDIF
  !    left_rcv_flag = 0

  !  ENDIF

  !  IF(chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) THEN

  !    IF (right_rcv_flag .EQ. 0) THEN
  !      CALL SHMEM_INT4_WAIT_UNTIL(right_rcv_flag, SHMEM_CMP_EQ, 1)
  !    ENDIF
  !    right_rcv_flag = 0

  !  ENDIF
  !ENDIF


  ! Unpack buffers in halo cells
  !IF(parallel%task.EQ.chunks(chunk)%task) THEN
  !  IF(use_fortran_kernels) THEN
  !    CALL unpack_left_right_buffers(chunks(chunk)%field%x_min,chunks(chunk)%field%x_max, &
  !                                   chunks(chunk)%field%y_min,chunks(chunk)%field%y_max, &
  !                                   chunks(chunk)%chunk_neighbours(chunk_left),          &
  !                                   chunks(chunk)%chunk_neighbours(chunk_right),         &
  !                                   external_face,                                       &
  !                                   x_inc,y_inc,depth,size,                              &
  !                                   field,left_rcv_buffer,right_rcv_buffer)
  !  ELSEIF(use_C_kernels)THEN
  !    CALL unpack_left_right_buffers_c(chunks(chunk)%field%x_min,chunks(chunk)%field%x_max, &
  !                                     chunks(chunk)%field%y_min,chunks(chunk)%field%y_max, &
  !                                     chunks(chunk)%chunk_neighbours(chunk_left),          &
  !                                     chunks(chunk)%chunk_neighbours(chunk_right),         &
  !                                     external_face,                                       &
  !                                     x_inc,y_inc,depth,size,                              &
  !                                     field,left_rcv_buffer,right_rcv_buffer)
  !  ENDIF

  !  
  !  IF(chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) THEN
  !    receiver=chunks(chunks(chunk)%chunk_neighbours(chunk_left))%task

  !    CALL SHMEM_PUT4_NB(right_write_flag, 1, 1, receiver)
  !  ENDIF

  !  IF(chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) THEN
  !    receiver=chunks(chunks(chunk)%chunk_neighbours(chunk_right))%task

  !    CALL SHMEM_PUT4_NB(left_write_flag, 1, 1, receiver)
  !  ENDIF

  !ENDIF


  ! Pack real data into buffers
  IF(parallel%task.EQ.chunks(chunk)%task) THEN
    size=(1+(chunks(chunk)%field%x_max+x_inc+depth)-(chunks(chunk)%field%x_min-depth))*depth

    ! Send/receive the data
    IF(chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) THEN
        receiver=chunks(chunks(chunk)%chunk_neighbours(chunk_bottom))%task

        WRITE(*,*) "Task: ", parallel%task, " sending down to: ", receiver, " size: ", size, &
                                    " source x: ", chunks(chunk)%field%x_min-depth,  &
                                    " source y: ", chunks(chunk)%field%y_min+y_inc, &
                                    " target x: ", chunks(chunk)%field%x_min-depth, &
                                    " target y: ", chunks(chunk)%field%y_max+y_inc+1 

        CALL SHMEM_PUT64_NB(field(chunks(chunk)%field%x_min-depth,chunks(chunk)%field%y_max+y_inc+1), &
                            field(chunks(chunk)%field%x_min-depth,chunks(chunk)%field%y_min+y_inc), &
                            size, receiver)

      !IF (bottom_write_flag .EQ. 0) THEN
      !  CALL SHMEM_INT4_WAIT_UNTIL(bottom_write_flag, SHMEM_CMP_EQ, 1)
      !ENDIF
      !CALL SHMEM_PUT64_NB(top_rcv_buffer, bottom_snd_buffer, size, receiver)

      !bottom_write_flag = 0
    ENDIF

    IF(chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) THEN
        receiver=chunks(chunks(chunk)%chunk_neighbours(chunk_top))%task

        WRITE(*,*) "Task: ", parallel%task, " sending up to: ", receiver, " size: ", size, &
                                    " source x: ", chunks(chunk)%field%x_min-depth, &
                                    " source y: ", chunks(chunk)%field%y_max+1-depth, &
                                    " target x: ", chunks(chunk)%field%x_min-depth, &
                                    " target y: ", chunks(chunk)%field%y_min-depth
                                    

        CALL SHMEM_PUT64_NB(field(chunks(chunk)%field%x_min-depth,chunks(chunk)%field%y_min-depth), &
                            field(chunks(chunk)%field%x_min-depth,chunks(chunk)%field%y_max+1-depth), &
                            size, receiver)


      !IF (top_write_flag .EQ. 0) THEN
      !  CALL SHMEM_INT4_WAIT_UNTIL(top_write_flag, SHMEM_CMP_EQ, 1)
      !ENDIF

      !CALL SHMEM_PUT64_NB(bottom_rcv_buffer, top_snd_buffer, size, receiver)

      !top_write_flag = 0 
    ENDIF

  ENDIF

  ! Wait for the messages
  !CALL SHMEM_QUIET
    CALL SHMEM_BARRIER_ALL


  !IF(parallel%task.EQ.chunks(chunk)%task) THEN

  !  ! Send/receive the data
  !  IF(chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) THEN
  !    receiver=chunks(chunks(chunk)%chunk_neighbours(chunk_bottom))%task

  !    CALL SHMEM_PUT4_NB(top_rcv_flag, 1, 1, receiver)
  !  ENDIF

  !  IF(chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) THEN
  !    receiver=chunks(chunks(chunk)%chunk_neighbours(chunk_top))%task

  !    CALL SHMEM_PUT4_NB(bottom_rcv_flag, 1, 1, receiver)
  !  ENDIF
  !ENDIF

  !IF(parallel%task.EQ.chunks(chunk)%task) THEN

  !  ! Send/receive the data
  !  IF(chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) THEN

  !    IF (bottom_rcv_flag .EQ. 0) THEN
  !      CALL SHMEM_INT4_WAIT_UNTIL(bottom_rcv_flag, SHMEM_CMP_EQ, 1)
  !    ENDIF

  !    bottom_rcv_flag = 0
  !  ENDIF

  !  IF(chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) THEN

  !    IF (top_rcv_flag .EQ. 0) THEN
  !      CALL SHMEM_INT4_WAIT_UNTIL(top_rcv_flag, SHMEM_CMP_EQ, 1)
  !    ENDIF

  !    top_rcv_flag = 0
  !  ENDIF
  !ENDIF


  ! Unpack buffers in halo cells
  !IF(parallel%task.EQ.chunks(chunk)%task) THEN
  !  IF(use_fortran_kernels) THEN
  !    CALL unpack_top_bottom_buffers(chunks(chunk)%field%x_min,chunks(chunk)%field%x_max, &
  !                                   chunks(chunk)%field%y_min,chunks(chunk)%field%y_max, &
  !                                   chunks(chunk)%chunk_neighbours(chunk_bottom),        &
  !                                   chunks(chunk)%chunk_neighbours(chunk_top),           &
  !                                   external_face,                                       &
  !                                   x_inc,y_inc,depth,size,                              &
  !                                   field,bottom_rcv_buffer,top_rcv_buffer)
  !  ELSEIF(use_C_kernels)THEN
  !    CALL unpack_top_bottom_buffers_c(chunks(chunk)%field%x_min,chunks(chunk)%field%x_max, &
  !                                     chunks(chunk)%field%y_min,chunks(chunk)%field%y_max, &
  !                                     chunks(chunk)%chunk_neighbours(chunk_bottom),        &
  !                                     chunks(chunk)%chunk_neighbours(chunk_top),           &
  !                                     external_face,                                       &
  !                                     x_inc,y_inc,depth,size,                              &
  !                                     field,bottom_rcv_buffer,top_rcv_buffer)
  !  ENDIF

  !  IF(chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) THEN
  !    receiver=chunks(chunks(chunk)%chunk_neighbours(chunk_bottom))%task

  !    CALL SHMEM_PUT4_NB(top_write_flag, 1, 1, receiver)
  !  ENDIF

  !  IF(chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) THEN
  !    receiver=chunks(chunks(chunk)%chunk_neighbours(chunk_top))%task

  !    CALL SHMEM_PUT4_NB(bottom_write_flag, 1, 1, receiver)
  !  ENDIF
  !ENDIF

END SUBROUTINE clover_exchange_message

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
