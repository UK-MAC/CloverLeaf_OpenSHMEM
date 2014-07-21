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

    INTEGER(KIND=4), VOLATILE :: right_pe_ready, left_pe_ready, right_pe_written, left_pe_written
    INTEGER(KIND=4), VOLATILE :: top_pe_ready, bottom_pe_ready, top_pe_written, bottom_pe_written 
    COMMON /ARRAY_FLAG/ right_pe_ready, left_pe_ready, top_pe_ready, bottom_pe_ready, & 
                        right_pe_written, left_pe_written, top_pe_written, bottom_pe_written 
    ! 1 = false 0 = true

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
    INTEGER :: SHMEM_MY_PE, SHMEM_N_PES

    right_pe_ready = 1
    left_pe_ready = 1
    top_pe_ready = 1
    bottom_pe_ready = 1

    right_pe_written = 1
    left_pe_written = 1
    top_pe_written = 1
    bottom_pe_written = 1

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

                IF(cx.EQ.1)       chunks(1)%chunk_neighbours(chunk_left)=external_face
                IF(cx.EQ.chunk_x) chunks(1)%chunk_neighbours(chunk_right)=external_face
                IF(cy.EQ.1)       chunks(1)%chunk_neighbours(chunk_bottom)=external_face
                IF(cy.EQ.chunk_y) chunks(1)%chunk_neighbours(chunk_top)=external_face
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

SUBROUTINE clover_exchange(fields,depth)

  IMPLICIT NONE

  INTEGER      :: fields(:),depth,location_of_tasks_chunk

  ! Assuming 1 patch per task, this will be changed
  ! Also, not packing all fields for each communication, doing one at a time

    location_of_tasks_chunk = 1

  IF(fields(FIELD_DENSITY0).EQ.1) THEN
    CALL clover_exchange_message(location_of_tasks_chunk,density0,      &
                                 depth,CELL_DATA)
  ENDIF

  IF(fields(FIELD_DENSITY1).EQ.1) THEN
    CALL clover_exchange_message(location_of_tasks_chunk,density1,      &
                                 depth,CELL_DATA)
  ENDIF

  IF(fields(FIELD_ENERGY0).EQ.1) THEN
    CALL clover_exchange_message(location_of_tasks_chunk,energy0,       &
                                 depth,CELL_DATA)
  ENDIF

  IF(fields(FIELD_ENERGY1).EQ.1) THEN
    CALL clover_exchange_message(location_of_tasks_chunk,energy1,       &
                                 depth,CELL_DATA)
  ENDIF

  IF(fields(FIELD_PRESSURE).EQ.1) THEN
    CALL clover_exchange_message(location_of_tasks_chunk,pressure,      &
                                 depth,CELL_DATA)
  ENDIF

  IF(fields(FIELD_VISCOSITY).EQ.1) THEN
    CALL clover_exchange_message(location_of_tasks_chunk,viscosity,     &
                                 depth,CELL_DATA)
  ENDIF

  IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
    CALL clover_exchange_message(location_of_tasks_chunk,soundspeed,    &
                                 depth,CELL_DATA)
  ENDIF

  IF(fields(FIELD_XVEL0).EQ.1) THEN
    CALL clover_exchange_message(location_of_tasks_chunk,xvel0,         &
                                 depth,VERTEX_DATA)
  ENDIF

  IF(fields(FIELD_XVEL1).EQ.1) THEN
    CALL clover_exchange_message(location_of_tasks_chunk,xvel1,         &
                                 depth,VERTEX_DATA)
  ENDIF

  IF(fields(FIELD_YVEL0).EQ.1) THEN
    CALL clover_exchange_message(location_of_tasks_chunk,yvel0,         &
                                 depth,VERTEX_DATA)
  ENDIF

  IF(fields(FIELD_YVEL1).EQ.1) THEN
    CALL clover_exchange_message(location_of_tasks_chunk,yvel1,         &
                                 depth,VERTEX_DATA)
  ENDIF

  IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
    CALL clover_exchange_message(location_of_tasks_chunk,vol_flux_x,    &
                                 depth,X_FACE_DATA)
  ENDIF

  IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
    CALL clover_exchange_message(location_of_tasks_chunk,vol_flux_y,    &
                                 depth,Y_FACE_DATA)
  ENDIF

  IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
    CALL clover_exchange_message(location_of_tasks_chunk,mass_flux_x,   &
                                 depth,X_FACE_DATA)
  ENDIF

  IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
    CALL clover_exchange_message(location_of_tasks_chunk,mass_flux_y,   &
                                 depth,Y_FACE_DATA)
  ENDIF


END SUBROUTINE clover_exchange
 
SUBROUTINE clover_exchange_left_put(depth, x_inc, y_inc, array_width, num_elements, field, chunk)

    IMPLICIT NONE

    INTEGER :: depth, x_inc, y_inc, array_width, num_elements, receiver, chunk
    REAL(KIND=8) :: field(-1:,-1:) ! This seems to work for any type of mesh data

    receiver=chunks(chunk)%chunk_neighbours(chunk_left) - 1

    IF (depth .EQ. 1) THEN

        !WRITE(*,*) "Task : ", parallel%task, " send to left to: ", receiver, " stride: ", array_width, " depth: ", depth, " num elements: ", num_elements, &
        !                                     " source x: ", chunks(chunk)%field%x_min+x_inc, &
        !                                     " source y: ", chunks(chunk)%field%y_min-depth, &
        !                                     " target x: ", chunks(chunk)%field%x_max+x_inc+1, &
        !                                     " target y: ", chunks(chunk)%field%y_min-depth 


        CALL SHMEM_IPUT64(field(chunks(chunk)%field%x_max+x_inc+1,chunks(chunk)%field%y_min-depth), &
                          field(chunks(chunk)%field%x_min+x_inc,chunks(chunk)%field%y_min-depth), &
                          array_width, array_width,                                  &
                          num_elements, receiver)

    ELSE

        CALL SHMEM_IPUT64(field(chunks(chunk)%field%x_max+x_inc+1,chunks(chunk)%field%y_min-depth), &
                          field(chunks(chunk)%field%x_min+x_inc,chunks(chunk)%field%y_min-depth), &
                          array_width, array_width,                                  &
                          num_elements, receiver)
        !WRITE(*,*) "Task: ", parallel%task, " first send to left to: ", receiver,  " stride: ", array_width, " depth: ", depth,  &
        !                                                       " source x: ", chunks(chunk)%field%x_min+x_inc,  &
        !                                                       " source y: ", chunks(chunk)%field%y_min-depth,  &
        !                                                       " target x: ", chunks(chunk)%field%x_max+x_inc+1,& 
        !                                                       " target y: ", chunks(chunk)%field%y_min-depth  

        CALL SHMEM_IPUT64(field(chunks(chunk)%field%x_max+x_inc+2,chunks(chunk)%field%y_min-depth), &
                          field(chunks(chunk)%field%x_min+x_inc+1,chunks(chunk)%field%y_min-depth), &
                          array_width, array_width,                                  &
                          num_elements, receiver)

        !WRITE(*,*) "Task: ", parallel%task, " second sec to left to: ", receiver, " stride: ", array_width, " depth: ", depth, &
        !                                                       " source x: ", chunks(chunk)%field%x_min+x_inc+1, &
        !                                                       " source y: ", chunks(chunk)%field%y_min-depth, &
        !                                                       " target x: ", chunks(chunk)%field%x_max+x_inc+2, &
        !                                                       " target y: ", chunks(chunk)%field%y_min-depth
    ENDIF

END SUBROUTINE clover_exchange_left_put

SUBROUTINE clover_exchange_right_put(depth, x_inc, y_inc, array_width, num_elements, field, chunk)

    IMPLICIT NONE

    INTEGER :: depth, x_inc, y_inc, array_width, num_elements, receiver, chunk
    REAL(KIND=8) :: field(-1:,-1:) ! This seems to work for any type of mesh data

    receiver=chunks(chunk)%chunk_neighbours(chunk_right) - 1

    IF (depth .EQ. 1) THEN

        !WRITE(*,*) "Task: ", parallel%task, " send right to: ", receiver, " depth: ", depth, "stride: ", array_width, " num elements: ", num_elements, &
        !                                    " source x: ", chunks(chunk)%field%x_max, &
        !                                    " source y: ", chunks(chunk)%field%y_min-depth, &
        !                                    " target x: ", chunks(chunk)%field%x_min-depth, &    
        !                                    " target y: ", chunks(chunk)%field%y_min-depth 

        CALL SHMEM_IPUT64(field(chunks(chunk)%field%x_min-depth,chunks(chunk)%field%y_min-depth), &
                          field(chunks(chunk)%field%x_max,chunks(chunk)%field%y_min-depth), &
                          array_width, array_width,                                  &
                          num_elements, receiver)

    ELSE
        CALL SHMEM_IPUT64(field(chunks(chunk)%field%x_min-1,chunks(chunk)%field%y_min-depth), &
                             field(chunks(chunk)%field%x_max,chunks(chunk)%field%y_min-depth), &
                             array_width, array_width,                                  &
                             num_elements, receiver)
        !WRITE(*,*) "Task: ", parallel%task, " first send to right to: ", receiver, " stride: ", array_width, " depth: ", depth, &
        !                                                " source x: ", chunks(chunk)%field%x_max, &
        !                                                " source y: ", chunks(chunk)%field%y_min-depth, &
        !                                                " target x: ", chunks(chunk)%field%x_min-1, &
        !                                                " target y: ", chunks(chunk)%field%y_min-depth 

        CALL SHMEM_IPUT64(field(chunks(chunk)%field%x_min-2,chunks(chunk)%field%y_min-depth), &
                          field(chunks(chunk)%field%x_max-1,chunks(chunk)%field%y_min-depth), &
                          array_width, array_width,                                  &
                          num_elements, receiver)
        !WRITE(*,*) "Task: ", parallel%task, " second send to right to: ", receiver, " stride: ", array_width, " depth: ", depth, &
        !                                                " source x: ", chunks(chunk)%field%x_max-1, &
        !                                                " source y: ", chunks(chunk)%field%y_min-depth, &
        !                                                " target x: ", chunks(chunk)%field%x_min-2, &
        !                                                " target y: ", chunks(chunk)%field%y_min-depth 
    ENDIF

END SUBROUTINE clover_exchange_right_put


SUBROUTINE clover_exchange_top_put(depth, x_inc, y_inc, size, field, chunk)

    IMPLICIT NONE

    INTEGER :: depth, x_inc, y_inc, chunk, size, receiver 
    REAL(KIND=8) :: field(-1:,-1:) ! This seems to work for any type of mesh data

    receiver=chunks(chunk)%chunk_neighbours(chunk_top) - 1

#ifdef CRAY_NONBLOCK
    CALL SHMEM_PUT64_NB( &
#else
    CALL SHMEM_PUT64( &
#endif
                     field(chunks(chunk)%field%x_min-depth,chunks(chunk)%field%y_min-depth), &
                     field(chunks(chunk)%field%x_min-depth,chunks(chunk)%field%y_max+1-depth), &
                     size, receiver)

END SUBROUTINE clover_exchange_top_put

SUBROUTINE clover_exchange_bottom_put(depth, x_inc, y_inc, size, field, chunk)

    IMPLICIT NONE

    INTEGER :: depth, x_inc, y_inc, chunk, size, receiver 
    REAL(KIND=8) :: field(-1:,-1:) ! This seems to work for any type of mesh data

    receiver=chunks(chunk)%chunk_neighbours(chunk_bottom) - 1

#ifdef CRAY_NONBLOCK
    CALL SHMEM_PUT64_NB( &
#else
    CALL SHMEM_PUT64( &
#endif
                     field(chunks(chunk)%field%x_min-depth,chunks(chunk)%field%y_max+y_inc+1), &
                     field(chunks(chunk)%field%x_min-depth,chunks(chunk)%field%y_min+y_inc), &
                     size, receiver)

END SUBROUTINE clover_exchange_bottom_put

SUBROUTINE clover_exchange_message(chunk,field,depth,field_type)
    
    IMPLICIT NONE

    REAL(KIND=8) :: field(-1:,-1:) ! This seems to work for any type of mesh data

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


    ! Pack real data into buffers
    IF(parallel%task.EQ.chunks(chunk)%task) THEN

        !tell left and right neighbours that this pe is ready to receive data
        IF(chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) THEN
            receiver=chunks(chunk)%chunk_neighbours(chunk_left) - 1

#ifdef CRAY_NONBLOCK
            CALL SHMEM_PUT4_NB(right_pe_ready, 0, 1, receiver) 
#else
            CALL SHMEM_PUT4(right_pe_ready, 0, 1, receiver) 
#endif
        ENDIF

        IF(chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) THEN
            receiver=chunks(chunk)%chunk_neighbours(chunk_right) - 1

#ifdef CRAY_NONBLOCK
            CALL SHMEM_PUT4_NB(left_pe_ready, 0, 1, receiver)
#else
            CALL SHMEM_PUT4(left_pe_ready, 0, 1, receiver)
#endif
        ENDIF


        num_elements = 1+(chunks(chunk)%field%y_max+y_inc+depth)-(chunks(chunk)%field%y_min-depth)
        !WRITE(*,*) "Task: ", parallel%task, " num elements comms: ", num_elements

        array_width = chunks(chunk)%field%x_max+4+x_inc
        !WRITE(*,*) "Task: ", parallel%task, " array width: ", array_width

        IF ( (chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) ) THEN

            DO WHILE ((left_pe_ready .EQ. 1) .AND. (right_pe_ready .EQ. 1)) 
            ENDDO
            IF (left_pe_ready .EQ. 0) THEN
                CALL clover_exchange_left_put(depth, x_inc, y_inc, array_width, num_elements, field, chunk)
                left_pe_ready = 1
                DO WHILE (right_pe_ready .EQ. 1) 
                ENDDO
                CALL clover_exchange_right_put(depth, x_inc, y_inc, array_width, num_elements, field, chunk)    
                right_pe_ready = 1
            ELSE
                CALL clover_exchange_right_put(depth, x_inc, y_inc, array_width, num_elements, field, chunk)
                right_pe_ready = 1
                DO WHILE (left_pe_ready .EQ. 1)
                ENDDO
                CALL clover_exchange_left_put(depth, x_inc, y_inc, array_width, num_elements, field, chunk)
                left_pe_ready = 1
            ENDIF

        ELSEIF (chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) THEN

            DO WHILE(left_pe_ready .EQ. 1) 
            ENDDO
            CALL clover_exchange_left_put(depth, x_inc, y_inc, array_width, num_elements, field, chunk)
            left_pe_ready = 1

        ELSEIF (chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) THEN

            DO WHILE(right_pe_ready .EQ. 1) 
            ENDDO
            CALL clover_exchange_right_put(depth, x_inc, y_inc, array_width, num_elements, field, chunk)
            right_pe_ready = 1

        ENDIF
    


#ifdef FENCE_NOT_QUIET
        CALL SHMEM_FENCE
#else
        CALL SHMEM_QUIET
#endif


        IF(chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) THEN
            receiver=chunks(chunk)%chunk_neighbours(chunk_left) - 1

#ifdef CRAY_NONBLOCK
            CALL SHMEM_PUT4_NB(right_pe_written, 0, 1, receiver)
#else
            CALL SHMEM_PUT4(right_pe_written, 0, 1, receiver)
#endif
        ENDIF

        IF(chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) THEN
            receiver=chunks(chunk)%chunk_neighbours(chunk_right) - 1

#ifdef CRAY_NONBLOCK
            CALL SHMEM_PUT4_NB(left_pe_written, 0, 1, receiver)
#else
            CALL SHMEM_PUT4(left_pe_written, 0, 1, receiver)
#endif
        ENDIF


        IF ((chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face)) THEN 
            DO WHILE ((left_pe_written .EQ. 1) .OR. (right_pe_written .EQ. 1))
            ENDDO
            left_pe_written = 1
            right_pe_written = 1
        ELSEIF (chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) THEN
            DO WHILE (left_pe_written .EQ. 1) 
            ENDDO
            left_pe_written = 1
        ELSEIF (chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) THEN
            DO WHILE (right_pe_written .EQ. 1) 
            ENDDO
            right_pe_written = 1  
        ENDIF



        !tell up and down neighbours that this pe is ready to receive 
        IF(chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) THEN
            receiver=chunks(chunk)%chunk_neighbours(chunk_bottom) - 1

#ifdef CRAY_NONBLOCK
            CALL SHMEM_PUT4_NB(top_pe_ready, 0, 1, receiver)
#else
            CALL SHMEM_PUT4(top_pe_ready, 0, 1, receiver)
#endif
        ENDIF
    
        IF(chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) THEN
            receiver=chunks(chunk)%chunk_neighbours(chunk_top) - 1

#ifdef CRAY_NONBLOCK
            CALL SHMEM_PUT4_NB(bottom_pe_ready, 0, 1, receiver)
#else
            CALL SHMEM_PUT4(bottom_pe_ready, 0, 1, receiver)
#endif
        ENDIF



        size=(1+(chunks(chunk)%field%x_max+x_inc+depth)-(chunks(chunk)%field%x_min-depth))*depth

        IF ((chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face)) THEN

            DO WHILE ((bottom_pe_ready .EQ. 1) .AND. (top_pe_ready .EQ. 1)) 
            ENDDO
            IF (bottom_pe_ready .EQ. 0) THEN
                CALL clover_exchange_bottom_put(depth, x_inc, y_inc, size, field, chunk)
                bottom_pe_ready = 1
                DO WHILE (top_pe_ready .EQ. 1)
                ENDDO
                CALL clover_exchange_top_put(depth, x_inc, y_inc, size, field, chunk)
                top_pe_ready = 1
            ELSEIF (top_pe_ready .EQ. 0) THEN
                CALL clover_exchange_top_put(depth, x_inc, y_inc, size, field, chunk)
                top_pe_ready = 1
                DO WHILE (bottom_pe_ready .EQ. 1)
                ENDDO 
                CALL clover_exchange_bottom_put(depth, x_inc, y_inc, size, field, chunk)
                bottom_pe_ready = 1
            ENDIF
        
        ELSEIF (chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) THEN
            DO WHILE (bottom_pe_ready .EQ. 1)
            ENDDO
            CALL clover_exchange_bottom_put(depth, x_inc, y_inc, size, field, chunk)
            bottom_pe_ready = 1
        ELSEIF (chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) THEN
            DO WHILE (top_pe_ready .EQ. 1)
            ENDDO
            CALL clover_exchange_top_put(depth, x_inc, y_inc, size, field, chunk)
            top_pe_ready = 1
        ENDIF


        ! Wait for the messages
#ifdef FENCE_NOT_QUIET
        CALL SHMEM_FENCE
#else
        CALL SHMEM_QUIET
#endif


        ! Send/receive the data
        IF(chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) THEN
            receiver=chunks(chunk)%chunk_neighbours(chunk_bottom) - 1

#ifdef CRAY_NONBLOCK
            CALL SHMEM_PUT4_NB(top_pe_written, 0, 1, receiver)
#else
            CALL SHMEM_PUT4(top_pe_written, 0, 1, receiver)
#endif
        ENDIF

        IF(chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) THEN
            receiver=chunks(chunk)%chunk_neighbours(chunk_top) - 1

#ifdef CRAY_NONBLOCK
            CALL SHMEM_PUT4_NB(bottom_pe_written, 0, 1, receiver)
#else
            CALL SHMEM_PUT4(bottom_pe_written, 0, 1, receiver)
#endif
        ENDIF



        IF ((chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face)) THEN
            DO WHILE ((bottom_pe_written .EQ. 1) .OR. (top_pe_written .EQ. 1))
            ENDDO
            bottom_pe_written = 1
            top_pe_written = 1
        ELSEIF (chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) THEN
            DO WHILE (bottom_pe_written .EQ. 1)
            ENDDO
            bottom_pe_written = 1
        ELSEIF (chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) THEN
            DO WHILE (top_pe_written .EQ. 1)
            ENDDO
            top_pe_written = 1
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
