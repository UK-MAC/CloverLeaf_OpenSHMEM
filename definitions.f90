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

!>  @brief Holds the high level Fortran data types
!>  @author Wayne Gaudin
!>  @details The high level data types used to store the mesh and field data
!>  are defined here.
!>
!>  Also the global variables used for defining the input and controlling the
!>  scheme are defined here.

MODULE definitions_module

   USE data_module
   
   IMPLICIT NONE

   TYPE state_type
      LOGICAL            :: defined

      REAL(KIND=8)       :: density          &
                           ,energy           &
                           ,xvel             &
                           ,yvel

      INTEGER            :: geometry

      REAL(KIND=8)       :: xmin               &
                           ,xmax               &
                           ,ymin               &
                           ,ymax               &
                           ,radius
   END TYPE state_type

   TYPE(state_type), ALLOCATABLE             :: states(:)
   INTEGER                                   :: number_of_states

   TYPE grid_type
     REAL(KIND=8)       :: xmin            &
                          ,ymin            &
                          ,xmax            &
                          ,ymax
                     
     INTEGER            :: x_cells              &
                          ,y_cells
   END TYPE grid_type

   INTEGER      :: step

   LOGICAL      :: advect_x

   INTEGER      :: error_condition

   INTEGER      :: test_problem
   LOGICAL      :: complete

   LOGICAL      :: use_fortran_kernels
   LOGICAL      :: use_C_kernels
   LOGICAL      :: use_OA_kernels

   LOGICAL      :: profiler_on ! Internal code profiler to make comparisons across systems easier

   TYPE profiler_type
     REAL(KIND=8)       :: timestep        &
                          ,acceleration    &
                          ,PdV             &
                          ,cell_advection  &
                          ,mom_advection   &
                          ,viscosity       &
                          ,ideal_gas       &
                          ,visit           &
                          ,summary         &
                          ,reset           &
                          ,revert          &
                          ,flux            &
                          ,halo_exchange
                     
   END TYPE profiler_type
   TYPE(profiler_type)  :: profiler

   REAL(KIND=8) :: end_time

   INTEGER      :: end_step

   REAL(KIND=8) :: dtold          &
                  ,dt             &
                  ,time           &
                  ,dtinit         &
                  ,dtmin          &
                  ,dtmax          &
                  ,dtrise         &
                  ,dtu_safe       &
                  ,dtv_safe       &
                  ,dtc_safe       &
                  ,dtdiv_safe     &
                  ,dtc            &
                  ,dtu            &
                  ,dtv            &
                  ,dtdiv

   INTEGER      :: visit_frequency   &
                  ,summary_frequency

   INTEGER         :: jdt,kdt

   TYPE field_type
     REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: density0,density1
     REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: energy0,energy1
     REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: pressure
     REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: viscosity
     REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: soundspeed
     REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: xvel0,xvel1
     REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: yvel0,yvel1
     REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: vol_flux_x,mass_flux_x
     REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: vol_flux_y,mass_flux_y
     REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: work_array1 !node_flux, stepbymass, volume_change, pre_vol
     REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: work_array2 !node_mass_post, post_vol
     REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: work_array3 !node_mass_pre,pre_mass
     REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: work_array4 !advec_vel, post_mass
     REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: work_array5 !mom_flux, advec_vol
     REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: work_array6 !pre_vol, post_ener
     REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: work_array7 !post_vol, ener_flux

     INTEGER         :: left            &
                       ,right           &
                       ,bottom          &
                       ,top             &
                       ,left_boundary   &
                       ,right_boundary  &
                       ,bottom_boundary &
                       ,top_boundary

     INTEGER         :: x_min  &
                       ,y_min  &
                       ,x_max  &
                       ,y_max

     REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: cellx    &
                                                 ,celly    &
                                                 ,vertexx  &
                                                 ,vertexy  &
                                                 ,celldx   &
                                                 ,celldy   &
                                                 ,vertexdx &
                                                 ,vertexdy

     REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: volume  &
                                                 ,xarea   &
                                                 ,yarea

   END TYPE field_type
   
  TYPE chunk_type

     INTEGER         :: task   !mpi task

     INTEGER         :: chunk_neighbours(8) ! Chunks, not tasks, so we can overload in the future

     ! Idealy, create an array to hold the buffers for each field so a commuincation only needs
     !  one send and one receive per face, rather than per field.
     ! If chunks are overloaded, i.e. more chunks than tasks, might need to pack for a task to task comm 
     !  rather than a chunk to chunk comm. See how performance is at high core counts before deciding

     !REAL(KIND=8),ALLOCATABLE:: left_rcv_buffer(:),right_rcv_buffer(:),bottom_rcv_buffer(:),top_rcv_buffer(:)
     !REAL(KIND=8),ALLOCATABLE:: left_snd_buffer(:),right_snd_buffer(:),bottom_snd_buffer(:),top_snd_buffer(:)

     TYPE(field_type):: field

  END TYPE chunk_type

  REAL(KIND=8) :: density0_left_snd_buffer(1)
  REAL(KIND=8) :: density0_right_snd_buffer(1)
  REAL(KIND=8) :: density0_bottom_snd_buffer(1)
  REAL(KIND=8) :: density0_top_snd_buffer(1)
  REAL(KIND=8) :: density0_left_rcv_buffer(1)
  REAL(KIND=8) :: density0_right_rcv_buffer(1)
  REAL(KIND=8) :: density0_bottom_rcv_buffer(1)
  REAL(KIND=8) :: density0_top_rcv_buffer(1)
  REAL(KIND=8) :: density0_left_top_snd_buffer(1)
  REAL(KIND=8) :: density0_right_top_snd_buffer(1)
  REAL(KIND=8) :: density0_left_bottom_snd_buffer(1)
  REAL(KIND=8) :: density0_right_bottom_snd_buffer(1)
  REAL(KIND=8) :: density0_left_top_rcv_buffer(1)
  REAL(KIND=8) :: density0_right_top_rcv_buffer(1)
  REAL(KIND=8) :: density0_left_bottom_rcv_buffer(1)
  REAL(KIND=8) :: density0_right_bottom_rcv_buffer(1)

  REAL(KIND=8) :: density1_left_snd_buffer(1)
  REAL(KIND=8) :: density1_right_snd_buffer(1)
  REAL(KIND=8) :: density1_bottom_snd_buffer(1)
  REAL(KIND=8) :: density1_top_snd_buffer(1)
  REAL(KIND=8) :: density1_left_rcv_buffer(1)
  REAL(KIND=8) :: density1_right_rcv_buffer(1)
  REAL(KIND=8) :: density1_bottom_rcv_buffer(1)
  REAL(KIND=8) :: density1_top_rcv_buffer(1)
  REAL(KIND=8) :: density1_left_top_snd_buffer(1)
  REAL(KIND=8) :: density1_right_top_snd_buffer(1)
  REAL(KIND=8) :: density1_left_bottom_snd_buffer(1)
  REAL(KIND=8) :: density1_right_bottom_snd_buffer(1)
  REAL(KIND=8) :: density1_left_top_rcv_buffer(1)
  REAL(KIND=8) :: density1_right_top_rcv_buffer(1)
  REAL(KIND=8) :: density1_left_bottom_rcv_buffer(1)
  REAL(KIND=8) :: density1_right_bottom_rcv_buffer(1)

  REAL(KIND=8) :: energy0_left_snd_buffer(1)
  REAL(KIND=8) :: energy0_right_snd_buffer(1)
  REAL(KIND=8) :: energy0_bottom_snd_buffer(1)
  REAL(KIND=8) :: energy0_top_snd_buffer(1)
  REAL(KIND=8) :: energy0_left_rcv_buffer(1)
  REAL(KIND=8) :: energy0_right_rcv_buffer(1)
  REAL(KIND=8) :: energy0_bottom_rcv_buffer(1)
  REAL(KIND=8) :: energy0_top_rcv_buffer(1)
  REAL(KIND=8) :: energy0_left_top_snd_buffer(1)
  REAL(KIND=8) :: energy0_right_top_snd_buffer(1)
  REAL(KIND=8) :: energy0_left_bottom_snd_buffer(1)
  REAL(KIND=8) :: energy0_right_bottom_snd_buffer(1)
  REAL(KIND=8) :: energy0_left_top_rcv_buffer(1)
  REAL(KIND=8) :: energy0_right_top_rcv_buffer(1)
  REAL(KIND=8) :: energy0_left_bottom_rcv_buffer(1)
  REAL(KIND=8) :: energy0_right_bottom_rcv_buffer(1)

  REAL(KIND=8) :: energy1_left_snd_buffer(1)
  REAL(KIND=8) :: energy1_right_snd_buffer(1)
  REAL(KIND=8) :: energy1_bottom_snd_buffer(1)
  REAL(KIND=8) :: energy1_top_snd_buffer(1)
  REAL(KIND=8) :: energy1_left_rcv_buffer(1)
  REAL(KIND=8) :: energy1_right_rcv_buffer(1)
  REAL(KIND=8) :: energy1_bottom_rcv_buffer(1)
  REAL(KIND=8) :: energy1_top_rcv_buffer(1)
  REAL(KIND=8) :: energy1_left_top_snd_buffer(1)
  REAL(KIND=8) :: energy1_right_top_snd_buffer(1)
  REAL(KIND=8) :: energy1_left_bottom_snd_buffer(1)
  REAL(KIND=8) :: energy1_right_bottom_snd_buffer(1)
  REAL(KIND=8) :: energy1_left_top_rcv_buffer(1)
  REAL(KIND=8) :: energy1_right_top_rcv_buffer(1)
  REAL(KIND=8) :: energy1_left_bottom_rcv_buffer(1)
  REAL(KIND=8) :: energy1_right_bottom_rcv_buffer(1)

  REAL(KIND=8) :: pressure_left_snd_buffer(1)
  REAL(KIND=8) :: pressure_right_snd_buffer(1)
  REAL(KIND=8) :: pressure_bottom_snd_buffer(1)
  REAL(KIND=8) :: pressure_top_snd_buffer(1)
  REAL(KIND=8) :: pressure_left_rcv_buffer(1)
  REAL(KIND=8) :: pressure_right_rcv_buffer(1)
  REAL(KIND=8) :: pressure_bottom_rcv_buffer(1)
  REAL(KIND=8) :: pressure_top_rcv_buffer(1)
  REAL(KIND=8) :: pressure_left_top_snd_buffer(1)
  REAL(KIND=8) :: pressure_right_top_snd_buffer(1)
  REAL(KIND=8) :: pressure_left_bottom_snd_buffer(1)
  REAL(KIND=8) :: pressure_right_bottom_snd_buffer(1)
  REAL(KIND=8) :: pressure_left_top_rcv_buffer(1)
  REAL(KIND=8) :: pressure_right_top_rcv_buffer(1)
  REAL(KIND=8) :: pressure_left_bottom_rcv_buffer(1)
  REAL(KIND=8) :: pressure_right_bottom_rcv_buffer(1)

  REAL(KIND=8) :: viscosity_left_snd_buffer(1)
  REAL(KIND=8) :: viscosity_right_snd_buffer(1)
  REAL(KIND=8) :: viscosity_bottom_snd_buffer(1)
  REAL(KIND=8) :: viscosity_top_snd_buffer(1)
  REAL(KIND=8) :: viscosity_left_rcv_buffer(1)
  REAL(KIND=8) :: viscosity_right_rcv_buffer(1)
  REAL(KIND=8) :: viscosity_bottom_rcv_buffer(1)
  REAL(KIND=8) :: viscosity_top_rcv_buffer(1)
  REAL(KIND=8) :: viscosity_left_top_snd_buffer(1)
  REAL(KIND=8) :: viscosity_right_top_snd_buffer(1)
  REAL(KIND=8) :: viscosity_left_bottom_snd_buffer(1)
  REAL(KIND=8) :: viscosity_right_bottom_snd_buffer(1)
  REAL(KIND=8) :: viscosity_left_top_rcv_buffer(1)
  REAL(KIND=8) :: viscosity_right_top_rcv_buffer(1)
  REAL(KIND=8) :: viscosity_left_bottom_rcv_buffer(1)
  REAL(KIND=8) :: viscosity_right_bottom_rcv_buffer(1)

  REAL(KIND=8) :: soundspeed_left_snd_buffer(1)
  REAL(KIND=8) :: soundspeed_right_snd_buffer(1)
  REAL(KIND=8) :: soundspeed_bottom_snd_buffer(1)
  REAL(KIND=8) :: soundspeed_top_snd_buffer(1)
  REAL(KIND=8) :: soundspeed_left_rcv_buffer(1)
  REAL(KIND=8) :: soundspeed_right_rcv_buffer(1)
  REAL(KIND=8) :: soundspeed_bottom_rcv_buffer(1)
  REAL(KIND=8) :: soundspeed_top_rcv_buffer(1)
  REAL(KIND=8) :: soundspeed_left_top_snd_buffer(1)
  REAL(KIND=8) :: soundspeed_right_top_snd_buffer(1)
  REAL(KIND=8) :: soundspeed_left_bottom_snd_buffer(1)
  REAL(KIND=8) :: soundspeed_right_bottom_snd_buffer(1)
  REAL(KIND=8) :: soundspeed_left_top_rcv_buffer(1)
  REAL(KIND=8) :: soundspeed_right_top_rcv_buffer(1)
  REAL(KIND=8) :: soundspeed_left_bottom_rcv_buffer(1)
  REAL(KIND=8) :: soundspeed_right_bottom_rcv_buffer(1)

  REAL(KIND=8) :: xvel0_left_snd_buffer(1)
  REAL(KIND=8) :: xvel0_right_snd_buffer(1)
  REAL(KIND=8) :: xvel0_bottom_snd_buffer(1)
  REAL(KIND=8) :: xvel0_top_snd_buffer(1)
  REAL(KIND=8) :: xvel0_left_rcv_buffer(1)
  REAL(KIND=8) :: xvel0_right_rcv_buffer(1)
  REAL(KIND=8) :: xvel0_bottom_rcv_buffer(1)
  REAL(KIND=8) :: xvel0_top_rcv_buffer(1)
  REAL(KIND=8) :: xvel0_left_top_snd_buffer(1)
  REAL(KIND=8) :: xvel0_right_top_snd_buffer(1)
  REAL(KIND=8) :: xvel0_left_bottom_snd_buffer(1)
  REAL(KIND=8) :: xvel0_right_bottom_snd_buffer(1)
  REAL(KIND=8) :: xvel0_left_top_rcv_buffer(1)
  REAL(KIND=8) :: xvel0_right_top_rcv_buffer(1)
  REAL(KIND=8) :: xvel0_left_bottom_rcv_buffer(1)
  REAL(KIND=8) :: xvel0_right_bottom_rcv_buffer(1)

  REAL(KIND=8) :: xvel1_left_snd_buffer(1)
  REAL(KIND=8) :: xvel1_right_snd_buffer(1)
  REAL(KIND=8) :: xvel1_bottom_snd_buffer(1)
  REAL(KIND=8) :: xvel1_top_snd_buffer(1)
  REAL(KIND=8) :: xvel1_left_rcv_buffer(1)
  REAL(KIND=8) :: xvel1_right_rcv_buffer(1)
  REAL(KIND=8) :: xvel1_bottom_rcv_buffer(1)
  REAL(KIND=8) :: xvel1_top_rcv_buffer(1)
  REAL(KIND=8) :: xvel1_left_top_snd_buffer(1)
  REAL(KIND=8) :: xvel1_right_top_snd_buffer(1)
  REAL(KIND=8) :: xvel1_left_bottom_snd_buffer(1)
  REAL(KIND=8) :: xvel1_right_bottom_snd_buffer(1)
  REAL(KIND=8) :: xvel1_left_top_rcv_buffer(1)
  REAL(KIND=8) :: xvel1_right_top_rcv_buffer(1)
  REAL(KIND=8) :: xvel1_left_bottom_rcv_buffer(1)
  REAL(KIND=8) :: xvel1_right_bottom_rcv_buffer(1)

  REAL(KIND=8) :: yvel0_left_snd_buffer(1)
  REAL(KIND=8) :: yvel0_right_snd_buffer(1)
  REAL(KIND=8) :: yvel0_bottom_snd_buffer(1)
  REAL(KIND=8) :: yvel0_top_snd_buffer(1)
  REAL(KIND=8) :: yvel0_left_rcv_buffer(1)
  REAL(KIND=8) :: yvel0_right_rcv_buffer(1)
  REAL(KIND=8) :: yvel0_bottom_rcv_buffer(1)
  REAL(KIND=8) :: yvel0_top_rcv_buffer(1)
  REAL(KIND=8) :: yvel0_left_top_snd_buffer(1)
  REAL(KIND=8) :: yvel0_right_top_snd_buffer(1)
  REAL(KIND=8) :: yvel0_left_bottom_snd_buffer(1)
  REAL(KIND=8) :: yvel0_right_bottom_snd_buffer(1)
  REAL(KIND=8) :: yvel0_left_top_rcv_buffer(1)
  REAL(KIND=8) :: yvel0_right_top_rcv_buffer(1)
  REAL(KIND=8) :: yvel0_left_bottom_rcv_buffer(1)
  REAL(KIND=8) :: yvel0_right_bottom_rcv_buffer(1)

  REAL(KIND=8) :: yvel1_left_snd_buffer(1)
  REAL(KIND=8) :: yvel1_right_snd_buffer(1)
  REAL(KIND=8) :: yvel1_bottom_snd_buffer(1)
  REAL(KIND=8) :: yvel1_top_snd_buffer(1)
  REAL(KIND=8) :: yvel1_left_rcv_buffer(1)
  REAL(KIND=8) :: yvel1_right_rcv_buffer(1)
  REAL(KIND=8) :: yvel1_bottom_rcv_buffer(1)
  REAL(KIND=8) :: yvel1_top_rcv_buffer(1)
  REAL(KIND=8) :: yvel1_left_top_snd_buffer(1)
  REAL(KIND=8) :: yvel1_right_top_snd_buffer(1)
  REAL(KIND=8) :: yvel1_left_bottom_snd_buffer(1)
  REAL(KIND=8) :: yvel1_right_bottom_snd_buffer(1)
  REAL(KIND=8) :: yvel1_left_top_rcv_buffer(1)
  REAL(KIND=8) :: yvel1_right_top_rcv_buffer(1)
  REAL(KIND=8) :: yvel1_left_bottom_rcv_buffer(1)
  REAL(KIND=8) :: yvel1_right_bottom_rcv_buffer(1)

  REAL(KIND=8) :: volflux_x_left_snd_buffer(1)
  REAL(KIND=8) :: volflux_x_right_snd_buffer(1)
  REAL(KIND=8) :: volflux_x_bottom_snd_buffer(1)
  REAL(KIND=8) :: volflux_x_top_snd_buffer(1)
  REAL(KIND=8) :: volflux_x_left_rcv_buffer(1)
  REAL(KIND=8) :: volflux_x_right_rcv_buffer(1)
  REAL(KIND=8) :: volflux_x_bottom_rcv_buffer(1)
  REAL(KIND=8) :: volflux_x_top_rcv_buffer(1)
  REAL(KIND=8) :: volflux_x_left_top_snd_buffer(1)
  REAL(KIND=8) :: volflux_x_right_top_snd_buffer(1)
  REAL(KIND=8) :: volflux_x_left_bottom_snd_buffer(1)
  REAL(KIND=8) :: volflux_x_right_bottom_snd_buffer(1)
  REAL(KIND=8) :: volflux_x_left_top_rcv_buffer(1)
  REAL(KIND=8) :: volflux_x_right_top_rcv_buffer(1)
  REAL(KIND=8) :: volflux_x_left_bottom_rcv_buffer(1)
  REAL(KIND=8) :: volflux_x_right_bottom_rcv_buffer(1)

  REAL(KIND=8) :: volflux_y_left_snd_buffer(1)
  REAL(KIND=8) :: volflux_y_right_snd_buffer(1)
  REAL(KIND=8) :: volflux_y_bottom_snd_buffer(1)
  REAL(KIND=8) :: volflux_y_top_snd_buffer(1)
  REAL(KIND=8) :: volflux_y_left_rcv_buffer(1)
  REAL(KIND=8) :: volflux_y_right_rcv_buffer(1)
  REAL(KIND=8) :: volflux_y_bottom_rcv_buffer(1)
  REAL(KIND=8) :: volflux_y_top_rcv_buffer(1)
  REAL(KIND=8) :: volflux_y_left_top_snd_buffer(1)
  REAL(KIND=8) :: volflux_y_right_top_snd_buffer(1)
  REAL(KIND=8) :: volflux_y_left_bottom_snd_buffer(1)
  REAL(KIND=8) :: volflux_y_right_bottom_snd_buffer(1)
  REAL(KIND=8) :: volflux_y_left_top_rcv_buffer(1)
  REAL(KIND=8) :: volflux_y_right_top_rcv_buffer(1)
  REAL(KIND=8) :: volflux_y_left_bottom_rcv_buffer(1)
  REAL(KIND=8) :: volflux_y_right_bottom_rcv_buffer(1)

  REAL(KIND=8) :: massflux_x_left_snd_buffer(1)
  REAL(KIND=8) :: massflux_x_right_snd_buffer(1)
  REAL(KIND=8) :: massflux_x_bottom_snd_buffer(1)
  REAL(KIND=8) :: massflux_x_top_snd_buffer(1)
  REAL(KIND=8) :: massflux_x_left_rcv_buffer(1)
  REAL(KIND=8) :: massflux_x_right_rcv_buffer(1)
  REAL(KIND=8) :: massflux_x_bottom_rcv_buffer(1)
  REAL(KIND=8) :: massflux_x_top_rcv_buffer(1)
  REAL(KIND=8) :: massflux_x_left_top_snd_buffer(1)
  REAL(KIND=8) :: massflux_x_right_top_snd_buffer(1)
  REAL(KIND=8) :: massflux_x_left_bottom_snd_buffer(1)
  REAL(KIND=8) :: massflux_x_right_bottom_snd_buffer(1)
  REAL(KIND=8) :: massflux_x_left_top_rcv_buffer(1)
  REAL(KIND=8) :: massflux_x_right_top_rcv_buffer(1)
  REAL(KIND=8) :: massflux_x_left_bottom_rcv_buffer(1)
  REAL(KIND=8) :: massflux_x_right_bottom_rcv_buffer(1)

  REAL(KIND=8) :: massflux_y_left_snd_buffer(1)
  REAL(KIND=8) :: massflux_y_right_snd_buffer(1)
  REAL(KIND=8) :: massflux_y_bottom_snd_buffer(1)
  REAL(KIND=8) :: massflux_y_top_snd_buffer(1)
  REAL(KIND=8) :: massflux_y_left_rcv_buffer(1)
  REAL(KIND=8) :: massflux_y_right_rcv_buffer(1)
  REAL(KIND=8) :: massflux_y_bottom_rcv_buffer(1)
  REAL(KIND=8) :: massflux_y_top_rcv_buffer(1)
  REAL(KIND=8) :: massflux_y_left_top_snd_buffer(1)
  REAL(KIND=8) :: massflux_y_right_top_snd_buffer(1)
  REAL(KIND=8) :: massflux_y_left_bottom_snd_buffer(1)
  REAL(KIND=8) :: massflux_y_right_bottom_snd_buffer(1)
  REAL(KIND=8) :: massflux_y_left_top_rcv_buffer(1)
  REAL(KIND=8) :: massflux_y_right_top_rcv_buffer(1)
  REAL(KIND=8) :: massflux_y_left_bottom_rcv_buffer(1)
  REAL(KIND=8) :: massflux_y_right_bottom_rcv_buffer(1)

  POINTER (density0_plr,  density0_left_rcv_buffer)
  POINTER (density0_prr,  density0_right_rcv_buffer)
  POINTER (density0_pbr,  density0_bottom_rcv_buffer)
  POINTER (density0_ptr,  density0_top_rcv_buffer)
  POINTER (density0_pls,  density0_left_snd_buffer)
  POINTER (density0_prs,  density0_right_snd_buffer)
  POINTER (density0_pbs,  density0_bottom_snd_buffer)
  POINTER (density0_pts,  density0_top_snd_buffer)
  POINTER (density0_pltr, density0_left_top_rcv_buffer)
  POINTER (density0_prtr, density0_right_top_rcv_buffer)
  POINTER (density0_plbr, density0_left_bottom_rcv_buffer)
  POINTER (density0_prbr, density0_right_bottom_rcv_buffer)
  POINTER (density0_plts, density0_left_top_snd_buffer)
  POINTER (density0_prts, density0_right_top_snd_buffer)
  POINTER (density0_plbs, density0_left_bottom_snd_buffer)
  POINTER (density0_prbs, density0_right_bottom_snd_buffer)

  POINTER (density1_plr,  density1_left_rcv_buffer)
  POINTER (density1_prr,  density1_right_rcv_buffer)
  POINTER (density1_pbr,  density1_bottom_rcv_buffer)
  POINTER (density1_ptr,  density1_top_rcv_buffer)
  POINTER (density1_pls,  density1_left_snd_buffer)
  POINTER (density1_prs,  density1_right_snd_buffer)
  POINTER (density1_pbs,  density1_bottom_snd_buffer)
  POINTER (density1_pts,  density1_top_snd_buffer)
  POINTER (density1_pltr, density1_left_top_rcv_buffer)
  POINTER (density1_prtr, density1_right_top_rcv_buffer)
  POINTER (density1_plbr, density1_left_bottom_rcv_buffer)
  POINTER (density1_prbr, density1_right_bottom_rcv_buffer)
  POINTER (density1_plts, density1_left_top_snd_buffer)
  POINTER (density1_prts, density1_right_top_snd_buffer)
  POINTER (density1_plbs, density1_left_bottom_snd_buffer)
  POINTER (density1_prbs, density1_right_bottom_snd_buffer)

  POINTER (energy0_plr,  energy0_left_rcv_buffer)
  POINTER (energy0_prr,  energy0_right_rcv_buffer)
  POINTER (energy0_pbr,  energy0_bottom_rcv_buffer)
  POINTER (energy0_ptr,  energy0_top_rcv_buffer)
  POINTER (energy0_pls,  energy0_left_snd_buffer)
  POINTER (energy0_prs,  energy0_right_snd_buffer)
  POINTER (energy0_pbs,  energy0_bottom_snd_buffer)
  POINTER (energy0_pts,  energy0_top_snd_buffer)
  POINTER (energy0_pltr, energy0_left_top_rcv_buffer)
  POINTER (energy0_prtr, energy0_right_top_rcv_buffer)
  POINTER (energy0_plbr, energy0_left_bottom_rcv_buffer)
  POINTER (energy0_prbr, energy0_right_bottom_rcv_buffer)
  POINTER (energy0_plts, energy0_left_top_snd_buffer)
  POINTER (energy0_prts, energy0_right_top_snd_buffer)
  POINTER (energy0_plbs, energy0_left_bottom_snd_buffer)
  POINTER (energy0_prbs, energy0_right_bottom_snd_buffer)

  POINTER (energy1_plr,  energy1_left_rcv_buffer)
  POINTER (energy1_prr,  energy1_right_rcv_buffer)
  POINTER (energy1_pbr,  energy1_bottom_rcv_buffer)
  POINTER (energy1_ptr,  energy1_top_rcv_buffer)
  POINTER (energy1_pls,  energy1_left_snd_buffer)
  POINTER (energy1_prs,  energy1_right_snd_buffer)
  POINTER (energy1_pbs,  energy1_bottom_snd_buffer)
  POINTER (energy1_pts,  energy1_top_snd_buffer)
  POINTER (energy1_pltr, energy1_left_top_rcv_buffer)
  POINTER (energy1_prtr, energy1_right_top_rcv_buffer)
  POINTER (energy1_plbr, energy1_left_bottom_rcv_buffer)
  POINTER (energy1_prbr, energy1_right_bottom_rcv_buffer)
  POINTER (energy1_plts, energy1_left_top_snd_buffer)
  POINTER (energy1_prts, energy1_right_top_snd_buffer)
  POINTER (energy1_plbs, energy1_left_bottom_snd_buffer)
  POINTER (energy1_prbs, energy1_right_bottom_snd_buffer)

  POINTER (pressure_plr,  pressure_left_rcv_buffer)
  POINTER (pressure_prr,  pressure_right_rcv_buffer)
  POINTER (pressure_pbr,  pressure_bottom_rcv_buffer)
  POINTER (pressure_ptr,  pressure_top_rcv_buffer)
  POINTER (pressure_pls,  pressure_left_snd_buffer)
  POINTER (pressure_prs,  pressure_right_snd_buffer)
  POINTER (pressure_pbs,  pressure_bottom_snd_buffer)
  POINTER (pressure_pts,  pressure_top_snd_buffer)
  POINTER (pressure_pltr, pressure_left_top_rcv_buffer)
  POINTER (pressure_prtr, pressure_right_top_rcv_buffer)
  POINTER (pressure_plbr, pressure_left_bottom_rcv_buffer)
  POINTER (pressure_prbr, pressure_right_bottom_rcv_buffer)
  POINTER (pressure_plts, pressure_left_top_snd_buffer)
  POINTER (pressure_prts, pressure_right_top_snd_buffer)
  POINTER (pressure_plbs, pressure_left_bottom_snd_buffer)
  POINTER (pressure_prbs, pressure_right_bottom_snd_buffer)

  POINTER (viscosity_plr,  viscosity_left_rcv_buffer)
  POINTER (viscosity_prr,  viscosity_right_rcv_buffer)
  POINTER (viscosity_pbr,  viscosity_bottom_rcv_buffer)
  POINTER (viscosity_ptr,  viscosity_top_rcv_buffer)
  POINTER (viscosity_pls,  viscosity_left_snd_buffer)
  POINTER (viscosity_prs,  viscosity_right_snd_buffer)
  POINTER (viscosity_pbs,  viscosity_bottom_snd_buffer)
  POINTER (viscosity_pts,  viscosity_top_snd_buffer)
  POINTER (viscosity_pltr, viscosity_left_top_rcv_buffer)
  POINTER (viscosity_prtr, viscosity_right_top_rcv_buffer)
  POINTER (viscosity_plbr, viscosity_left_bottom_rcv_buffer)
  POINTER (viscosity_prbr, viscosity_right_bottom_rcv_buffer)
  POINTER (viscosity_plts, viscosity_left_top_snd_buffer)
  POINTER (viscosity_prts, viscosity_right_top_snd_buffer)
  POINTER (viscosity_plbs, viscosity_left_bottom_snd_buffer)
  POINTER (viscosity_prbs, viscosity_right_bottom_snd_buffer)

  POINTER (soundspeed_plr,  soundspeed_left_rcv_buffer)
  POINTER (soundspeed_prr,  soundspeed_right_rcv_buffer)
  POINTER (soundspeed_pbr,  soundspeed_bottom_rcv_buffer)
  POINTER (soundspeed_ptr,  soundspeed_top_rcv_buffer)
  POINTER (soundspeed_pls,  soundspeed_left_snd_buffer)
  POINTER (soundspeed_prs,  soundspeed_right_snd_buffer)
  POINTER (soundspeed_pbs,  soundspeed_bottom_snd_buffer)
  POINTER (soundspeed_pts,  soundspeed_top_snd_buffer)
  POINTER (soundspeed_pltr, soundspeed_left_top_rcv_buffer)
  POINTER (soundspeed_prtr, soundspeed_right_top_rcv_buffer)
  POINTER (soundspeed_plbr, soundspeed_left_bottom_rcv_buffer)
  POINTER (soundspeed_prbr, soundspeed_right_bottom_rcv_buffer)
  POINTER (soundspeed_plts, soundspeed_left_top_snd_buffer)
  POINTER (soundspeed_prts, soundspeed_right_top_snd_buffer)
  POINTER (soundspeed_plbs, soundspeed_left_bottom_snd_buffer)
  POINTER (soundspeed_prbs, soundspeed_right_bottom_snd_buffer)

  POINTER (xvel0_plr,  xvel0_left_rcv_buffer)
  POINTER (xvel0_prr,  xvel0_right_rcv_buffer)
  POINTER (xvel0_pbr,  xvel0_bottom_rcv_buffer)
  POINTER (xvel0_ptr,  xvel0_top_rcv_buffer)
  POINTER (xvel0_pls,  xvel0_left_snd_buffer)
  POINTER (xvel0_prs,  xvel0_right_snd_buffer)
  POINTER (xvel0_pbs,  xvel0_bottom_snd_buffer)
  POINTER (xvel0_pts,  xvel0_top_snd_buffer)
  POINTER (xvel0_pltr, xvel0_left_top_rcv_buffer)
  POINTER (xvel0_prtr, xvel0_right_top_rcv_buffer)
  POINTER (xvel0_plbr, xvel0_left_bottom_rcv_buffer)
  POINTER (xvel0_prbr, xvel0_right_bottom_rcv_buffer)
  POINTER (xvel0_plts, xvel0_left_top_snd_buffer)
  POINTER (xvel0_prts, xvel0_right_top_snd_buffer)
  POINTER (xvel0_plbs, xvel0_left_bottom_snd_buffer)
  POINTER (xvel0_prbs, xvel0_right_bottom_snd_buffer)

  POINTER (xvel1_plr,  xvel1_left_rcv_buffer)
  POINTER (xvel1_prr,  xvel1_right_rcv_buffer)
  POINTER (xvel1_pbr,  xvel1_bottom_rcv_buffer)
  POINTER (xvel1_ptr,  xvel1_top_rcv_buffer)
  POINTER (xvel1_pls,  xvel1_left_snd_buffer)
  POINTER (xvel1_prs,  xvel1_right_snd_buffer)
  POINTER (xvel1_pbs,  xvel1_bottom_snd_buffer)
  POINTER (xvel1_pts,  xvel1_top_snd_buffer)
  POINTER (xvel1_pltr, xvel1_left_top_rcv_buffer)
  POINTER (xvel1_prtr, xvel1_right_top_rcv_buffer)
  POINTER (xvel1_plbr, xvel1_left_bottom_rcv_buffer)
  POINTER (xvel1_prbr, xvel1_right_bottom_rcv_buffer)
  POINTER (xvel1_plts, xvel1_left_top_snd_buffer)
  POINTER (xvel1_prts, xvel1_right_top_snd_buffer)
  POINTER (xvel1_plbs, xvel1_left_bottom_snd_buffer)
  POINTER (xvel1_prbs, xvel1_right_bottom_snd_buffer)

  POINTER (yvel0_plr,  yvel0_left_rcv_buffer)
  POINTER (yvel0_prr,  yvel0_right_rcv_buffer)
  POINTER (yvel0_pbr,  yvel0_bottom_rcv_buffer)
  POINTER (yvel0_ptr,  yvel0_top_rcv_buffer)
  POINTER (yvel0_pls,  yvel0_left_snd_buffer)
  POINTER (yvel0_prs,  yvel0_right_snd_buffer)
  POINTER (yvel0_pbs,  yvel0_bottom_snd_buffer)
  POINTER (yvel0_pts,  yvel0_top_snd_buffer)
  POINTER (yvel0_pltr, yvel0_left_top_rcv_buffer)
  POINTER (yvel0_prtr, yvel0_right_top_rcv_buffer)
  POINTER (yvel0_plbr, yvel0_left_bottom_rcv_buffer)
  POINTER (yvel0_prbr, yvel0_right_bottom_rcv_buffer)
  POINTER (yvel0_plts, yvel0_left_top_snd_buffer)
  POINTER (yvel0_prts, yvel0_right_top_snd_buffer)
  POINTER (yvel0_plbs, yvel0_left_bottom_snd_buffer)
  POINTER (yvel0_prbs, yvel0_right_bottom_snd_buffer)

  POINTER (yvel1_plr,  yvel1_left_rcv_buffer)
  POINTER (yvel1_prr,  yvel1_right_rcv_buffer)
  POINTER (yvel1_pbr,  yvel1_bottom_rcv_buffer)
  POINTER (yvel1_ptr,  yvel1_top_rcv_buffer)
  POINTER (yvel1_pls,  yvel1_left_snd_buffer)
  POINTER (yvel1_prs,  yvel1_right_snd_buffer)
  POINTER (yvel1_pbs,  yvel1_bottom_snd_buffer)
  POINTER (yvel1_pts,  yvel1_top_snd_buffer)
  POINTER (yvel1_pltr, yvel1_left_top_rcv_buffer)
  POINTER (yvel1_prtr, yvel1_right_top_rcv_buffer)
  POINTER (yvel1_plbr, yvel1_left_bottom_rcv_buffer)
  POINTER (yvel1_prbr, yvel1_right_bottom_rcv_buffer)
  POINTER (yvel1_plts, yvel1_left_top_snd_buffer)
  POINTER (yvel1_prts, yvel1_right_top_snd_buffer)
  POINTER (yvel1_plbs, yvel1_left_bottom_snd_buffer)
  POINTER (yvel1_prbs, yvel1_right_bottom_snd_buffer)

  POINTER (volflux_x_plr,  volflux_x_left_rcv_buffer)
  POINTER (volflux_x_prr,  volflux_x_right_rcv_buffer)
  POINTER (volflux_x_pbr,  volflux_x_bottom_rcv_buffer)
  POINTER (volflux_x_ptr,  volflux_x_top_rcv_buffer)
  POINTER (volflux_x_pls,  volflux_x_left_snd_buffer)
  POINTER (volflux_x_prs,  volflux_x_right_snd_buffer)
  POINTER (volflux_x_pbs,  volflux_x_bottom_snd_buffer)
  POINTER (volflux_x_pts,  volflux_x_top_snd_buffer)
  POINTER (volflux_x_pltr, volflux_x_left_top_rcv_buffer)
  POINTER (volflux_x_prtr, volflux_x_right_top_rcv_buffer)
  POINTER (volflux_x_plbr, volflux_x_left_bottom_rcv_buffer)
  POINTER (volflux_x_prbr, volflux_x_right_bottom_rcv_buffer)
  POINTER (volflux_x_plts, volflux_x_left_top_snd_buffer)
  POINTER (volflux_x_prts, volflux_x_right_top_snd_buffer)
  POINTER (volflux_x_plbs, volflux_x_left_bottom_snd_buffer)
  POINTER (volflux_x_prbs, volflux_x_right_bottom_snd_buffer)

  POINTER (volflux_y_plr,  volflux_y_left_rcv_buffer)
  POINTER (volflux_y_prr,  volflux_y_right_rcv_buffer)
  POINTER (volflux_y_pbr,  volflux_y_bottom_rcv_buffer)
  POINTER (volflux_y_ptr,  volflux_y_top_rcv_buffer)
  POINTER (volflux_y_pls,  volflux_y_left_snd_buffer)
  POINTER (volflux_y_prs,  volflux_y_right_snd_buffer)
  POINTER (volflux_y_pbs,  volflux_y_bottom_snd_buffer)
  POINTER (volflux_y_pts,  volflux_y_top_snd_buffer)
  POINTER (volflux_y_pltr, volflux_y_left_top_rcv_buffer)
  POINTER (volflux_y_prtr, volflux_y_right_top_rcv_buffer)
  POINTER (volflux_y_plbr, volflux_y_left_bottom_rcv_buffer)
  POINTER (volflux_y_prbr, volflux_y_right_bottom_rcv_buffer)
  POINTER (volflux_y_plts, volflux_y_left_top_snd_buffer)
  POINTER (volflux_y_prts, volflux_y_right_top_snd_buffer)
  POINTER (volflux_y_plbs, volflux_y_left_bottom_snd_buffer)
  POINTER (volflux_y_prbs, volflux_y_right_bottom_snd_buffer)

  POINTER (massflux_x_plr,  massflux_x_left_rcv_buffer)
  POINTER (massflux_x_prr,  massflux_x_right_rcv_buffer)
  POINTER (massflux_x_pbr,  massflux_x_bottom_rcv_buffer)
  POINTER (massflux_x_ptr,  massflux_x_top_rcv_buffer)
  POINTER (massflux_x_pls,  massflux_x_left_snd_buffer)
  POINTER (massflux_x_prs,  massflux_x_right_snd_buffer)
  POINTER (massflux_x_pbs,  massflux_x_bottom_snd_buffer)
  POINTER (massflux_x_pts,  massflux_x_top_snd_buffer)
  POINTER (massflux_x_pltr, massflux_x_left_top_rcv_buffer)
  POINTER (massflux_x_prtr, massflux_x_right_top_rcv_buffer)
  POINTER (massflux_x_plbr, massflux_x_left_bottom_rcv_buffer)
  POINTER (massflux_x_prbr, massflux_x_right_bottom_rcv_buffer)
  POINTER (massflux_x_plts, massflux_x_left_top_snd_buffer)
  POINTER (massflux_x_prts, massflux_x_right_top_snd_buffer)
  POINTER (massflux_x_plbs, massflux_x_left_bottom_snd_buffer)
  POINTER (massflux_x_prbs, massflux_x_right_bottom_snd_buffer)

  POINTER (massflux_y_plr,  massflux_y_left_rcv_buffer)
  POINTER (massflux_y_prr,  massflux_y_right_rcv_buffer)
  POINTER (massflux_y_pbr,  massflux_y_bottom_rcv_buffer)
  POINTER (massflux_y_ptr,  massflux_y_top_rcv_buffer)
  POINTER (massflux_y_pls,  massflux_y_left_snd_buffer)
  POINTER (massflux_y_prs,  massflux_y_right_snd_buffer)
  POINTER (massflux_y_pbs,  massflux_y_bottom_snd_buffer)
  POINTER (massflux_y_pts,  massflux_y_top_snd_buffer)
  POINTER (massflux_y_pltr, massflux_y_left_top_rcv_buffer)
  POINTER (massflux_y_prtr, massflux_y_right_top_rcv_buffer)
  POINTER (massflux_y_plbr, massflux_y_left_bottom_rcv_buffer)
  POINTER (massflux_y_prbr, massflux_y_right_bottom_rcv_buffer)
  POINTER (massflux_y_plts, massflux_y_left_top_snd_buffer)
  POINTER (massflux_y_prts, massflux_y_right_top_snd_buffer)
  POINTER (massflux_y_plbs, massflux_y_left_bottom_snd_buffer)
  POINTER (massflux_y_prbs, massflux_y_right_bottom_snd_buffer)



  TYPE(chunk_type),  ALLOCATABLE       :: chunks(:)
  INTEGER                              :: number_of_chunks
  INTEGER                              :: num_chunks_x, num_chunks_y
  INTEGER                              :: num_neighbours 

  TYPE(grid_type)                      :: grid

END MODULE definitions_module
