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

   LOGICAL      :: use_vector_loops ! Some loops work better in serial depending on the hardware

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
     !REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: density0,density1
     !REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: energy0,energy1
     !REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: pressure
     !REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: viscosity
     !REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: soundspeed
     !REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: xvel0,xvel1
     !REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: yvel0,yvel1
     !REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: vol_flux_x,mass_flux_x
     !REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: vol_flux_y,mass_flux_y
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

     INTEGER         :: chunk_neighbours(4) ! Chunks, not tasks, so we can overload in the future

     ! Idealy, create an array to hold the buffers for each field so a commuincation only needs
     !  one send and one receive per face, rather than per field.
     ! If chunks are overloaded, i.e. more chunks than tasks, might need to pack for a task to task comm 
     !  rather than a chunk to chunk comm. See how performance is at high core counts before deciding

     TYPE(field_type):: field

  END TYPE chunk_type

    !REAL(KIND=8),    DIMENSION(1,1) :: density0,density1
    !REAL(KIND=8),    DIMENSION(1,1) :: energy0,energy1
    !REAL(KIND=8),    DIMENSION(1,1) :: pressure
    !REAL(KIND=8),    DIMENSION(1,1) :: viscosity
    !REAL(KIND=8),    DIMENSION(1,1) :: soundspeed
    !REAL(KIND=8),    DIMENSION(1,1) :: xvel0,xvel1
    !REAL(KIND=8),    DIMENSION(1,1) :: yvel0,yvel1
    !REAL(KIND=8),    DIMENSION(1,1) :: vol_flux_x,mass_flux_x
    !REAL(KIND=8),    DIMENSION(1,1) :: vol_flux_y,mass_flux_y

    !REAL(KIND=8) :: density0(-1:482,-1:962)
    !REAL(KIND=8) :: density1(-1:482,-1:962)
    !REAL(KIND=8) :: energy0(-1:482,-1:962)
    !REAL(KIND=8) :: energy1(-1:482,-1:962)
    !REAL(KIND=8) :: pressure(-1:482,-1:962)
    !REAL(KIND=8) :: viscosity(-1:482,-1:962)
    !REAL(KIND=8) :: soundspeed(-1:482,-1:962)
    !REAL(KIND=8) :: xvel0(-1:483,-1:963)
    !REAL(KIND=8) :: xvel1(-1:483,-1:963)
    !REAL(KIND=8) :: yvel0(-1:483,-1:963)
    !REAL(KIND=8) :: yvel1(-1:483,-1:963)
    !REAL(KIND=8) :: vol_flux_x(-1:483,-1:962)
    !REAL(KIND=8) :: mass_flux_x(-1:483,-1:962)
    !REAL(KIND=8) :: vol_flux_y(-1:482,-1:963)
    !REAL(KIND=8) :: mass_flux_y(-1:482,-1:963)

    !REAL(KIND=8) :: density0(-1,-1)
    !REAL(KIND=8) :: density1(-1,-1)
    !REAL(KIND=8) :: energy0(-1,-1)
    !REAL(KIND=8) :: energy1(-1,-1)
    !REAL(KIND=8) :: pressure(-1,-1)
    !REAL(KIND=8) :: viscosity(-1,-1)
    !REAL(KIND=8) :: soundspeed(-1,-1)
    !REAL(KIND=8) :: xvel0(-1,-1)
    !REAL(KIND=8) :: xvel1(-1,-1)
    !REAL(KIND=8) :: yvel0(-1,-1)
    !REAL(KIND=8) :: yvel1(-1,-1)
    !REAL(KIND=8) :: vol_flux_x(-1,-1)
    !REAL(KIND=8) :: mass_flux_x(-1,-1)
    !REAL(KIND=8) :: vol_flux_y(-1,-1)
    !REAL(KIND=8) :: mass_flux_y(-1,-1)

    INTEGER, PARAMETER :: xmaxplustwo   = 122
    INTEGER, PARAMETER :: xmaxplusthree = 123
    INTEGER, PARAMETER :: ymaxplustwo   = 82
    INTEGER, PARAMETER :: ymaxplsuthree = 83

    REAL(KIND=8) :: density0(-1:xmaxplustwo,-1:ymaxplustwo)
    REAL(KIND=8) :: density1(-1:xmaxplustwo,-1:ymaxplustwo)
    REAL(KIND=8) :: energy0(-1:xmaxplustwo,-1:ymaxplustwo)
    REAL(KIND=8) :: energy1(-1:xmaxplustwo,-1:ymaxplustwo)
    REAL(KIND=8) :: pressure(-1:xmaxplustwo,-1:ymaxplustwo)
    REAL(KIND=8) :: viscosity(-1:xmaxplustwo,-1:ymaxplustwo)
    REAL(KIND=8) :: soundspeed(-1:xmaxplustwo,-1:ymaxplustwo)
    REAL(KIND=8) :: xvel0(-1:xmaxplusthree,-1:ymaxplsuthree)
    REAL(KIND=8) :: xvel1(-1:xmaxplusthree,-1:ymaxplsuthree)
    REAL(KIND=8) :: yvel0(-1:xmaxplusthree,-1:ymaxplsuthree)
    REAL(KIND=8) :: yvel1(-1:xmaxplusthree,-1:ymaxplsuthree)
    REAL(KIND=8) :: vol_flux_x(-1:xmaxplusthree,-1:ymaxplustwo)
    REAL(KIND=8) :: mass_flux_x(-1:xmaxplusthree,-1:ymaxplustwo)
    REAL(KIND=8) :: vol_flux_y(-1:xmaxplustwo,-1:ymaxplsuthree)
    REAL(KIND=8) :: mass_flux_y(-1:xmaxplustwo,-1:ymaxplsuthree)

    POINTER(pDensity0,   density0)
    POINTER(pDensity1,   density1)
    POINTER(pEnergy0,    energy0)
    POINTER(pEnergy1,    energy1)
    POINTER(pPressure,   pressure)
    POINTER(pViscosity,  viscosity)
    POINTER(pSoundspeed, soundspeed)
    POINTER(pXvel0,      xvel0)
    POINTER(pXvel1,      xvel1)
    POINTER(pYvel0,      yvel0)
    POINTER(pYvel1,      yvel1)
    POINTER(pVolFluxX,   vol_flux_x)
    POINTER(pMassFluxX,  mass_flux_x)
    POINTER(pVolFluxY,   vol_flux_y)
    POINTER(pMassFluxY,  mass_flux_y)

  TYPE(chunk_type),  ALLOCATABLE       :: chunks(:)
  INTEGER                              :: number_of_chunks
  INTEGER                              :: num_chunks_x, num_chunks_y

  TYPE(grid_type)                      :: grid

END MODULE definitions_module
