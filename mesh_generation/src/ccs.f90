    ! Conformal Cubed Sphere mesh generator
    program CCS
      use parameters_mod
      use mesh_mod
      use output_mod
      implicit none
      
      integer i,j,iPatch
      
      call initParameters
      
      call initMesh
      
      call generateConformalCubedSphere
      
      call history_init
    
    end program CCS
