
    MODULE mesh_mod
      use constants_mod
      use parameters_mod
      use math_mod
      use jimcc_m
      use ghost_mod
      implicit none
      
      ! coordinate
      type mesh_info
        real, dimension(:,:,:    ), allocatable :: xi       ! central angle on x direction for cells on each patch, unit: radian
        real, dimension(:,:,:    ), allocatable :: eta      ! central angle on y direction for cells on each patch, unit: radian
        real, dimension(:,:,:    ), allocatable :: x        ! x cartesian coordinate on a unit sphere
        real, dimension(:,:,:    ), allocatable :: y        ! y cartesian coordinate on a unit sphere
        real, dimension(:,:,:    ), allocatable :: z        ! z cartesian coordinate on a unit sphere
        real, dimension(:,:,:    ), allocatable :: lon      ! longitude on cells
        real, dimension(:,:,:    ), allocatable :: lat      ! latitude on cells
        real, dimension(:,:,:    ), allocatable :: sqrtG    ! jacobian of Transformation, sqrt(G)
        real, dimension(:,:,:,:,:), allocatable :: matrixG  ! horizontal metric Tensor, which transform covariant vectors to contravariant vectors
        real, dimension(:,:,:,:,:), allocatable :: matrixIG ! horizontal metric Tensor, which transform contravariant vectors to covariant vectors
        real, dimension(:,:,:,:,:), allocatable :: matrixA  ! horizontal metric Tensor, which transform 
        real, dimension(:,:,:,:,:), allocatable :: matrixIA ! horizontal metric Tensor, which transform
        real, dimension(:,:,:,:,:), allocatable :: jab      ! jacobian matrix
        
        real, dimension(:,:,:    ), allocatable :: f       ! Coriolis parameter
        
        real, dimension(:,:,:    ), allocatable :: sinlon  ! sin(longitude)
        real, dimension(:,:,:    ), allocatable :: coslon  ! cos(longitude)
        real, dimension(:,:,:    ), allocatable :: sinlat  ! sin(latitude)
        real, dimension(:,:,:    ), allocatable :: coslat  ! cos(latitude)
        
      end type mesh_info
      
      type(mesh_info), target :: mesh
      
      contains
      
      subroutine initMesh
        integer :: iPV, jPV, iCell, jCell, iPatch, iVar, iDOF, jDOF
        integer :: iPVs, iPVe, jPVs, jPVe
        
        ! Allocate arrays in structures
        allocate( mesh%xi       (      ids:ide, jds:jde, ifs:ife) )
        allocate( mesh%eta      (      ids:ide, jds:jde, ifs:ife) )
        allocate( mesh%x        (      ids:ide, jds:jde, ifs:ife) )
        allocate( mesh%y        (      ids:ide, jds:jde, ifs:ife) )
        allocate( mesh%z        (      ids:ide, jds:jde, ifs:ife) )
        allocate( mesh%lon      (      ids:ide, jds:jde, ifs:ife) )
        allocate( mesh%lat      (      ids:ide, jds:jde, ifs:ife) )
        allocate( mesh%sqrtG    (      ids:ide, ids:ide, ifs:ife) )
        allocate( mesh%matrixG  (2, 2, ids:ide, ids:ide, ifs:ife) )
        allocate( mesh%matrixIG (2, 2, ids:ide, ids:ide, ifs:ife) )
        allocate( mesh%matrixA  (2, 2, ids:ide, ids:ide, ifs:ife) )
        allocate( mesh%matrixIA (2, 2, ids:ide, ids:ide, ifs:ife) )
        allocate( mesh%jab      (3, 2, ids:ide, ids:ide, ifs:ife) )
        allocate( mesh%f        (      ids:ide, ids:ide, ifs:ife) )
        allocate( mesh%sinlon   (      ids:ide, ids:ide, ifs:ife) )
        allocate( mesh%coslon   (      ids:ide, ids:ide, ifs:ife) )
        allocate( mesh%sinlat   (      ids:ide, ids:ide, ifs:ife) )
        allocate( mesh%coslat   (      ids:ide, ids:ide, ifs:ife) )
        
        ! Calculate mesh infomation on VIA
        do iPatch = ifs, ife
          do jCell = jds, jde
            do iCell = ids, ide
              !mesh%x(iCell, jCell, iPatch) = (iCell - 0.5) * dx * R2D * D2M + x_min * D2M
              !mesh%y(iCell, jCell, iPatch) = (jCell - 0.5) * dy * R2D * D2M + y_min * D2M
              mesh%xi (iCell, jCell, iPatch) = (iCell - 1) * dx * R2D * D2M / 2. + x_min * D2M
              mesh%eta(iCell, jCell, iPatch) = (jCell - 1) * dy * R2D * D2M / 2. + y_min * D2M
            end do
          end do
        end do
      end subroutine initMesh
      
      subroutine generateConformalCubedSphere
      
        real lon,lat
        real xx,yy,xc(3)
        real dlondx,dlatdx
        real dlondy,dlatdy
        real dlondz,dlatdz
        real sph2cart_matrix(2,3)
        
        real dxdxi ,dydxi ,dzdxi
        real dxdeta,dydeta,dzdeta
        real dxdr  ,dydr  ,dzdr
        
        integer, parameter :: stencil_width = 2 * xhalo+1
        real fc(stencil_width)
        real dh
        
        real, dimension(:,:,:), allocatable :: x_tmp
        real, dimension(:,:,:), allocatable :: y_tmp
        real, dimension(:,:,:), allocatable :: z_tmp
        real, dimension(:,:,:), allocatable :: lon_tmp
        real, dimension(:,:,:), allocatable :: lat_tmp
        
        real, dimension(2,2) :: IA1,IA2
        
        integer i,j,k
        integer iPatch
        integer stat
      
        call inrot()
        
        ! Generating grids
        do k = ifs,ife
          !$OMP PARALLEL DO PRIVATE(i,xx,yy,xc,lon,lat)
          do j = jds, jde
            do i = ids, ide
              xx = mesh%xi (i,j,k)
              yy = mesh%eta(i,j,k)
              
              call mtoc(xx, yy, k, xc)
              call cart2sph(lon,lat,xc(1),xc(2),xc(3))
              
              ! reset longitude to [0,2*pi)
              if(lon<0    ) lon = 2. * pi + lon
              if(lon>2.*pi) lon = lon - 2. * pi
              
              if(k == 1) iPatch = 4
              if(k == 2) iPatch = 1
              if(k == 3) iPatch = 2
              if(k == 4) iPatch = 3
              if(k == 5) iPatch = 5
              if(k == 6) iPatch = 6
              !iPatch = k
              
              mesh%x  (i,j,iPatch) = xc(1)
              mesh%y  (i,j,iPatch) = xc(2)
              mesh%z  (i,j,iPatch) = xc(3)
              mesh%lon(i,j,iPatch) = lon
              mesh%lat(i,j,iPatch) = lat
            enddo
          enddo
          !$OMP END PARALLEL DO
        enddo
        
        ! rotate relative position in panel 5 and 6
        allocate(x_tmp  (ips:ipe,jps:jpe,ifs:ife))
        allocate(y_tmp  (ips:ipe,jps:jpe,ifs:ife))
        allocate(z_tmp  (ips:ipe,jps:jpe,ifs:ife))
        allocate(lon_tmp(ips:ipe,jps:jpe,ifs:ife))
        allocate(lat_tmp(ips:ipe,jps:jpe,ifs:ife))
        
        x_tmp  (ids:ide,jds:jde,ifs:ife) = mesh%x
        y_tmp  (ids:ide,jds:jde,ifs:ife) = mesh%y
        z_tmp  (ids:ide,jds:jde,ifs:ife) = mesh%z
        lon_tmp(ids:ide,jds:jde,ifs:ife) = mesh%lon
        lat_tmp(ids:ide,jds:jde,ifs:ife) = mesh%lat
        
        k = 5
        do i = ids,ide
          mesh%x  (ids:ide,i,k) = x_tmp  (ide-i+1,jds:jde,k)
          mesh%y  (ids:ide,i,k) = y_tmp  (ide-i+1,jds:jde,k)
          mesh%z  (ids:ide,i,k) = z_tmp  (ide-i+1,jds:jde,k)
          mesh%lon(ids:ide,i,k) = lon_tmp(ide-i+1,jds:jde,k)
          mesh%lat(ids:ide,i,k) = lat_tmp(ide-i+1,jds:jde,k)
        enddo
        
        k = 6
        do i = ids,ide
          mesh%x  (ids:ide,ide-i+1,k) = x_tmp  (ide-i+1,ide:ids:-1,k)
          mesh%y  (ids:ide,ide-i+1,k) = y_tmp  (ide-i+1,ide:ids:-1,k)
          mesh%z  (ids:ide,ide-i+1,k) = z_tmp  (ide-i+1,ide:ids:-1,k)
          mesh%lon(ids:ide,ide-i+1,k) = lon_tmp(ide-i+1,ide:ids:-1,k)
          mesh%lat(ids:ide,ide-i+1,k) = lat_tmp(ide-i+1,ide:ids:-1,k)
        enddo
        
        x_tmp  (ids:ide,jds:jde,ifs:ife) = mesh%x
        y_tmp  (ids:ide,jds:jde,ifs:ife) = mesh%y
        z_tmp  (ids:ide,jds:jde,ifs:ife) = mesh%z
        lon_tmp(ids:ide,jds:jde,ifs:ife) = mesh%lon
        lat_tmp(ids:ide,jds:jde,ifs:ife) = mesh%lat
        
        call CubedSphereFillHalo(x_tmp)
        call CubedSphereFillHalo(y_tmp)
        call CubedSphereFillHalo(z_tmp)
        
        ! Calculate matrices
        dh = dx / 2.
        do k = ifs,ife
          !$OMP PARALLEL DO PRIVATE(i,lon,lat,dlondx,dlondy,dlondz,sph2cart_matrix,fc,stat,IA1,IA2)
          do j = jds, jde
            do i = ids, ide
              lon = mesh%lon(i,j,k)
              lat = mesh%lat(i,j,k)
              
              dlondx = -sin(lon); dlatdx = -cos(lon)*sin(lat)
              dlondy =  cos(lon); dlatdy = -sin(lon)*sin(lat)
              dlondz = 0.       ; dlatdz =  cos(lat)
              
              sph2cart_matrix(1,1) = dlondx
              sph2cart_matrix(1,2) = dlondy
              sph2cart_matrix(1,3) = dlondz
              sph2cart_matrix(2,1) = dlatdx
              sph2cart_matrix(2,2) = dlatdy
              sph2cart_matrix(2,3) = dlatdz
              
              !dxdxi
              fc(1:stencil_width) = x_tmp(i-xhalo:i+xhalo,j,k)
              mesh%jab(1,1,i,j,k) = center_difference(fc,dh)
              !dydxi
              fc(1:stencil_width) = y_tmp(i-xhalo:i+xhalo,j,k)
              mesh%jab(2,1,i,j,k) = center_difference(fc,dh)
              !dzdxi
              fc(1:stencil_width) = z_tmp(i-xhalo:i+xhalo,j,k)
              mesh%jab(3,1,i,j,k) = center_difference(fc,dh)
              !dxdeta
              fc(1:stencil_width) = x_tmp(i,j-xhalo:j+xhalo,k)
              mesh%jab(1,2,i,j,k) = center_difference(fc,dh)
              !dydeta
              fc(1:stencil_width) = y_tmp(i,j-xhalo:j+xhalo,k)
              mesh%jab(2,2,i,j,k) = center_difference(fc,dh)
              !dzdeta
              fc(1:stencil_width) = z_tmp(i,j-xhalo:j+xhalo,k)
              mesh%jab(3,2,i,j,k) = center_difference(fc,dh)
              
              mesh%matrixG(:,:,i,j,k) = matmul( transpose(mesh%jab(:,:,i,j,k)), mesh%jab(:,:,i,j,k) )
              
              call BRINV(2,mesh%matrixG(:,:,i,j,k),mesh%matrixIG(:,:,i,j,k),stat)
              if(stat==0)stop 'Fail to calculate inverse of matrixG'
              
              mesh%sqrtG(i,j,k) = sqrt(det(mesh%matrixG(:,:,i,j,k)))
              
              mesh%matrixA(:,:,i,j,k) = matmul(sph2cart_matrix,mesh%jab(:,:,i,j,k))
              
              if(abs(lat*R2D)/=90)then
                call BRINV(2,mesh%matrixA(:,:,i,j,k),mesh%matrixIA(:,:,i,j,k),stat)
                if(stat==0)then
                  print*,'Fail to calculate inverse of matrixA at i,j,k,det :',i,j,k,det(mesh%matrixA(:,:,i,j,k))
                  stop
                endif
              elseif(lat*R2D==90)then
                ! Replace by analytical approximation on north pole
                IA1(1,1) = cos(lon)
                IA1(1,2) = -sin(lon)
                IA1(2,1) = sin(lon)
                IA1(2,2) = cos(lon)
                
                IA2(1,1) = cos(pi/4.)**2
                IA2(1,2) = 0.
                IA2(2,1) = 0.
                IA2(2,2) = cos(pi/4.)**2
                
                mesh%matrixIA(:,:,i,j,k) = matmul(IA2,IA1)
                !print*,mesh%matrixIA(:,:,i,j,k)
              elseif(lat*R2D==-90)then
                ! Replace by analytical approximation on south pole
                IA1(1,1) = -cos(lon)
                IA1(1,2) = sin(lon)
                IA1(2,1) = sin(lon)
                IA1(2,2) = cos(lon)
                
                IA2(1,1) = cos(pi/4.)**2
                IA2(1,2) = 0.
                IA2(2,1) = 0.
                IA2(2,2) = cos(pi/4.)**2
                
                mesh%matrixIA(:,:,i,j,k) = matmul(IA2,IA1)
                !print*,mesh%matrixIA(:,:,i,j,k)
              endif
              
            enddo
          enddo
          !$OMP END PARALLEL DO
        enddo
        
        print*,'Total spherical area',sum(mesh%sqrtG(ids+1:ide-1:2,jds+1:jde-1:2,:)) * dx**2
        
      end subroutine generateConformalCubedSphere
      
    END MODULE mesh_mod

