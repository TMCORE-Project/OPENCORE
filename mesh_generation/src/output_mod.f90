module output_mod
  use netcdf
  use constants_mod
  use parameters_mod
  use mesh_mod
  implicit none
    
    character(13) :: ncFile = 'ccs_output.nc'
    
    contains
    subroutine history_init
      
      integer status
      integer ncid
      integer lon_dim_id,lat_dim_id
      integer two_dim_id,three_dim_id
      integer patch_dim_id
      integer x_id,y_id,z_id
      integer lon_id,lat_id
      integer xi_id,eta_id
      integer jab_id
      integer sqrtG_id
      integer matrixG_id
      integer matrixIG_id
      integer matrixA_id
      integer matrixIA_id
      
      status = nf90_create(ncFile, NF90_CLOBBER + NF90_NETCDF4 , ncid)
      if(status/=nf90_noerr) call handle_err(status)
      
      status = nf90_def_dim(ncid,'lon'   ,Nx,lon_dim_id   )
      status = nf90_def_dim(ncid,'lat'   ,Ny,lat_dim_id   )
      status = nf90_def_dim(ncid,'two'   ,2 ,two_dim_id   )
      status = nf90_def_dim(ncid,'three' ,3 ,three_dim_id )
      status = nf90_def_dim(ncid,'nPatch',Nf,patch_dim_id )
      if(status/=nf90_noerr) call handle_err(status)
      
      status = nf90_def_var(ncid,'x'  ,NF90_DOUBLE,(/lon_dim_id,lat_dim_id,patch_dim_id/),x_id       )
      status = nf90_def_var(ncid,'y'  ,NF90_DOUBLE,(/lon_dim_id,lat_dim_id,patch_dim_id/),y_id       )
      status = nf90_def_var(ncid,'z'  ,NF90_DOUBLE,(/lon_dim_id,lat_dim_id,patch_dim_id/),z_id       )
      status = nf90_def_var(ncid,'lon',NF90_DOUBLE,(/lon_dim_id,lat_dim_id,patch_dim_id/),lon_id     )
      status = nf90_def_var(ncid,'lat',NF90_DOUBLE,(/lon_dim_id,lat_dim_id,patch_dim_id/),lat_id     )
      status = nf90_def_var(ncid,'xi' ,NF90_DOUBLE,(/lon_dim_id,lat_dim_id,patch_dim_id/),xi_id      )
      status = nf90_def_var(ncid,'eta',NF90_DOUBLE,(/lon_dim_id,lat_dim_id,patch_dim_id/),eta_id     )
      
      status = nf90_def_var(ncid,'jab'     ,NF90_DOUBLE,(/three_dim_id,two_dim_id,lon_dim_id,lat_dim_id,patch_dim_id/),jab_id     )
      status = nf90_def_var(ncid,'sqrtG'   ,NF90_DOUBLE,(/                        lon_dim_id,lat_dim_id,patch_dim_id/),sqrtG_id   )
      status = nf90_def_var(ncid,'matrixG' ,NF90_DOUBLE,(/two_dim_id  ,two_dim_id,lon_dim_id,lat_dim_id,patch_dim_id/),matrixG_id )
      status = nf90_def_var(ncid,'matrixIG',NF90_DOUBLE,(/two_dim_id  ,two_dim_id,lon_dim_id,lat_dim_id,patch_dim_id/),matrixIG_id)
      status = nf90_def_var(ncid,'matrixA' ,NF90_DOUBLE,(/two_dim_id  ,two_dim_id,lon_dim_id,lat_dim_id,patch_dim_id/),matrixA_id )
      status = nf90_def_var(ncid,'matrixIA',NF90_DOUBLE,(/two_dim_id  ,two_dim_id,lon_dim_id,lat_dim_id,patch_dim_id/),matrixIA_id)
      if(status/=nf90_noerr) call handle_err(status)
      
      !print*,'nf90_put_att'
      status = nf90_put_att(ncid,nf90_global       ,'dx'       ,dx*R2D)
      status = nf90_put_att(ncid,nf90_global       ,'dy'       ,dy*R2D)
      status = nf90_put_att(ncid,nf90_global       ,'ids'      ,ids)
      status = nf90_put_att(ncid,nf90_global       ,'ide'      ,ide)
      status = nf90_put_att(ncid,nf90_global       ,'jds'      ,jds)
      status = nf90_put_att(ncid,nf90_global       ,'jde'      ,jde)
      status = nf90_put_att(ncid,nf90_global       ,'ifs'      ,ifs)
      status = nf90_put_att(ncid,nf90_global       ,'ife'      ,ife)
      
      status = nf90_put_att(ncid,lon_id            ,'units'    ,'degree_east' )
      status = nf90_put_att(ncid,lat_id            ,'units'    ,'degree_north')
      
      status = nf90_put_att(ncid,xi_id             ,'long_name','x coordinate in local panel' )
      status = nf90_put_att(ncid,eta_id            ,'long_name','y coordinate in local panel'  )
      status = nf90_put_att(ncid,lon_id            ,'long_name','longitude on sphere coordinate for Cells' )
      status = nf90_put_att(ncid,lat_id            ,'long_name','latitude on sphere coordinate for Cells'  )
      status = nf90_put_att(ncid,jab_id            ,'long_name','jacobian matrix between cartesian coordinate and local panel'  )
      status = nf90_put_att(ncid,sqrtG_id          ,'long_name','mertric tensor'  )
      status = nf90_put_att(ncid,matrixG_id        ,'long_name','G matrix for converting contravariant to covariant'  )
      status = nf90_put_att(ncid,matrixIG_id       ,'long_name','IG matrix for converting covariant to contravariant'  )
      status = nf90_put_att(ncid,matrixA_id        ,'long_name','A matrix for converting contravariant to spherical vector'  )
      status = nf90_put_att(ncid,matrixIA_id       ,'long_name','IA matrix for converting spherical vector to contravariant'  )
      
      ! Define coordinates
      status = nf90_put_att(ncid, xi_id             ,'_CoordinateAxisTypes','lon lat nPatch')
      status = nf90_put_att(ncid, eta_id            ,'_CoordinateAxisTypes','lon lat nPatch')
      status = nf90_put_att(ncid, lon_id            ,'_CoordinateAxisTypes','lon lat nPatch')
      status = nf90_put_att(ncid, lat_id            ,'_CoordinateAxisTypes','lon lat nPatch')
      status = nf90_put_att(ncid, x_id              ,'_CoordinateAxisTypes','lon lat nPatch')
      status = nf90_put_att(ncid, y_id              ,'_CoordinateAxisTypes','lon lat nPatch')
      status = nf90_put_att(ncid, z_id              ,'_CoordinateAxisTypes','lon lat nPatch')
      status = nf90_put_att(ncid, sqrtG_id          ,'_CoordinateAxisTypes','lon lat nPatch')
      if(status/=nf90_noerr) call handle_err(status)
      
      status = nf90_enddef(ncid)
      if(status/=nf90_noerr) call handle_err(status)
      
      status = nf90_put_var(ncid,xi_id      , mesh%xi  / D2M)
      status = nf90_put_var(ncid,eta_id     , mesh%eta / D2M)
      status = nf90_put_var(ncid,lon_id     , mesh%lon * R2D)
      status = nf90_put_var(ncid,lat_id     , mesh%lat * R2D)
      status = nf90_put_var(ncid,x_id       , mesh%x        )
      status = nf90_put_var(ncid,y_id       , mesh%y        )
      status = nf90_put_var(ncid,z_id       , mesh%z        )
      
      status = nf90_put_var(ncid,jab_id     , mesh%jab     )
      status = nf90_put_var(ncid,sqrtG_id   , mesh%sqrtG   )
      status = nf90_put_var(ncid,matrixG_id , mesh%matrixG )
      status = nf90_put_var(ncid,matrixIG_id, mesh%matrixIG)
      status = nf90_put_var(ncid,matrixA_id , mesh%matrixA )
      status = nf90_put_var(ncid,matrixIA_id, mesh%matrixIA)
      if(status/=nf90_noerr) call handle_err(status)
      
      status = nf90_close(ncid)
      if(status/=nf90_noerr) call handle_err(status)
      
    end subroutine history_init
    
    subroutine handle_err(status)
      implicit none
      integer,intent(in)::status
            
      if(status/=nf90_noerr)then
          print*, trim(nf90_strerror(status))
          stop "Stopped by netCDF"
      endif  
    endsubroutine handle_err
end module output_mod
    