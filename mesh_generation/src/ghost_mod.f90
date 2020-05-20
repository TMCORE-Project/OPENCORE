MODULE ghost_mod
  use parameters_mod
  
  implicit none
  contains
  
  !------------------------------------------------------------------------------
  ! SUBROUTINE CubedSphereFillHalo
  !
  ! Description:
  !   Recompute the cubed sphere data storage array, with the addition of a
  !   halo region around the specified panel.
  !------------------------------------------------------------------------------
  SUBROUTINE CubedSphereFillHalo(field)

    IMPLICIT NONE

    REAL, DIMENSION(ips:ipe,jps:jpe,ifs:ife), INTENT(INOUT) :: field
    
    REAL, DIMENSION(ids:ide,jds:jde,ifs:ife) :: field_inner

    ! Local variables
    INTEGER :: i
    
    !zarg = 0.0 !DBG
    field_inner = field(ids:ide,jds:jde,ifs:ife)

    ! Equatorial panels
    !IF (np==1) THEN
       DO i=1,nPVHalo
          field(1-i    ,jds:jde,1) = field_inner(ide-i  ,jds:jde,4)  !exchange left
          field(ide+i  ,jds:jde,1) = field_inner(i+1    ,jds:jde,2)  !exchange right
          field(ids:ide,1-i    ,1) = field_inner(ids:ide,jde-i  ,6)  !exchange below
          field(ids:ide,jde+i  ,1) = field_inner(ids:ide,i+1    ,5)  !exchange over
       ENDDO
    !ELSE IF (np==2) THEN
       DO i=1,nPVHalo
          field(1-i    ,jds:jde,2) = field_inner(ide-i  ,jds:jde   ,1)  !exchange left
          field(ide+i  ,jds:jde,2) = field_inner(i+1    ,jds:jde   ,3)  !exchange right
          field(ids:ide,1-i    ,2) = field_inner(ide-i  ,jde:jds:-1,6)  !exchange below
          field(ids:ide,jde+i  ,2) = field_inner(ide-i  ,jds:jde   ,5)  !exchange over
       ENDDO
    !ELSE IF (np==3) THEN
       DO i=1,nPVHalo
          field(1-i    ,jds:jde,3) = field_inner(ide-i     ,jds:jde,2)  !exchange left
          field(ide+i  ,jds:jde,3) = field_inner(i+1       ,jds:jde,4)  !exchange right
          field(ids:ide,1-i    ,3) = field_inner(ide:ids:-1,i+1    ,6)  !exchange below
          field(ids:ide,jde+i  ,3) = field_inner(ide:ids:-1,jde-i  ,5)  !exchange over
       ENDDO
    !ELSE IF (np==4) THEN
       DO i=1,nPVHalo
          field(1-i    ,jds:jde,4) = field_inner(ide-i     ,jds:jde   ,3) !exchange left
          field(ide+i  ,jds:jde,4) = field_inner(i+1       ,jds:jde   ,1) !exchange right
          field(ids:ide,1-i    ,4) = field_inner(i+1       ,jds:jde   ,6) !exchange below
          field(ids:ide,jde+i  ,4) = field_inner(i+1       ,jde:jds:-1,5) !exchange over
       ENDDO
    ! Top panel
    !ELSE IF (np==5) THEN
       DO i=1,nPVHalo
          field(1-i    ,jds:jde,5) = field_inner(ide:ids:-1,jde-i,4) !exchange left
          field(ide+i  ,jds:jde,5) = field_inner(ids:ide   ,jde-i,2) !exchange right
          field(ids:ide,1-i    ,5) = field_inner(ids:ide   ,jde-i,1) !exchange below
          field(ids:ide,jde+i  ,5) = field_inner(ide:ids:-1,jde-i,3) !exchange over
       ENDDO
    ! Bottom panel
    !ELSE IF (np==6) THEN
       DO i=1,nPVHalo
          field(1-i    ,jds:jde,6) = field_inner(ids:ide   ,i+1  ,4) !exchange left
          field(ide+i  ,jds:jde,6) = field_inner(ide:ids:-1,i+1  ,2) !exchange right
          field(ids:ide,1-i    ,6) = field_inner(ide:ids:-1,i+1  ,3) !exchange below
          field(ids:ide,jde+i  ,6) = field_inner(ids:ide   ,i+1  ,1) !exchange over
       ENDDO
    !ELSE
    !   WRITE (*,*) 'Fatal error: In CubedSphereFillHalo'
    !   WRITE (*,*) 'Invalid panel id ', np
    !   STOP
    !ENDIF
  END SUBROUTINE CubedSphereFillHalo
  
END MODULE ghost_mod