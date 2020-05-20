! Conformal Cubic Atmospheric Model
    
! Copyright 2015 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
! This file is part of the Conformal Cubic Atmospheric Model (CCAM)
!
! CCAM is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! CCAM is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with CCAM.  If not, see <http://www.gnu.org/licenses/>.

!------------------------------------------------------------------------------
    
module jimcc_m

  use precis_m
  implicit none

  public :: jimcc, ctom  ! ctom needs to be public for newdepts.
  !private :: inrot, rgrid, inhedra, mtoc, mtocd, flip, stoc, cttom, &
  !           stoe, tay, tayd, taydd, lufm, invmm, vmtoc, vmtocd, vtay, vtaydd

! From cstoc (internal Purser)
  complex(kind=rx), save, private :: ci, cip4, cip3oss
  real(kind=rx), parameter, private, dimension(30) :: a =                   &
      (/ 1.47713062600964_rx, -0.38183510510174_rx, -0.05573058001191_rx,   &
        -0.00895883606818_rx, -0.00791315785221_rx, -0.00486625437708_rx,   &
        -0.00329251751279_rx, -0.00235481488325_rx, -0.00175870527475_rx,   &
        -0.00135681133278_rx, -0.00107459847699_rx, -0.00086944475948_rx,   &
        -0.00071607115121_rx, -0.00059867100093_rx, -0.00050699063239_rx,   &
        -0.00043415191279_rx, -0.00037541003286_rx, -0.00032741060100_rx,   &
        -0.00028773091482_rx, -0.00025458777519_rx, -0.00022664642371_rx,   &
        -0.00020289261022_rx, -0.00018254510830_rx, -0.00016499474460_rx,   &
        -0.00014976117167_rx, -0.00013646173947_rx, -0.00012478875822_rx,   &
        -0.00011449267279_rx, -0.00010536946150_rx, -0.00009725109376_rx /)
  real(kind=rx), parameter, private, dimension(30) :: b =                   &
      (/ 0.67698819751739_rx, 0.11847293456554_rx, 0.05317178134668_rx,     &
         0.02965810434052_rx, 0.01912447304028_rx, 0.01342565621117_rx,     &
         0.00998873323180_rx, 0.00774868996406_rx, 0.00620346979888_rx,     &
         0.00509010874883_rx, 0.00425981184328_rx, 0.00362308956077_rx,     &
         0.00312341468940_rx, 0.00272360948942_rx, 0.00239838086555_rx,     &
         0.00213001905118_rx, 0.00190581316131_rx, 0.00171644156404_rx,     &
         0.00155493768255_rx, 0.00141600715207_rx, 0.00129556597754_rx,     &
         0.00119042140226_rx, 0.00109804711790_rx, 0.00101642216628_rx,     &
         0.00094391366522_rx, 0.00087919021224_rx, 0.00082115710311_rx,     &
         0.00076890728775_rx, 0.00072168382969_rx, 0.00067885087750_rx /)

  real(kind=rx), private, save :: ss, third, fourth

! Initialize some of the arrays, particularly those associated
! with the implicit representation of the group of symmetries of the cube.

! From chedra (internal Purser)
  integer, parameter, private, dimension(0:47) :: kda =                     &
              (/ 0,1,1,2,2,0, 0,1,1,5,5,0, 3,1,1,2,2,3, 3,1,1,5,5,3,        &
                 0,4,4,2,2,0, 0,4,4,5,5,0, 3,4,4,2,2,3, 3,4,4,5,5,3 /)
  integer, parameter, private, dimension(0:47) :: kdb =                     &
              (/ 3,10,4,11,5,9, 6,7,1,11,5,0, 3,1,7,8,2,9,  6,4,10,8,2,0,   &
                 0,10,4,2,8,6, 9,7,1,2,8,3,  0,1,7,5,11,6, 9,4,10,5,11,3 /)
  integer, parameter, private, dimension(0:47) :: kdc =                     &
              (/ 5,11,3,9,4,10, 5,11,6,0,1,7, 2,8,3,9,7,1,  2,8,6,0,10,4,   &
                 8,2,0,6,4,10, 8,2,9,3,1,7,  11,5,0,6,7,1, 11,5,9,3,10,4 /)
  integer, parameter, private, dimension(0:47) :: kna =                     &
              (/ 12,25,26,9,10,17,   18,31,32,3,4,23,                       &
                  0,37,38,21,22,5,   6,43,44,15,16,11,                      &
                 36,1,2,33,34,41,   42,7,8,27,28,47,                        &
                 24,13,14,45,46,29, 30,19,20,39,40,35 /)
  integer, parameter, private, dimension(0:47) :: knb =                     &
              (/  5,2,1,4,3,0,       11,8,7,10,9,6,                         &
                 17,14,13,16,15,12, 23,20,19,22,21,18,                      &
                 29,26,25,28,27,24, 35,32,31,34,33,30,                      &
                 41,38,37,40,39,36, 47,44,43,46,45,42 /)
  integer, parameter, private, dimension(0:47) :: knc =                     &
              (/  1,0,3,2,5,4,       7,6,9,8,11,10,                         &
                 13,12,15,14,17,16, 19,18,21,20,23,22,                      &
                 25,24,27,26,29,28, 31,30,33,32,35,34,                      &
                 37,36,39,38,41,40, 43,42,45,44,47,46 /)

  integer, private, save, dimension(0:47) :: ipofig, kgofig
  integer, private, save, dimension(8,6) :: igofkg
  integer, save, private :: ig = 0
  real(kind=rx), parameter, private, dimension(2,2,8) :: flip8 =            &
           reshape ( (/ 1.0,0.0,0.0,1.0,  -1.0,0.0,0.0,1.0,                 &
                        1.0,0.0,0.0,-1.0,  -1.0,0.0,0.0,-1.0,               &
                        0.0,1.0,1.0,0.0,  0.0,-1.0,1.0,0.0,                 &
                        0.0,1.0,-1.0,0.0,  0.0,-1.0,-1.0,0.0 /),            &
                     (/ 2, 2, 8 /) )

  real(kind=rx), private, save, dimension(3,0:5) :: f
  real(kind=rx), private, save, dimension(3,0:11) :: e
  real(kind=rx), private, save, dimension(3,3) :: txe
  real(kind=rx), private, save, dimension(3,3,0:47) :: rotg

contains

subroutine jimcc ( np, xx4, yy4, em4, ax4, ay4, az4 )

!     like jim6.f but without stretch option
!     hedra1 data is hardwired
!     xx-->xx4, yy-->yy4, fm-->em4, dxa-->ax4, dxb-->ay4, dxc-->az4
!     dya, dyb, dyc commented out (not needed)

  integer, intent(in) :: np  ! = iquad
! Actual dimensions(iquad,iqiad) = (np,np)
  real(kind=rx), intent(out), dimension(:,:) :: xx4, yy4, em4, ax4, ay4, az4

  integer, parameter :: ipanel=2, ngr=1

  real(kind=rx), dimension(np,np) :: xa, xb, xc
  integer :: i, j

  call inrot()

  call rgrid ( xa, xb, xc, ax4, ay4, az4, em4, np, ipanel, ngr )

! Impose 8 fold symmetry (jlm)
  do i=1,(np+1)/2
     xc(i,1) = max(-1.0_rx,xc(i,1)) ! for rounding errors
  end do
  do j=1,(np+1)/2
     do i=j+1,(np+1/2)  ! rest of bottom LH corner
        xa(j,i) = xa(i,j)
        xb(j,i) = xc(i,j)
        xc(j,i) = xb(i,j)
     end do
     do i=1,np/2        ! all of bottom RH corner
        xa(np+1-i,j) = xa(i,j)
        xb(np+1-i,j) = -xb(i,j)
        xc(np+1-i,j) = xc(i,j)
     end do
     do i=1,np          ! all of top half
        xa(i,np+1-j) = xa(i,j)
        xb(i,np+1-j) = xb(i,j)
        xc(i,np+1-j) = -xc(i,j)
     end do
  end do

! Now convert these to just x,y values on the cube
  xx4 = xb/xa
  yy4 = xc/xa

end subroutine jimcc

!-----------------------------------------------------------------------------
!   R.J.Purser, National Meteorological Center, Washington D.C.  1994
!                   SUBROUTINE  INROT
!
!   Initialize the rotation matrix needed to set the orientation
!   of the cube to user's specification, then call general initialization
!   routine INHEDRA to make the transformation code ready for use.
!
!  ROT is the user-defined orthogonal matrix determining the orientation of
!      the mapping-cube with respect to the standard earth-centered
!      coordinate basis defined in the comments below.
!  IG1  is an array of 6 group indices. For designated map-panel IPANEL, the
!           element IG1(IPANEL) is associated with the octant of this map
!           in contact with the map-origin and that abuts the x-axis
!   Note: it is up to the user to ensure that ROT is indeed orthogonal and
!   that the 6 elements specified in IG1 do indeed belong to the 6 distinct
!   faces of the cube in accordance with the numbering convention adopted
!   for the triangulation into "group elements".
!----------------------------------------------------------------------------
subroutine inrot()
  real(kind=rx), dimension(3,3) :: rot
  integer, parameter, dimension(6) :: ig1 = (/ 22, 13, 40, 43, 9, 45 /)
  real(kind=rx) :: r2, r3, r6

!  SPECIFY THE ORIENTATION OF THE TRIANGULATED MAPPING-CUBE IN TERMS OF 3
!  ORTHOGONAL UNIT-VECTORS, REFERRED TO HERE AS p, q, r, DEFINED AS FOLLOWS:
!  VECTOR q POINTS TOWARDS THE MID-POINT OF EDGE-3 WHERE TRIANGULAR ELEMENTS
!  24, 29, 41 AND 36 MEET; VECTOR r POINTS TOWARDS VERTEX-0 WHERE ELEMENTS
!  0, 1, 2, 3, 4 AND 5 MEET; VECTOR p POINTS TOWARDS A POSITION ON THE
!  COMMON BOUNDARY OF ELEMENTS 14 AND 15 (IN FACE-3) SUCH THAT IT IS
!  ORTHOGONAL TO BOTH q AND r, OR EQUIVALENTLY, p IS DEFINED AS THE
!  (RIGHT-HANDED) CROSS-PRODUCT, q*r. THE BASIS-VECTORS USED TO EXPRESS
!  p, q, AND r ARE (1): THE UNIT VECTOR POINTING TO LAT,LONG=0,0; (2): THE
!  VECTOR POINTING TO LAT,LONG=0,90E; (3) THE UNIT VECTOR POINTING TO THE
!  NORTH POLE.

  r2 = sqrt(2.0_rx)
  r3 = sqrt(3.0_rx)
  r6 = r2*r3

! DEFINE VECTOR p AND MAKE IT THE FIRST COLUMN OF MATRIX ROT:
  rot(1,1) =  r6/6.0_rx
  rot(2,1) = -r6/6.0_rx
  rot(3,1) = -r6/3.0_rx

! DEFINE VECTOR q AND MAKE IT THE SECOND COLUMN OF MATRIX ROT:
  rot(1,2) = r2/2.0_rx
  rot(2,2) = r2/2.0_rx
  rot(3,2) = 0.0_rx

! DEFINE VECTOR r AND MAKE IT THE THIRD COLUMN OF MATRIX ROT:
  rot(1,3) =  r3/3.0_rx
  rot(2,3) = -r3/3.0_rx
  rot(3,3) =  r3/3.0_rx

!  CUSTOMIZATION OF THE MAPPING TRANSFORMATION IS COMPLETED BY SPECIFYING,
!  FOR EACH NUMBERED MAP-PANEL (FROM 1 TO 6) THE TRIANGULAR ELEMENT THAT
!  CORRESPONDS TO THE 3 RESTRICTIONS IN THE LOCAL COORDINATES x,y:
!          a: x < .5;
!          b: 0. < y;
!          c: y < x.
!  FOR EACH MAP-PANEL, IP, THIS BASIC ELEMENT IS PRESCRIBED IN IG1(IP).
!  THESE, TOGETHER WITH THE ORTHOGONAL MATRIX ROT MADE UP OF p,q,r, ARE
!  PASSED TO THE GENERAL INITIALIZATION ROUTINE INHEDRA, AFTER WHICH, THE
!  MAP-TRANSFORMATION ROUTINES ARE READY FOR USE.

  call inhedra(rot,ig1)

end subroutine inrot

!-------------------------------------------------------------------------
!   r.j.purser, national meteorological center, washington d.c.  1994
!                   subroutine  rgrid
!
!   set up array of standard earth-centered coordinates of a chosen panel
!   of the rancic map. convention used for map coordinates here is that
!   each origin is the pole (north or south as appropriate) and the (x,y)
!   map coordinates in each panel form a right-handed pair, with x and y both
!   between 0 and 1. to choose another panel-corner for origin, or to alter
!   chiral convention, just rearrange the table igofkg. for more radical
!   change in map-coordinate convention, rewrite this routine!
!
!  <-- xe,ye,ze     earth-centered coordinate of regular map-grid
!  --> np           number of points along each grid line (edges included)
!  --> ipanel       map-panel index [0,5]
!----------------------------------------------------------------------------
subroutine rgrid ( xe, ye, ze, dxa, dxb, dxc, em4, np, ipanel, ngr )

!     ngr = 1  at unstaggered positions
!         = 2  at i+.5,j+.5   positions
!         = 3  at i+.5,j      positions
!         = 4  at i,j+.5      positions
!     for staggered grids, actually only use 1,n part of the array
! set up earth-centered coordinates for points of chosen panel of the
!  rancic map

  integer, intent(in) :: np, ipanel, ngr
! Actually dimension(np,np)
  real(kind=rx), intent(out), dimension(:,:) :: xe, ye, ze,                 &
                                                dxa, dxb, dxc,              &
                                                em4

  real(kind=rx), parameter :: stretch  = 1.0_rx
  real(kind=rx), parameter :: stretchm = 1.0_rx - stretch

! Local variables
  real(kind=rx), dimension(3,np)   :: ex
  real(kind=rx), dimension(3,3,np) :: xc
  real(kind=rx), dimension(2,np)   :: dfdx
  integer :: n, j, jp, i, ip
  real(kind=rx) :: d, xadd, yadd, y, x, den
  real(kind=rx), dimension(np) :: xvec

  n = np - 1
  d = 1.0_rx/real(n,kind=rx)
  xadd = 0.0
  yadd = 0.0
  if ( ngr==2 .or. ngr==3 ) then
     xadd = 0.5_rx*d
  endif
  if ( ngr==2 .or. ngr==4 ) then
     yadd = 0.5_rx*d
  endif
  do j = 0, n
     jp = j + 1
     y = j*d + yadd                          ! jlm allows staggered v
     y = 0.5_rx + (y - 0.5_rx)*(stretch + stretchm*(2.0*y - 1.0)**2)
     do i = 0, n
        ip = i + 1
        x = i*d + xadd                       ! jlm allows staggered u
        xvec(ip) = 0.5_rx + (x - 0.5_rx)*(stretch + stretchm*(2.0*x - 1.0)**2)
     end do
     call vmtoc (xvec, y, ipanel, ex, np)
     do ip=1,np
        xe(ip,jp) = ex(1,ip)
        ye(ip,jp) = ex(2,ip)
        ze(ip,jp) = ex(3,ip)
!       print*, ' XYZ ', ip, jp, xe(ip,jp), ye(ip,jp), ze(ip,jp)
     end do
     call vmtocd (xvec, y, ipanel, xc, em4(:,jp), dfdx, np)
!    return dxa etc as unit vectors
     do ip=1,np
        den = sqrt(xc(1,1,ip)**2+xc(2,1,ip)**2+xc(3,1,ip)**2)
        if (den < 1.0e-6_rx) then
           den = 1.0
        endif
        dxa(ip,jp) = xc(1,1,ip)/den ! the three components of a vector along dx
        dxb(ip,jp) = xc(2,1,ip)/den
        dxc(ip,jp) = xc(3,1,ip)/den
     end do
  end do
end subroutine rgrid

!                                               *****************
!                                               *   HED6A.FOR   *
!                                               *  PURSER 1994  *
!                                               *****************
!
!  SUBROUTINES NEEDED TO PERFORM RANCIC-CUBE TRANSFORMATIONS
!

!----------------------------------------------------------------------------
!   R.J.Purser, National Meteorological Center, Washington D.C.  1994
!                   SUBROUTINE INHEDRA
!   Initialize variables needed to perform the Rancic-transformation and
!   its inverse
!----------------------------------------------------------------------------
subroutine inhedra(rot,ig1)
!  real, intent(in), dimension(3,3) :: rot
!  integer, intent(in), dimension(6) :: ig1
  real(kind=rx), intent(in), dimension(:,:) :: rot
  integer, intent(in), dimension(:) :: ig1

  integer, dimension(8) :: igk
  integer :: ip, kg, i, j, k, lg
  real(kind=rx) :: r2, r3, r6, r2o2, r3o2, r3o3, r6o3, r6o6, r3o6

!  SET UP CROSS-REFERENCE TABLES CONNECTING GROUP-ELEMENTS IG TO
!  USER-DEFINED NUMBERING AND ORIENTATIONS OF PANELS:
  do ip = 1,6
     igk(1) = ig1(ip)
     igk(2) = kna(igk(1))
     igk(5) = knc(igk(1))
     igk(6) = knc(igk(2))
     igk(7) = kna(igk(5))
     igk(8) = kna(igk(6))
     igk(3) = knc(igk(7))
     igk(4) = knc(igk(8))
     do kg = 1,8
        ig = igk(kg)
        igofkg(kg,ip) = ig
        ipofig(ig) = ip
        kgofig(ig) = kg
     end do
  end do
  r2 = sqrt(2.0_rx)
  r3 = sqrt(3.0_rx)
  r6 = sqrt(6.0_rx)
  r2o2 = r2/2.0_rx
  r3o2 = r3/2.0_rx
  r3o3 = r3/3.0_rx
  r6o3 = r6/3.0_rx
  r6o6 = r6/6.0_rx
  r3o6 = r3/6.0_rx
  ss = r2
  third = 1.0_rx/3.0_rx
  fourth = 1.0_rx/4.0_rx
  ci = cmplx(0.0,1.0,kind=rx)
  cip4 = ci**fourth
  cip3oss = ci**third/ss
  f(1,0) = -r6o3
  f(2,0) = 0.0
  f(3,0) = r3o3
  f(1,1) = r6o6
  f(2,1) = -r2o2
  f(3,1) = r3o3
  f(1,2) = r6o6
  f(2,2) = r2o2
  f(3,2) = r3o3
  e(1,0) = r3o3
  e(2,0) = 0.0
  e(3,0) = r6o3
  e(1,1) = -r3o6
  e(2,1) = 0.5_rx
  e(3,1) = r6o3
  e(1,2) = -r3o6
  e(2,2) = -0.5_rx
  e(3,2) = r6o3
  e(1,3) = 0.0
  e(2,3) = 1.0
  e(3,3) = 0.0
  e(1,4) = -r3o2
  e(2,4) = -0.5_rx
  e(3,4) = 0.0
  e(1,5) = r3o2
  e(2,5) = -0.5_rx
  e(3,5) = 0.0
  do j = 0,2
     k = j+3
     do i = 1,3
        f(i,k) = -f(i,j)
     enddo
  enddo
  do j = 0,5
     k = j+6
     do i = 1,3
        e(i,k) = -e(i,j)
     enddo
  enddo
  do i = 1,3
     txe(1,i) = f(i,0)
     txe(2,i) = e(i,3)
     txe(3,i) = e(i,5)
  enddo

!  call invmm(txe,txe,3,3,3)
  call invmm(txe)
! ROTATE THE 6 FACE-VECTORS, F, TO USER-DEFINED ORIENTATION:
  f = matmul ( rot, f )

! ROTATE THE 12 EDGE-VECTORS, E, TO USER-DEFINED ORIENTATION:
  e = matmul ( rot, e )

!  BASED ON THE PRESCRIBED ORIENTATION (DEFINED BY "ROT"),
!  CONSTRUCT THE ROTATION MATRIX ASSOCIATED WITH EACH GROUP ELEMENT LG:

  do lg = 0,47
     do j = 1,3
        do i = 1,3
           rotg(i,j,lg) = f(i,kda(lg))*txe(j,1) + e(i,kdb(lg))*txe(j,2) +    &
                          e(i,kdc(lg))*txe(j,3)
        enddo
     enddo
  enddo

end subroutine inhedra

!----------------------------------------------------------------------------
!   R.J.Purser, National Meteorological Center, Washington D.C.  1994
!                   SUBROUTINE MTOC
!   Transform from map-coordinates to standard earth-centered cartesians
!
!  -->    XX,YY  map-coordinates
!  -->    IPANEL map-panel index
!  <--    XC     standard earth-centered cartesian coordinates
!----------------------------------------------------------------------------

subroutine mtoc(xx, yy, ipanel, xc)

  integer, intent(in) :: ipanel
  real(kind=rx), intent(in) :: xx
  real(kind=rx), intent(in) :: yy
! dimension(3)
  real(kind=rx), dimension(:), intent(out) :: xc

  integer :: kg
  real(kind=rx), dimension(3) :: xv
  real(kind=rx) :: x, y, t, xw, yw, h
  complex(kind=rx) :: w, z, arg
!-----------------------------------------------
  x = xx
  y = yy
  kg = 1
  if ( x > 0.5_rx ) then
     kg = kg + 1
     x = 1.0 - x
  endif
  if (y > 0.5_rx) then
     kg = kg + 2
     y = 1.0 - y
  endif
  if (y > x) then
     kg = kg + 4
     t = x
     x = y
     y = t
  endif
! z=cmplx(x,y)**4
  z = cmplx(x,y,kind=rx)*cmplx(x,y,kind=rx)
  z = z*z
  call tay (z, a, 30, w)
  arg = -ci*w                                ! mrd
  if (abs(arg) == 0.0) then                 ! mrd
     w = (0.0,0.0)                           ! mrd
  else                                       ! mrd
     w = cip3oss*(-ci*w)**third              ! mrd
  endif                                      ! mrd
  xw = real(w)
  yw = aimag(w)
  h = 2.0/(1.0 + xw*xw + yw*yw)
  xv(1) = xw*h
  xv(2) = yw*h
  xv(3) = h - 1.0
  ig = igofkg(kg,ipanel)
  xc = matmul ( rotg(:,:,ig), xv )

end subroutine mtoc

 subroutine vmtoc(xx, yy, ipanel, xc, n)

  integer, intent(in) :: ipanel
  integer, intent(in) :: n
  real(kind=rx), intent(in), dimension(:) :: xx
  real(kind=rx), intent(in) :: yy
! dimension(3,n)
  real(kind=rx), dimension(:,:), intent(out) :: xc

  integer :: i
  integer, dimension(n) :: kg
  real(kind=rx), dimension(3) :: xv
  real(kind=rx), dimension(n) :: x, y, t, xw, yw, h
  complex(kind=rx), dimension(n) :: w, z, arg
!-----------------------------------------------
  x = xx
  y = yy
  kg = 1
  where ( x > 0.5_rx ) 
     kg = kg + 1
     x = 1.0 - x
  endwhere
  where (y > 0.5_rx) 
     kg = kg + 2
     y = 1.0 - y
  endwhere
  where (y > x) 
     kg = kg + 4
     t = x
     x = y
     y = t
  endwhere

! z=cmplx(x,y)**4
  z = cmplx(x,y,kind=rx)*cmplx(x,y,kind=rx)
  z = z*z
  call vtay (z, a, 30, w)
  arg = -ci*w                                ! mrd
  where ( abs(arg) == 0.0 )                  ! mrd
     w = (0.0,0.0)                           ! mrd
  elsewhere                                  ! mrd
     w = cip3oss*(-ci*w)**third              ! mrd
  endwhere                                   ! mrd
  xw = real(w)
  yw = aimag(w)
  h = 2.0/(1.0 + xw*xw + yw*yw)
  do i=1,n
     xv(1) = xw(i)*h(i)
     xv(2) = yw(i)*h(i)
     xv(3) = h(i) - 1.0
     ig = igofkg(kg(i),ipanel)
     xc(:,i) = matmul ( rotg(:,:,ig), xv )
  end do


end subroutine vmtoc




!---------------------------------------------------------------------------
!   R.J.Purser, National Meteorological Center, Washington D.C.  1994
!                   SUBROUTINE MTOCD
!   Transform from map-coordinates to standard earth-centered cartesians
!   and simultaneously accumulate the derivative of the transformation in
!   order to provide information on map-scaling factor and relative
!   orientation
!
!  --> XX,YY  map-coordinates (from corner IG, each panel a unit square)
!  --> IPANEL map-panel index
!  <-- XC   augmented jacobian matrix: first two columns represent the
!           derivative with respect to X and Y of the earth-centered
!           cartesian coordinates of the point corresponding to map-image
!           X,Y. These cartesian coordinates themselves are inserted into
!           column-3 of XDC.
!  <-- em4   map-factor at this point
!  <-- DFDX x- and y-derivatives of map-factor here
!---------------------------------------------------------------------------

subroutine mtocd(xx, yy, ipanel, xc, em4, dfdx)

  integer, intent(in) :: ipanel
  real(kind=rx), intent(in) :: xx
  real(kind=rx), intent(in) :: yy
  real(kind=rx), intent(out) :: em4
! Actually dimension(3,3)
  real(kind=rx), intent(out), dimension(:,:)  :: xc
! dimension(2)
  real(kind=rx), intent(out), dimension(:)    :: dfdx

  integer :: kg
  real(kind=rx), dimension(3,3) :: xdc
  real(kind=rx), dimension(3,2) :: xd
  real(kind=rx), dimension(2) :: v1
  real(kind=rx) :: x, y, t, xw, yw, xwxw, xwyw, ywyw, h, hh, rd, qd, rdd, qdd, s, &
         dsdx, dsdy, dhdx, dhdy
  complex(kind=rx) :: w, z, zu, wu, cd, cdd, arg
!-----------------------------------------------
  x = xx
  y = yy
  kg = 1
  if ( x > 0.5_rx ) then
     kg = kg + 1
     x = 1.0 - x
  endif
  if ( y > 0.5_rx ) then
     kg = kg + 2
     y = 1.0 - y
  endif
  if ( y > x ) then
     kg = kg + 4
     t = x
     x = y
     y = t
  endif
  zu = cmplx(x,y,kind=rx)
  z = zu**4
  call taydd (z, a, 30, w, cd, cdd)
  arg = -ci*w                                ! mrd
  if ( abs(arg) == 0.0 ) then               ! mrd
     wu = (0.0,0.0)                          ! mrd
  else                                       ! mrd
     wu = cip3oss*(-ci*w)**third             ! mrd
  endif                                      ! mrd
  xw = real(wu,kind=rx)
  yw = aimag(wu)
  xwxw = xw*xw
  xwyw = xw*yw
  ywyw = yw*yw
  h = 2.0/(1.0 + xwxw + ywyw)
  hh = h*h
  xdc(1,3) = xw*h
  xdc(2,3) = yw*h
  xdc(3,3) = h - 1.0
  xdc(1,1) = h - hh*xwxw
  xdc(2,1) = -hh*xwyw
  xdc(3,1) = -hh*xw
  xdc(1,2) = xdc(2,1)
  xdc(2,2) = h - hh*ywyw
  xdc(3,2) = -hh*yw
  if ( abs(z) == 0.0 ) then
     cd    = 0.0
     cdd   = 0.0
     em4   = 0.0
     v1(1) = 0.0
     v1(2) = 0.0
  else
     cd = 4.0*wu*cd*z/(3.0*w*zu)
     cdd = 3.0*cd/zu - 2.0*cd*cd/wu + 16.0*wu*z**2*cdd/(3.0*w*zu**2)
     rd = real(cd,kind=rx)
     qd = aimag(cd)
     rdd = real(cdd,kind=rx)
     qdd = aimag(cdd)
     s = sqrt(rd*rd + qd*qd)
     dsdx = (rdd*rd + qdd*qd)/s
     dsdy = (rdd*qd - qdd*rd)/s
     dhdx = -hh*(xw*rd + yw*qd)
     dhdy = -hh*((-xw*qd) + yw*rd)
     em4 = h*s
     v1(1) = dhdx*s + h*dsdx
     v1(2) = dhdy*s + h*dsdy
  endif
  rd = real(cd,kind=rx)
  qd = aimag(cd)
  xd(:,1) = xdc(:,1)*rd + xdc(:,2)*qd
  xd(:,2) = (-xdc(:,1)*qd) + xdc(:,2)*rd
  ig = igofkg(kg,ipanel)
!  call mulmm (v1, flip8(1,1,kg), dfdx, 1, 2, 2, 1, 2, 1)
  dfdx = matmul ( transpose(flip8(:,:,kg)), v1 )
  xdc(:,1:2) = matmul ( xd, flip8(:,:,kg) )
  xc = matmul ( rotg(:,:,ig), xdc )

end subroutine mtocd

subroutine vmtocd(xx, yy, ipanel, xc, em4, dfdx, n)

  integer, intent(in) :: ipanel
  integer, intent(in) :: n
  real(kind=rx), intent(in), dimension(:) :: xx
  real(kind=rx), intent(in) :: yy
  real(kind=rx), intent(out), dimension(:) :: em4
! Actually dimension(3,3,n)
  real(kind=rx), intent(out), dimension(:,:,:)  :: xc
! dimension(2,n)
  real(kind=rx), intent(out), dimension(:,:)    :: dfdx

  integer, dimension(n) :: kg
  integer :: i
  real(kind=rx), dimension(3,3,n) :: xdc
  real(kind=rx), dimension(3,2,n) :: xd
  real(kind=rx), dimension(2,n) :: v1
  real(kind=rx), dimension(n) :: x, y, t, xw, yw, xwxw, xwyw, ywyw, h, hh, rd, qd, &
                        rdd, qdd, s, dsdx, dsdy, dhdx, dhdy
  complex(kind=rx), dimension(n) :: w, z, zu, wu, cd, cdd, arg
!-----------------------------------------------
  x = xx
  y = yy
  kg = 1
  where ( x > 0.5_rx ) 
     kg = kg + 1
     x = 1.0 - x
  endwhere
  where ( y > 0.5_rx ) 
     kg = kg + 2
     y = 1.0 - y
  endwhere
  where ( y > x ) 
     kg = kg + 4
     t = x
     x = y
     y = t
  endwhere
  zu = cmplx(x,y,kind=rx)
  z = zu**4
  call vtaydd (z, a, 30, w, cd, cdd)
  arg = -ci*w                                ! mrd
  where ( abs(arg) == 0.0 )                  ! mrd
     wu = (0.0,0.0)                          ! mrd
  elsewhere                                  ! mrd
     wu = cip3oss*(-ci*w)**third             ! mrd
  endwhere                                   ! mrd
  xw = real(wu,kind=rx)
  yw = aimag(wu)
  xwxw = xw*xw
  xwyw = xw*yw
  ywyw = yw*yw
  h = 2.0/(1.0 + xwxw + ywyw)
  hh = h*h
  xdc(1,3,:) = xw*h
  xdc(2,3,:) = yw*h
  xdc(3,3,:) = h - 1.0
  xdc(1,1,:) = h - hh*xwxw
  xdc(2,1,:) = -hh*xwyw
  xdc(3,1,:) = -hh*xw
  xdc(1,2,:) = xdc(2,1,:)
  xdc(2,2,:) = h - hh*ywyw
  xdc(3,2,:) = -hh*yw
  where ( abs(z) == 0.0 ) 
     cd    = 0.0
     cdd   = 0.0
     em4   = 0.0
     v1(1,:) = 0.0
     v1(2,:) = 0.0
  elsewhere
     cd = 4.0*wu*cd*z/(3.0*w*zu)
     cdd = 3.0*cd/zu - 2.0*cd*cd/wu + 16.0*wu*z**2*cdd/(3.0*w*zu**2)
     rd = real(cd,kind=rx)
     qd = aimag(cd)
     rdd = real(cdd,kind=rx)
     qdd = aimag(cdd)
     s = sqrt(rd*rd + qd*qd)
     dsdx = (rdd*rd + qdd*qd)/s
     dsdy = (rdd*qd - qdd*rd)/s
     dhdx = -hh*(xw*rd + yw*qd)
     dhdy = -hh*((-xw*qd) + yw*rd)
     em4 = h*s
     v1(1,:) = dhdx*s + h*dsdx
     v1(2,:) = dhdy*s + h*dsdy
  endwhere
  rd = real(cd,kind=rx)
  qd = aimag(cd)
  do i=1,3
     xd(i,1,:) = xdc(i,1,:)*rd + xdc(i,2,:)*qd
     xd(i,2,:) = (-xdc(i,1,:)*qd) + xdc(i,2,:)*rd
  end do
  do i=1,n
     ig = igofkg(kg(i),ipanel)
!  call mulmm (v1, flip8(1,1,kg), dfdx, 1, 2, 2, 1, 2, 1)
     dfdx(:,i) = matmul ( transpose(flip8(:,:,kg(i))), v1(:,i) )
     xdc(:,1:2,i) = matmul ( xd(:,:,i), flip8(:,:,kg(i)) )
     xc(:,:,i) = matmul ( rotg(:,:,ig), xdc(:,:,i) )
  end do

end subroutine vmtocd


!-----------------------------------------------------------------------------
!   R.J.Purser, National Meteorological Center, Washington D.C.  1994
!                   SUBROUTINE FLIP
!   Use standard earth-centered cartesian coordinates of a point to determine
!   the group-element IG of the orthogonal transformation needed to transform
!   this point into the standard triangular wedge, and the new cartesian
!   coordinates of the point once it has indergone this transformation. The
!   transformed position puts the point close to the pole and with small non-
!   negative components XC(1),XC(2), to ensure that it lies well inside the
!   circle of convergence of the Taylor series and with appropriate azimuth
!   to ensure intended solution is obtained when fractional power is taken.
!  --> XE standard earth-centered cartesian coordinates [3]
!  <-- XC transformed earth-centered cartesian coordinates [3]
!-----------------------------------------------------------------------------
subroutine flip(xe,xc)
!  real, dimension(3) :: xe, xc
  real(kind=rx), intent(in),  dimension(:) :: xe
  real(kind=rx), intent(out), dimension(:) :: xc
  real(kind=rx) :: dax, dbx, dcx
  integer :: i

! VALIDITY OF CONDITIONS A & B & C UNCERTAIN:
  dax = dot_product ( f(:,kda(ig)), xe)
  dbx = dot_product ( e(:,kdb(ig)), xe)
  if ( dax < 0.0 ) then
     dax = -dax
     ig = kna(ig)
  endif
  if ( dbx < 0.0 ) then
     dbx = -dbx
     ig = knb(ig)
  endif

! VALIDITY OF CONDITIONS A & B TRUE, C UNCERTAIN:
  do 
     dcx = dot_product ( e(:,kdc(ig)), xe )
     if ( dcx < 0.0 ) then
        dcx = -dcx
        ig = knc(ig)

!       VALIDITY OF CONDITIONS A & B UNCERTAIN, C TRUE:
        dax = dot_product ( f(:,kda(ig)), xe)
        dbx = dot_product ( e(:,kdb(ig)), xe)
        if ( dax < 0.0 ) then
           dax = -dax
           ig = kna(ig)
           if ( dbx < 0.0 ) then
              dbx = -dbx
              ig = knb(ig)
           endif
        else
           if ( dbx >= 0.0 ) then
              exit
           endif
           dbx = -dbx
           ig = knb(ig)
        endif
        cycle
     endif
     exit
  end do
!  VALIDITY OF CONDITIONS A & B & C TRUE:
   do i=1,3
      xc(i) = txe(i,1)*dax+txe(i,2)*dbx+txe(i,3)*dcx
   enddo
 end subroutine flip

!----------------------------------------------------------------------------
!   R.J.Purser, National Meteorological Center, Washington D.C.  1994
!                   SUBROUTINE STOC
!   Transform from latitude,longitude, to corner-origin coordinates of one
!   of the maps.
!   -->   DLAT, DLON    Latitude and longitude (degrees)
!   <--   X,Y           map coordinates scaled to the unit square
!   <--   IPANEL        map-panel index
!----------------------------------------------------------------------------
subroutine stoc(dlat,dlon,x,y,ipanel)
  real(kind=rx), intent(in)     :: dlat
  real(kind=rx), intent(in)     :: dlon
  real(kind=rx), intent(out)    :: x
  real(kind=rx), intent(out)    :: y
  integer, intent(out) :: ipanel

  real(kind=rx), dimension(3) :: xe

  call stoe(dlat,dlon,xe)
  call ctom(xe,x,y,ipanel)

end subroutine stoc


!--------------------------------------------------------------------------
!   R.J.Purser, National Meteorological Center, Washington D.C.  1994
!                   SUBROUTINE CTOM
!   Transform from group-oriented cartesians to map-coordinates.
!   Use CTTOM to transform first to corner-coordinates. Then use group
!   element IG to determine map-panel and perform the conversion from
!   corner-coordinates (0<y<x<.5) to full map-panel coordinates (0<x<1,
!  -->  XE      earth-centered cartesians
!  <--  X,Y     map-coordinates
!  <--  IPANEL  map-panel index
!--------------------------------------------------------------------------

subroutine ctom(xe, x, y, ipanel)
! dimension(3)
  real(kind=rx), dimension(:), intent(in) :: xe
  integer, intent(out) :: ipanel
  real(kind=rx), intent(out) :: x
  real(kind=rx), intent(out) :: y

  integer :: kg
  real(kind=rx), dimension(3) :: xc
  real(kind=rx) :: t
!-----------------------------------------------
  call flip (xe, xc)
  ipanel = ipofig(ig)
  kg = kgofig(ig)
  call cttom (xc, x, y)
  if (kg > 4) then
     t = x
     x = y
     y = t
     kg = kg - 4
  endif
  if (kg > 2) then
     y = 1.0 - y
     kg = kg - 2
  endif
  if (kg > 1) then
     x = 1.0 - x
  endif
end subroutine ctom



!---------------------------------------------------------------------------
!   R.J.Purser, National Meteorological Center, Washington D.C.  1994
!                   SUBROUTINE CTTOM
!   Transform from earth-centered cartesians to map-coordinates
!   Point must lie in standard wedge near the pole. Transform first to
!   complex-position in the stereographic plane. Then apply fractional-
!   power of Taylor-series.
!  -->  XE  earth-centered cartesians
!  <--  X,Y corner-origin map-coordinates oriented with respect to principal
!           edge.
!---------------------------------------------------------------------------
subroutine cttom(xe, x, y)
!  dimension(3)
   real(kind=rx), intent(in), dimension(:) :: xe
   real(kind=rx), intent(out) :: x
   real(kind=rx), intent(out) :: y

   real(kind=rx) :: hi
   complex(kind=rx) :: w, z
!-----------------------------------------------
!  SS=SQRT(2.) AND REPRESENTS THE SCALE FACTOR NEEDED
!  TO PLACE THE NEIGHBORING 3 VERTICES OF THE NOMINATED POLE OF THE
!      6-HEDRON ON THE UNIT CIRCLE OF THE RESCALED STEREOGRAPHIC MAP
!  CIP4 IS THE COMPLEX 4th-ROOT OF i (i.e., 8th-ROOT OF -1)

   hi = ss/(1.0 + xe(3))
   w = (cmplx(xe(1),xe(2),kind=rx)*hi)**3
   call tay (w, b, 30, z)

!  ROTATE AWAY FROM THE BRANCH-CUT THAT MARKS THE NEGATIVE-REAL AXIS:
!  TAKE 4th-ROOT AND ROTATE BACK AGAIN

   z = cip4*(-ci*z)**fourth
   x = real(z,kind=rx)
   y = aimag(z)

end subroutine cttom



!---------------------------------------------------------------------------
!   R.J.Purser, National Meteorological Center, Washington D.C.  1994
!                   SUBROUTINE STOE
!   Transform latitude and longitude (degrees) to earth-centered cartesian
!   coordinates.
!  --> DLAT     latitude
!  --> DLON     longitude
!  <-- XE       three cartesian components.
!---------------------------------------------------------------------------

subroutine stoe(dlat, dlon, xe)

   use parm_m, only : dtor
   real(kind=rx), intent(in) :: dlat
   real(kind=rx), intent(in) :: dlon
   !  dimension(3)
   real(kind=rx), intent(out), dimension(:) :: xe

   real(kind=rx) :: rlat, rlon, sla, cla, slo, clo

   rlat = dtor*dlat
   rlon = dtor*dlon
   sla = sin(rlat)
   cla = cos(rlat)
   slo = sin(rlon)
   clo = cos(rlon)
   xe(1) = cla*clo
   xe(2) = cla*slo
   xe(3) = sla

end subroutine stoe



!--------------------------------------------------------------------------
!   R.J.Purser, National Meteorological Center, Washington D.C.  1994
!                   SUBROUTINE  TAY
!  Evaluate the complex function W of Z whose real
!  Taylor series coefficients are RA.
!
!  --> Z    function argument (complex)
!  --> RA   Taylor coefficients (real)
!  --> N    number of coefficients (starting with the linear term)
!  <-- W    Taylor-series approximation of the function (complex)
!--------------------------------------------------------------------------
subroutine tay(z, ra, n, w)
   integer, intent(in) :: n
   complex(kind=rx), intent(in) :: z
   complex(kind=rx), intent(out) :: w
   real(kind=rx), intent(in), dimension(:) :: ra

   integer :: i

   w = 0.0
   do i = n, 1, -1
      w = (w + ra(i))*z
   end do

end subroutine tay

subroutine vtay(z, ra, n, w)
   integer, intent(in) :: n
   complex(kind=rx), intent(in), dimension(:) :: z
   complex(kind=rx), intent(out), dimension(:) :: w
   real(kind=rx), intent(in), dimension(:) :: ra

   integer :: i

   w = 0.0
   do i = n, 1, -1
      w = (w + ra(i))*z
   end do

end subroutine vtay


!----------------------------------------------------------------------------
!   R.J.Purser, National Meteorological Center, Washington D.C.  1994
!                   SUBROUTINE  TAYD
!  Evaluate the complex function W of Z whose real
!  Taylor series coefficients are RA, together with its derivative WD.
!
!  --> Z    function argument (complex)
!  --> RA   Taylor coefficients (real)
!  --> N    number of coefficients (starting with the linear term)
!  <-- W    Taylor-series approximation of the function (complex)
!  <-- WD   Taylor series approximation of the derivative of the function W
!----------------------------------------------------------------------------

subroutine tayd(z, ra, n, w, wd)
   integer, intent(in) :: n
   complex(kind=rx), intent(in) :: z
   complex(kind=rx), intent(out) :: w
   complex(kind=rx), intent(out) :: wd
   real(kind=rx), intent(in), dimension(:) :: ra

   integer :: i

   w  = 0.0
   wd = 0.0
   do i = n, 1, -1
      w = (w + ra(i))*z
      wd = z*wd + i*ra(i)
   end do

end subroutine tayd



!---------------------------------------------------------------------------
!   R.J.Purser, National Meteorological Center, Washington D.C.  1994
!                   SUBROUTINE  TAYDD
!  Evaluate the complex function W of Z whose real
!  Taylor series coefficients are RA, together with its derivative WD.
!  and second derivative WDD
!
!  --> Z    function argument (complex)
!  --> RA   Taylor coefficients (real)
!  --> N    number of coefficients (starting with the linear term)
!  <-- W    Taylor-series approximation of the function (complex)
!  <-- WD   Taylor series approximation of the derivative of the function W
!  <-- WDD  Taylor series approximation of the derivative of the function WD
!---------------------------------------------------------------------------

subroutine taydd(z, ra, n, w, wd, wdd)
   integer, intent(in) :: n
   complex(kind=rx), intent(in) :: z
   complex(kind=rx), intent(out) :: w
   complex(kind=rx), intent(out) :: wd
   complex(kind=rx), intent(out) :: wdd
   real(kind=rx), intent(in), dimension(:) :: ra

   integer :: i

   w    = 0.0
   wd   = 0.0
   wdd  = 0.0
   do i = n, 1, -1
      w = (w + ra(i))*z
      wd = z*wd + i*ra(i)
   end do
   do i = n, 2, -1
      wdd = z*wdd + i*(i - 1)*ra(i)
   end do

end subroutine taydd

subroutine vtaydd(z, ra, n, w, wd, wdd)
   integer, intent(in) :: n
   complex(kind=rx), intent(in), dimension(:) :: z
   complex(kind=rx), intent(out), dimension(:)  :: w
   complex(kind=rx), intent(out), dimension(:)  :: wd
   complex(kind=rx), intent(out), dimension(:)  :: wdd
   real(kind=rx), intent(in), dimension(:) :: ra

   integer :: i

   w    = 0.0
   wd   = 0.0
   wdd  = 0.0
   do i = n, 1, -1
      w = (w + ra(i))*z
      wd = z*wd + i*ra(i)
   end do
   do i = n, 2, -1
      wdd = z*wdd + i*(i - 1)*ra(i)
   end do

end subroutine vtaydd

!--------------------------------------------------------------------------
!   R.J.Purser, National Meteorological Center, Washington D.C.  1993
!                   SUBROUTINE LUFM
!  perform l-u decomposition of square matrix a in place with
!  partial pivoting
!  For DOUBLE PRECISION version see DLUFM
!
!  --> a    square matrix to be factorized
!  <-- ipiv array encoding the pivoting sequence
!  <-- d    indicator for possible sign change of determinant
!  --> m    degree of (active part of) a
!  --> na   first fortran dimension of a
!
!--------------------------------------------------------------------------
subroutine lufm(a, ipiv, d)
   real(kind=rx), intent(inout), dimension(:,:) :: a
   integer, intent(out), dimension(:) :: ipiv
   real(kind=rx), intent(out) :: d

   integer :: m
   integer :: j, jp, iquad, i, k, jm
   real(kind=rx)    :: abig, aa, t, ajj, ajji, aij

   m = size(a,1)
   d = 1.0
   ipiv(m) = m
   do j = 1, m - 1
      jp = j + 1
      abig = abs(a(j,j))
      iquad = j
      do i = jp, m
         aa = abs(a(i,j))
         if (aa <= abig) then
            cycle
         endif
         iquad = i
         abig = aa
      end do
      !  swap rows, recording changed sign of determinant
      ipiv(j) = iquad
      if (iquad /= j) then
         d = -d
         do k = 1, m
            t = a(j,k)
            a(j,k) = a(iquad,k)
            a(iquad,k) = t
         end do
      endif
      ajj = a(j,j)
      if (ajj == 0.0) then
         jm = j - 1
         print *, "failure in lufact: matrix singular, rank=", jm
         stop
      endif
      ajji = 1.0/ajj
      do i = jp, m
         aij = ajji*a(i,j)
         a(i,j) = aij
         a(i,jp:m) = a(i,jp:m) - aij*a(j,jp:m)
      end do
   end do
   return
end subroutine lufm



!---------------------------------------------------------------------------
!   R.J.Purser, National Meteorological Center, Washington D.C.  1993
!                   SUBROUTINE INVMM
!  invert matrix in place using the l-u decomposition method
!  For DOUBLE PRECISION version see DINVMM
!
!  <--> a   matrix
!  --> m    degree of (active part of) b and a
!  --> nb   first fortran dimension of b
!  --> na   first fortran dimension of a
!
!   LIMITATION:
!    ipiv is an index array, internal to this array, encoding the
!    pivoting sequence used. It is given a fortran dimension of NN=500
!    in the parameter statement below. If the order of the linear system
!    exceeds 500, increase this parameter accordingly
!
!----------------------------------------------------------------------------
subroutine invmm(a)
   real(kind=rx), intent(inout), dimension(:,:) :: a

   integer, dimension(size(a,1)) :: ipiv
   integer :: m
   integer :: j, i, l
   real(kind=rx) :: d, s, t
!-----------------------------------------------

!  Check it's a square matrix
   if ( size(a, 1) /= size(a, 2) ) then
      print*, " Can't calculate inverse of non-square matrix "
      print*, " Shape is ", shape(a)
      stop
   end if

   m = size(a, 1)
   call lufm (a, ipiv, d)
!  invert u in place:
   do i = 1, m
      a(i,i) = 1.0/a(i,i)
   end do
   do i = 1, m - 1
      do j = i + 1, m
         s = 0.0
         s = -sum(a(i,i:j-1)*a(i:j-1,j))
         a(i,j) = a(j,j)*s
      end do
   end do
!  invert l in place assuming implicitly diagonal elements of unity
   do j = 1, m - 1
      do i = j + 1, m
         s = -a(i,j)
         s = s - sum(a(i,j+1:i-1)*a(j+1:i-1,j))
         a(i,j) = s
      end do
   end do
!  form the product of u**-1 and l**-1 in place
   do j = 1, m - 1
      do i = 1, j
         s = a(i,j)
         s = s + sum(a(i,j+1:m)*a(j+1:m,j))
         a(i,j) = s
      end do
      do i = j + 1, m
         s = 0.0
         s = sum(a(i,i:m)*a(i:m,j))
         a(i,j) = s
      end do
   end do
!  permute columns according to ipiv
   do j = m - 1, 1, -1
      l = ipiv(j)
      do i = 1, m
         t = a(i,j)
         a(i,j) = a(i,l)
         a(i,l) = t
      end do
   end do
end subroutine invmm

end module jimcc_m
