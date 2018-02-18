
subroutine C3D8(mesh, node, stiff)
  use util
  implicit none
  type(meshdef) :: mesh
  integer(kint) :: i, n, ndof
  integer(kint) :: ir1, ir2, ir3
  integer(kint) :: ngsp
  real(kdouble) :: node(3,8), stiff(24,24)
  real(kdouble) :: r(3), gsp(2), wg, det
  real(kdouble) :: B(6,24), D(6,6), dndx(8,3)

  ngsp   = 2
  gsp(1) = -0.577350269189626d0
  gsp(2) =  0.577350269189626d0
  wg     = 1.0d0
  stiff  = 0.0d0

  do ir1=1,ngsp
    do ir2=1,ngsp
      do ir3=1,ngsp
        r(1) = gsp(ir1)
        r(2) = gsp(ir2)
        r(3) = gsp(ir3)
        call get_global_deriv(node, r, dndx, det)
        call Bmat(B, dndx)
        call Dmat(mesh, D)
        call Kmat(stiff, D, B, wg, det)
      enddo
    enddo
  enddo

end subroutine C3D8

subroutine Shapefunc_C3D8(local, func)
  use util
  implicit none
  real(kdouble) ::  local(3), func(8)

  func(1) = 0.125d0*(1.d0-local(1))*(1.d0-local(2))*(1.d0-local(3))
  func(2) = 0.125d0*(1.d0+local(1))*(1.d0-local(2))*(1.d0-local(3))
  func(3) = 0.125d0*(1.d0+local(1))*(1.d0+local(2))*(1.d0-local(3))
  func(4) = 0.125d0*(1.d0-local(1))*(1.d0+local(2))*(1.d0-local(3))
  func(5) = 0.125d0*(1.d0-local(1))*(1.d0-local(2))*(1.d0+local(3))
  func(6) = 0.125d0*(1.d0+local(1))*(1.d0-local(2))*(1.d0+local(3))
  func(7) = 0.125d0*(1.d0+local(1))*(1.d0+local(2))*(1.d0+local(3))
  func(8) = 0.125d0*(1.d0-local(1))*(1.d0+local(2))*(1.d0+local(3))
end subroutine Shapefunc_C3D8

subroutine Shapefunc_deriv_C3D8(local, func)
  use util
  implicit none
  real(kdouble) :: local(3), func(8,3)

  func(1,1) = -0.125d0*(1.d0-local(2))*(1.d0-local(3))
  func(2,1) =  0.125d0*(1.d0-local(2))*(1.d0-local(3))
  func(3,1) =  0.125d0*(1.d0+local(2))*(1.d0-local(3))
  func(4,1) = -0.125d0*(1.d0+local(2))*(1.d0-local(3))
  func(5,1) = -0.125d0*(1.d0-local(2))*(1.d0+local(3))
  func(6,1) =  0.125d0*(1.d0-local(2))*(1.d0+local(3))
  func(7,1) =  0.125d0*(1.d0+local(2))*(1.d0+local(3))
  func(8,1) = -0.125d0*(1.d0+local(2))*(1.d0+local(3))

  func(1,2) = -0.125d0*(1.d0-local(1))*(1.d0-local(3))
  func(2,2) = -0.125d0*(1.d0+local(1))*(1.d0-local(3))
  func(3,2) =  0.125d0*(1.d0+local(1))*(1.d0-local(3))
  func(4,2) =  0.125d0*(1.d0-local(1))*(1.d0-local(3))
  func(5,2) = -0.125d0*(1.d0-local(1))*(1.d0+local(3))
  func(6,2) = -0.125d0*(1.d0+local(1))*(1.d0+local(3))
  func(7,2) =  0.125d0*(1.d0+local(1))*(1.d0+local(3))
  func(8,2) =  0.125d0*(1.d0-local(1))*(1.d0+local(3))

  func(1,3) = -0.125d0*(1.d0-local(1))*(1.d0-local(2))
  func(2,3) = -0.125d0*(1.d0+local(1))*(1.d0-local(2))
  func(3,3) = -0.125d0*(1.d0+local(1))*(1.d0+local(2))
  func(4,3) = -0.125d0*(1.d0-local(1))*(1.d0+local(2))
  func(5,3) =  0.125d0*(1.d0-local(1))*(1.d0-local(2))
  func(6,3) =  0.125d0*(1.d0+local(1))*(1.d0-local(2))
  func(7,3) =  0.125d0*(1.d0+local(1))*(1.d0+local(2))
  func(8,3) =  0.125d0*(1.d0-local(1))*(1.d0+local(2))
end subroutine Shapefunc_deriv_C3D8

subroutine get_global_deriv(node, r, dndx, det)
  use util
  implicit none
  real(kdouble) :: node(3,8), r(3), dndx(8,3), deriv(8,3), xj(3,3), inv(3,3), det, detinv

!write(*,*)"node"
!write(*,"(3e12.5)")node

  call Shapefunc_deriv_C3D8(r, deriv)

  xj = matmul( node, deriv )

  det = xj(1,1) * xj(2,2) * xj(3,3) &
      + xj(2,1) * xj(3,2) * xj(1,3) &
      + xj(3,1) * xj(1,2) * xj(2,3) &
      - xj(3,1) * xj(2,2) * xj(1,3) &
      - xj(2,1) * xj(1,2) * xj(3,3) &
      - xj(1,1) * xj(3,2) * xj(2,3)

  if( det==0.d0 ) stop "determinant = 0.0"

  detinv = 1.d0/det
  inv(1,1) = detinv * ( xj(2,2)*xj(3,3) - xj(3,2)*xj(2,3))
  inv(1,2) = detinv * (-xj(1,2)*xj(3,3) + xj(3,2)*xj(1,3))
  inv(1,3) = detinv * ( xj(1,2)*xj(2,3) - xj(2,2)*xj(1,3))
  inv(2,1) = detinv * (-xj(2,1)*xj(3,3) + xj(3,1)*xj(2,3))
  inv(2,2) = detinv * ( xj(1,1)*xj(3,3) - xj(3,1)*xj(1,3))
  inv(2,3) = detinv * (-xj(1,1)*xj(2,3) + xj(2,1)*xj(1,3))
  inv(3,1) = detinv * ( xj(2,1)*xj(3,2) - xj(3,1)*xj(2,2))
  inv(3,2) = detinv * (-xj(1,1)*xj(3,2) + xj(3,1)*xj(1,2))
  inv(3,3) = detinv * ( xj(1,1)*xj(2,2) - xj(2,1)*xj(1,2))

  dndx = matmul( deriv, inv )
end subroutine get_global_deriv

subroutine Bmat(B, dndx)
  use util
  implicit none
  type(meshdef) :: mesh
  integer(kint) :: i, i1, i2, i3
  real(kdouble) :: B(6,24), dndx(8,3)

  B = 0.0d0

  do i = 1,8
    i1 = 3*i-2
    i2 = 3*i-1
    i3 = 3*i
    B(1,i1) = dndx(i,1)
    B(2,i2) = dndx(i,2)
    B(3,i3) = dndx(i,3)
    B(4,i1) = dndx(i,2)
    B(4,i2) = dndx(i,1)
    B(5,i2) = dndx(i,3)
    B(5,i3) = dndx(i,2)
    B(6,i1) = dndx(i,3)
    B(6,i3) = dndx(i,1)
  enddo
end subroutine Bmat

subroutine Dmat(mesh, D)
  use util
  implicit none
  type(meshdef) :: mesh
  real(kdouble) :: D(6,6), E, mu, g

  E = mesh%E
  mu= mesh%mu

  D = 0.0d0
  g = E / ((1.0d0+mu) * (1.0d0-2.0d0*mu))

  D(1,1) = g*(1.0d0-mu)
  D(1,2) = g*mu
  D(1,3) = g*mu
  D(2,1) = g*mu
  D(2,2) = g*(1.0d0-mu)
  D(2,3) = g*mu
  D(3,1) = g*mu
  D(3,2) = g*mu
  D(3,3) = g*(1.0d0-mu)
  D(4,4) = g/(2.0d0-mu)
  D(5,5) = g/(2.0d0-mu)
  D(6,6) = g/(2.0d0-mu)
end subroutine Dmat

subroutine Kmat(stiff, D, B, wg, det)
  use util
  implicit none
  integer(kint) :: i, j
  real(kdouble) :: stiff(24,24), D(6,6), B(6,24), DB(6,24), wg, det

  DB = MATMUL(D, B)
  forall(i=1:24, j=1:24)
    stiff(i,j) = stiff(i,j) + DOT_PRODUCT(B(:,i), DB(:,j)) * wg * det
  end forall
end subroutine Kmat
