
subroutine C3D8_stiff(mesh, node, stiff)
  use util
  implicit none
  type(meshdef) :: mesh
  integer(kint) :: i, n, ndof
  integer(kint) :: ir1, ir2, ir3
  real(kdouble) :: node(3,8), stiff(24,24)
  real(kdouble) :: r(3), wg, det
  real(kdouble) :: B(6,24), D(6,6), dndx(8,3)
  real(kdouble) :: stress(6)

  wg     = 1.0d0
  stiff  = 0.0d0

  do i=1,8
    call C3D8_integral_point(i, r)
    call C3D8_get_global_deriv(node, r, dndx, det)
    call C3D8_Bmat(node, B, dndx)
    call C3D8_Dmat(mesh, D)
    call C3D8_Kmat(mesh, stress, D, B, wg, det, stiff)
  enddo

end subroutine C3D8_stiff

subroutine C3D8_integral_point(i, r)
  use util
  implicit none
  integer(kint) :: i
  real(kdouble) :: gsp(3,8), r(3)

  gsp(1,1) = -0.577350269189626d0; gsp(2,1) = -0.577350269189626d0; gsp(3,1) = -0.577350269189626d0
  gsp(1,2) =  0.577350269189626d0; gsp(2,2) = -0.577350269189626d0; gsp(3,2) = -0.577350269189626d0
  gsp(1,3) = -0.577350269189626d0; gsp(2,3) =  0.577350269189626d0; gsp(3,3) = -0.577350269189626d0
  gsp(1,4) =  0.577350269189626d0; gsp(2,4) =  0.577350269189626d0; gsp(3,4) = -0.577350269189626d0
  gsp(1,5) = -0.577350269189626d0; gsp(2,5) = -0.577350269189626d0; gsp(3,5) =  0.577350269189626d0
  gsp(1,6) =  0.577350269189626d0; gsp(2,6) = -0.577350269189626d0; gsp(3,6) =  0.577350269189626d0
  gsp(1,7) = -0.577350269189626d0; gsp(2,7) =  0.577350269189626d0; gsp(3,7) =  0.577350269189626d0
  gsp(1,8) =  0.577350269189626d0; gsp(2,8) =  0.577350269189626d0; gsp(3,8) =  0.577350269189626d0

  r(1) = gsp(1,i)
  r(2) = gsp(2,i)
  r(3) = gsp(3,i)
end subroutine C3D8_integral_point

subroutine C3D8_shapefunc(local, func)
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
end subroutine C3D8_shapefunc

subroutine C3D8_shapefunc_deriv(local, func)
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
end subroutine C3D8_shapefunc_deriv

subroutine C3D8_get_inverse_matrix(xj, inv, det)
  use util
  implicit none
  real(kdouble) :: xj(3,3), inv(3,3), det, detinv

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
end subroutine C3D8_get_inverse_matrix

subroutine C3D8_get_global_deriv(node, r, dndx, det)
  use util
  implicit none
  real(kdouble) :: node(3,8), r(3), dndx(8,3), deriv(8,3)
  real(kdouble) :: xj(3,3), inv(3,3), det, detinv

!write(*,*)"node"
!write(*,"(3e12.5)")node

  call C3D8_shapefunc_deriv(r, deriv)
  xj = matmul( node, deriv )
  call C3D8_get_inverse_matrix(xj, inv, det)
  dndx = matmul( deriv, inv )

end subroutine C3D8_get_global_deriv

subroutine C3D8_Bmat(node, B, dndx)
  use util
  implicit none
  type(meshdef) :: mesh
  integer(kint) :: i, i1, i2, i3
  real(kdouble) :: node(3,8), B(6,24), dndx(8,3), dudx(3,3)

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

  if(isNLGeom)then
    dudx = 0.0d0
    dudx = matmul(node, dndx)
    do i = 1, 8
      i1 = 3*i-2
      i2 = 3*i-1
      i3 = 3*i
      B(1,i1) = B(1,i1) + dudx(1,1)*dndx(i,1)
      B(1,i2) = B(1,i2) + dudx(2,1)*dndx(i,1)
      B(1,i3) = B(1,i3) + dudx(3,1)*dndx(i,1)
      B(2,i1) = B(2,i1) + dudx(1,2)*dndx(i,2)
      B(2,i2) = B(2,i2) + dudx(2,2)*dndx(i,2)
      B(2,i3) = B(2,i3) + dudx(3,2)*dndx(i,2)
      B(3,i1) = B(3,i1) + dudx(1,3)*dndx(i,3)
      B(3,i2) = B(3,i2) + dudx(2,3)*dndx(i,3)
      B(3,i3) = B(3,i3) + dudx(3,3)*dndx(i,3)
      B(4,i1) = B(4,i1) + dudx(1,2)*dndx(i,1) + dudx(1,1)*dndx(i,2)
      B(4,i2) = B(4,i2) + dudx(2,2)*dndx(i,1) + dudx(2,1)*dndx(i,2)
      B(4,i3) = B(4,i3) + dudx(3,2)*dndx(i,1) + dudx(3,1)*dndx(i,2)
      B(5,i1) = B(5,i1) + dudx(1,2)*dndx(i,3) + dudx(1,3)*dndx(i,2)
      B(5,i2) = B(5,i2) + dudx(2,2)*dndx(i,3) + dudx(2,3)*dndx(i,2)
      B(5,i3) = B(5,i3) + dudx(3,2)*dndx(i,3) + dudx(3,3)*dndx(i,2)
      B(6,i1) = B(6,i1) + dudx(1,3)*dndx(i,1) + dudx(1,1)*dndx(i,3)
      B(6,i2) = B(6,i2) + dudx(2,3)*dndx(i,1) + dudx(2,1)*dndx(i,3)
      B(6,i3) = B(6,i3) + dudx(3,3)*dndx(i,1) + dudx(3,1)*dndx(i,3)
    enddo
  endif
end subroutine C3D8_Bmat

subroutine C3D8_Dmat(mesh, D)
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
end subroutine C3D8_Dmat

subroutine C3D8_Kmat(mesh, stress, D, B, wg, det, stiff)
  use util
  implicit none
  type(meshdef) :: mesh
  integer(kint) :: i, j
  real(kdouble) :: stiff(24,24), D(6,6), B(6,24), DB(6,24), wg, det
  real(kdouble) :: stress(6), S(9,9), BN(9,24), SBN(9,24)

  DB = matmul(D, B)
  forall(i=1:24, j=1:24)
    stiff(i,j) = stiff(i,j) + DOT_PRODUCT(B(:,i), DB(:,j)) * wg * det
  end forall

  !if(isNLGeom)then
  !  stress(1:6) = gausses(LX)%stress
  !  do j = 1, nn+3
  !    BN(1, 3*j-2) = gderiv(j, 1)
  !    BN(2, 3*j-1) = gderiv(j, 1)
  !    BN(3, 3*j  ) = gderiv(j, 1)
  !    BN(4, 3*j-2) = gderiv(j, 2)
  !    BN(5, 3*j-1) = gderiv(j, 2)
  !    BN(6, 3*j  ) = gderiv(j, 2)
  !    BN(7, 3*j-2) = gderiv(j, 3)
  !    BN(8, 3*j-1) = gderiv(j, 3)
  !    BN(9, 3*j  ) = gderiv(j, 3)
  !  end do
  !  Smat(:, :) = 0.0D0
  !  do j = 1, 3
  !    Smat(j  , j  ) = stress(1)
  !    Smat(j  , j+3) = stress(4)
  !    Smat(j  , j+6) = stress(6)
  !    Smat(j+3, j  ) = stress(4)
  !    Smat(j+3, j+3) = stress(2)
  !    Smat(j+3, j+6) = stress(5)
  !    Smat(j+6, j  ) = stress(6)
  !    Smat(j+6, j+3) = stress(5)
  !    Smat(j+6, j+6) = stress(3)
  !  end do
  !  SBN(1:9, 1:(nn+3)*ndof) = matmul( Smat(1:9, 1:9), BN(1:9, 1:(nn+3)*ndof) )
  !  forall( i=1:(nn+3)*ndof, j=1:(nn+3)*ndof )
  !    tmpstiff(i, j) = tmpstiff(i, j)+dot_product( BN(:, i), SBN(:, j) )*wg
  !  end forall
  !endif
end subroutine C3D8_Kmat

subroutine C3D8_update(mesh, icel, node, u)
  use util
  implicit none
  type(meshdef) :: mesh
  integer(kint) :: i, j, icel
  real(kdouble) :: node(3,8), u(3,8), r(3), deriv(8,3), xj(3,3), D(6,6)
  real(kdouble) :: strain(6), stress(6), det

  do i=1,8
    call C3D8_integral_point(i, r)
    call C3D8_get_global_deriv(node, r, deriv, det)
    call C3D8_Dmat(mesh, D)
    xj = matmul(u, deriv)
    strain(1) = xj(1,1)
    strain(2) = xj(2,2)
    strain(3) = xj(3,3)
    strain(4) =(xj(1,2) + xj(2,1))
    strain(5) =(xj(2,3) + xj(3,2))
    strain(6) =(xj(3,1) + xj(1,3))
    stress = matmul(D, strain)
    mesh%gauss(i,icel)%strain = strain
    mesh%gauss(i,icel)%stress = stress
  enddo
end subroutine C3D8_update

subroutine  C3D8_get_nodal_values(mesh, icel, inv, nstrain, nstress, estrain, estress)
  use util
  implicit none
  type(meshdef) :: mesh
  integer(kint) :: i, j, k, icel
  real(kdouble) :: inv(8,8)
  real(kdouble) :: nstrain(8,6), nstress(8,6)
  real(kdouble) :: estrain(6), estress(6)

  nstrain  = 0.0d0
  nstress  = 0.0d0
  estrain  = 0.0d0
  estress  = 0.0d0

  do i=1,8
    do j=1,8
      do k=1,6
        nstrain(i,k) = nstrain(i,k) + inv(i,j) * mesh%gauss(j,icel)%strain(k)
        nstress(i,k) = nstress(i,k) + inv(i,j) * mesh%gauss(j,icel)%stress(k)
      enddo
    enddo
  enddo

  do i=1,8
    do j=1,6
      estrain(j) = estrain(j) + mesh%gauss(i,icel)%strain(j)
      estress(j) = estress(j) + mesh%gauss(i,icel)%stress(j)
    enddo
  enddo
  estrain = estrain/8.0d0
  estress = estress/8.0d0

end subroutine C3D8_get_nodal_values

