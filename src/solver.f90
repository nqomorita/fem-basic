
subroutine solver(mesh, monolis)
  use util
  use mod_monolis_util
  use mod_monolis_solve
  implicit none
  type(meshdef) :: mesh
  type(monolis_structure) :: monolis

  monolis%MAT%N = mesh%nnode
  monolis%MAT%NP = mesh%nnode
  monolis%MAT%NZ = mesh%index(mesh%nnode)
  monolis%MAT%NDOF = 3
  monolis%MAT%index => mesh%index
  monolis%MAT%item => mesh%item
  monolis%MAT%A => mesh%A
  monolis%MAT%B => mesh%B
  monolis%MAT%X => mesh%X

  call monolis_solve(monolis%PRM, monolis%COM, monolis%MAT)
  !call LU(ndof*mesh%nnode, mesh%A, mesh%B, mesh%X)
  !call residual(ndof*mesh%nnode, mesh%A, mesh%B, mesh%X)
write(*,"(1p3e12.4)")mesh%X

  nullify(monolis%MAT%index)
  nullify(monolis%MAT%item)
  nullify(monolis%MAT%A)
  nullify(monolis%MAT%B)
  nullify(monolis%MAT%X)
end subroutine solver

subroutine LU(nndof, Ain, B, X)
  use util
  implicit none
  integer(kint) :: i, j, k, nndof
  real(kdouble) :: Ain(nndof, nndof), B(nndof), X(nndof)
  real(kdouble), allocatable :: A(:,:)
  real(kdouble) :: t

  allocate(A(nndof, nndof))
  do i=1, nndof
    do j=1, nndof
      A(j,i) = Ain(j,i)
    enddo
  enddo

  do i=1, nndof
    A(i,i) = 1.0d0/A(i,i)
    do j=i+1, nndof
      A(i,j) = A(i,j) * A(i,i)
      do k=i+1, nndof
        A(k,j) = A(k,j) - A(i,j) * A(k,i)
      enddo
    enddo
  enddo

  X = B

  do i=1,nndof
    t = 0.0d0
    do j=1,i-1
      t = t + A(j,i) * X(j)
    enddo
    X(i) = X(i) - t
  enddo

  do i=nndof, 1, -1
    t = 0.0d0
    do j=nndof, i+1, -1
      t = t + A(j,i) * X(j)
    enddo
    X(i) = X(i) - t
    X(i) = X(i) * A(i,i)
  enddo

  deallocate(A)
end subroutine LU

subroutine residual(nndof, A, B, X)
  use util
  implicit none
  integer(kint) :: i, j, nndof
  real(kdouble) :: A(nndof, nndof), B(nndof), X(nndof), res, tmp
  real(kdouble), allocatable :: Y(:)

  allocate(Y(nndof))

  do i=1, nndof
    Y(i) = 0.0d0
    do j=1, nndof
      Y(i) = Y(i) + A(j,i) * X(j)
    enddo
  enddo

  res = 0.0d0
  do i=1, nndof
    tmp = B(i) - Y(i)
    res = res + tmp * tmp
  enddo
  res = dsqrt(res)

  write(*,"(a, 1pe12.5)")"  ** solver residual: ", res

  deallocate(Y)
end subroutine residual
