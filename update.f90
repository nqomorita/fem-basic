
subroutine delta_u_update(mesh)
  use util
  implicit none
  type(meshdef) :: mesh
  integer(kint) :: i

  do i=1,mesh%nnode
    mesh%du(3*i-2) = mesh%du(3*i-2) + mesh%X(3*i-2)
    mesh%du(3*i-1) = mesh%du(3*i-1) + mesh%X(3*i-1)
    mesh%du(3*i  ) = mesh%du(3*i  ) + mesh%X(3*i  )
  enddo
end subroutine delta_u_update

subroutine u_update(mesh)
  use util
  implicit none
  type(meshdef) :: mesh
  integer(kint) :: i

  do i=1,mesh%nnode
    mesh%u(3*i-2) = mesh%u(3*i-2) + mesh%du(3*i-2)
    mesh%u(3*i-1) = mesh%u(3*i-1) + mesh%du(3*i-1)
    mesh%u(3*i  ) = mesh%u(3*i  ) + mesh%du(3*i  )
  enddo
end subroutine u_update

subroutine stress_update(mesh)
  use util
  implicit none
  type(meshdef) :: mesh
  integer(kint) :: i, j, in, icel, n, ndof, elem(8)
  real(kdouble) :: node(3,8), u(3,8), r(3), tmp
  real(kdouble) :: func(8,8), inv(8,8)
  real(kdouble) :: nstrain(8,6), nstress(8,6)
  real(kdouble) :: estrain(6),   estress(6)
  real(kdouble) :: q(24)
  integer(kint), allocatable :: inode(:)

  do i=1,mesh%nnode
    do j=1,6
      mesh%nstrain(j,i) = 0.0d0
      mesh%nstress(j,i) = 0.0d0
    enddo
  enddo

  allocate(inode(mesh%nnode))
  inode = 0

  do i=1,8
    call C3D8_integral_point(i, r)
    call C3D8_shapefunc(r, func(i,1:8))
  enddo
  call get_inverse_matrix(8, func, inv)

  mesh%q = 0.0d0

  do icel=1,mesh%nelem
    call C3D8_update(mesh, icel, q)
    call C3D8_get_nodal_values(mesh, icel, inv, nstrain, nstress, estrain, estress)

    do i=1,8
      in = mesh%elem(i,icel)
      inode(in) = inode(in) + 1
      do j=1,6
        mesh%nstrain(j,in) = mesh%nstrain(j,in) + nstrain(i,j)
        mesh%nstress(j,in) = mesh%nstress(j,in) + nstress(i,j)
      enddo
      mesh%q(3*in-2) = mesh%q(3*in-2) + q(3*i-2)
      mesh%q(3*in-1) = mesh%q(3*in-1) + q(3*i-1)
      mesh%q(3*in  ) = mesh%q(3*in  ) + q(3*i  )
    enddo

    do j=1,6
      mesh%estrain(j,icel) = estrain(j)
      mesh%estress(j,icel) = estress(j)
    enddo
  enddo

  do i=1,mesh%nnode
    tmp = 1.0d0/dble(inode(i))
    do j=1,6
      mesh%nstrain(j,i) = mesh%nstrain(j,i) * tmp
      mesh%nstress(j,i) = mesh%nstress(j,i) * tmp
    enddo
  enddo

  do i=1,mesh%nnode
    call get_mises(mesh%nstress(1:6,i), mesh%nmises(i))
  enddo
  do i=1,mesh%nelem
    call get_mises(mesh%estress(1:6,i), mesh%emises(i))
  enddo

  deallocate(inode)
end subroutine stress_update

subroutine get_mises(s, mises)
  use util
  implicit none
  real(kdouble) :: mises, s(6)
  real(kdouble) :: s11, s22, s33, s12, s23, s13, ps, smises

  s11 = s(1)
  s22 = s(2)
  s33 = s(3)
  s12 = s(4)
  s23 = s(5)
  s13 = s(6)
  ps = (s11 + s22 + s33) / 3.0d0
  smises = 0.5d0 * ((s11-ps)**2 + (s22-ps)**2 + (s33-ps)**2) + s12**2 + s23**2 + s13**2
  mises  = dsqrt( 3.0d0 * smises )
end subroutine get_mises

