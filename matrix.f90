
subroutine stiff_matrix(mesh)
  use util
  implicit none
  type(meshdef) :: mesh
  integer(kint) :: i, in, j
  integer(kint) :: icel
  integer(kint) :: elem(8)
  real(kdouble) :: node(3,8)
  real(kdouble) :: stiff(24,24)

  do icel=1,mesh%nelem
    do i=1,8
      elem(i) = mesh%elem(i,icel)
      in = elem(i)
      node(1,i) = mesh%node(1,in)
      node(2,i) = mesh%node(2,in)
      node(3,i) = mesh%node(3,in)
    enddo
    call C3D8_stiff(mesh, node, stiff)
    call merge(mesh, elem, stiff)
  enddo
end subroutine stiff_matrix

subroutine merge(mesh, elem, stiff)
  use util
  implicit none
  type(meshdef) :: mesh
  integer(kint) :: i, in, j, jn, n
  integer(kint) :: elem(8)
  real(kdouble) :: stiff(24,24)

  do i=1, 8
    in = elem(i)
    do j=1, 8
      jn = elem(j)
      mesh%A(3*jn-2, 3*in-2) = mesh%A(3*jn-2, 3*in-2) + stiff(3*j-2, 3*i-2)
      mesh%A(3*jn-2, 3*in-1) = mesh%A(3*jn-2, 3*in-1) + stiff(3*j-2, 3*i-1)
      mesh%A(3*jn-2, 3*in  ) = mesh%A(3*jn-2, 3*in  ) + stiff(3*j-2, 3*i  )
      mesh%A(3*jn-1, 3*in-2) = mesh%A(3*jn-1, 3*in-2) + stiff(3*j-1, 3*i-2)
      mesh%A(3*jn-1, 3*in-1) = mesh%A(3*jn-1, 3*in-1) + stiff(3*j-1, 3*i-1)
      mesh%A(3*jn-1, 3*in  ) = mesh%A(3*jn-1, 3*in  ) + stiff(3*j-1, 3*i  )
      mesh%A(3*jn  , 3*in-2) = mesh%A(3*jn  , 3*in-2) + stiff(3*j  , 3*i-2)
      mesh%A(3*jn  , 3*in-1) = mesh%A(3*jn  , 3*in-1) + stiff(3*j  , 3*i-1)
      mesh%A(3*jn  , 3*in  ) = mesh%A(3*jn  , 3*in  ) + stiff(3*j  , 3*i  )
    enddo
  enddo
end subroutine merge

subroutine load_condition(mesh)
  use util
  implicit none
  type(meshdef) :: mesh
  integer(kint) :: i, in, dof, j, jn
  real(kdouble) :: val

  !cload
  do i=1,mesh%ncload
    in  = mesh%icload(1, i)
    dof = mesh%icload(2, i)
    val = mesh%cload(i)
    if(3 < dof) stop "*** error: 3 < dof"
    mesh%f(3*in-3+dof) = val
  enddo
end subroutine load_condition

subroutine bound_condition(mesh)
  use util
  implicit none
  type(meshdef) :: mesh
  integer(kint) :: i, in, dof, j, jn
  real(kdouble) :: val

  !boundary
  do i=1,mesh%nbound
    in  = mesh%ibound(1, i)
    dof = mesh%ibound(2, i)
    if(3 < dof) stop "*** error: 3 < dof"
    !val = mesh%bound(i)
    !mesh%B(3*in-3+dof) = val
    jn  = 3*in-3+dof
    do j=1,3*mesh%nnode
      mesh%A(jn, j) = 0.0d0
      mesh%A(j, jn) = 0.0d0
    enddo
    mesh%A(jn, jn) = 1.0d0

    mesh%B(jn) = 0.0d0
  enddo
end subroutine bound_condition

subroutine stress_update(mesh)
  use util
  implicit none
  type(meshdef) :: mesh
  integer(kint) :: i, j, in, icel, n, ndof, elem(8)
  real(kdouble) :: node(3,8), u(3,8), r(3), det
  real(kdouble) :: func(8,8), inv(8,8)
  real(kdouble) :: nstrain(8,6), nstress(8,6)
  real(kdouble) :: estrain(6),   estress(6)
  real(kdouble) :: q_elem(24)
  integer(kint), allocatable :: inode(:)

  allocate(inode(mesh%nnode))
  inode = 0

  do i=1,8
    call C3D8_integral_point(i, r)
    call C3D8_shapefunc(r, func(i,1:8))
  enddo
  call get_inverse_matrix(8, func, inv)

  do icel=1,mesh%nelem
    do i=1,8
      in = mesh%elem(i,icel)
      u(1,i) = mesh%u(3*i-2) + mesh%du(3*i-2) + mesh%X(3*in-2)
      u(2,i) = mesh%u(3*i-1) + mesh%du(3*i-1) + mesh%X(3*in-1)
      u(3,i) = mesh%u(3*i  ) + mesh%du(3*i  ) + mesh%X(3*in  )
      node(1,i) = mesh%node(1,in)
      node(2,i) = mesh%node(2,in)
      node(3,i) = mesh%node(3,in)
    enddo
    call C3D8_update(mesh, icel, node, u, q_elem)
    call C3D8_get_nodal_values(mesh, icel, inv, nstrain, nstress, estrain, estress)

    do i=1,8
      in = mesh%elem(i,icel)
      inode(in) = inode(in) + 1
      do j=1,6
        mesh%nstrain(j,in) = mesh%nstrain(j,in) + estrain(j)
        mesh%nstress(j,in) = mesh%nstress(j,in) + estress(j)
      enddo
      mesh%q(3*in-2) = q_elem(3*i-2)
      mesh%q(3*in-1) = q_elem(3*i-1)
      mesh%q(3*in  ) = q_elem(3*i  )
    enddo

    do j=1,6
      mesh%estrain(j,icel) = estrain(j)
      mesh%estress(j,icel) = estress(j)
    enddo
  enddo

  do i=1,mesh%nnode
    det = 1.0d0/dble(inode(i))
    do j=1,6
      mesh%nstrain(j,i) = mesh%nstrain(j,i) * det
      mesh%nstress(j,i) = mesh%nstress(j,i) * det
    enddo
  enddo

  do i=1,mesh%nnode
    call get_mises(mesh%nstress(1:6,i),  mesh%nmises(i))
  enddo
  do i=1,mesh%nelem
    call get_mises(mesh%estress(1:6,i),  mesh%emises(i))
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
