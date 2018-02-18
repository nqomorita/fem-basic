
subroutine stiff_matrix(mesh)
  use util
  implicit none
  type(meshdef) :: mesh
  integer(kint) :: i, in, j, n, ndof, nndof
  integer(kint) :: icel
  integer(kint) :: elem(8)
  real(kdouble) :: node(3,8)
  real(kdouble) :: stiff(24,24)

  mesh%ndof = 3

  n    = mesh%nnode
  ndof = mesh%ndof
  nndof= n*ndof
  allocate(mesh%A(nndof,nndof))
  mesh%A = 0.0d0

  do icel=1,mesh%nelem
    do i=1,8
      elem(i) = mesh%elem(i,icel)
      in = elem(i)
      node(1,i) = mesh%node(1,in)
      node(2,i) = mesh%node(2,in)
      node(3,i) = mesh%node(3,in)
    enddo
    call C3D8(mesh, node, stiff)
    call merge(mesh, elem, stiff)
  enddo

end subroutine stiff_matrix

subroutine merge(mesh, elem, stiff)
  use util
  implicit none
  type(meshdef) :: mesh
  integer(kint) :: i, in, j, jn, n, ndof
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

subroutine bound_condition(mesh)
  use util
  implicit none
  type(meshdef) :: mesh
  integer(kint) :: i, in, dof, j, jn, n, ndof, nndof
  real(kdouble) :: val

  n    = mesh%nnode
  ndof = mesh%ndof
  nndof= n*ndof
  allocate(mesh%B(n*ndof))
  mesh%B = 0.0d0

  !boundary
  do i=1,mesh%nbound
    in  = mesh%ibound(1, i)
    dof = mesh%ibound(2, i)
    !val = mesh%bound(i)
    !mesh%B(3*in-3+dof) = val
    jn  = 3*in-3+dof
    do j=1,nndof
      mesh%A(jn, j) = 0.0d0
      mesh%A(j, jn) = 0.0d0
    enddo
    mesh%A(jn, jn) = 1.0d0
  enddo

  !cload
  do i=1,mesh%ncload
    in  = mesh%icload(1, i)
    dof = mesh%icload(2, i)
    val = mesh%cload(i)
    mesh%B(3*in-3+dof) = val
  enddo
end subroutine bound_condition

subroutine stress_update(mesh)
  use util
  implicit none
  type(meshdef) :: mesh
  integer(kint) :: i, n, ndof

end subroutine stress_update
