
subroutine stiff_matrix(mesh)
  use util
  implicit none
  type(meshdef) :: mesh
  integer(kint) :: i, icel
  integer(kint) :: elem(8)
  real(kdouble) :: stiff(24,24)

  do icel=1,mesh%nelem
    do i=1,8
      elem(i) = mesh%elem(i, icel)
    enddo
    call C3D8_stiff(mesh, icel, elem, stiff)
    call merge(mesh, elem, stiff)
  enddo
end subroutine stiff_matrix

subroutine merge(mesh, elem, stiff)
  use util
  implicit none
  type(meshdef) :: mesh
  integer(kint) :: i, in, j, jn
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
  integer(kint) :: i, in, dof
  real(kdouble) :: val

  mesh%f = 0.0d0

  !cload
  do i=1,mesh%ncload
    in  = mesh%icload(1, i)
    dof = mesh%icload(2, i)
    val = mesh%cload(i)
    if(3 < dof) stop "*** error: 3 < dof"
    mesh%f(3*in-3+dof) = val
  enddo
end subroutine load_condition

subroutine get_RHS(mesh)
  use util
  implicit none
  type(meshdef) :: mesh
  integer(kint) :: i

  do i=1,mesh%nnode
    mesh%B(3*i-2) = mesh%f(3*i-2) - mesh%q(3*i-2)
    mesh%B(3*i-1) = mesh%f(3*i-1) - mesh%q(3*i-1)
    mesh%B(3*i  ) = mesh%f(3*i  ) - mesh%q(3*i  )
  enddo
end subroutine get_RHS

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

