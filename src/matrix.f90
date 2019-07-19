
subroutine stiff_matrix(mesh)
  use util
  implicit none
  type(meshdef) :: mesh
  integer(kint) :: i, icel
  integer(kint) :: elem(8)
  real(kdouble) :: stiff(24,24)

  mesh%A = 0.0d0
  do icel = 1,mesh%nelem
    do i = 1, 8
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
  integer(kint) :: i, j, k, in, jn, jS, jE
  integer(kint) :: elem(8)
  real(kdouble) :: stiff(24,24)

  do i = 1, 8
    in = elem(i)
    jS = mesh%index(in-1) + 1
    jE = mesh%index(in)
    aa:do j = 1, 8
      do k = jS, jE
        jn = mesh%item(k)
        if(jn == elem(j))then
          mesh%A(9*k-8) = mesh%A(9*k-8) + stiff(3*j-2,3*i-2)
          mesh%A(9*k-7) = mesh%A(9*k-7) + stiff(3*j-2,3*i-1)
          mesh%A(9*k-6) = mesh%A(9*k-6) + stiff(3*j-2,3*i  )
          mesh%A(9*k-5) = mesh%A(9*k-5) + stiff(3*j-1,3*i-2)
          mesh%A(9*k-4) = mesh%A(9*k-4) + stiff(3*j-1,3*i-1)
          mesh%A(9*k-3) = mesh%A(9*k-3) + stiff(3*j-1,3*i  )
          mesh%A(9*k-2) = mesh%A(9*k-2) + stiff(3*j  ,3*i-2)
          mesh%A(9*k-1) = mesh%A(9*k-1) + stiff(3*j  ,3*i-1)
          mesh%A(9*k  ) = mesh%A(9*k  ) + stiff(3*j  ,3*i  )
          cycle aa
        endif
      enddo
    enddo aa
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
  do i = 1, mesh%ncload
    in  = mesh%icload(1, i)
    dof = mesh%icload(2, i)
    val = mesh%cload(i)
    if(ndof < dof) stop "*** error: 3 < dof"
    mesh%f(ndof*(in-1) + dof) = val
  enddo
end subroutine load_condition

subroutine get_RHS(mesh)
  use util
  implicit none
  type(meshdef) :: mesh
  integer(kint) :: i

  do i = 1, mesh%nnode
    mesh%B(3*i-2) = mesh%f(3*i-2) - mesh%q(3*i-2)
    mesh%B(3*i-1) = mesh%f(3*i-1) - mesh%q(3*i-1)
    mesh%B(3*i  ) = mesh%f(3*i  ) - mesh%q(3*i  )
  enddo
end subroutine get_RHS

subroutine bound_condition(mesh)
  use util
  implicit none
  type(meshdef) :: mesh
  integer(kint) :: i, in, dof, j, jS, jE, k, jn, nb
  real(kdouble) :: val

  !boundary
  do nb = 1, mesh%nbound
    in  = mesh%ibound(1, nb)
    dof = mesh%ibound(2, nb)
    if(ndof < dof) stop "*** error: 3 < dof"
    !val = mesh%bound(i)
    !mesh%B(3*in-3+dof) = val
    do i = 1, mesh%nnode
      jS = mesh%index(i-1) + 1
      jE = mesh%index(i)
      bb:do j = jS, jE
        jn = mesh%item(j)
        if(in < jn) cycle bb
        if(in == jn)then
          if(dof == 1)then
            mesh%A(9*j-8) = 0.0d0
            mesh%A(9*j-5) = 0.0d0
            mesh%A(9*j-2) = 0.0d0
          elseif(dof == 2)then
            mesh%A(9*j-7) = 0.0d0
            mesh%A(9*j-4) = 0.0d0
            mesh%A(9*j-1) = 0.0d0
          elseif(dof == 3)then
            mesh%A(9*j-6) = 0.0d0
            mesh%A(9*j-3) = 0.0d0
            mesh%A(9*j  ) = 0.0d0
          endif
        endif
      enddo bb
    enddo

    jS = mesh%index(in-1) + 1
    jE = mesh%index(in)
    do j = jS, jE
      if(dof == 1)then
        mesh%A(9*j-8) = 0.0d0
        mesh%A(9*j-7) = 0.0d0
        mesh%A(9*j-6) = 0.0d0
      elseif(dof == 2)then
        mesh%A(9*j-5) = 0.0d0
        mesh%A(9*j-4) = 0.0d0
        mesh%A(9*j-3) = 0.0d0
      elseif(dof == 3)then
        mesh%A(9*j-2) = 0.0d0
        mesh%A(9*j-1) = 0.0d0
        mesh%A(9*j  ) = 0.0d0
      endif
      jn = mesh%item(j)
      if(jn == in)then
        if(dof == 1)then
          mesh%A(9*j-8) = 1.0d0
        elseif(dof == 2)then
          mesh%A(9*j-4) = 1.0d0
        elseif(dof == 3)then
          mesh%A(9*j  ) = 1.0d0
        endif
      endif
    enddo
  enddo
end subroutine bound_condition

