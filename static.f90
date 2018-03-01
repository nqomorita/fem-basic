
subroutine static(mesh)
  use util
  implicit none
  type(meshdef) :: mesh
  integer(kint) :: i

  call stiff_matrix(mesh)
  call bound_condition(mesh)

  do i=1,mesh%nnode
    mesh%B(3*i-2) = mesh%f(3*i-2)
    mesh%B(3*i-1) = mesh%f(3*i-1)
    mesh%B(3*i  ) = mesh%f(3*i  )
  enddo

  call solver(mesh)
  call stress_update(mesh)

  do i=1,mesh%nnode
    mesh%u(3*i-2) = mesh%X(3*i-2)
    mesh%u(3*i-1) = mesh%X(3*i-1)
    mesh%u(3*i  ) = mesh%X(3*i  )
  enddo
end subroutine static

subroutine nonlinear_static(mesh)
  use util
  implicit none
  type(meshdef) :: mesh
  integer(kint) :: i, NRiter
  real(kdouble) :: b0nrm, bnrm, rnrm, rnrmmax
  real(kdouble), allocatable :: q(:), f(:)

  mesh%u  = 0.0d0
  mesh%du = 0.0d0
  mesh%cur_tstep = 1

  do NRiter=1,4
    mesh%cur_nrstep = NRiter
    write(*,"(a)")""
    write(*,"(a,i8)")"** NRiter:", NRiter

    call stiff_matrix(mesh)
    call load_condition(mesh)

    !write(*,"(a)")"mesh%f"
    !write(*,"(3e12.5)")mesh%f
    !write(*,"(a)")"mesh%q"
    !write(*,"(3e12.5)")mesh%q

    do i=1,mesh%nnode
      mesh%B(3*i-2) = mesh%f(3*i-2) - mesh%q(3*i-2)
      mesh%B(3*i-1) = mesh%f(3*i-1) - mesh%q(3*i-1)
      mesh%B(3*i  ) = mesh%f(3*i  ) - mesh%q(3*i  )
    enddo

    call bound_condition(mesh)

    !write(*,"(a)")"mesh%b"
    !write(*,"(3e12.5)")mesh%b

    bnrm = 0.0d0
    do i=1,3*mesh%nnode
      bnrm = bnrm + mesh%B(i)*mesh%B(i)
    enddo
    if(NRiter == 1)then
      b0nrm = bnrm
    else
      rnrm    = dsqrt(bnrm/b0nrm)
      rnrmmax = dabs(maxval(mesh%B))
      write(*,"(a,1pe12.5,a,1pe12.5)")"  ** NR     residual: ", rnrm, ", ", rnrmmax
    endif

    call solver(mesh)

    !write(*,"(a)")"mesh%x"
    !write(*,"(3e12.5)")mesh%x

    do i=1,mesh%nnode
      mesh%du(3*i-2) = mesh%du(3*i-2) + mesh%X(3*i-2)
      mesh%du(3*i-1) = mesh%du(3*i-1) + mesh%X(3*i-1)
      mesh%du(3*i  ) = mesh%du(3*i  ) + mesh%X(3*i  )
    enddo

    call stress_update(mesh)
  enddo

  do i=1,mesh%nnode
    mesh%u(3*i-2) = mesh%u(3*i-2) + mesh%du(3*i-2)
    mesh%u(3*i-1) = mesh%u(3*i-1) + mesh%du(3*i-1)
    mesh%u(3*i  ) = mesh%u(3*i  ) + mesh%du(3*i  )
  enddo
end subroutine nonlinear_static
