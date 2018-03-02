
subroutine static(mesh)
  use util
  implicit none
  type(meshdef) :: mesh
  integer(kint) :: i

  call stiff_matrix(mesh)
  call load_condition(mesh)
  call get_RHS(mesh)
  call bound_condition(mesh)
  call solver(mesh)
  call stress_update(mesh)
  call delta_u_update(mesh)
  call u_update(mesh)
end subroutine static

subroutine nonlinear_static(mesh)
  use util
  implicit none
  type(meshdef) :: mesh
  integer(kint) :: i, NRiter
  real(kdouble), allocatable :: q(:), f(:)

  mesh%u  = 0.0d0
  mesh%du = 0.0d0
  mesh%cur_tstep = 1

  do NRiter=1,5
    mesh%cur_nrstep = NRiter
    write(*,"(a)")""
    write(*,"(a,i8)")"** NRiter:", NRiter

    call stiff_matrix(mesh)
    call load_condition(mesh)
    call get_RHS(mesh)
    call bound_condition(mesh)
    call is_convergence(mesh)
    call solver(mesh)
    call stress_update(mesh)
    call delta_u_update(mesh)
  enddo

  call u_update(mesh)
end subroutine nonlinear_static
