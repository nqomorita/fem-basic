
subroutine static(mesh)
  use util
  implicit none
  type(meshdef) :: mesh
  integer(kint) :: i, id

  call stiff_matrix(mesh)
  call bound_condition(mesh)
  call solver(mesh)
  !call disp_update(mesh)
  call stress_update(mesh)
end subroutine static
