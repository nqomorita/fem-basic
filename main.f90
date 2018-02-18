module util
  integer(4), parameter :: kint    = 4
  integer(4), parameter :: kdouble = 8

  type meshdef
    integer(kint) :: nnode
    integer(kint) :: nelem
    integer(kint) :: ndof
    integer(kint), pointer :: elem(:,:)
    integer(kint), pointer :: ibound(:,:)
    integer(kint), pointer :: icload(:,:)
    integer(kint) :: nbound, ncload

    real(kdouble), pointer :: A(:,:)
    real(kdouble), pointer :: B(:)
    real(kdouble), pointer :: X(:)
    real(kdouble), pointer :: node(:,:)
    real(kdouble), pointer :: bound(:)
    real(kdouble), pointer :: cload(:)
    real(kdouble) :: E, mu, rho
  end type meshdef
end module util

program main
  use util
  implicit none
  type(meshdef) :: mesh

  call input_mesh(mesh)
  call static(mesh)
  call outout_res(mesh)

end program main
