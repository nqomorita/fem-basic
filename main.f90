module util
  integer(4), parameter :: kint    = 4
  integer(4), parameter :: kdouble = 8

  logical, save :: isNLGeom = .false.

  type gaussdef
    real(kdouble) :: strain(6) = 0.0d0
    real(kdouble) :: stress(6) = 0.0d0
  end type gaussdef

  type meshdef
    integer(kint) :: cur_tstep
    integer(kint) :: cur_nrstep

    integer(kint) :: nnode
    real(kdouble), pointer :: node(:,:)

    integer(kint) :: nelem
    integer(kint), pointer :: elem(:,:)

    integer(kint) :: nbound
    integer(kint), pointer :: ibound(:,:)
    real(kdouble), pointer :: bound(:)

    integer(kint) :: ncload
    integer(kint), pointer :: icload(:,:)
    real(kdouble), pointer :: cload(:)

    real(kdouble) :: E, mu, rho

    real(kdouble), pointer :: u(:)
    real(kdouble), pointer :: du(:)
    real(kdouble), pointer :: q(:)
    real(kdouble), pointer :: f(:)

    real(kdouble), pointer :: A(:,:)
    real(kdouble), pointer :: B(:)
    real(kdouble), pointer :: X(:)

    type(gaussdef), pointer :: gauss(:,:)
    real(kdouble), pointer :: nstrain(:,:)
    real(kdouble), pointer :: nstress(:,:)
    real(kdouble), pointer :: nmises(:)
    real(kdouble), pointer :: estrain(:,:)
    real(kdouble), pointer :: estress(:,:)
    real(kdouble), pointer :: emises(:)
  end type meshdef
end module util

program main
  use util
  implicit none
  type(meshdef) :: mesh

  call input_mesh(mesh)
  call init_mesh(mesh)

  if(isNLGeom)then
    call nonlinear_static(mesh)
  else
    !call static(mesh)
    call nonlinear_static(mesh)
  endif

  call outout_res(mesh)
  call finalize_mesh(mesh)
end program main
