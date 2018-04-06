
subroutine C3D8_integral_point(i, r)
  use util
  implicit none
  integer(kint) :: i
  real(kdouble) :: gsp(3,8), r(3)

  gsp(1,1) = -0.577350269189626d0; gsp(2,1) = -0.577350269189626d0; gsp(3,1) = -0.577350269189626d0
  gsp(1,2) =  0.577350269189626d0; gsp(2,2) = -0.577350269189626d0; gsp(3,2) = -0.577350269189626d0
  gsp(1,3) = -0.577350269189626d0; gsp(2,3) =  0.577350269189626d0; gsp(3,3) = -0.577350269189626d0
  gsp(1,4) =  0.577350269189626d0; gsp(2,4) =  0.577350269189626d0; gsp(3,4) = -0.577350269189626d0
  gsp(1,5) = -0.577350269189626d0; gsp(2,5) = -0.577350269189626d0; gsp(3,5) =  0.577350269189626d0
  gsp(1,6) =  0.577350269189626d0; gsp(2,6) = -0.577350269189626d0; gsp(3,6) =  0.577350269189626d0
  gsp(1,7) = -0.577350269189626d0; gsp(2,7) =  0.577350269189626d0; gsp(3,7) =  0.577350269189626d0
  gsp(1,8) =  0.577350269189626d0; gsp(2,8) =  0.577350269189626d0; gsp(3,8) =  0.577350269189626d0

  r(1) = gsp(1,i)
  r(2) = gsp(2,i)
  r(3) = gsp(3,i)
end subroutine C3D8_integral_point

subroutine C3D8_shapefunc(local, func)
  use util
  implicit none
  real(kdouble) ::  local(3), func(8)

  func(1) = 0.125d0*(1.0d0-local(1))*(1.0d0-local(2))*(1.0d0-local(3))
  func(2) = 0.125d0*(1.0d0+local(1))*(1.0d0-local(2))*(1.0d0-local(3))
  func(3) = 0.125d0*(1.0d0+local(1))*(1.0d0+local(2))*(1.0d0-local(3))
  func(4) = 0.125d0*(1.0d0-local(1))*(1.0d0+local(2))*(1.0d0-local(3))
  func(5) = 0.125d0*(1.0d0-local(1))*(1.0d0-local(2))*(1.0d0+local(3))
  func(6) = 0.125d0*(1.0d0+local(1))*(1.0d0-local(2))*(1.0d0+local(3))
  func(7) = 0.125d0*(1.0d0+local(1))*(1.0d0+local(2))*(1.0d0+local(3))
  func(8) = 0.125d0*(1.0d0-local(1))*(1.0d0+local(2))*(1.0d0+local(3))
end subroutine C3D8_shapefunc

subroutine C3D8_shapefunc_deriv(local, func)
  use util
  implicit none
  real(kdouble) :: local(3), func(8,3)

  func(1,1) = -0.125d0*(1.0d0-local(2))*(1.0d0-local(3))
  func(2,1) =  0.125d0*(1.0d0-local(2))*(1.0d0-local(3))
  func(3,1) =  0.125d0*(1.0d0+local(2))*(1.0d0-local(3))
  func(4,1) = -0.125d0*(1.0d0+local(2))*(1.0d0-local(3))
  func(5,1) = -0.125d0*(1.0d0-local(2))*(1.0d0+local(3))
  func(6,1) =  0.125d0*(1.0d0-local(2))*(1.0d0+local(3))
  func(7,1) =  0.125d0*(1.0d0+local(2))*(1.0d0+local(3))
  func(8,1) = -0.125d0*(1.0d0+local(2))*(1.0d0+local(3))

  func(1,2) = -0.125d0*(1.0d0-local(1))*(1.0d0-local(3))
  func(2,2) = -0.125d0*(1.0d0+local(1))*(1.0d0-local(3))
  func(3,2) =  0.125d0*(1.0d0+local(1))*(1.0d0-local(3))
  func(4,2) =  0.125d0*(1.0d0-local(1))*(1.0d0-local(3))
  func(5,2) = -0.125d0*(1.0d0-local(1))*(1.0d0+local(3))
  func(6,2) = -0.125d0*(1.0d0+local(1))*(1.0d0+local(3))
  func(7,2) =  0.125d0*(1.0d0+local(1))*(1.0d0+local(3))
  func(8,2) =  0.125d0*(1.0d0-local(1))*(1.0d0+local(3))

  func(1,3) = -0.125d0*(1.0d0-local(1))*(1.0d0-local(2))
  func(2,3) = -0.125d0*(1.0d0+local(1))*(1.0d0-local(2))
  func(3,3) = -0.125d0*(1.0d0+local(1))*(1.0d0+local(2))
  func(4,3) = -0.125d0*(1.0d0-local(1))*(1.0d0+local(2))
  func(5,3) =  0.125d0*(1.0d0-local(1))*(1.0d0-local(2))
  func(6,3) =  0.125d0*(1.0d0+local(1))*(1.0d0-local(2))
  func(7,3) =  0.125d0*(1.0d0+local(1))*(1.0d0+local(2))
  func(8,3) =  0.125d0*(1.0d0-local(1))*(1.0d0+local(2))
end subroutine C3D8_shapefunc_deriv
