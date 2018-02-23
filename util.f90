
subroutine get_inverse_matrix(n, a, inv)
  use util
  implicit none
  integer(kint) :: n, i, j, k
  real(kdouble) :: a(n,n), inv(n,n), tmp

  inv = 0.0d0
  do i=1,8
    inv(i,i) = 1.0d0
  enddo

  do i=1,n
    tmp = 1.0d0/a(i,i)
    do j=1,n
        a(j,i) =   a(j,i) * tmp
      inv(j,i) = inv(j,i) * tmp
    enddo
    do j=1,n
      if(i /= j) then
        tmp = a(i,j)
        do k=1,n
            a(k,j) =   a(k,j) -   a(k,i) * tmp
          inv(k,j) = inv(k,j) - inv(k,i) * tmp
        enddo
      endif
    enddo
  enddo
end subroutine get_inverse_matrix
