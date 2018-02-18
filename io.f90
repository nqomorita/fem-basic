
subroutine input_mesh(mesh)
  use util
  implicit none
  type(meshdef) :: mesh
  integer(kint) :: i, id

  open(10, file="input.dat", status='old')
    read(10,*) mesh%nnode
    allocate(mesh%node(3, mesh%nnode))
    do i=1, mesh%nnode
      read(10,*) id, mesh%node(1,i), mesh%node(2,i), mesh%node(3,i)
    enddo

    read(10,*) mesh%nelem
    allocate(mesh%elem(8, mesh%nelem))
    do i=1, mesh%nelem
      read(10,*)id, mesh%elem(1,i), mesh%elem(2,i), mesh%elem(3,i), mesh%elem(4,i), &
                  & mesh%elem(5,i), mesh%elem(6,i), mesh%elem(7,i), mesh%elem(8,i)
    enddo

    read(10,*)mesh%nbound
    allocate(mesh%ibound(2, mesh%nbound))
    allocate(mesh%bound (mesh%nbound))
    do i=1,mesh%nbound
      read(10,*) mesh%ibound(1,i), mesh%ibound(2,i), mesh%bound(i)
    enddo

    read(10,*)mesh%ncload
    allocate(mesh%icload(2, mesh%ncload))
    allocate(mesh%cload (mesh%ncload))
    do i=1,mesh%ncload
      read(10,*) mesh%icload(1,i), mesh%icload(2,i), mesh%cload(i)
    enddo

    read(10,*) mesh%E
    read(10,*) mesh%mu
    read(10,*) mesh%rho
  close(10)
end subroutine input_mesh

subroutine outout_res(mesh)
  use util
  implicit none
  type(meshdef) :: mesh
  integer(kint) :: i, id

  open(10, file='result.inp', status='replace')
    write(10,"(5i3)")mesh%nnode, mesh%nelem, 4, 0, 0

    do i=1,mesh%nnode
      write(10,"(i8, 3e12.4)")i, mesh%node(1,i), mesh%node(2,i), mesh%node(3,i)
    end do

    do i=1,mesh%nelem
      write(10,"(2i8, a5, 8i8)")i, 1, " hex ", mesh%elem(1,i), mesh%elem(2,i), mesh%elem(3,i), mesh%elem(4,i), &
                                             & mesh%elem(5,i), mesh%elem(6,i), mesh%elem(7,i), mesh%elem(8,i)
    end do

    write(10,"(3i8)")2, 3, 1
    write(10,"(a)")"disp, mm"
    write(10,"(a)")"stress, MPa"
    do i=1,mesh%nnode
      write(10,"(i8, 4e12.4)")i, mesh%X(3*i-2), mesh%X(3*i-1), mesh%X(3*i), 0.0
    end do
  close(10)
end subroutine outout_res
