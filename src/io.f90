
subroutine input_mesh(mesh)
  use util
  implicit none
  type(meshdef) :: mesh
  integer(kint) :: i, id

  open(10, file="input.dat", status='old')
    read(10,*) i
    if(i == 1) isNLGeom = .true.
    read(10,*) mesh%max_nrstep
    read(10,*) mesh%nnode
    allocate(mesh%node(ndof, mesh%nnode))
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

subroutine init_mesh(mesh)
  use util
  implicit none
  type(meshdef) :: mesh

  !allocate section
  allocate(mesh%gauss(8, mesh%nelem))
  allocate(mesh%nstrain(6, mesh%nnode))
  allocate(mesh%nstress(6, mesh%nnode))
  allocate(mesh%nmises (mesh%nnode))
  allocate(mesh%estrain(6, mesh%nelem))
  allocate(mesh%estress(6, mesh%nelem))
  allocate(mesh%emises (mesh%nelem))
  allocate(mesh%u (ndof*mesh%nnode))
  allocate(mesh%du(ndof*mesh%nnode))
  allocate(mesh%q (ndof*mesh%nnode))
  allocate(mesh%f (ndof*mesh%nnode))
  allocate(mesh%A (ndof*mesh%nnode, ndof*mesh%nnode))
  allocate(mesh%x (ndof*mesh%nnode))
  allocate(mesh%b (ndof*mesh%nnode))
  mesh%nstrain = 0.0d0
  mesh%nstress = 0.0d0
  mesh%nmises  = 0.0d0
  mesh%estrain = 0.0d0
  mesh%estress = 0.0d0
  mesh%emises  = 0.0d0
  mesh%u       = 0.0d0
  mesh%du      = 0.0d0
  mesh%q       = 0.0d0
  mesh%f       = 0.0d0
  mesh%A       = 0.0d0
  mesh%x       = 0.0d0
  mesh%b       = 0.0d0
end subroutine init_mesh

subroutine outout_res(mesh)
  use util
  implicit none
  type(meshdef) :: mesh
  integer(kint) :: i, id

  open(10, file='result.inp', status='replace')
    write(10,"(5i3)")mesh%nnode, mesh%nelem, 16, 13, 0

    do i=1,mesh%nnode
      write(10,"(i8, 1p3e12.4)")i, mesh%node(1,i), mesh%node(2,i), mesh%node(3,i)
    enddo

    do i=1,mesh%nelem
      write(10,"(2i8, a5, 8i8)")i, 1, " hex ", mesh%elem(1,i), mesh%elem(2,i), mesh%elem(3,i), mesh%elem(4,i), &
                                             & mesh%elem(5,i), mesh%elem(6,i), mesh%elem(7,i), mesh%elem(8,i)
    enddo

    write(10,"(5i8)")4, 3, 6, 6, 1
    write(10,"(a)")"disp, mm"
    write(10,"(a)")"strain, "
    write(10,"(a)")"stress, MPa"
    write(10,"(a)")"mises, MPa"
    do i=1,mesh%nnode
      write(10,"(i8, 1p16e12.4)")i, mesh%u(3*i-2), mesh%u(3*i-1), mesh%u(3*i), &
      & mesh%nstrain(1,i), mesh%nstrain(2,i), mesh%nstrain(3,i), mesh%nstrain(4,i), mesh%nstrain(5,i), mesh%nstrain(6,i), &
      & mesh%nstress(1,i), mesh%nstress(2,i), mesh%nstress(3,i), mesh%nstress(4,i), mesh%nstress(5,i), mesh%nstress(6,i), &
      & mesh%nmises(i)
    enddo

    write(10,"(4i8)")3, 6, 6, 1
    write(10,"(a)")"strain, "
    write(10,"(a)")"stress, MPa"
    write(10,"(a)")"mises, MPa"
    do i=1,mesh%nelem
      write(10,"(i8, 1p13e12.4)")i, &
      & mesh%estrain(1,i), mesh%estrain(2,i), mesh%estrain(3,i), mesh%estrain(4,i), mesh%estrain(5,i), mesh%estrain(6,i), &
      & mesh%estress(1,i), mesh%estress(2,i), mesh%estress(3,i), mesh%estress(4,i), mesh%estress(5,i), mesh%estress(6,i), &
      & mesh%emises(i)
    enddo
  close(10)
end subroutine outout_res

subroutine finalize_mesh(mesh)
  use util
  implicit none
  type(meshdef) :: mesh

  deallocate(mesh%gauss)
  deallocate(mesh%nstrain)
  deallocate(mesh%nstress)
  deallocate(mesh%nmises)
  deallocate(mesh%estrain)
  deallocate(mesh%estress)
  deallocate(mesh%emises)
  deallocate(mesh%u)
  deallocate(mesh%du)
  deallocate(mesh%q)
  deallocate(mesh%f)
  deallocate(mesh%A)
  deallocate(mesh%x)
  deallocate(mesh%b)
end subroutine finalize_mesh

