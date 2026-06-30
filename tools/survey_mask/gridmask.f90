! mask onto a grid by using mocks
program gridmask
    implicit none

    !in
    integer :: n
    real, allocatable :: x(:), y(:), z(:)
    
    !timing
    integer*8 :: t1, t2, trate, tmax

    !local
    logical :: fixed_box 
    integer :: ip,  i, j, k
    integer :: i2, j2, k2
    integer :: nx
    integer*1, allocatable :: mask(:,:,:)
    real :: dx, adj, shift
    real :: boxlim(6), side, lower, upper

    !out
    real, allocatable :: vx(:), vy(:), vz(:)

    ! Read DATA
    call system_clock(t1,trate,tmax)
    open(unit=10, file='../STORAGE/rands.dat', form = 'unformatted')
    read(10) n
    WRITE(*,*) "Number of mock/mask points:", n
    allocate(x(n), y(n), z(n))
    read(10) x
    read(10) y
    read(10) z
    close(10)

    call system_clock(t2,trate,tmax)
    WRITE(*,*) "Time taken to read points:", float(t2 - t1)/float(trate), "seconds"
    WRITE(*,*)

    !data already shifted and properly in -L/2, L/2
    fixed_box = .true.

    if (fixed_box .eqv. .true.) then
        write(*,*) 'Introduce side of box:'
        read(*,*) side
        write(*,*) 'side:', side
        lower = -side/2
        upper = side/2

    !shifts to centre and -L/2, L/2 with rand 
    else 
        !Set upper and lower limits
        boxlim(1) = minval(x)
        boxlim(2) = maxval(x)
        boxlim(3) = minval(y)
        boxlim(4) = maxval(y)
        boxlim(5) = minval(z)
        boxlim(6) = maxval(z)
        WRITE(*,*) 'X', boxlim(1), boxlim(2)
        WRITE(*,*) 'Y', boxlim(3), boxlim(4)
        WRITE(*,*) 'Z', boxlim(5), boxlim(6)

        lower = boxlim(1)
        upper = boxlim(2)
        if (boxlim(3) < lower) lower = boxlim(3)
        if (boxlim(4) > upper) upper = boxlim(4)
        if (boxlim(5) < lower) lower = boxlim(5)
        if (boxlim(6) > upper) upper = boxlim(6)
        shift = (upper + lower) / 2.
    
        !shift to (-L/2, L/2)
        x = x - shift
        y = y - shift
        z = z - shift

        boxlim(1) = minval(x)
        boxlim(2) = maxval(x)
        boxlim(3) = minval(y)
        boxlim(4) = maxval(y)
        boxlim(5) = minval(z)
        boxlim(6) = maxval(z)
        WRITE(*,*) 'X (shifted)', boxlim(1), boxlim(2)
        WRITE(*,*) 'Y (shifted)', boxlim(3), boxlim(4)
        WRITE(*,*) 'Z (shifted)', boxlim(5), boxlim(6)

        lower = boxlim(1)
        upper = boxlim(2)
        if (boxlim(3) < lower) lower = boxlim(3)
        if (boxlim(4) > upper) upper = boxlim(4)
        if (boxlim(5) < lower) lower = boxlim(5)
        if (boxlim(6) > upper) upper = boxlim(6)
        shift = (upper + lower) / 2. 

        !Side and buffer
        side = upper - lower
        adj = 0.01*side
        lower = lower - adj
        upper = upper + adj
        side = upper - lower

        WRITE(*,*) '----------------------------------------------------------------'
        WRITE(*,*) 'lower, upper, side, centre', lower, upper, side, shift
        WRITE(*,*) '----------------------------------------------------------------'
        WRITE(*,*)
    endif

    !grid and mask fill -> TSC
    nx = 500
    write(*,*) 'nx:', nx
    dx = side / real(nx)
    allocate(mask(nx,nx,nx))
    mask = 0
    !$OMP PARALLEL SHARED(mask,x,y,z,lower,dx,n,nx), &
    !$OMP PRIVATE(ip,i,j,k), &
    !$OMP DEFAULT(NONE)
    !$OMP DO REDUCTION(+:mask)
    do ip=1,n
        i = int((x(ip) - lower) / dx) + 1
        j = int((y(ip) - lower) / dx) + 1
        k = int((z(ip) - lower) / dx) + 1
        DO k2=k-1,k+1
        IF (k2.GT.nx .OR. k2.LT.1) CYCLE
        DO j2=j-1,j+1
        IF (j2.GT.nx .OR. j2.LT.1) CYCLE
        DO i2=i-1,i+1
        IF (i2.GT.nx .OR. i2.LT.1) CYCLE
            mask(i2,j2,k2) = 1
        enddo
        enddo
        enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL
    WRITE(*,*) 'Filled cells, %', COUNT(mask>0), &
                            real(COUNT(mask>0)) / real(nx**3)*100
    WRITE(*,*)
    
    !to 1 and 0
    where(mask > 0) mask = 1

    open(unit=1, file='mask.dat', form = 'unformatted', access='stream', action='write', status='replace')
    write(1) nx
    write(1) lower, upper, side
    write(1) mask
    close(1)

end program
