!#######################################################
module cosmokdtree
!#######################################################
    implicit none
    private  
    public :: build_kdtree, deallocate_kdtree, & 
              knn_search, ball_search, &
              KDTreeNode
    !+++++++++++++++++++++++++++++++
    !++++ Dimensionality (default 3D)
    !+++++++++++++++++++++++++++++++   
#if defined(dimen)
    integer, parameter :: ndim = dimen
#else 
    integer, parameter :: ndim = 3
#endif

    !+++++++++++++++++++++++++++++++
    !++++ Precision and integer kind
    !+++++++++++++++++++++++++++++++
#if doubleprecision == 1 
    integer, parameter :: prec = 8 
#else
    integer, parameter :: prec = 4
#endif

#if longint == 1
    integer, parameter :: intkind = 8
#else
    integer, parameter :: intkind = 4
#endif

    !+++++++++++++++++++++++++++++++
    !++++ Periodic boundary conditions
    !+++++++++++++++++++++++++++++++
#if periodic == 1
    real(kind=prec) :: L(ndim) ! Will be initialized in build_kdtree
#endif

    !+++++++++++++++++++++++++++++++
    !++++ Type definitions
    !+++++++++++++++++++++++++++++++
    type :: KDTreeNode
        !basic ---------------------
        real(kind=prec) :: point !splitting point coordinates
        integer :: axis !splitting axis (1 for x, 2 for y, 3 for z, 4 for w, ...)
        type(KDTreeNode), pointer :: left => null()  !left child node
        type(KDTreeNode), pointer :: right => null() !right child node
        !leaf variables
        logical :: is_leaf  !indicates if the node is a leaf
        integer(kind=intkind) :: start_idx
        integer(kind=intkind) :: end_idx
    end type KDTreeNode
    !+++++++++++++++++++++++++++++++

contains

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Initialize kd-tree construction (not in-place for thread-safety)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#if periodic == 1
    function build_kdtree(points, idx, L_in, leaf) result(tree)
#else
    function build_kdtree(points, idx, leaf) result(tree)
#endif
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        use omp_lib
        implicit none
        !in
        real(kind=prec), intent(inout) :: points(:, :)
        integer(kind=intkind), intent(inout) :: idx(:)
        integer, intent(in), optional :: leaf 
#if periodic == 1
        real(kind=prec), intent(in) :: L_in(:)
        integer :: flag_stop
        integer :: ndim_par
#endif
        !local
        integer :: j
        integer(kind=intkind) :: n, i
        integer :: depth, max_depth, nproc, leafsize
        real(kind=prec) :: bounds(2*ndim)
        
        !out
        type(KDTreeNode), pointer :: tree

        !!!! nested parallelism!!!!!!!
        !-- deprecated
        ! call omp_set_nested(.true.) 
        !-- new -> nested parallelism limit (30-40 are already extremely high values)
        !     for security, we set it to 50. For a leaf of 16 would mean npart=10^16 
        !     far beyond the actual manageable limit of 10^10-10^11 points.
        call omp_set_max_active_levels(50)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! Check dimensionality of input
        if (size(points, 1) /= ndim) then
            STOP 'Input points must have the same dimensionality as the tree!'
        end if

        ! Number of points
        n = size(points, 2, kind=intkind)

#if periodic == 1
        ! Check dimensionality of input
        if (size(L_in, 1) /= ndim) then
            STOP 'Input L must have the same dimensionality as the tree!'
        end if
        
        !(xmin, xmax, ymin, ymax, zmin, zmax, ...)
        L = L_in

        !points within (-L/2, L/2)
        flag_stop = 0
        ndim_par = ndim
        !$OMP PARALLEL SHARED(points,n,L,ndim_par), &
        !$OMP PRIVATE(i,j)
        !$OMP DO REDUCTION(+:flag_stop)
        do i=1,n
            do j=1,ndim_par
                if (points(j,i) .lt. -L(j)/2 .or. points(j,i) .gt. L(j)/2) flag_stop = 1
            enddo
        enddo
        !$OMP ENDDO
        !$OMP END PARALLEL
        
        if (flag_stop .gt. 0) STOP 'Points outside (-L/2, L/2) range !!'
#endif

        !depth for parallelism control
        depth = 0

        !$OMP PARALLEL
        !$OMP SINGLE
        nproc = omp_get_num_threads()
        !$OMP END SINGLE
        !$OMP END PARALLEL

        max_depth = compute_max_depth(omp_get_max_threads())

        !define leaf
        if ( present(leaf) ) then 
            leafsize = leaf
        else
            leafsize = 64
        endif

        !init point bounds -> check periodic
#if periodic == 1
        bounds(1) = -L(1)/2.
        bounds(2) =  L(1)/2.
        bounds(3) = -L(2)/2.
        bounds(4) =  L(2)/2.
        bounds(5) = -L(3)/2.
        bounds(6) =  L(3)/2.
#else
        !bounds
        do j = 1, ndim
            bounds(2*j-1) = HUGE(0.0)
            bounds(2*j)   = -HUGE(0.0)
        end do

        !find all min/max without array temporaries
        do i = 1, n
            do j = 1, ndim
                if (points(j, i) < bounds(2*j-1)) bounds(2*j-1) = points(j, i)
                if (points(j, i) > bounds(2*j))   bounds(2*j)   = points(j, i)
            end do
        end do
#endif

        !recursive build call
        tree => build_kdtree_recursive(points, idx, 1_intkind, n, depth, max_depth, leafsize, bounds)

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end function build_kdtree
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Deallocation of the kd-tree
    ! This subroutine recursively deallocates the kd-tree nodes and their associated data
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    recursive subroutine deallocate_kdtree(node)
        implicit none
        type(KDTreeNode), pointer, intent(inout) :: node
        if (.not. associated(node)) return

        !children
        call deallocate_kdtree(node%left)
        call deallocate_kdtree(node%right)

        !self
        deallocate(node)

    end subroutine deallocate_kdtree
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !Compute maximum depth of the tree for parallelism (Ncores)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function compute_max_depth(nproc) result(max_depth)
        implicit none
        integer, intent(in) :: nproc
        integer :: max_depth
        !nproc: number of processes/threads available.
        !max_depth Maximum depth which guarantees that there is at least one idle process.
        max_depth = int (log(dble(nproc)+0.1) / log(2.))
    end function compute_max_depth
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Recursive kd-tree building function according to ---> sliding midpoint rule <----
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    recursive function build_kdtree_recursive(points,idx,p_start,p_end,depth,max_depth,leafsize,bounds) result(node)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    use omp_lib
    implicit none
    !in
    real(kind=prec), intent(inout) :: points(:, :) !(ndim, n)
    integer(kind=intkind), intent(inout) :: idx(:)
    integer(kind=intkind), intent(in) :: p_start, p_end !range of points indices to partition
    integer, intent(in) :: depth
    integer, intent(in) :: max_depth
    integer, intent(in) :: leafsize
    real(kind=prec), intent(in) :: bounds(:)
    !local
    real(kind=prec) :: maxside
    real(kind=prec) :: bounds_left(ndim*2), bounds_right(ndim*2) !bounds for left and right child nodes
    real(kind=prec) :: minvals(ndim), maxvals(ndim), spreads(ndim), side(ndim) !for splitting axis selection
    real(kind=prec) :: midpoint
    integer :: axis, j
    integer(kind=intkind) :: i, size_points, n_left !data partitioning variables
    type(KDTreeNode), pointer :: node 
    !edge cases (sliding)
    real(kind=prec) :: minx, maxx
    integer(kind=intkind) :: minind, maxind
    real(kind=prec) :: temp_point(ndim) !for swapping
    integer(kind=intkind) :: temp_idx !for swapping

    allocate(node)

    size_points = p_end - p_start + 1

    if (size_points == 0) then
        !should never happen, return a null node
        node => null()
        return
    end if

    if (size_points <= leafsize) then
        node%is_leaf = .true.
        ! Store start and end indices of the points in the leaf node
        node%start_idx = p_start
        node%end_idx = p_end
        node%left => null()
        node%right => null()
        return
    else
        node%is_leaf = .false.
    end if

    ! Compute necessary statistics to choose splitting axis and point
    side = 0.
    do j = 1, ndim
        minvals(j) = HUGE(0.0)
        maxvals(j) = -HUGE(0.0)
    end do

    do i = p_start, p_end
        do j = 1, ndim
            if (points(j, i) < minvals(j)) minvals(j) = points(j, i)
            if (points(j, i) > maxvals(j)) maxvals(j) = points(j, i)
        end do
    end do

    do j = 1, ndim
        side(j)    = bounds(2*j) - bounds(2*j-1)
        spreads(j) = maxvals(j) - minvals(j)
    end do

    maxside = side(1)
    axis = 1
    do j = 2, ndim
        if (side(j) > maxside) then
            maxside = side(j)
            axis = j
        end if
        if (side(j) == maxside .and. spreads(j) > spreads(axis)) then
            axis = j
        end if
    end do
    
    midpoint = 0.5 * (bounds(2*axis-1) + bounds(2*axis))

    !check sliding
    ! minx = minvals(axis)
    ! maxx = maxvals(axis)
    ! minind = minloc(points(axis, :), dim = 1)
    ! maxind = maxloc(points(axis, :), dim = 1)
    minx = HUGE(0.0)
    maxx = -HUGE(0.0)
    do i = p_start, p_end
        if (points(axis, i) < minx) then
            minx = points(axis, i)
            minind = i
        end if
        if (points(axis, i) > maxx) then
            maxx = points(axis, i)
            maxind = i
        end if
    end do

    !slide the hyperplane to the closest point 
    if (minx >= midpoint) then
        midpoint = minx

        !swaping mindind to the first position
        temp_point = points(:, p_start)
        points(:, p_start) = points(:, minind)
        points(:, minind) = temp_point

        temp_idx = idx(p_start)
        idx(p_start) = idx(minind)
        idx(minind) = temp_idx

        n_left = 1 !just one point on the left

    !slide the hyperplane to the closest point 
    else if (maxx < midpoint) then
        midpoint = maxx
        
        !swaping maxind to the last position
        temp_point = points(:, p_end)
        points(:, p_end) = points(:, maxind)
        points(:, maxind) = temp_point

        temp_idx = idx(p_end)
        idx(p_end) = idx(maxind)
        idx(maxind) = temp_idx

        n_left = size_points - 1 !all points on the left except one

    !do not slide, default: midpoint splitting rule
    else
        
        call partition(points, idx, p_start, p_end, midpoint, axis, n_left)
        
    end if
    ! all points with points(axis, :) < midpoint are in the left side
    ! all points with points(axis, :) >= midpoint are in the right side
    
    !splitting plane
    node%point = midpoint
    node%axis = axis

    ! bounds for left and right child nodes
    bounds_left = bounds
    bounds_right = bounds
    bounds_left(2*axis) = midpoint
    bounds_right(2*axis-1) = midpoint

    ! Subtree construction (parallel at the top levels)
    if (depth < max_depth) then
    !$OMP PARALLEL NUM_THREADS(2)
    !$OMP SINGLE

    !$OMP TASK
    node%left => build_kdtree_recursive(points, idx, p_start, p_start + n_left - 1, &
                                        depth+1, max_depth, leafsize, bounds_left)
    !$OMP END TASK

    !$OMP TASK
    node%right => build_kdtree_recursive(points, idx, p_start + n_left, p_end, &
                                        depth+1, max_depth, leafsize, bounds_right)
    !$OMP END TASK

    !$OMP END SINGLE
    !$OMP END PARALLEL
    else
    node%left => build_kdtree_recursive(points, idx, p_start, p_start + n_left - 1, &
                                        depth+1, max_depth, leafsize, bounds_left)
    node%right => build_kdtree_recursive(points, idx, p_start + n_left, p_end, &
                                        depth+1, max_depth, leafsize, bounds_right)
    end if

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end function build_kdtree_recursive
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Partition function to split points around a pivot value along a given axis
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subroutine partition(points, idx, p_start, p_end, pivot_value, axis, n_left)
        implicit none
        real(kind=prec), intent(inout) :: points(:, :)
        integer(kind=intkind), intent(inout) :: idx(:)
        integer(kind=intkind), intent(in) :: p_start, p_end
        integer, intent(in) :: axis
        real(kind=prec):: pivot_value
        integer(kind=intkind) :: n_left
        integer(kind=intkind) :: i, j
        real(kind=prec):: temp_point(ndim)
        integer(kind=intkind) :: temp_idx
        
        n_left = 0
        do i = p_start, p_end
            if (points(axis, i) < pivot_value) then
                j = p_start + n_left
                if (i /= j) then
                    temp_point(:) = points(:, i)
                    points(:, i)  = points(:, j)
                    points(:, j)  = temp_point(:)

                    temp_idx = idx(i)
                    idx(i) = idx(j)
                    idx(j) = temp_idx
                end if
                n_left = n_left + 1
            end if
        end do
    end subroutine partition
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! k-nearest neighbor search main function
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subroutine knn_search(node, targett, k, points, q_idx, q_dist, sorted)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        implicit none
        !in
        integer, intent(in) :: k
        real(kind=prec), intent(in) :: points(:, :)
        type(KDTreeNode), pointer, intent(in) :: node
        real(kind=prec), intent(in) :: targett(ndim)
        logical, intent(in), optional :: sorted
        !out
        real(kind=prec), intent(inout):: q_dist(k)
        integer(kind=intkind), intent(inout) :: q_idx(k)

        q_dist = HUGE(0.0)
        q_idx = -1

        call knn_search_recursive_hyperp(node, targett, points, q_idx, q_dist, k)

        !Perform full sort with quicksort
        ! by default, sorts
        if ( present(sorted) ) then 
            if (sorted) then
                call quicksort(q_dist, q_idx, k)
            end if
        else
            call quicksort(q_dist, q_idx, k)
        endif

        q_dist = sqrt(q_dist)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end subroutine knn_search
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Standard hyperplane method to find k-nearest neighbors (SIMPLEST)
    ! Just slightly slower than queue bbox version at large k
    ! Can be faster, depending on the case, for small k
    ! (periodicity available)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    recursive subroutine knn_search_recursive_hyperp(node, targett, points, q_idx, q_dist, k)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        implicit none
        !in
        integer :: k ! Number of nearest neighbors to find
        type(KDTreeNode), pointer, intent(in) :: node ! root node
        real(kind=prec), intent(in) :: points(:, :)
        real(kind=prec), intent(in) :: targett(ndim) ! Target point (k-D)  
        !local
        real(kind=prec), intent(inout) :: q_dist(k)  
        integer(kind=intkind), intent(inout) :: q_idx(k) 
        integer(kind=intkind) :: i ! Index running over leaf points
        real(kind=prec):: dist_current, dist_furthest, d1d
        integer :: axis
        integer(kind=intkind) :: startidx, endidx
        logical :: look_opposite
        real(kind=prec):: temp_point(ndim)

        if (.not. associated(node)) return

        ! check if it is leaf
        if (node%is_leaf .eqv. .true.) then
            dist_furthest = q_dist(1)
            startidx = node%start_idx
            endidx = node%end_idx
            ! brute-force
            do i = startidx, endidx
                temp_point(1:3) = points(1:3, i)
                dist_current = distance(temp_point, targett)
                !add new (better) candidates
                if (dist_current < dist_furthest) then
                    q_dist(1) = dist_current
                    q_idx(1) = i
                    call max_heap_insert(q_dist, q_idx, k)
                    dist_furthest = q_dist(1)
                end if
            end do

        else 

            axis = node%axis
            ! 1D distance from target to the splitting plane
            d1d = targett(axis) - node%point

            ! Recursively search the subtree that contains the target
            if (d1d < 0) then
                call knn_search_recursive_hyperp(node%left, targett, points, q_idx, q_dist, k)
                dist_furthest = q_dist(1)
                !check if we need to search the other subtree
                look_opposite = .false.
                if (d1d**2 < dist_furthest) look_opposite = .true.
#if periodic == 1
                if (targett(axis) - sqrt(dist_furthest) <= -L(axis) / 2. ) look_opposite = .true.
#endif
                if (look_opposite .eqv. .true.) then
                    call knn_search_recursive_hyperp(node%right, targett, points, q_idx, q_dist, k)
                end if
            else
                call knn_search_recursive_hyperp(node%right, targett, points, q_idx, q_dist, k)
                dist_furthest = q_dist(1)
                !check if we need to search the other subtree
                look_opposite = .false.
                if (d1d**2 < dist_furthest) look_opposite = .true.
#if periodic == 1
                if (targett(axis) + sqrt(dist_furthest) >= L(axis) / 2. ) look_opposite = .true.
#endif
                if (look_opposite .eqv. .true.) then
                    call knn_search_recursive_hyperp(node%left, targett, points, q_idx, q_dist, k)
                end if
            end if


        end if !node%is_leaf
       
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end subroutine knn_search_recursive_hyperp
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Max-heap routine to insert new nearest neighbour candidates from the top
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subroutine max_heap_insert(q_dist, q_idx, k)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        implicit none
        ! Inputs
        integer, intent(in) :: k
        real(kind=prec), intent(inout) :: q_dist(k)
        integer(kind=intkind), intent(inout) :: q_idx(k)
        ! Local
        integer :: i, largest, left, right
        real(kind=prec) :: tmp_dist
        integer(kind=intkind) :: tmp_idx
    
        ! The new value has replaced the root value (furthest) at index 1
        ! dist(1), idx(1)

        ! Heapify Down from root to restore max-heap property (every parent bigger than its children)
        i = 1
        do
            left = 2 * i
            right = 2 * i + 1
            largest = i
    
            if (left <= k) then
                if (q_dist(left) > q_dist(largest)) largest = left
            end if
            if (right <= k) then
                if (q_dist(right) > q_dist(largest)) largest = right
            end if
    
            if (largest /= i) then
                !swap to correct position
                tmp_dist = q_dist(i)
                q_dist(i) = q_dist(largest)
                q_dist(largest) = tmp_dist
    
                tmp_idx = q_idx(i)
                q_idx(i) = q_idx(largest)
                q_idx(largest) = tmp_idx
    
                i = largest

            else !already fulfills max-heap property
                exit
            end if
        end do
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end subroutine max_heap_insert  
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Search for points within a given radius using bounding box pruning
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subroutine ball_search(node, targett, radius, points, max_size, q_idx, q_dist, count, sorted)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    implicit none
    !in
    real(kind=prec), intent(in) :: radius !ball radius
    type(KDTreeNode), pointer, intent(in) :: node
    real(kind=prec), intent(in) :: targett(ndim)
    real(kind=prec), intent(in) :: points(:, :)
    logical, intent(in), optional :: sorted
    !inout
    integer, intent(inout) :: max_size
    real(kind=prec), allocatable, intent(inout) :: q_dist(:)
    integer(kind=intkind), allocatable, intent(inout) :: q_idx(:)
    !out
    integer, intent(out) :: count
    
    !result buffers allocated
    if (.not. allocated(q_dist)) allocate(q_dist(max_size))
    if (.not. allocated(q_idx)) allocate(q_idx(max_size))
        
    count = 0

    call ball_search_recursive_hyperp(node, targett, points, q_idx, q_dist, count, radius)

    if (count > 0) then
        if ( present(sorted) ) then 
            if (sorted) call quicksort(q_dist, q_idx, count)
        endif
        q_dist(1:count) = sqrt(q_dist(1:count))
    end if

    max_size = size(q_dist)

    !do not resize, do not touch the buffer

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end subroutine ball_search
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
    recursive subroutine ball_search_recursive_hyperp(node, targett, points, q_idx, q_dist, count, radius)
    ! This versions uses hyperplanes to prune the search
    ! For moderate R is as fast or just a bit faster than the bbox version
    ! For large R is x1.5-2 slower than the bbox version
    ! works also with periodicity
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
        implicit none
        !in
        real(kind=prec):: radius 
        type(KDTreeNode), pointer, intent(in) :: node ! Starting node (usually the root)
        real(kind=prec), intent(in) :: targett(ndim) 
        real(kind=prec), intent(in) :: points(:, :)
        !out
        integer(kind=intkind), allocatable, intent(inout) :: q_idx(:) 
        real(kind=prec), allocatable, intent(inout) :: q_dist(:)
        integer, intent(inout) :: count
        !local
        integer(kind=intkind) :: i
        real(kind=prec):: dist_current, d1d
        integer :: axis
        integer(kind=intkind) :: startidx, endidx
        logical :: look_opposite 
        real(kind=prec) :: temp_point(ndim)

        if (.not. associated(node)) then
            return
        end if

        ! check if it is leaf
        if (node%is_leaf .eqv. .true.) then
            ! brute-force
            startidx = node%start_idx
            endidx = node%end_idx
            do i = startidx, endidx
                temp_point(1:3) = points(1:3, i)
                dist_current = distance(temp_point, targett)
                if ( dist_current <= radius**2) then
                    ! append new neighbors
                    call add_to_list(q_idx, q_dist, i, dist_current, count)
                end if
            end do

        else
            axis = node%axis
            ! 1D distance from target to the splitting plane
            d1d = targett(axis) - node%point
            ! Recursively search the subtree that contains the target
            if (d1d < 0) then
                call ball_search_recursive_hyperp(node%left, targett, points, q_idx, q_dist, count, radius)
                !check if we need to search the other subtree
                look_opposite = .false.
                if (abs(d1d) <= radius) look_opposite = .true.
#if periodic == 1
                if (targett(axis) - radius <= -L(axis) / 2. ) look_opposite = .true.
#endif
                if (look_opposite .eqv. .true.) then
                    call ball_search_recursive_hyperp(node%right, targett, points, q_idx, q_dist, count, radius)
                end if
            else
                call ball_search_recursive_hyperp(node%right, targett, points, q_idx, q_dist, count, radius)
                look_opposite = .false.
                !check if we need to search the other subtree
                if (abs(d1d) <= radius) look_opposite = .true.
#if periodic == 1
                if (targett(axis) + radius >= L(axis) / 2. ) look_opposite = .true.
#endif
                if (look_opposite .eqv. .true.) then
                    call ball_search_recursive_hyperp(node%left, targett, points, q_idx, q_dist, count, radius)
                end if
            end if

        endif !node%is_leaf
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end subroutine ball_search_recursive_hyperp
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !subroutines to append an element to an array, called by BALL_SEARCH and BOX_SEARCH
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subroutine add_to_list(q_idx, q_dist, idx_new, dist_new, count)
        implicit none
        integer(kind=intkind), allocatable, intent(inout) :: q_idx(:)
        real(kind=prec), allocatable, intent(inout) :: q_dist(:)
        integer(kind=intkind), intent(in) :: idx_new
        real(kind=prec), intent(in) :: dist_new
        integer, intent(inout) :: count
        integer(kind=intkind), allocatable :: temp_idx(:)
        real(kind=prec), allocatable :: temp_dist(:)
        integer :: n

        if (count >= size(q_idx)) then
            ! Resize the array
            n = size(q_idx)
            allocate(temp_idx(2 * n))
            allocate(temp_dist(2 * n))
            temp_idx(1:n) = q_idx
            temp_dist(1:n) = q_dist
            call move_alloc(temp_idx, q_idx)
            call move_alloc(temp_dist, q_dist)
        end if

        count = count + 1
        q_idx(count) = idx_new
        q_dist(count) = dist_new
        
    end subroutine add_to_list
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Euclidean distance between two points (Minkowski p=2)
    ! Takes into account periodicity if required
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function distance(p1, p2) result(dist)
        implicit none
        real(kind=prec), intent(in) :: p1(ndim), p2(ndim)
        real(kind=prec):: dist
        real(kind=prec):: dx(ndim)
        !local
        integer :: i

        do i=1,ndim
            dx(i) = p1(i) - p2(i)
        end do

#if periodic == 1
        do i=1,ndim
            dx(i) = min(abs(dx(i)), L(i) - abs(dx(i)))
        end do
#endif

        dist = 0.
        do i=1,ndim
            dist = dist + dx(i)*dx(i)
        end do
        ! dist = sqrt(dist)

    end function distance
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 1D (absolute) distance with periodicity if required
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function distance_1D(x1, x2, dim) result(dx)
        implicit none
        real(kind=prec), intent(in) :: x1, x2
        real(kind=prec):: dx
        integer, intent(in) :: dim

        dx = x1 - x2

#if periodic == 1
        dx = min(abs(dx), L(dim) - abs(dx))
#endif

        dx = abs(dx)
    end function distance_1D
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !Quicksort algorithm for sorting FINAL RESULTS
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subroutine quicksort(dist, idx, n)
        implicit none
        !in/out
        real(kind=prec), intent(inout) :: dist(n)   
        integer(kind=intkind), intent(inout) :: idx(n)
        integer, intent(in) :: n  
        !local
        integer :: low, high
        
        low = 1
        high = n

        call quicksort_recursive(dist, idx, low, high)

    end subroutine quicksort    
    
    recursive subroutine quicksort_recursive(dist, idx, low, high)
        implicit none
        !in/out
        real(kind=prec), intent(inout) :: dist(:)
        integer(kind=intkind), intent(inout) :: idx(:)
        integer, intent(in) :: low, high
        !local
        integer :: pivot_index

        if (low < high) then
            ! Partition the array and get the pivot index
            call partition2(dist, idx, low, high, pivot_index)

            ! Recursively sort the subarrays
            call quicksort_recursive(dist, idx, low, pivot_index - 1)
            call quicksort_recursive(dist, idx, pivot_index + 1, high)
        end if
 
    end subroutine quicksort_recursive

    subroutine partition2(dist, idx, low, high, pivot_index)
        implicit none
        !in/out
        real(kind=prec), intent(inout) :: dist(:)
        integer(kind=intkind), intent(inout) :: idx(:)
        integer, intent(in) :: low, high
        integer, intent(out) :: pivot_index
        !local
        integer :: mid
        real(kind=prec):: pivot_value
        integer :: i, j

        ! pivot: median of low, mid, high
        mid = (low + high) / 2

        ! find median
        if (dist(low) > dist(mid)) call swap2(dist, idx, low, mid)
        if (dist(low) > dist(high)) call swap2(dist, idx, low, high)
        if (dist(mid) > dist(high)) call swap2(dist, idx, mid, high)
        ! now dist(low) <= dist(mid) <= dist(high)
        pivot_value = dist(mid)
        ! move pivot to end
        call swap2(dist, idx, mid, high)
        i = low - 1
        ! do partitioning
        do j = low, high - 1
            if (dist(j) <= pivot_value) then
                i = i + 1
                call swap2(dist, idx, i, j)
            end if
        end do

        ! pivot to correct position
        i = i + 1
        call swap2(dist, idx, i, high)
        pivot_index = i

    end subroutine partition2

    subroutine swap2(dist, idx, i, j)
        implicit none
        real(kind=prec), intent(inout) :: dist(:)
        integer(kind=intkind), intent(inout) :: idx(:)
        integer, intent(in) :: i, j
        real(kind=prec):: temp_dist
        integer(kind=intkind) :: temp_idx

        temp_dist = dist(i)
        dist(i) = dist(j)
        dist(j) = temp_dist

        temp_idx = idx(i)
        idx(i) = idx(j)
        idx(j) = temp_idx
    end subroutine swap2
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!#######################################################
end module cosmokdtree
!#######################################################