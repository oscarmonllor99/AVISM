!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! cosmokdtree module
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! GRUP DE COSMOLOGIA COMPUTACIONAL (GCC) UNIVERSITAT DE VALÈNCIA
! Authors: Óscar Monllor Berbegal and David Vallés Pérez
! Date: Genuary 30th 2025
! Last update: December 5th 2025
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Brief description:
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! - This module implements an optimized parallel k-d tree construction and 
!   fast search for k-nearest neighbors and range queries.
!
! - Sliding midpoint splitting rule is used. Boundind boxes are stored
!   for faster search (faster than hyperplane distance checks).
!
! - The tree is built in parallel using OpenMP tasks due 
!   to Divide and Conquer nature of the algorithm. Nested parallelism
!   must be enabled in OpenMP to allow parallelism at the top levels.
!
! - The search for k-nearest neighbors uses a Max-Heap storage
!   to efficiently keep track of the k nearest points and a Min-Heap together
!   with a priority queue for best-first traversal of the tree.
!   
! - Ball search leverages dynamic arrays to store results and bounding box
!   checks for efficient traversal.
! 
! - Quicksort is used to sort the distances and indices of    
!   points in the final knn_search max_heap and range queries results.
!   We have seen that even having an already existing max_heap structure,
!   quicksort outperforms heapsort by a 10-20% factor.
!
! - For knn (ball) queries, the results are sorted (unsorted) by default.
!
! - INT*8 indices and REAL*8 points are allowed through conditional compilation
!
! - Similarly, periodic boundary conditions are allowed
!
! - Dimensionality is set to 3 by default, but can be changed to an
!   arbitrary value (e.g. 6D phase space).
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Pending improvements:
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! - Better parallelism (scaling) for tree building (if possible)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!#######################################################
module cosmokdtree
!#######################################################
    implicit none
    private  
    public :: build_kdtree, deallocate_kdtree, & 
              knn_search, ball_search, box_search, &
              KDTreeNode, KDTreeResult

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
        real(kind=prec) :: point(ndim) ! Splitting point coordinates
        integer :: axis ! Splitting axis (1 for x, 2 for y, 3 for z, 4 for w, ...)
        type(KDTreeNode), pointer :: left => null()  ! Left child node
        type(KDTreeNode), pointer :: right => null() ! Right child node
        real(kind=prec) :: maxbounds(ndim)
        real(kind=prec) :: minbounds(ndim)
        !leaf, for faster search and building
        logical :: is_leaf  ! Flag to indicate if the node is a leaf
        real(kind=prec), pointer :: leaf_points(:, :) => null()  ! Points in the leaf (for leaf nodes)
        integer(kind=intkind), pointer :: leaf_indices(:) => null()  ! Indices of points in the leaf
    end type KDTreeNode

    type :: KDTreeResult
        integer(kind=intkind), allocatable :: idx(:)
        real(kind=prec), allocatable :: dist(:)
    end type KDTreeResult

    ! type: array of nodes for priority queue
    type :: node_queue
        type(KDTreeNode), pointer :: node 
    end type node_queue

    !+++++++++++++++++++++++++++++++

contains

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Initialize kd-tree construction
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#if periodic == 1
    function build_kdtree(points_in, L_in, leaf) result(tree)
#else
    function build_kdtree(points_in, leaf) result(tree)
#endif
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        use omp_lib
        implicit none
        !in
        real(kind=prec), intent(in) :: points_in(:, :)
        integer, intent(in), optional :: leaf 
#if periodic == 1
        real(kind=prec), intent(in) :: L_in(:)
        integer :: flag_stop
        integer :: j, ndim_par
#endif
        !local
        real(kind=prec), allocatable :: points(:, :)
        integer(kind=intkind), allocatable :: indices(:)
        integer(kind=intkind) :: n, i
        integer :: depth, max_depth, nproc, leafsize
        real(kind=prec) :: bounds(2*ndim)
        
        !out
        type(KDTreeNode), pointer :: tree

        ! Enable nested parallelism
        call omp_set_nested(.true.) 

        ! Check dimensionality of input
        if (size(points_in, 2) /= ndim) then
            STOP 'Input points must have the same dimensionality as the tree!'
        end if

        ! Number of points
        n = size(points_in, 1, kind=intkind)

#if periodic == 1
        ! Check dimensionality of input
        if (size(L_in, 1) /= ndim) then
            STOP 'Input L must have the same dimensionality as the tree!'
        end if
        
        ! Set the box size: (xmin, xmax, ymin, ymax, zmin, zmax, ...)
        L = L_in

        ! Check all points are within (-L/2, L/2)
        flag_stop = 0
        ndim_par = ndim
        !$OMP PARALLEL SHARED(points,n,L,ndim_par), &
        !$OMP PRIVATE(i,j)
        !$OMP DO REDUCTION(+:flag_stop)
        do j=1,ndim_par
            do i=1,n
                if (points_in(i,j) .lt. -L(j)/2 .or. points_in(i,j) .gt. L(j)/2) flag_stop = 1
            enddo
        enddo
        !$OMP ENDDO
        !$OMP END PARALLEL
        
        if (flag_stop .gt. 0) STOP 'Points outside (-L/2, L/2) range !!'
#endif

        ! Initialize global indices and points
        allocate(points(n, ndim))
        points = points_in
        allocate(indices(n))
        indices = [(i, i=1, n)]
        
        ! Init depth = 0
        depth = 0

        !$OMP PARALLEL
        !$OMP SINGLE
        nproc = omp_get_num_threads()
        !$OMP END SINGLE
        !$OMP END PARALLEL

        ! Build KD-tree
        max_depth = compute_max_depth(omp_get_max_threads())

        ! Leafsize
        if ( present(leaf) ) then 
            leafsize = leaf
        else
            leafsize = 16
        endif

        ! Initial bounds for the tree:
#if periodic == 1
        bounds(1) = -L(1)/2.
        bounds(2) =  L(1)/2.
        bounds(3) = -L(2)/2.
        bounds(4) =  L(2)/2.
        bounds(5) = -L(3)/2.
        bounds(6) =  L(3)/2.
#else
        do i = 1, ndim
            bounds(2*i-1) = minval(points(:, i))  ! Min value for dimension i
            bounds(2*i) = maxval(points(:, i))   ! Max value for dimension i
        end do
#endif

        tree => build_kdtree_recursive(points, indices, depth, max_depth, leafsize, bounds)

        deallocate(points, indices)
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

        ! Deallocate left and right subtrees
        call deallocate_kdtree(node%left)
        call deallocate_kdtree(node%right)

        ! Deallocate leaf points and indices if leaf node
        if (node%is_leaf .eqv. .true.) then
            if (associated(node%leaf_points)) then
                deallocate(node%leaf_points)
            end if
            if (associated(node%leaf_indices)) then
                deallocate(node%leaf_indices)
            end if
        end if

        ! Deallocate the node itself
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
    ! Recursive kd-tree building function according to sliding midpoint rule
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    recursive function build_kdtree_recursive(points,indices,depth,max_depth,leafsize,bounds) result(node)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    use omp_lib
    implicit none
    !in
    real(kind=prec), intent(inout) :: points(:, :) ! 2D array of points
    integer(kind=intkind), intent(inout) :: indices(:) ! 1D array of indices
    integer, intent(in) :: depth ! Current depth in the tree
    integer, intent(in) :: max_depth ! max_depth to allow parallelism
    integer, intent(in) :: leafsize ! size of leaf nodes
    real(kind=prec), intent(in) :: bounds(:)
    !local
    real(kind=prec) :: maxside, minx, maxx
    real(kind=prec) :: bounds_left(ndim*2), bounds_right(ndim*2)
    real(kind=prec) :: minvals(ndim), maxvals(ndim), spreads(ndim), side(ndim)
    real :: midpoint
    integer :: axis, j !axis to split points
    integer(kind=intkind) :: i, size_points, n_left ! Size of the data to be split
    type(KDTreeNode), pointer :: node   ! New node to be created

    ! Allocate node
    allocate(node)

    !Due to the building algorithm, no empty nodes can be created.
    size_points = size(points, 1)

    !Assign bounding box for this node
    do j=1,ndim
        node%minbounds(j) = bounds(2*j-1)
        node%maxbounds(j) = bounds(2*j)
    enddo

    ! Check if this is a leaf node
    if (size_points <= leafsize) then
        node%is_leaf = .true.
        ! Store all points and indices in the leaf node
        ! Leaf points are stored (ndim, n) for contiguous memory access in query
        allocate(node%leaf_points(ndim, size(points, 1, kind=intkind)))
        allocate(node%leaf_indices(size(points, 1, kind=intkind)))
        do i = 1, size_points
            do j = 1, ndim
                node%leaf_points(j,i) = points(i,j)
            enddo
            node%leaf_indices(i) = indices(i)
        enddo
        node%left => null()
        node%right => null()
        return
    else
        node%is_leaf = .false.
    end if

    ! Compute necessary statistics to choose splitting axis and point
    side = 0.
    do j = 1, ndim
        minvals(j) = minval(points(:, j))
        maxvals(j) = maxval(points(:, j))
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
    
    ! Midpoint: split bounding box in half across the longest side
    midpoint = 0.5 * (bounds(2*axis-1) + bounds(2*axis))

    ! Check if data lies in both sides of the hyperplane
    minx = minvals(axis)
    maxx = maxvals(axis)

    !slide the hyperplane to the closest point 
    if (minx > midpoint) then
        midpoint = minx
        call partition(points, indices, midpoint, axis, size_points, n_left)
        n_left = 1 !at least the smallest one is on the left side
    else if (maxx < midpoint) then
        midpoint = maxx
        call partition(points, indices, midpoint, axis, size_points, n_left)
        n_left = size_points - 1 !at least the biggest one is on the right side
    else
        !default: midpoint splitting rule
        call partition(points, indices, midpoint, axis, size_points, n_left)
    end if
    ! all points with points(:,axis) < midpoint are in the left side
    ! all points with points(:,axis) >= midpoint are in the right side
    
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
    node%left => build_kdtree_recursive(points(1:n_left,:),indices(1:n_left),depth+1,max_depth,leafsize, &
                                        bounds_left)
    !$OMP END TASK

    !$OMP TASK
    node%right => build_kdtree_recursive(points(n_left+1:size_points,:), &
                                        indices(n_left+1:size_points), depth+1, max_depth, leafsize, &
                                        bounds_right)
    !$OMP END TASK

    !$OMP END SINGLE
    !$OMP END PARALLEL
    else
    node%left => build_kdtree_recursive(points(1:n_left,:),indices(1:n_left),depth+1,max_depth,leafsize, &
                                        bounds_left)
    node%right => build_kdtree_recursive(points(n_left+1:size_points,:), &
                                        indices(n_left+1:size_points), depth+1, max_depth, leafsize, &
                                        bounds_right)
    end if

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end function build_kdtree_recursive
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Partition function to split points around a pivot value along a given axis
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subroutine partition(points, indices, pivot_value, axis, size_points, n_left) 
        implicit none
        real(kind=prec), intent(inout) :: points(:, :)  ! 2D array of points
        integer(kind=intkind), intent(inout) :: indices(:) ! 1D array of indices
        integer, intent(in) :: axis
        integer(kind=intkind), intent(in) :: size_points ! Size of the data to be partitioned
        integer(kind=intkind) :: i, j
        real(kind=prec):: pivot_value
        integer(kind=intkind) :: n_left
        
        n_left = 0
        do i = 1, size_points
            if (points(i, axis) < pivot_value) then
                j = n_left + 1
                if (i /= j) then
                    call swap(points, indices, i, j)
                end if
                n_left = n_left + 1
            end if
        end do

    end subroutine partition

    subroutine swap(points, indices, i, j)
        implicit none
        real(kind=prec), intent(inout) :: points(:, :)
        integer(kind=intkind), intent(inout) :: indices(:)
        integer(kind=intkind), intent(in) :: i, j
        real(kind=prec):: temp_point(ndim)
        integer(kind=intkind) :: temp_index
    
        ! Swap points
        temp_point = points(i, :)
        points(i, :) = points(j, :)
        points(j, :) = temp_point
    
        ! Swap indices
        temp_index = indices(i)
        indices(i) = indices(j)
        indices(j) = temp_index
    end subroutine swap
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! (USED) k-nearest neighbor search main function
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function knn_search(node, targett, k, sorted) result(query)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        implicit none
        !in
        integer, intent(in) :: k ! Number of nearest neighbors to find
        type(KDTreeNode), pointer, intent(in) :: node
        real(kind=prec), intent(in) :: targett(ndim) ! Target point (k-D)
        logical, intent(in), optional :: sorted ! If true, sort the points by distance to the target
        !out
        real(kind=prec):: dist(k)
        integer(kind=intkind) :: idx(k)
        type(KDTreeResult) :: query

        !Initialize 
        dist = HUGE(0.0)
        idx = -1

        call knn_search_queue(node, targett, idx, dist, k)

        !Perform full sort with quicksort
        ! by default, sorts
        if ( present(sorted) ) then 
            if (sorted) then
                call quicksort(dist, idx, k)
            end if
        else
            call quicksort(dist, idx, k)
        endif

        query%idx = idx
        query%dist = sqrt(dist)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    endfunction knn_search
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! k-nearest neighbor search with priority queue and bounding box checks
    ! Leverages a combined max-heap and min-heap for efficient search
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subroutine knn_search_queue(root_node, targett, idx, dist, k)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    implicit none
    !in
    integer :: k ! Number of nearest neighbors to find
    type(KDTreeNode), pointer, intent(in) :: root_node ! Starting node (usually the root)
    real(kind=prec), intent(inout) :: dist(k)  
    integer(kind=intkind), intent(inout) :: idx(k) 
    real(kind=prec), intent(in) :: targett(ndim) ! Target point (k-D)  
    !local
    integer :: i
    integer :: axis
    !min-heap priority queue
    integer :: Nqueue ! Number of nodes in the queue
    integer :: Nstack, Nstack2
    type(KDTreeNode), pointer :: node
    type(node_queue), allocatable :: queue(:), tmp_queue(:)
    real(kind=prec), allocatable :: nodes_mindist(:), tmp_nodes_mindist(:) 
    real(kind=prec) :: mindistkd, dist_furthest, mindistkd_right, mindistkd_left
    real(kind=prec) :: mindist_save(ndim), maxdist_save(ndim)
    real(kind=prec):: dist_current
    ! Temporary point for contiguous memory access
    real(kind=prec):: temp_point(ndim)


    ! Allocate queue. We estimate that, at least, "k" leaf nodes have to be visited.
    ! Initialize to null.
    Nstack = 100
    allocate(queue(Nstack))
    allocate(nodes_mindist(Nstack))
    nodes_mindist(:) = HUGE(0.0)

    ! Initialize queue with the root node
    Nqueue = 1
    queue(1)%node => root_node
    call bbox_mindist_kd(root_node%minbounds, root_node%maxbounds, targett, mindistkd)
    nodes_mindist(1) = mindistkd

    !########################################
    do while (Nqueue > 0) ! Priority search
    !########################################

        ! Stack resize if priority queue is almost full (and thus children may overflow it)
        if (Nqueue > Nstack - 1) then
            ! Reallocate properly
            Nstack2 = 10 * Nstack
            allocate(tmp_queue(Nstack2))
            allocate(tmp_nodes_mindist(Nstack2))
            tmp_nodes_mindist(1:Nstack) = nodes_mindist(1:Nstack)
            do i = 1, Nstack
                tmp_queue(i)%node => queue(i)%node
            end do
            call move_alloc(tmp_queue, queue)
            call move_alloc(tmp_nodes_mindist, nodes_mindist)
            Nstack = Nstack2
        end if

        ! Current node to process
        node => queue(1)%node
        mindistkd = nodes_mindist(1)

        ! Pop this node from the list -> move last to the first position
        queue(1)%node => queue(Nqueue)%node
        nodes_mindist(1) = nodes_mindist(Nqueue)
        Nqueue = Nqueue - 1

        ! Restore min-heap from the top
        if (Nqueue > 0) then
            call min_heap_queue_down(queue(1:Nqueue), nodes_mindist(1:Nqueue), Nqueue)
        end if

        if (.not. associated(node)) then
            cycle
        end if
    
        ! Furthest "nearest-neighbor" distance found so far
        dist_furthest = dist(1)

        ! Early exit of this branch if cannot improve current nearest neighbors
        if (mindistkd > dist_furthest) then
            cycle
        end if


        ! If it is not a leaf but was closer than the furthest in the queue,
        ! check if left and right children can improve the nearest neighbors
        if (node%is_leaf .eqv. .false.) then

            ! Process children nodes:
            axis = node%axis

            call bbox_mindist_kd(node%left%minbounds, node%left%maxbounds, &
                                    targett, mindistkd_left)
            call bbox_mindist_kd(node%right%minbounds, node%right%maxbounds, &
                                    targett, mindistkd_right)
            
            ! Push LEFT child 
            if (mindistkd_left <= dist_furthest) then
                Nqueue = Nqueue + 1
                queue(Nqueue)%node => node%left
                nodes_mindist(Nqueue) = mindistkd_left
                ! Restore min-heap property from the bottom
                call min_heap_queue_up(queue(1:Nqueue), nodes_mindist(1:Nqueue), Nqueue)
            end if

            ! Push RIGHT child
            if (mindistkd_right <= dist_furthest) then
                Nqueue = Nqueue + 1
                queue(Nqueue)%node => node%right
                nodes_mindist(Nqueue) = mindistkd_right
                ! Restore min-heap property from the bottom
                call min_heap_queue_up(queue(1:Nqueue), nodes_mindist(1:Nqueue), Nqueue)
            end if

        ! If it is a leaf, then search all points and use max_heap
        ! to keep track of candidates
        else if (node%is_leaf .eqv. .true.) then

            ! Check all points in the leaf with brute force
            do i = 1, size(node%leaf_indices)
                temp_point = node%leaf_points(:, i)
                dist_current = distance(temp_point, targett)
                ! If the current point is closer than the furthest, replace it
                if (dist_current < dist_furthest) then
                    dist(1) = dist_current
                    idx(1) = node%leaf_indices(i)
                    call max_heap_insert(dist, idx, k)
                    dist_furthest = dist(1)
                end if
            end do

        end if !node%is_leaf

    !##################################
    end do ! Priority search
    !##################################

    deallocate(queue, nodes_mindist)

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end subroutine knn_search_queue
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! k-nearest neighbor search (slower version with hyperplane distance checks and no
    ! priority queue)
    ! DEPRECATED: use knn_search
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function knn_search_old(node, targett, k, sorted) result(query)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        implicit none
        !in
        integer, intent(in) :: k ! Number of nearest neighbors to find
        type(KDTreeNode), pointer, intent(in) :: node
        real(kind=prec), intent(in) :: targett(ndim) ! Target point (k-D)
        logical, intent(in), optional :: sorted ! If true, sort the points by distance to the target
        !out
        real(kind=prec):: dist(k)
        integer(kind=intkind) :: idx(k)
        type(KDTreeResult) :: query
        !local
        integer :: i
        real(kind=prec) :: mindist(ndim), maxdist(ndim)

        !Initialize 
        dist = HUGE(0.0)
        idx = -1

        !mindist -> minimum distance between targett and node bbox along each dimension
        !maxdist -> maximum distance "
        mindist = 0.
        maxdist = 0.
        do i = 1, ndim
            call bbox_distance_1D(node%minbounds, node%maxbounds, targett, mindist, maxdist, i)
        enddo

        call knn_search_recursive_hybrid(node, targett, dist, idx, mindist, maxdist, k)
        !call knn_search_recursive_hyperp(node, targett, dist, idx, k)

        !Perform full sort with quicksort
        ! by default, sorts
        if ( present(sorted) ) then 
            if (sorted) then
                call quicksort(dist, idx, k)
            end if
        else
            call quicksort(dist, idx, k)
        endif

        query%idx = idx
        query%dist = sqrt(dist)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    endfunction knn_search_old
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Hyperplane method with early bbox prunning (just a bit faster than pure hyperplane)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    recursive subroutine knn_search_recursive_hybrid(node, targett, dist, idx, mindist, maxdist, k)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        implicit none
        !in
        integer :: k ! Number of nearest neighbors to find
        type(KDTreeNode), pointer, intent(in) :: node ! Starting node (usually the root)
        real(kind=prec), intent(in) :: targett(ndim) ! Target point (k-D)  
        real(kind=prec), intent(inout) :: mindist(ndim) ! Minimum distance to bbox in each dimension
        real(kind=prec), intent(inout) :: maxdist(ndim) ! Maximum distance to bbox
        !local
        real(kind=prec), intent(inout) :: dist(k)  
        integer(kind=intkind), intent(inout) :: idx(k) 
        integer :: i ! Index running over leaf points
        real(kind=prec):: dist_current, dist_furthest, d1d
        integer :: axis
        logical :: look_opposite
        real(kind=prec) :: mindist3D, maxdist3D
        real(kind=prec) :: mindist_save(ndim), maxdist_save(ndim)
        ! Temporary point for contiguous memory access
        real(kind=prec):: temp_point(ndim)

        if (.not. associated(node)) then
            return
        end if

        ! Early exit
        dist_furthest = dist(1)
        mindist3D = sum(mindist**2)
        if (mindist3D > dist_furthest) then
            return
        end if

        ! First, check if it is a leaf node
        if (node%is_leaf .eqv. .true.) then
            dist_furthest = dist(1)
            ! Check all points in the leaf with brute force
            do i = 1, size(node%leaf_indices)
                temp_point = node%leaf_points(:, i)
                dist_current = distance(temp_point, targett)
                ! If the current point is closer than the furthest, replace it
                if (dist_current < dist_furthest) then
                    dist(1) = dist_current
                    idx(1) = node%leaf_indices(i)
                    call max_heap_insert(dist, idx, k)
                    dist_furthest = dist(1)
                end if
            end do

        else 

            axis = node%axis
            ! 1D distance from target to the splitting plane
            d1d = targett(axis) - node%point(axis)

            ! Recursively search the subtree that contains the target
            if (d1d < 0) then
                mindist_save = mindist
                maxdist_save = maxdist
                call bbox_distance_1D(node%left%minbounds, node%left%maxbounds, &
                                        targett, mindist, maxdist, axis)
                call knn_search_recursive_hybrid(node%left, targett, dist, idx, &
                                                mindist, maxdist, k)
                dist_furthest = dist(1)
                !Check if we need to search the right subtree 
                look_opposite = .false.
                if (d1d**2 < dist_furthest) look_opposite = .true.
#if periodic == 1
                if (targett(axis) - dist_furthest <= -L(axis) / 2. ) look_opposite = .true.
#endif
                if (look_opposite .eqv. .true.) then
                    mindist = mindist_save
                    maxdist = maxdist_save
                    call bbox_distance_1D(node%right%minbounds, node%right%maxbounds, &
                                            targett, mindist, maxdist, axis)
                    call knn_search_recursive_hybrid(node%right, targett, dist, idx, &
                                                    mindist, maxdist, k)
                end if
            else
                mindist_save = mindist
                maxdist_save = maxdist
                call bbox_distance_1D(node%right%minbounds, node%right%maxbounds, &
                                        targett, mindist, maxdist, axis)
                call knn_search_recursive_hybrid(node%right, targett, dist, idx, &
                                                mindist, maxdist, k)
                dist_furthest = dist(1)
                !Check if we need to search the left subtree
                look_opposite = .false.
                if (d1d**2 < dist_furthest) look_opposite = .true.
#if periodic == 1
                if (targett(axis) + dist_furthest >= L(axis) / 2. ) look_opposite = .true.
#endif
                if (look_opposite .eqv. .true.) then
                    mindist = mindist_save
                    maxdist = maxdist_save
                    call bbox_distance_1D(node%left%minbounds, node%left%maxbounds, &
                                            targett, mindist, maxdist, axis)
                    call knn_search_recursive_hybrid(node%left, targett, dist, idx, &
                                                    mindist, maxdist, k)
                end if
            end if


        end if !node%is_leaf
       
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end subroutine knn_search_recursive_hybrid
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Standard hyperplane method to find k-nearest neighbors (SIMPLEST)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    recursive subroutine knn_search_recursive_hyperp(node, targett, dist, idx, k)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        implicit none
        !in
        integer :: k ! Number of nearest neighbors to find
        type(KDTreeNode), pointer, intent(in) :: node ! Starting node (usually the root)
        real(kind=prec), intent(in) :: targett(ndim) ! Target point (k-D)  
        !local
        real(kind=prec), intent(inout) :: dist(k)  
        integer(kind=intkind), intent(inout) :: idx(k) 
        integer :: i ! Index running over leaf points
        real(kind=prec):: dist_current, dist_furthest, d1d
        integer :: axis
        logical :: look_opposite
        ! Temporary point for contiguous memory access
        real(kind=prec):: temp_point(ndim)

        if (.not. associated(node)) return

        ! First, check if it is a leaf node
        if (node%is_leaf .eqv. .true.) then

            dist_furthest = dist(1)
            ! Check all points in the leaf with brute force
            do i = 1, size(node%leaf_indices)
                temp_point = node%leaf_points(:, i)
                dist_current = distance(temp_point, targett)
                ! If the current point is closer than the furthest, replace it
                if (dist_current < dist_furthest) then
                    dist(1) = dist_current
                    idx(1) = node%leaf_indices(i)
                    call max_heap_insert(dist, idx, k)
                    dist_furthest = dist(1)
                end if
            end do

        else 

            axis = node%axis
            ! 1D distance from target to the splitting plane
            d1d = targett(axis) - node%point(axis)

            ! Recursively search the subtree that contains the target
            if (d1d < 0) then
                call knn_search_recursive_hyperp(node%left, targett, dist, idx, k)
                dist_furthest = dist(1)
                !Check if we need to search the right subtree 
                look_opposite = .false.
                if (d1d**2 < dist_furthest) look_opposite = .true.
#if periodic == 1
                if (targett(axis) - dist_furthest <= -L(axis) / 2. ) look_opposite = .true.
#endif
                if (look_opposite .eqv. .true.) then
                    call knn_search_recursive_hyperp(node%right, targett, dist, idx, k)
                end if
            else
                call knn_search_recursive_hyperp(node%right, targett, dist, idx, k)
                dist_furthest = dist(1)
                !Check if we need to search the left subtree
                look_opposite = .false.
                if (d1d**2 < dist_furthest) look_opposite = .true.
#if periodic == 1
                if (targett(axis) + dist_furthest >= L(axis) / 2. ) look_opposite = .true.
#endif
                if (look_opposite .eqv. .true.) then
                    call knn_search_recursive_hyperp(node%left, targett, dist, idx, k)
                end if
            end if


        end if !node%is_leaf
       
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end subroutine knn_search_recursive_hyperp
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Max-heap routine to insert new nearest neighbour candidates from the top
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subroutine max_heap_insert(dist, idx, k)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        implicit none
        ! Inputs
        integer, intent(in) :: k
        real(kind=prec), intent(inout) :: dist(k)
        integer(kind=intkind), intent(inout) :: idx(k)
        ! Local
        integer :: i, largest, left, right
        real(kind=prec) :: tmp_dist
        integer(kind=intkind) :: tmp_idx
    
        ! The new value has replaced the root value (furthest) at index 1
        ! dist(1), idx(1)

        ! Heapify down from root to restore max-heap property (every parent bigger than its children)
        i = 1
        do
            left = 2 * i
            right = 2 * i + 1
            largest = i
    
            if (left <= k) then
                if (dist(left) > dist(largest)) largest = left
            end if
            if (right <= k) then
                if (dist(right) > dist(largest)) largest = right
            end if
    
            if (largest /= i) then
                ! Swap i-th element with largest
                tmp_dist = dist(i)
                dist(i) = dist(largest)
                dist(largest) = tmp_dist
    
                tmp_idx = idx(i)
                idx(i) = idx(largest)
                idx(largest) = tmp_idx
    
                i = largest

            else !already fulfills max-heap property
                exit
            end if
        end do
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end subroutine max_heap_insert  
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Min-heap routine to restore the heap property from the top according to nodes_mindist
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subroutine min_heap_queue_down(queue, nodes_mindist, k)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        implicit none
        ! Inputs
        integer, intent(in) :: k
        real(kind=prec), intent(inout) :: nodes_mindist(k)
        type(node_queue), intent(inout) :: queue(k)
        ! Local
        integer :: i, smallest, left, right
        real(kind=prec) :: tmp_dist
        type(KDTreeNode), pointer :: tmp_ptr
    
        ! The new value has replaced the root value (min) at index 1
        ! queue(1), nodes_mindist(1)

        ! Heapify down from root to restore min-heap property (every parent smaller than its children)
        i = 1
        do
            left = 2 * i
            right = 2 * i + 1
            smallest = i
    
            if (left <= k) then
                if (nodes_mindist(left) < nodes_mindist(smallest)) smallest = left
            end if
            if (right <= k) then
                if (nodes_mindist(right) < nodes_mindist(smallest)) smallest = right
            end if
    
            if (smallest /= i) then
                ! Swap i-th element with smallest
                tmp_dist = nodes_mindist(i)
                nodes_mindist(i) = nodes_mindist(smallest)
                nodes_mindist(smallest) = tmp_dist
                
                tmp_ptr => queue(i)%node
                queue(i)%node => queue(smallest)%node
                queue(smallest)%node => tmp_ptr
    
                i = smallest

            else !already fulfills min-heap property
                exit
            end if
        end do
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end subroutine min_heap_queue_down
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Min-heap routine to restore the heap property from the bottom according to nodes_mindist
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subroutine min_heap_queue_up(queue, nodes_mindist, k)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        implicit none
        ! Inputs
        integer, intent(in) :: k
        real(kind=prec), intent(inout) :: nodes_mindist(k)
        type(node_queue), intent(inout) :: queue(k)
        ! Local
        integer :: i, father, largest
        real(kind=prec) :: tmp_dist
        type(KDTreeNode), pointer :: tmp_ptr
    
        ! The new value has replaced the last value at index k
        ! queue(k), nodes_mindist(k)

        ! Heapify up from last to root to restore min-heap property (every parent smaller than its children)
        i = k
        do
            father = i / 2
            largest = i
            
            if (father >= 1) then
                if (nodes_mindist(father) > nodes_mindist(largest)) largest = father
            end if

            if (largest /= i) then
                ! Swap i-th element with smallest
                tmp_dist = nodes_mindist(i)
                nodes_mindist(i) = nodes_mindist(largest)
                nodes_mindist(largest) = tmp_dist
                
                tmp_ptr => queue(i)%node
                queue(i)%node => queue(largest)%node
                queue(largest)%node => tmp_ptr
    
                i = largest

            else !already fulfills min-heap property
                exit
            end if
            
        end do
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end subroutine min_heap_queue_up
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Search for points within a given radius using bounding box pruning
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function ball_search(node, targett, radius, sorted) result(query)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    implicit none
    !in
    real(kind=prec):: radius ! Radius of the ball
    type(KDTreeNode), pointer, intent(in) :: node
    real(kind=prec), intent(in) :: targett(ndim) ! Target point (k-D)
    logical, intent(in), optional :: sorted ! If true, sort the points by distance to the target
    !local
    integer :: i, max_size
    integer :: count
    integer(kind=intkind), allocatable :: temp_idx(:)
    real(kind=prec), allocatable :: temp_dist(:)
    real(kind=prec) :: mindist(ndim), maxdist(ndim)
    !out
    real(kind=prec), allocatable :: dist(:) ! Distance of the points within the radius
    integer(kind=intkind), allocatable :: idx(:) !index of the points within the radius
    type(KDTreeResult) :: query

    max_size = 100
    ! Later dist and idx will be reallocated to the correct size
    allocate(dist(max_size)) 
    allocate(idx(max_size))
    dist = 0
    idx = 0
    count = 0

    !mindist -> minimum distance between targett and node bbox along each dimension
    !maxdist -> maximum distance "
    mindist = 0.
    maxdist = 0.
    do i = 1, ndim
        call bbox_distance_1D(node%minbounds, node%maxbounds, targett, mindist, maxdist, i)
    enddo

    call ball_search_recursive_bbox(node, targett, dist, idx, count, mindist, maxdist, radius)
  
    if (count == 0) then
        ! No points found
        temp_dist = dist(1:0)
        call move_alloc(temp_dist, dist)
        temp_idx = idx(1:0)
        call move_alloc(temp_idx, idx)

    else
        ! Reallocation to the correct size
        temp_dist = dist(1:count)
        call move_alloc(temp_dist, dist)
        temp_idx = idx(1:count)
        call move_alloc(temp_idx, idx)
        ! Last step, sort the distances (slightly decreases performance)
        if ( present(sorted) ) then 
            if (sorted) then
                call quicksort(dist, idx, count)
            end if
        endif
    end if

    query%idx = idx
    query%dist = sqrt(dist)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end function ball_search
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
    recursive subroutine ball_search_recursive_bbox(node, targett, dist, idx, count, mindist, maxdist, radius)
    ! This version uses bounding boxes to prune the search
    ! For large R is faster than the hyperplane version
    ! For small-moderate R is similar to the hyperplane version or just a bit slower
    ! Periodicity is taken into account inside the bbox distance calculations
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
        implicit none
        !in
        real(kind=prec):: radius ! Radius of the ball
        type(KDTreeNode), pointer, intent(in) :: node ! Starting node (usually the root)
        real(kind=prec), intent(in) :: targett(ndim) ! Target point (k-D)
        real(kind=prec) :: mindist(ndim), maxdist(ndim)
        !out
        integer(kind=intkind), allocatable, intent(inout) :: idx(:)  ! Index of the points within the radius
        real(kind=prec), allocatable, intent(inout) :: dist(:)
        integer, intent(inout) :: count
        !local
        real(kind=prec) :: mindist_save(ndim), maxdist_save(ndim)
        real(kind=prec) :: mindist3D, maxdist3D
        integer :: i
        integer :: axis
        real(kind=prec):: dist_current
        ! Temporary point for contiguous memory access
        real(kind=prec) :: temp_point(ndim)

        if (.not. associated(node)) then
            return
        end if

        !Calculate 3D mindist and maxdist
        mindist3D = sum(mindist**2)
        maxdist3D = sum(maxdist**2)

        !Early prunning
        if (mindist3D > radius**2) then
            return
            
        !Completely inside the ball
        else if (maxdist3D < radius**2) then
            call ball_reduction_recursive(node, targett, dist, idx, count)
            return
        end if

        ! If just intersecting, proceed as usual
        if (node%is_leaf) then
            ! Check all points in the leaf with brute force
            do i = 1, size(node%leaf_indices)
                temp_point = node%leaf_points(:, i)
                dist_current = distance(temp_point, targett)
                if ( dist_current <= radius**2) then
                    ! Append the index to the list
                    call add_to_list(idx, dist, node%leaf_indices(i), dist_current, count)
                end if
            end do

        ! Then check which childs to visit
        else
            axis = node%axis

            !Save mindist, maxdist for right child
            mindist_save = mindist
            maxdist_save = maxdist

            call bbox_distance_1D(node%left%minbounds, node%left%maxbounds, targett, & 
                                    mindist, maxdist, axis)
            call ball_search_recursive_bbox(node%left, targett, dist, idx, count, mindist, maxdist, radius)

            mindist = mindist_save
            maxdist = maxdist_save
            
            call bbox_distance_1D(node%right%minbounds, node%right%maxbounds, targett, & 
                                    mindist, maxdist, axis)
            call ball_search_recursive_bbox(node%right, targett, dist, idx, count, mindist, maxdist, radius)

        endif

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end subroutine ball_search_recursive_bbox
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
    recursive subroutine ball_search_recursive_hyperp(node, targett, dist, idx, radius, count)
    ! This versions uses hyperplanes to prune the search
    ! For moderate R is as fast or just a bit faster than the bbox version
    ! For large R is x2 slower than the bbox version
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
        implicit none
        !in
        real(kind=prec):: radius ! Radius of the ball
        type(KDTreeNode), pointer, intent(in) :: node ! Starting node (usually the root)
        real(kind=prec), intent(in) :: targett(ndim) ! Target point (k-D)
        !out
        integer(kind=intkind), allocatable, intent(inout) :: idx(:)  ! Index of the points within the radius
        real(kind=prec), allocatable, intent(inout) :: dist(:)
        integer, intent(inout) :: count
        !local
        integer :: i
        real(kind=prec):: dist_current, d1d
        integer :: axis
        logical :: look_opposite 
        ! Temporary point for contiguous memory access
        real(kind=prec) :: temp_point(ndim)

        if (.not. associated(node)) then
            return
        end if

        ! If just intersecting, proceed as usual
        if (node%is_leaf .eqv. .true.) then
            ! Check all points in the leaf with brute force
            do i = 1, size(node%leaf_indices)
                temp_point = node%leaf_points(:, i)
                dist_current = distance(temp_point, targett)
                if ( dist_current <= radius**2) then
                    ! Append the index to the list
                    call add_to_list(idx, dist, node%leaf_indices(i), dist_current, count)
                end if
            end do

        ! Then check which childs to visit
        else
            axis = node%axis
            ! 1D distance from target to the splitting plane
            d1d = targett(axis) - node%point(axis)
            ! Recursively search the subtree that contains the target
            if (d1d < 0) then
                call ball_search_recursive(node%left, targett, dist, idx, radius, count)
                ! Check if we need to search the other subtree
                look_opposite = .false.
                if (abs(d1d) <= radius) look_opposite = .true.
#if periodic == 1
                if (targett(axis) - radius <= -L(axis) / 2. ) look_opposite = .true.
#endif
                if (look_opposite .eqv. .true.) then
                    call ball_search_recursive(node%right, targett, dist, idx, radius, count)
                end if
            else
                call ball_search_recursive(node%right, targett, dist, idx, radius, count)
                look_opposite = .false.
                if (abs(d1d) <= radius) look_opposite = .true.
#if periodic == 1
                if (targett(axis) + radius >= L(axis) / 2. ) look_opposite = .true.
#endif
                if (look_opposite .eqv. .true.) then
                    call ball_search_recursive(node%left, targett, dist, idx, radius, count)
                end if
            end if

        endif
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end subroutine ball_search_recursive_hyperp
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Search for points within a given box (cuboid) aligned with cartesian axes
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    function box_search(node, box) result(query)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        implicit none
        !in
        real(kind=prec) :: box(2*ndim) ! Box limits (xmin, xmax, ymin, ymax, zmin, zmax,...)
        type(KDTreeNode), pointer, intent(in) :: node
        !local
        integer :: init_depth = 0
        integer :: count_idx ! Counters for the number of elements in idx and dist
        integer(kind=intkind), allocatable :: temp_idx(:)
        !out
        integer(kind=intkind), allocatable :: idx(:) !index of the points within the radius
        type(KDTreeResult) :: query
    
        allocate(idx(1000))
        idx = -1
        count_idx = 0
    
        call box_search_recursive(node, init_depth, idx, box, count_idx)
    
        if (.not. allocated(idx)) STOP 'idx array is not allocated!'
            
        if (count_idx == 0) then
            ! No points found
            temp_idx = idx(1:0)
            call move_alloc(temp_idx, idx)
        else
            ! Reallocation to the correct size
            temp_idx = idx(1:count_idx)
            call move_alloc(temp_idx, idx)
            
            !For the moment, we do not sort the points by distance to the target

        end if

        query%idx = idx

        ! Distances are not calculated in this function, hence dist is not allocated
        deallocate(idx)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end function box_search
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
    recursive subroutine box_search_recursive(node, depth, idx, box, count_idx)
    ! Since only 1D comparisons are needed, this routine is extremely fast
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
        implicit none
        !in
        real(kind=prec) :: box(2*ndim) ! Box limits (xmin, xmax, ymin, ymax, zmin, zmax,...)
        type(KDTreeNode), pointer, intent(in) :: node ! Starting node (usually the root)
        integer, intent(in) :: depth     
        !out
        integer(kind=intkind), allocatable, intent(inout) :: idx(:)  ! Index of the points within the box
        integer, intent(inout) :: count_idx ! Counters for the number of elements in idx 
        !local
        real(kind=prec) :: split_value(3)
        integer :: i
        logical :: in_box ! tells if a point is inside the query box
        integer :: axis
        ! Temporary point for contiguous memory access
        real(kind=prec):: temp_point(ndim)

        if (.not. associated(node)) then
            return
        end if

        ! First, check if it is a leaf node
        if (node%is_leaf .eqv. .true.) then
            
            ! Check all points in the leaf with brute force
            do i = 1, size(node%leaf_indices)
                temp_point = node%leaf_points(:, i)
                call is_in_box(temp_point, box, in_box)
                if ( in_box .eqv. .true. ) then
                    ! Append the index to the list
                    call int_add_to_list(idx, node%leaf_indices(i), count_idx)
                end if
            end do

        else

            split_value = node%point
            axis = node%axis

            ! Recursively search the subtrees intersecting the box
            ! Periodicity does not makes sense as the query box is assumed to be contained
            ! inside the (periodic) bounding box of all points
            if (box(2*axis-1) < split_value(axis)) then
                call box_search_recursive(node%left, depth + 1, idx, box, count_idx)
            endif

            if (box(2*axis) > split_value(axis)) then
                call box_search_recursive(node%right, depth + 1, idx, box, count_idx)
            endif

        endif !node%is_leaf

    contains

        !subroutine to check if a point is inside the query box
        subroutine is_in_box(point, box, in_box)
            implicit none
            real(kind=prec), intent(in) :: point(ndim) ! Point to check
            real(kind=prec), intent(in) :: box(2*ndim) ! Box limits (xmin, xmax, ymin, ymax, zmin, zmax,...)
            logical, intent(out) :: in_box ! Result: true if the point is inside the box
            !local
            integer :: j

            in_box = .true.
            do j=1,ndim
                if (point(j) < box(2*j-1) .or. point(j) > box(2*j)) then
                    in_box = .false.
                    return
                end if
            end do

        end subroutine is_in_box

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end subroutine box_search_recursive
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
    recursive subroutine ball_reduction_recursive(node, targett, dist, idx, count)
    !Function that returns all indices and distances of points within a given node
    !Used when ball_search checks that a node and all its children are fully inside the search ball
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
        implicit none
        !in
        type(KDTreeNode), pointer, intent(in) :: node ! Starting node (usually the root)
        real(kind=prec), intent(in) :: targett(ndim) ! Target point (k-D)
        !out
        integer(kind=intkind), allocatable, intent(inout) :: idx(:)  ! Index of the points within the radius
        real(kind=prec), allocatable, intent(inout) :: dist(:)
        integer, intent(inout) :: count
        !local
        real(kind=prec):: dist_current
        integer :: i
        ! Temporary point for contiguous memory access
        real(kind=prec):: temp_point(ndim)

        if (.not. associated(node)) then
            return
        end if

        if (node%is_leaf .eqv. .true.) then
            do i = 1, size(node%leaf_indices)
                temp_point = node%leaf_points(:, i)
                dist_current = distance(temp_point, targett)
                call add_to_list(idx, dist, node%leaf_indices(i), dist_current, count)
            end do

        else
            call ball_reduction_recursive(node%left, targett, dist, idx, count)
            call ball_reduction_recursive(node%right, targett, dist, idx, count)
        endif 

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    end subroutine ball_reduction_recursive
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !subroutines to append an element to an array, called by BALL_SEARCH and BOX_SEARCH
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subroutine int_add_to_list(idx, new_value, count)
        implicit none
        integer(kind=intkind), allocatable, intent(inout) :: idx(:)
        integer(kind=intkind), intent(in) :: new_value
        integer, intent(inout) :: count 
        integer(kind=intkind), allocatable :: temp(:)
        integer :: n

        if (count >= size(idx)) then
            ! Resize the array
            n = size(idx)
            allocate(temp(10 * n))
            temp(1:n) = idx
            call move_alloc(temp, idx)
        end if

        count = count + 1
        idx(count) = new_value 

    end subroutine int_add_to_list


    subroutine add_to_list(idx, dist, idx_new, dist_new, count)
        implicit none
        integer(kind=intkind), allocatable, intent(inout) :: idx(:)
        integer(kind=intkind), intent(in) :: idx_new
        real(kind=prec), allocatable, intent(inout) :: dist(:)
        real(kind=prec), intent(in) :: dist_new
        integer, intent(inout) :: count
        integer(kind=intkind), allocatable :: temp_idx(:)
        real(kind=prec), allocatable :: temp_dist(:)
        integer :: n

        if (count >= size(idx)) then
            ! Resize the array
            n = size(idx)
            allocate(temp_idx(10 * n))
            allocate(temp_dist(10 * n))
            temp_idx(1:n) = idx
            temp_dist(1:n) = dist
            call move_alloc(temp_idx, idx)
            call move_alloc(temp_dist, dist)
        end if

        count = count + 1
        idx(count) = idx_new
        dist(count) = dist_new
        
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
    ! Minimum and maximum distance between box and point
    ! Uses distance_1D
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subroutine bbox_distance_1D(minbounds, maxbounds, point, mindist, maxdist, dim)
        implicit none
        real(kind=prec), intent(in) :: minbounds(ndim), maxbounds(ndim)
        real(kind=prec), intent(in) :: point(ndim)
        real(kind=prec), intent(inout) :: mindist(ndim), maxdist(ndim)
        integer, intent(in) :: dim
        real(kind=prec):: dmin, dmax

        ! Minimum distance
        if (point(dim) < minbounds(dim)) then
            dmin = distance_1D(point(dim), minbounds(dim), dim)
        else if (point(dim) > maxbounds(dim)) then
            dmin = distance_1D(point(dim), maxbounds(dim), dim)
        else
            dmin = 0.
        end if

        ! Maximum distance
        dmax = max( distance_1D(point(dim), minbounds(dim), dim), &
                    distance_1D(point(dim), maxbounds(dim), dim) )

        mindist(dim) = dmin
        maxdist(dim) = dmax
        
    end subroutine bbox_distance_1D
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Minimum k-d distance between box and point
    ! Uses distance_1D
    ! Used in knn search for the priority queue and prunnning
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subroutine bbox_mindist_kd(minbounds, maxbounds, point, mindistkd)
        implicit none
        real(kind=prec), intent(in) :: minbounds(ndim), maxbounds(ndim)
        real(kind=prec), intent(in) :: point(ndim)
        real(kind=prec), intent(inout) :: mindistkd
        integer :: dim
        real(kind=prec):: dmin

        mindistkd = 0.
        do dim = 1, ndim
            ! Minimum distance along this dimension
            if (point(dim) < minbounds(dim)) then
                dmin = distance_1D(point(dim), minbounds(dim), dim)
            else if (point(dim) > maxbounds(dim)) then
                dmin = distance_1D(point(dim), maxbounds(dim), dim)
            else
                dmin = 0.
            end if
            mindistkd = mindistkd + dmin**2
        end do

    end subroutine bbox_mindist_kd
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !Quicksort algorithm for sorting FINAL RESULTS
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subroutine quicksort(dist, idx, n)
        implicit none
        !in/out
        real(kind=prec), intent(inout) :: dist(n)    ! Distance array to be sorted
        integer(kind=intkind), intent(inout) :: idx(n) ! Corresponding indices
        integer, intent(in) :: n          ! Number of elements to sort
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

        ! Choose the pivot (median of low, mid, high)
        mid = (low + high) / 2

        ! Find median of low, mid, high
        if (dist(low) > dist(mid)) call swap2(dist, idx, low, mid)
        if (dist(low) > dist(high)) call swap2(dist, idx, low, high)
        if (dist(mid) > dist(high)) call swap2(dist, idx, mid, high)
        ! Now dist(low) <= dist(mid) <= dist(high)
        pivot_value = dist(mid)
        ! Move pivot to the end
        call swap2(dist, idx, mid, high)

        i = low - 1

        ! Partition the array
        do j = low, high - 1
            if (dist(j) <= pivot_value) then
                i = i + 1
                call swap2(dist, idx, i, j)
            end if
        end do

        ! Place the pivot in its correct position
        i = i + 1
        call swap2(dist, idx, i, high)

        ! Return the pivot index
        pivot_index = i
    end subroutine partition2

    subroutine swap2(dist, idx, i, j)
        implicit none
        real(kind=prec), intent(inout) :: dist(:)
        integer(kind=intkind), intent(inout) :: idx(:)
        integer, intent(in) :: i, j
        real(kind=prec):: temp_dist
        integer(kind=intkind) :: temp_idx

        ! Swap dist
        temp_dist = dist(i)
        dist(i) = dist(j)
        dist(j) = temp_dist

        ! Swap idx
        temp_idx = idx(i)
        idx(i) = idx(j)
        idx(j) = temp_idx
    end subroutine swap2
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!#######################################################
end module cosmokdtree
!#######################################################