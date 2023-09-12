subroutine mcNormSpatialHet(SH, array, n_permutes)
    implicit none
    integer :: y1, y2, x1, x2
    integer subset_y, subset_x
    integer :: i, j, k, count
    real(kind=8), intent(out) :: SH ! note that kind=8 for double precision is specific to gfortran!!
    integer, intent(in) :: n_permutes
    real(kind=8), dimension(:, :), intent(in) :: array
    real(kind=8), dimension(:, :), allocatable :: hstack_array
    real(kind=8), dimension(:, :), allocatable :: v_hstack_array
    integer, dimension(:, :), allocatable :: permutation_array
    real(kind=8) :: array_mean, subarray_mean
    real(kind=8) :: G, N
    integer, dimension(2) :: dim
    integer :: total_permutes
    real(kind=8) :: num, denom
    integer :: d1, d2
    G = 0
    
    dim = shape(array)
    d1 = dim(1)
    d2 = dim(2)

    ! via normalizedSpatialHet (nsh.f90), first term is # iterations over loop 1,
    ! second term is # iterations over loop 2 
    total_permutes = d1*d2*(d1-1)*(d2-1) + (d2*d2-d2)

    if (n_permutes > total_permutes) then
        print *, "Number of permutations is greater than total possible permutations. Exiting."
        return
    end if

    allocate(hstack_array(d1, 2*d2))
    allocate(v_hstack_array(2*d1, 2*d2))
    allocate(permutation_array(4, n_permutes))
    
    array_mean = sum(array)/max(1, size(array))
    ! reshape does horizontal stacking, so vertical stack requires transpose, horiz stack, transpose
    hstack_array = reshape([array, array], [d1, 2*d2])
    v_hstack_array = transpose(reshape([transpose(hstack_array), transpose(hstack_array)], [2*d1, 2*d2]))
    N = array_mean*(1.5*d1*d2*(d1-1)*(d2-1) + d1*(d1-1) + d2*(d2-1))/((d1*(d1-1)+1)*(d2*(d2-1)+1))

    do i = 1, n_permutes
        call rectangleBounds(x1, x2, 1, d1)
        call rectangleBounds(y1, y2, 1, d2)

        permutation_array(1, i) = x1 
        permutation_array(2, i) = x2 
        permutation_array(3, i) = y1 
        permutation_array(4, i) = y2 

        ! decided to skip including code for running without replacement since 
        ! its very computational expensive and not all that beneficial

        subset_x = x2  - x1
        subset_y = y2  - y1
        subarray_mean = 0.0
        count = 0
        do k = y1, y2-1
            do j = x1, x2-1
                subarray_mean = subarray_mean + v_hstack_array(k, j)
                count = count + 1
            end do
        end do
        subarray_mean = subarray_mean / count
        G = G + abs(subarray_mean - array_mean)
    end do

    num = total_permutes*G 
    denom = real(n_permutes)*(N*(d1*(d1-1)+1)*(d2*(d2-1)+1))
    SH = num / denom
    deallocate(hstack_array)
    deallocate(v_hstack_array)

end subroutine mcNormSpatialHet

subroutine rectangleBounds(l_bound, r_bound, min, max)
    implicit none 
    integer :: l_bound, r_bound, min, max, temp
    call randint(l_bound, min, max+1)
    call randint(r_bound, l_bound+1, l_bound+max)
    if (l_bound == r_bound) then
            r_bound = max
    end if
    if (r_bound < l_bound) then ! swap the values of l and r bounds
        temp = r_bound 
        r_bound = l_bound 
        l_bound = temp
    end if

end subroutine rectangleBounds

subroutine randint(j, min, max)
    implicit none
    real :: r
    integer, intent(out) :: j 
    integer, intent(in) :: min, max
    call random_number(r) ! float value between 0 and 1
    j = min + FLOOR((max+1-min)*r) ! scale value by min max range and round to int
end subroutine randint
    