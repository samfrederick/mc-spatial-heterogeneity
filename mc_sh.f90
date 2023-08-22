subroutine mcSpatialHet(SH, array, n_permutes, replacement)
    implicit none
    integer :: y1, y2, x1, x2
    integer subset_y, subset_x
    integer :: i, j, k, l, count
    real(kind=8), intent(out) :: SH ! note that kind=8 for double precision is specific to gfortran!!
    integer, intent(in) :: n_permutes
    real(kind=8), dimension(:, :), intent(in) :: array
    logical, intent(in) :: replacement
    real(kind=8), dimension(:, :), allocatable :: hstack_array
    real(kind=8), dimension(:, :), allocatable :: v_hstack_array
    integer, dimension(:, :), allocatable :: permutation_array
    real(kind=8) :: array_mean, subarray_mean
    real(kind=8) :: G
    integer, dimension(2) :: dim
    integer :: n, m
    integer :: total_permutes
    real(kind=8) :: num, denom
    logical :: loop, duplicates
    G = 0

    dim = shape(array)
    n = dim(1)
    m = dim(2)
    total_permutes = n*m*(n+1)*(m+1)

    if (n_permutes > total_permutes) then
        print *, "Number of permutations is greater than total possible permutations. Exiting."
        return
    end if

    allocate(hstack_array(n, 2*m))
    allocate(v_hstack_array(2*n, 2*m))
    allocate(permutation_array(4, n_permutes))
    
    array_mean = sum(array)/max(1, size(array))
    ! reshape does horizontal stacking, so vertical stack requires transpose, horiz stack, transpose
    hstack_array = reshape([array, array], [n, 2*m])
    v_hstack_array = transpose(reshape([transpose(hstack_array), transpose(hstack_array)], [2*n, 2*m]))

    do i = 1, n_permutes
        call rectangleBounds(x1, x2, 1, n)
        call rectangleBounds(y1, y2, 1, m)

        permutation_array(1, i) = x1 
        permutation_array(2, i) = x2 
        permutation_array(3, i) = y1 
        permutation_array(4, i) = y2 

        if (.not. replacement) then
            if (i > 1) then
                loop = .true.
                duplicates = .false.
                do while (loop)
                    ! check the permutations for duplicates, replace any if found
                    duplicates = .false.
                    do l=1, i-1
                        if(all(permutation_array(:, i) == permutation_array(:, l))) then
                            ! replace the values of permutation array at entry i
                            call rectangleBounds(x1, x2, 1, n)
                            call rectangleBounds(y1, y2, 1, m)

                            permutation_array(1, i) = x1 
                            permutation_array(2, i) = x2 
                            permutation_array(3, i) = y1 
                            permutation_array(4, i) = y2 
                            duplicates = .true.
                            exit
                        end if
                    end do

                    ! if duplicates were found, check the replacements for duplicates, if none
                    ! found then can exit while loop
                    if (duplicates) then
                        do l=1, i-1
                            if(all(permutation_array(:, i) == permutation_array(:, l))) then
                                exit
                            end if
                        end do
                    else
                        loop = .false.
                    end if 
                end do
            end if
        end if

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
    denom = real(n_permutes)*n*n*m*m
    SH = num / denom
    deallocate(hstack_array)
    deallocate(v_hstack_array)

end subroutine mcSpatialHet

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
    