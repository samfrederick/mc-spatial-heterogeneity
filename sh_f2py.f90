subroutine spatialHet(SH, array)
    implicit none
    integer :: y1, y2, x1, x2
    integer subset_y, subset_x
    integer :: k, j, count
    real(kind=8), intent(out) :: SH ! note that kind=8 for double precision is specific to gfortran!!
    real(kind=8), dimension(:, :), intent(in) :: array
    real(kind=8), dimension(:, :), allocatable :: hstack_array
    real(kind=8), dimension(:, :), allocatable :: v_hstack_array
    real(kind=8) :: array_mean, subarray_mean
    real(kind=8) :: G
    integer, dimension(2) :: dim
    integer :: n, m
    G = 0

    dim = shape(array)
    n = dim(1)
    m = dim(2)

    allocate(hstack_array(n, 2*m))
    allocate(v_hstack_array(2*n, 2*m))
    
    array_mean = sum(array)/max(1, size(array))
    ! reshape does horizontal stacking, so vertical stack requires transpose, horiz stack, transpose
    hstack_array = reshape([array, array], [n, 2*m])
    v_hstack_array = transpose(reshape([transpose(hstack_array), transpose(hstack_array)], [2*n, 2*m]))

    do y1 = 1, n+1
        do y2 = y1+1, y1+n
            do x1 = 1, m+1
                do x2 = x1+1, x1+m
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
                    !print *, 'mean = ', subarray_mean, 'G = ', G
                    !print *, G, x1,x2,y1,y2
                end do
            end do
        end do
    end do
    
    SH = G / (n*n*m*m)
    !G = 0
    deallocate(hstack_array)
    deallocate(v_hstack_array)
end subroutine spatialHet
