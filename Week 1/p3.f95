program test_performance
    implicit none
    integer,parameter      :: dp = kind (1.d0)
    integer(4)             :: n, m, l, i
    real (dp), allocatable :: mat1(:,:), mat2(:,:), mat_prod(:,:), mat_result(:,:), mat_result_2(:,:)
    real (dp)              :: finish, start, time_matmul, time_matxmat_1, time_matxmat_2
    
    open(11, file = 'results_p3_Ofast.dat', status = 'old', action = 'write')
    write(11, 120)
    120 format ('#' 25x, 'Results of problem 3: Time performance multiplication of nxm and mxl',/, &
                '#', 6x,'n', 6x,'m', 6x,'l',15x,'matmul',18x,  'loop (order 1)',18x, 'loop (order 2)' )
    
    n = 100
    l = 100

    do m = 100,1000, 100
        
        allocate(mat1(1:n,1:m), mat2(1:m, 1:l), mat_prod(1:n, 1:l), mat_result(1:n, 1:l), mat_result_2(1:n, 1:l))
        
        call random_number(mat1)
        call random_number(mat2)

        
            
        call cpu_time(start)
        mat_prod = matmul(mat1, mat2)
        call cpu_time(finish)
        time_matmul = finish-start
        
        call cpu_time(start)
        call matxmat_1(mat1, mat2, n, m, l, mat_result)
        call cpu_time(finish)
        time_matxmat_1 = finish-start
     
        call cpu_time(start)    
        call matxmat_2(mat1, mat2, n, m, l, mat_result_2)
        call cpu_time(finish)
        time_matxmat_2 = finish-start
!       
        write(11,*) n, m, l, time_matmul, time_matxmat_1, time_matxmat_2 
!                    100 format ( 4x, I12, 3x, I12, 3x, I12, 6x, E24.16, 2x,E24.16, 8x, E24.16)
         
        
        deallocate(mat1, mat2, mat_result, mat_result_2, mat_prod)
    
    end do

close(11)    
end program

subroutine matxmat_1(mat_1, mat_2, nn, ll, mm, mat_p)
    implicit none
    integer,parameter           :: dp = kind (1.d0)
    integer(4) :: i, j, k
    integer(4), intent(in) ::  nn, ll, mm
    real (dp), intent(in)       :: mat_1(1:nn, 1:ll), mat_2(1:ll, 1:mm)
    real (dp), intent(out) :: mat_p(1:nn,1:mm)  
    
    mat_p(i,j) = 0.0
    do i = 1, nn
        do j=1,mm
            do k=1,ll
                mat_p(i,j)=mat_p(i,j)+mat_1(i,k)*mat_2(k,j)
            end do
        end do
    end do

    
end subroutine
 


subroutine matxmat_2(mat_1, mat_2, nn, ll, mm, mat_p)
    implicit none
    integer,parameter           :: dp = kind (1.d0)
    integer(4) :: i, j, k
    integer(4), intent(in) ::  nn, ll, mm
    real (dp), intent(in)       :: mat_1(1:nn, 1:ll), mat_2(1:ll, 1:mm)
    real (dp), intent(out) :: mat_p(1:nn,1:mm)  
    
    mat_p(i,j) = 0.0
    do i = 1, mm
        do j=1,nn
            do k=1,ll
                mat_p(i,j)=mat_p(i,j)+mat_1(i,k)*mat_2(k,j)
            end do
        end do
    end do

    
end subroutine

