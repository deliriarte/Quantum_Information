program test_performance
    !Documentation: Compare  performance of  matrix multiplication by using the matmul function and two define matrices that perfoms the multiplication changing the rows and the colums. Results are saved in a file.
    
     !Arguments: 
     ! - Dimensions of the two matrices (n x k and lx m - rows and columns)
     ! - debug : logical value. If .true. then it will check for bugs. 
     ! - message : character for the bash message.
    use precision, dp=>dp
    use debugging
    implicit none
    integer(dp)             :: n, m, l, k, i, tmp
    character(:), allocatable :: message
    real (dp), allocatable :: mat1(:,:), mat2(:,:), mat_prod(:,:), mat_result(:,:), mat_result_2(:,:)
    real (dp)              :: finish, start, time_matmul, time_matxmat_1, time_matxmat_2
    logical :: debugg 
    debugg = .True.
    
  
    open(11, file = 'yata.dat', status = 'new', action = 'write')
    write(11, 120)
    120 format ('#' 25x, 'Results of problem 3: Time performance multiplication of nxm and mxl',/, &
                '#', 6x,'n', 6x,'m', 6x,'l',15x,'matmul',18x,  'loop (  ijk  )',18x, 'loop (  jik  )' )
 
    n = 100
    m = 100

    do tmp = 100,3000, 100
        write(*,*) tmp
        l = tmp
        !First check: negative dimension. 
        !write(*,*) 'Matrix 1 setup'
        !call check_dimension(debugg,n, l)
        !write(*,*) 'Matrix 2 setup'
        !call check_dimension(debugg, l, m)
        
        !Second check: inner dimension of the matrices.
        !call product_versality(debugg,  n, l, l, m)
        
        allocate(mat1(1:n,1:l), mat2(1:l, 1:m), mat_prod(1:n, 1:m), mat_result(1:n, 1:m), mat_result_2(1:n, 1:m))
            
        call random_number(mat1)
        call random_number(mat2)
        
        !third check: parameters of the matrices
        !call matrix_parameter(mat1)
        
        call cpu_time(start)
        mat_prod = matmul(mat1, mat2)
        call cpu_time(finish)
        time_matmul = finish-start
        
        call cpu_time(start)
        call matxmat_1(mat1, mat2, n, l,m, mat_result)
        call cpu_time(finish)
        time_matxmat_1 = finish-start
        
        call cpu_time(start)    
        call matxmat_2(mat1, mat2, n, l,m, mat_result_2)
        call cpu_time(finish)
        time_matxmat_2 = finish-start
        
        write(11,*) n, l,m, time_matmul, time_matxmat_1, time_matxmat_2 
        
        deallocate(mat1, mat2, mat_result, mat_result_2, mat_prod)
        !write(*,*) '-------------------------------------------------------------------------------------'
    end do
close(11)    
end program

subroutine matxmat_1(mat_1, mat_2, nn, ll, mm, mat_p)
     !Documentation: Matrix multiplication by looping first rows and then by columns.
     !Arguments: 
     ! - Dimensions of the two matrices (n x l and lx m - rows and columns) !The INNER are the SAME! and matrix elements. 
     ! - mat_p (nn x mm) : matrix whcih corresponds on the result of the product of the two previous matrices.
     !Conditions:
     ! PRE: two real matrices with dimension 2.
     ! POST: one real matrix (the correspondent result) with dimension 2
    use precision, dp=>dp
    implicit none
    integer(dp) :: i, j, k
    integer(dp), intent(in) ::  nn, ll, mm
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
     !Documentation: Matrix multiplication by looping first columns and then by rows.
     !Arguments: 
     ! - Dimensions of the two matrices (n x l and lx m - rows and columns) !The INNER are the SAME! and matrix elements. 
     ! - mat_p (nn x mm) : matrix whcih corresponds on the result of the product of the two previous matrices.
     !Conditions:
     ! PRE: two real matrices with dimension 2.
     ! POST: one real matrix (the correspondent result) with dimension 2
    use precision, dp=>dp
    implicit none
    integer(dp) :: i, j, k
    integer(dp), intent(in) ::  nn, ll, mm
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
