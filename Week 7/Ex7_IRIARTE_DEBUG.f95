MODULE debugging 
    !Module precision holds the precision of the data that will be involve. This will easily convert to double precision to single precision in all the call by just tokenizing the correct one.
    use precision, dp=>dp
    implicit none
    CONTAINS 
 
    SUBROUTINE check_dimension(debug, nrow, ncol)
     !Documentation: Check if the dimension of the matrix inserted are correct. If one of the rows or columns are negative then it will call in the bash to insert new values. 
     !Arguments: 
     ! - nrow: number of rows. 
     ! - ncol : number of columns
     ! - debug : logical value. If .true. then it will check for bugs. 
     ! - msg : character for the bash message.
     !Conditions:
     ! PRE: debug is a logical value. nrow and ncol integers
     ! POST: nrow and ncol integers.
    IMPLICIT NONE
    integer(dp), intent(inout) :: nrow, ncol
    logical, intent(in) :: debug
    character(:), allocatable :: msg

    msg = 'Checking if the dimensions are positive:'
    call check_debug(debug, msg)
    
    if (debug.eqv..true.) then 
        
        do while (nrow <= 0 .or. ncol <= 0) 
            write(*,*) 'Dimension of the matrix not valid. Insert again:'
            write(*,*) '------------------------------------------------'
            write(*,*) "Number of rows desire:"
            read(*,*) nrow   
            
            write(*,*) "Number of cols desire:"
            read(*,*) ncol
        end do
        write(*,*) 'All perfect'
    end if 
        
    END SUBROUTINE 
    
    subroutine product_versality( debug, nn, kk, ll, mm)
    !Documentation: Check if the dimension of the matrix inserted are correct in order to make a multiplication. 
    !Arguments: 
    !takes 2 matrices (number of rows x ncols) of nn x kk and ll x mm. 
    ! - debug : logical value. If .true. then it will check for bugs. 
     ! - msg : character for the bash message.
     !Conditions:
     ! PRE: debug is a logical value. nrows and ncols integers for the two matrices.
     ! POST:  nrows and ncols integers for the two matrices.
    implicit none
    integer(dp), intent(inout) ::  nn, ll, mm, kk
    integer(dp) ::inner
    logical, intent(in) :: debug
    character(:), allocatable :: msg
    
    
    msg = 'Checking if the matrices dimension for product purposes are ok:'
    call check_debug(debug, msg)
    
    if (debug.eqv..true.) then 
    
        if (ll .ne. kk) then
            write(*,*) 'Do you really think that you can do a matrix product with this dimensions?'
            
            write(*,*) '------------------------------------------------'
            write(*,*) "Inner dimension of the matrices:"
            read(*,*) inner
            
            ll = inner
    
        end if

    write(*,*) 'All perfect'
    
    end if 
    
    end subroutine
    
    subroutine check_square_matrix(debug, mat)
    !Documentation: Check if the matrix is square or rectangular
    !Arguments: 
    ! - debug : logical value. If .true. then it will check for bugs.
    ! - mat : real matriz value.
    ! - msg : character for the bash message.
     !Conditions:
     ! PRE: debug is a logical value and matrix.
     ! POST:  a msg which say if its square or rectangular
    implicit none
    real (dp), intent(in) :: mat(:,:)
    logical, intent(in) :: debug
    character(:), allocatable :: msg
    
    msg = 'Square or rectangular matrix'
    call check_debug(debug, msg)
    
    if (size(mat,dim=1)==size(mat,dim=2)) then
        write(*,*) "Matrix is square"
    else 
        write(*,*) "Matrix is rectangular"
    end if
    
    end subroutine
 
    subroutine check_debug(debug, msg, input)
    !Documentation: prints a debug message and optionally prints the type of the input variable.
    !Arguments: 
    ! - debug : logical value. If .true. then it will check for bugs.
    ! - input: input variable which if given we will give the kind.
    ! - msg : character for the bash message.
     !Conditions:
     ! PRE: debug is a logical value and matrix. msg is optional and is a character. input is optional as well and can be any class.
     ! POST:  messages with bugs.
        implicit none
        logical, intent(in) :: debug 
        character(:), allocatable, optional :: msg
        class(*), intent (in), optional :: input
        
        if (debug.eqv..true.) then 
            write(*,*) '-----------'
            write(*,*) msg
            
            if (present(input)) then
                write(*,*) 'type of the varible:'
                
                select type(input)
                            
                    type is (integer(2))
                    write(*,*) "Got an integer of kind 2"
                    
                    type is (integer(4))
                    write(*,*) "Got an integer of kind 4"
                    
                    type is (integer(8))
                    write(*,*) "Got an integer of kind 8"
                    
                    type is (real(4))
                    write(*,*) "Got a real of kind 4"
                    
                    type is (real(8))
                    write(*,*) "Got a real of kind 8"
                    
                    type is (logical)
                    write(*,*) "Got a logic"
                    
                    type is (complex)
                    write(*,*) "Got a complex"
                end select
        
            end if
        end if 
    
    end SUBROUTINE
    
    subroutine matrix_parameter(mtx)
    !Documentation: Prints the most important attributes of a matrix
    !Arguments: 
    ! - matrix : real 2 dimensional array.
     !Conditions:
     ! PRE: real matrix with 2 dimension.
     ! POST:  messages with attributes such as: rows, columns, kind and elements.
        implicit none
        real(dp), intent(in) :: mtx(:,:)
        write(*,*) 'Matrix parameters'
        write(*,*) '-----------------'
        write(*,*) 'Rows of the matrix:'
        write(*,*) size(mtx,1)
        write(*,*) 'columns of the matrix'
        write(*,*) size(mtx,2) 
        write(*,*)"Kind of matrix", kind(mtx)
        write(*,*) 'Matrix elements'
        call print_matrix(mtx)
    end subroutine 
    
    subroutine print_matrix(mat)
    !Documentation: Prints the elements of a matrix
    !Arguments: 
    ! - matrix : real 2 dimensional array.
     !Conditions:
     ! PRE: real matrix with 2 dimension.
     ! POST:  Elements of the matrix in the bash.
        implicit none
        real(dp), intent(in) :: mat(:,:)
        integer(dp) :: i, j, m, n
        m = size(mat,1)
        n = size(mat,2) 
        do, i=1,m
            write(*,*) ( mat(i,j), j=1,n )
        enddo
    end subroutine
end MODULE  
        
        
        
        
        
        
        
        
