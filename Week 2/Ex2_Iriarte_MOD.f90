MODULE matmodule
    use precision, pr=>dp
    implicit none
    type matrix
        double complex, allocatable :: mat_elements(: , :)
        integer(pr) :: n_row, n_col
        double complex :: mat_trace
        double complex :: mat_det
    end type matrix
        
    interface operator (.Adj.)
    module procedure adjoin
     end interface

     interface operator (.Tr.)
     module procedure diagonal_sum
     end interface
CONTAINS 
    
    SUBROUTINE initialization (mat, nrow, ncol)
        use precision, pr=>dp
        implicit none
        integer(pr) :: ii, jj
        integer(pr), intent(in) :: nrow, ncol
        type (matrix), intent(inout) :: mat
        mat%n_row = nrow
        mat%n_col = ncol
        
        allocate (mat%mat_elements(nrow, ncol))
        
        do ii = 1, nrow
            do jj = 1, ncol
            mat%mat_elements(ii, jj) = CMPLX(1d0*ii,1_pr *jj) 
            end do
        end do
       
       100 format (*('('sf6.2xspf6.2x'i)':x))
       do ii = 1, nrow
            write(*,100) (mat%mat_elements(ii, jj), jj = 1, ncol)
        end do
        mat%mat_trace = (0._pr, 0._pr)
        mat%mat_det = (0._pr, 0._pr)
        
    END SUBROUTINE 
    
    function diagonal_sum(matt)
        use precision, pr=>dp
        implicit none
        type(matrix), intent(in) :: matt
        integer(pr) :: kk
        double complex diagonal_sum
        diagonal_sum = (0._pr, 0._pr)
        
        !the trace is defined for square matrices only
        if (abs(matt%n_row - matt%n_col) < 0.5) then
            do kk = 1, matt%n_row 
                diagonal_sum = diagonal_sum + matt%mat_elements(kk,kk)
            end do
        end if
        
    end function
    


    function adjoin(matrix_a) 
        use precision, pr=>dp
        implicit none
        type(matrix), intent(in) :: matrix_a
        type(matrix) :: adjoin
        integer (pr) :: xx, yy
        
        adjoin%n_col = matrix_a%n_row
        adjoin%n_row = matrix_a%n_col
        
        allocate(adjoin%mat_elements(adjoin%n_row,  adjoin%n_col))
        
        do xx = 1, matrix_a%n_row
            do yy = 1, matrix_a%n_col
                adjoin%mat_elements(xx,yy) = dconjg(matrix_a%mat_elements(yy,xx))
            end do
        end do
        
       100 format (*('('sf6.2xspf6.2x'i)':x))
       do xx = 1, adjoin%n_row
            write(*,100) (adjoin%mat_elements(xx, yy), yy = 1, adjoin%n_col)
        end do
       
    end function
        
    SUBROUTINE writing_file(txt_name, mat)
        
        use precision, pr=>dp
        implicit none
        character(:), allocatable :: txt_name
        type (matrix), intent(in) :: mat
        integer(pr) ::ii, jj
        
        
        open(11, file = txt_name, status = 'unknown', action = 'write')
        
        write(11, 120)
        120 format ('#' 25x, 'Matrix values ',/)
        
        100 format (*('('sf6.2xspf6.2x'i)':x))
        do ii = 1, mat%n_row
             write(11,100) (mat%mat_elements(ii, jj), jj = 1, mat%n_col)
        end do
        
        write(11, *) 'Number of rows:', mat%n_row
        write(11, *) 'Number of cols:', mat%n_col        
        write(11, *) 'Trace:', mat%mat_trace
        
    
    end SUBROUTINE
    
        
end MODULE

        
        
