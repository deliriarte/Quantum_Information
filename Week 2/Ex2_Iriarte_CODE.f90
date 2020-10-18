program main
    use matmodule
    IMPLICIT NONE
    
    type(matrix) :: matriz 
    type(matrix) :: matriz_1 
    integer(pr) :: row, col, i, j 
    double complex :: trace
    character(:), allocatable :: file_1, file_2
    
    file_1= "matrix.txt"
    file_2= "matrix_adjacente.txt"
    
    write(*,*) "number of rows desire:"
    read(*,*) row
    
    write(*,*) "number of cols desire:"
    read(*,*) col
    
    call initialization (matriz, row, col)
    
   
    matriz%mat_trace =  .Tr.(matriz)
      
    102 format  (*('('sf6.2xspf6.2x'i)':x))
    write(*,102)  matriz%mat_trace
    
    
    matriz_1 = .Adj.(matriz) 
    matriz_1%mat_trace= .Tr.(matriz_1)
    
    write(*,102)  matriz_1%mat_trace
    
    call writing_file(file_1, matriz)
    
    call writing_file(file_2, matriz_1)
    
END PROGRAM 
