program eigenproblem
    use precision, pr=>dp
    use debugging
    use utilities
    implicit none
    complex (pr), allocatable :: matrix(:,:), work(:) 
    real (pr), allocatable    :: eig(:), rwork(:), space(:), s_i(:, :), rat(:)
    integer(pr)               :: n, i, j, lwork, info, nrep, rep,  k
    real(pr)                  :: xx, yy, avg_space, resuls
    integer , allocatable :: seed(:)
    integer(4), dimension(:), allocatable ::counts
    integer(8), dimension(7) ::spacings
    logical :: diag
    CHARACTER(pr) :: file_name
    
    !matrix dimension
    n = 2500
    
    !number of iterations
    nrep = 150
     
    !Set diagonal matrix ON
    diag = .True.
    
    !Load files where results of eigenvalue going to be saved
    if (diag.eqv..True.) then 
        open(12, file = 'eigen_space_diag.txt', status = 'old', action = 'write')
    else 
        open(12, file = 'eigenvalues_spacing.txt', status = 'old', action = 'write')
    end if
    
    !Load files where results of ratio going to be saved
    if (diag.eqv..True.) then 
        open(13, file = 'ratio_diag.txt', status = 'old', action = 'write')
    else 
        open(13, file = 'ratio.txt', status = 'old', action = 'write')
    end if
    
    !local spacing - set to be fixed
    spacings = [n/1500, n/1000, n/800, n/500, n/100, n/50, n/10]
    

    !Files to write "local spacing" results
    k = 40
    if (diag.eqv..True.) then 
        do i=1,SIZE(spacings)
           
            write(file_name,"(I4.4)") spacings(i)
            open(k, file="Results/local_space_diag"//TRIM("_"//file_name)//".txt", status='new', action='write')
            k = k +1
            
        end do
   else 
        do i=1,SIZE(spacings)
            write(file_name,"(I4.4)") spacings(i)
            open(k, file="local_spacing"//TRIM("_"//file_name)//".txt", status='old', action='write')
            k = k +1
        end do
    end if
    

    !Set random seed
   ! allocate(seed(1234))
  !  call random_seed(put = seed)
    
    allocate(matrix(1:n, 1:n ), eig(1:n), space(1:n-1), s_i(size(spacings), n-1))
    
    rep = 1
    do while (rep <= nrep)
    
        !Initialize matrix
        if (diag.eqv..True.) then
            call diagonalize_hermitian(matrix, n)
        else
            call random_hermitian(matrix, n)
        end if
        
        !crescent order
        info = 0 
!         100 format (*('('sf6.2xspf6.2x'i)':x))
!         do i = 1, n
!              write(*,100) (matrix(i,j), j = 1, n)
!         end do
        !Compute eigenvalues
        call eigenvalues(matrix, n, work, eig, rwork,  lwork, info)
        
!         do i=1,n
!              write(*,*)eig(i)
!          enddo
        
        !Compute spacing
        do i = 2 , n
            space(i-1) =  eig(i) - eig(i-1)
        end do
        
        !Compute ratio
        rat = ratio(space, n-1)        
        
        !Saving ratio results
        do i=1,n-1
            write(13,*)rat(i)
        enddo
        
        !Normalize spacing
        avg_space = sum(space)/(n-1)
        space =  space/avg_space
        
        !Saving normalize spacing results
        do i=1,n-1
            write(12,*)space(i)
        enddo
        
        !Calculate local spacing
        do i = 1, size(spacings)
            do j = 1, n-1
            call calculate_spacing(spacings, i, eig, j, n, resuls)
            s_i(i,j) = resuls
            end do
        end do
        
        !Saving local spacing results
        k = 40 
        do i=1,size(spacings)
            do j=1,n-1
                write(k,"(F13.7)")s_i(i,j)
            enddo
            k = k+1
        end do
        
        rep = rep +1

    end do

end program
