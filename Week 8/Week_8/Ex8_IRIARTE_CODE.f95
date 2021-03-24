program Exercise_8
    !This programs compute the density matrix of pure (separable and not separable) states and the partial trace of the generic matrix case.
    
    use precision, pr=>dp
    use utilities
    
    implicit none
    complex(pr), allocatable :: random_states(:), density_matrix(:,:), density_matrix_square(:,:), reduce_density_matrix(:,:)
    complex(pr), allocatable :: separable_states(:), density_matrix_separable(:,:), density_matrix_square_separable(:,:)
    complex(pr),allocatable  :: spin(:,:), red_spin(:,:), separable_states_dn(:), separable_states_dnn(:)
    integer(pr)              :: num_part, ham_dim, i, j,k , reduce_dimension, xx, yy
    complex(pr)              :: tra
    integer(pr)              :: left_particle
    logical                  :: size_memory, hermitian
    character(:),allocatable :: file_1, file_2, file_3, file_4, file_5, file_6, file_7, file_8
    integer(pr)              :: ind_s, ind_t
   
   !allocate size memory available   
   !Flag for compute this or not - The bad thing is that it breaks the whole program
   size_memory = .False.
   if (size_memory.EQV..True.) then
        num_part  = 30
        ham_dim = 2
        do i = 2, ham_dim
            do j = 1, num_part
                call random_pure_state(i,random_states,j)
                deallocate(random_states)
                write(*,*) 'dimension:', i, 'number of partcle:', j
            enddo
        enddo
       !memory only allows up to 25 number of particules with d = 2. 
    end if
    
   !We will only focus on 2 particles in d = 2
   num_part  = 2_pr
   ham_dim = 2_pr
   
   !------------------------------------------- general/non separable - --------------------------------------
   !generate the random states
   call random_pure_state(ham_dim,random_states,num_part)
   
    !print states in the bash
    100 format (*('('sf6.2xspf6.2x'i)':x))
    write(*,*) 'Printing states pure (generic) ='
    do i = 1, ham_dim**num_part
        write(*,100) (random_states(i))
    enddo
   
   
   !generate density matrix 
   allocate(density_matrix(1:ham_dim ** num_part, 1:ham_dim ** num_part))
   allocate(density_matrix_square(1:ham_dim ** num_part, 1:ham_dim ** num_part))
   
   do i = 1, ham_dim ** num_part
       do j = 1, ham_dim ** num_part
            density_matrix(i, j) =  random_states(i) * conjg(random_states(j))
       enddo
   enddo


   
   !Checking if the matrix is hermitian
   hermitian = .True.
   call check_hermitian(density_matrix, hermitian)
   
   !write in bash
   write(*,*) 'Generic density matrix='
   call print_complex_matrix(density_matrix)
   
   !saving results in file
   file_1= "generic_rho.txt"
   call writing_file(file_1, density_matrix)
   
   !rho²
   density_matrix_square = matmul (density_matrix, density_matrix)
   
   !print matrix in bash
   write(*,*) 'Generic density sqaure matrix='
   call print_complex_matrix(density_matrix_square)
   
   !compute the trace of rho²
   tra = trace(density_matrix_square, ham_dim ** num_part, ham_dim ** num_part)
   
   write(*,*) 'Trace of rho square= ', tra 
   
   !saving results in file
   file_2 = "generic_rho_square.txt"
   call writing_file(file_2, density_matrix)
   
   ! ---------------------------------------- separable case ------------------------------------------------------------
   
   call random_pure_separable_state(ham_dim,separable_states,num_part)
   
    !print states in the bash 
    write(*,*) 'Printing separable states = '
    do i = 1, ham_dim*num_part
        write(*,100) (separable_states(i))
    enddo
   
   
   !Density matrix for n = 2 
   allocate(density_matrix_separable(ham_dim** num_part, ham_dim ** num_part))
   allocate(density_matrix_square_separable(ham_dim ** num_part, ham_dim ** num_part))
 
   
   !We need to map this to a d**n state, this is done by
   	call transform_separable(separable_states,separable_states_dn, ham_dim, num_part)
   
   do i = 1, ham_dim ** num_part
       do j = 1, ham_dim ** num_part
            density_matrix_separable(i, j) =  conjg(separable_states_dn(i)) * separable_states_dn(j)
       enddo
   enddo
   
   !Checking if the matrix is hermitian
   hermitian = .True.
   call check_hermitian(density_matrix_separable, hermitian)
   
   !Printing in bash the separable density matrix
   write(*,*) 'separable density matrix='
   call print_complex_matrix(density_matrix_separable)
   
   !saving results in file
   file_3 = "separable_rho.txt"
   call writing_file(file_3, density_matrix_separable)
   
   !Computing rho²
   density_matrix_square_separable = matmul (density_matrix_separable, density_matrix_separable)
   
   write(*,*) 'separable density square matrix='
   call print_complex_matrix(density_matrix_square_separable)
   
   !saving results in file
   file_4 = "separable_rho_square.txt"
   call writing_file(file_4, density_matrix_square_separable)
   
   !trace of rho sqaure
   tra = trace(density_matrix_square_separable, ham_dim * num_part, ham_dim * num_part)
   write(*,*) 'Trace of rho square of separable states= ', tra 
   
   
   ! -------------------------------- COMPUTATION OF REDUCE MATRIX OF GENERIC CASE ------------------------------------
   
    reduce_dimension = ham_dim **(num_part -1)
    allocate(reduce_density_matrix(1:reduce_dimension,1:reduce_dimension))
   
    !--------------------------------------- left particle ------------------------------
   
    !Compute the reducre matrix for the left particle
    left_particle = 1._pr
    reduce_density_matrix = partial_trace(density_matrix, num_part, ham_dim, left_particle)
    
    !Printing results for left particles
    write(*,*) "reduce generic matrix"
    call print_complex_matrix(reduce_density_matrix)
    
   !Saving file for left particle
   file_5 = "reduce_matrix_left.txt"
   call writing_file(file_5, reduce_density_matrix)
   
    !--------------------------------------- right particle ------------------------------
    !Compute the reducre matrix for the left particle
    left_particle = 2._pr
    reduce_density_matrix = partial_trace(density_matrix, num_part, ham_dim, left_particle)
    
    !Printing results for left particles
    write(*,*) "reduce generic matrix"
    call print_complex_matrix(reduce_density_matrix)
    
    !Saving file for left particle
    file_6 = "reduce_matrix_right.txt"
    call writing_file(file_6, reduce_density_matrix)
   
   
    ! Generating the most common density matrix  in the all physics:
    
    allocate(spin(1:4, 1:4), red_spin(1:2, 1:2))
    spin = complex(0._pr, 0._pr)
    spin(1, 1) = real(0.5_pr)
    spin(4, 4) = real(0.5_pr)
    spin(1, 4) = real(-0.5_pr)
    spin(4, 1) = real(-0.5_pr)
    
    call print_complex_matrix(spin)
    
    !Saving file for left particle
    file_7 = "spins.txt"
    call writing_file(file_7, spin)
    
    !Computing reduce matrix for left particle (is the same for right and left)
    red_spin = partial_trace(spin, 2_pr, 2_pr, 1_pr)
    call print_complex_matrix(red_spin)
    
    !Saving file for left particle
    file_8 = "spins_reduce.txt"
    call writing_file(file_8, red_spin)
   
   end program

   
