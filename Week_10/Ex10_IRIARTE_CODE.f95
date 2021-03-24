program Exercise_8

!This programs computes the eigenvalues of the Ising model of size d^n x d^n using the lapack modules
    
    use precision, pr=>dp
    use utilities
    use ex10
    
    implicit none
    real(pr)                                :: eig_val_results(21, 5) !array with results
    !Pauli matrices
    complex(pr), dimension(2, 2), parameter :: sigma_z = reshape((/&
                                                cmplx(1._pr, 0._pr), cmplx(0._pr, 0._pr), &
                                                cmplx(0._pr, 0._pr), cmplx(-1._pr, 0._pr)/), shape(sigma_z))
    complex(pr), dimension(2, 2), parameter :: I = reshape((/&
                                                cmplx(1._pr, 0._pr), cmplx(0._pr, 0._pr), &
                                                cmplx(0._pr, 0._pr), cmplx(1._pr, 0._pr)/), shape(I))
    complex(pr), dimension(2, 2), parameter :: sigma_y = reshape((/&
                                                cmplx(0._pr, 0._pr), cmplx(0._pr, 1._pr), &
                                                cmplx(0._pr, -1._pr), cmplx(0._pr, 0._pr)/), shape(sigma_y))
    complex(pr), dimension(2, 2), parameter :: sigma_x = reshape((/&
                                                cmplx(0._pr, 0._pr), cmplx(1._pr, 0._pr), &
                                                cmplx(1._pr, 0._pr), cmplx(0._pr, 0._pr)/), shape(sigma_x))
    
    complex(pr), allocatable                :: hamiltonian (:,:), hamiltonian_magentic(:,:), hamiltonian_J(:,:) 
    complex(pr), allocatable                :: h_l(:,:), h_r(:,:), I_n(:,:)
    integer(pr)                             :: n_part, ham_dim, step !number of particles, dimension of the hamiltonian
    integer(pr), allocatable                :: hamiltonian_diagonal(:) !diagonal term of the hamiltonian
    integer(pr)                             :: ii, jj, tmpvar !support variables
    real(pr)                                :: lambda !interaction term
    complex (8), allocatable                :: work(:) !to call zheev
    real (pr), allocatable                  :: rwork(:) !to call zheev
    integer(pr)                             :: lwork, info !to call zheev
    real(pr), allocatable                   :: eigenvalue(:)
    complex(pr), allocatable                :: hamiltonian_2N (:,:), h_2n(:,:) !Save the eigenvalues 
    character(4)                            :: NN !to save files with differents names
    logical                                 :: debug !logical value to debug
    
    
    debug = .True.
    !reading variables 
    open(16, file = "n_particles.txt", status = 'old', action= 'read')
    read(16, *) n_part
    close(16)
    !in particular the maximum obtained correspond to N = 11
    
    open(16, file = "n_particles.txt", status = 'old', action= 'read')
    read(16, *) NN
    close(16)
    

    open(16, file = "hamiltonian_dimension.txt", status = 'old', action= 'read')
    read(16, *) ham_dim
    close(16)
    
    !allocating variables 
    allocate(hamiltonian(ham_dim**n_part,ham_dim**n_part))
    
    allocate(h_l(ham_dim**n_part,ham_dim**n_part))
    allocate(h_r(ham_dim**n_part,ham_dim**n_part))
    
    !support variables 
    allocate(hamiltonian_2N(ham_dim**(2*n_part),ham_dim**(2*n_part)))
    allocate(eigenvalue(ham_dim**(2*n_part)))
    allocate(h_2n(ham_dim**(2*n_part), ham_dim**(2*n_part)))
    
    
    !building  magnetic field hamiltonian part
    call ising_model_h(hamiltonian_magentic, ham_dim, n_part)
    call check_hermitian(hamiltonian_magentic,debug)
    
    
    !building interactive hamiltonian
    call ising_model_J(hamiltonian_J, ham_dim, n_part)
    call check_hermitian(hamiltonian_J,debug)
        
    allocate(I_n(ham_dim**n_part,ham_dim**n_part))

    !building identity matrix of dimension n x n
    I_n = 0
    do ii = 1, ham_dim**n_part
    I_n(ii, ii) = 1
    enddo
    
    
    hamiltonian =0
    
    do ii = 0, 10
        
        lambda = 3._pr/(10._pr)*ii
        
        !Building total hamiltonian
        hamiltonian = 0
        hamiltonian(:,:) = hamiltonian_J(:,:) + lambda *hamiltonian_magentic(:,:)
        
        if (ii == 0) then 
            call print_complex_matrix(hamiltonian)
            write(*,*) lambda
        endif
        call check_hermitian(hamiltonian,debug)
        
        
        call intializing_h_AB(h_l, h_r, ham_dim, n_part)
        
        if (ii == 0) then 
            call print_complex_matrix(h_l)
        endif
        
        call check_hermitian(h_l,debug)
        call check_hermitian(h_r,debug)
        
        !Real Space Renormalization Group
        
        do step = 1, 50
        
            hamiltonian_2N = tensor_product(hamiltonian, I_n) + tensor_product(I_n, hamiltonian) + tensor_product(h_l, h_r)
            !if step = 1 call print_complex_matrix(hamiltonian_2N)
            !call check_hermitian(hamiltonian_2N,debug)

            h_2n = hamiltonian_2N

            
            info =0
            
            !compute eigenvalues calling subroutine
            
            call eigenvalues(h_2n, ham_dim**(2*n_part), work, eigenvalue, rwork,  lwork, info)
            
            if(info.NE.0)then
                print*,'DIAGONALIZATION FAILED'
                stop
            end if
            
            hamiltonian = matmul(transpose(conjg(h_2n(:,1:ham_dim**n_part))),matmul(hamiltonian_2N, h_2n(:,1:ham_dim**n_part))) 
            
!             if (step == 1) then 
!                 call print_complex_matrix(hamiltonian)
!             endif
            
            ! Updating hamilt_l and hamilt_r (first N eigenvectors)
			h_l = matmul( transpose(conjg(h_2n(:,1:ham_dim**n_part))), matmul(tensor_product(h_l,I_n),h_2n(:,1:ham_dim**n_part)) )
			h_r = matmul( transpose(conjg(h_2n(:,1:ham_dim**n_part))), matmul(tensor_product(I_n,h_r),h_2n(:,1:ham_dim**n_part)) )
            
         enddo
         
            ! Storing the first four eigenvalues for every lambda
            eig_val_results(ii+1,1) = lambda
            eig_val_results(ii+1, 2) = eigenvalue(1)/(n_part*2**50)
        
        write (*,*) '========================='
        write(*,*) "Ground state energy:", eig_val_results(ii+1, 2)
        write(*,*) "lambda= ",lambda
   enddo
  
      ! Saving eigenvalues to file
    OPEN(12, file="results"//trim("_"//NN)//".txt", status='REPLACE', action='WRITE')

    DO ii=1,11
        WRITE(12,"(F12.8)",advance="No") eig_val_results(ii,1)

        DO jj=2,5-1
            WRITE(12,"(F18.8)",advance="No") eig_val_results(ii,jj)/abs(eig_val_results(1,2))
        END DO
        
        WRITE(12,"(F18.8)",advance="Yes") eig_val_results(ii,5)/abs(eig_val_results(1,2))
    END DO

     CLOSE(12)
  

   
    deallocate(hamiltonian,hamiltonian_magentic,hamiltonian_J)
    
    deallocate(eigenvalue)
   
   
end program
