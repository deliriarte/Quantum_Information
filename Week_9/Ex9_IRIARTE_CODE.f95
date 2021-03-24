program Exercise_8
    !This programs computes the eigenvalues of the Ising model of size d^n x d^n using the lapack modules
    
    use precision, pr=>dp
    use utilities
    
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
    integer(pr)                             :: n_part, ham_dim !number of particles, dimension of the hamiltonian
    integer(pr), allocatable                :: hamiltonian_diagonal(:) !diagonal term of the hamiltonian
    complex(pr), allocatable                :: temp(:,:) !support matrix
    integer(pr)                             :: ii, jj, tmpvar !support variables
    real(pr)                                :: lambda !interaction term
    complex (8), allocatable                :: work(:) !to call zheev
    real (pr), allocatable                  :: rwork(:) !to call zheev
    integer(pr)                             :: lwork, info !to call zheev
    real(pr), allocatable                   :: eigen(:) !Save the eigenvalues 
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
    allocate(hamiltonian_magentic(ham_dim**n_part,ham_dim**n_part))
    allocate(hamiltonian_J(ham_dim**n_part,ham_dim**n_part))
    
    !support variables 
    allocate(temp(ham_dim**n_part,ham_dim**n_part))
    allocate(hamiltonian_diagonal(ham_dim**n_part))
    allocate(eigen(ham_dim**n_part))
   
   
    ! Building diagonal term 
    hamiltonian_diagonal = 0
    do ii=1,n_part
        tmpvar = ham_dim**(n_part+1-ii)
        do jj=0,ham_dim**n_part-1,tmpvar
            hamiltonian_diagonal((jj+1):(jj+tmpvar/2)) = hamiltonian_diagonal((jj+1):(jj+tmpvar/2))+1
            hamiltonian_diagonal((jj+tmpvar/2+1):(jj+tmpvar)) = hamiltonian_diagonal((jj+tmpvar/2+1):(jj+tmpvar))-1
        enddo
    enddo

 
   !Printing diagonal term
    do ii = 1, size(hamiltonian_diagonal)
        write(*,*) (hamiltonian_diagonal(ii))
    end do
    hamiltonian =0
   hamiltonian_magentic = 0
   do ii = 1, ham_dim**n_part
        hamiltonian_magentic(ii,ii)=hamiltonian_diagonal(ii)
   enddo
   
    
    call check_hermitian(hamiltonian_magentic,debug)

   !initializing matrix for interaction term 

   hamiltonian_J = 0
    do ii=1,n_part-1
        temp = 0
        do jj=1,n_part
        
            if (jj==1) then
            
                if (((n_part+1-jj).eq.ii).or.((n_part+1-jj).eq.(ii+1))) then
                
                    temp(1:ham_dim**(jj),1:ham_dim**(jj)) = sigma_x(:,:) !sigma x
                
                else
                
                    temp(1:ham_dim**(jj),1:ham_dim**(jj)) = I(:,:) !identity
                
                endif
                
            else

                if (((n_part+1-jj).eq.ii).or.((n_part+1-jj).eq.(ii+1))) then
                    !perform tensor product
                    temp(1:ham_dim**(jj),1:ham_dim**(jj)) = tensor_product(sigma_x, temp( 1:ham_dim**(jj-1),1:ham_dim**(jj-1) ))
                
                else 
                    temp(1:ham_dim**(jj),1:ham_dim**(jj)) = tensor_product(I, temp( 1:ham_dim**(jj-1),1:ham_dim**(jj-1) ))
                    
                endif
        
            endif
        
        enddo
        
        hamiltonian_J(:,:) = hamiltonian_J(:,:) + temp(:,:)
        
    enddo
    
    call check_hermitian(hamiltonian_J,debug)

       
    do ii = 0, 20
        
        lambda = 3._pr/(20._pr)*ii
        
        
        hamiltonian = 0
        hamiltonian(:,:) = hamiltonian_J(:,:) + lambda *hamiltonian_magentic(:,:)
        
        call check_hermitian(hamiltonian,debug)
        
        info =0
        !compute eigenvalues calling subroutine
        call eigenvalues(hamiltonian, ham_dim**n_part, work, eigen, rwork,  lwork, info)
        
        if(info.NE.0)then
            print*,'DIAGONALIZATION FAILED'
            stop
        end if
        
        ! Storing the first four eigenvalues for every lambda
        eig_val_results(ii+1,1) = lambda
        eig_val_results(ii+1,2:5) = eigen(1:4)

    enddo


    ! Saving eigenvalues to file
    OPEN(12, file="results"//trim("_"//NN)//".txt", status='REPLACE', action='WRITE')

    DO ii=1,21
        WRITE(12,"(F12.8)",advance="No") eig_val_results(ii,1)

        DO jj=2,5-1
            WRITE(12,"(F15.8)",advance="No") eig_val_results(ii,jj)
        END DO
        
        WRITE(12,"(F15.8)",advance="Yes") eig_val_results(ii,5)
    END DO

    CLOSE(12)

   
    deallocate(hamiltonian,hamiltonian_magentic,hamiltonian_J)
    
    deallocate(temp,hamiltonian_diagonal,eigen)
   
   
end program
