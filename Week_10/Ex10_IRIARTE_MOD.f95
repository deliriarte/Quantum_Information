module ex10

    use precision, pr=>dp
    use utilities
    contains 
    
    subroutine ising_model_h(h_m, d, n_particles)
        !Initialize the magnetic part of the Ising model
        implicit none
        integer(pr), intent(in)                  :: d, n_particles !Dimension of the pure vector, number of particles
        complex (pr), allocatable, intent(inout) :: h_m(:,:) !Pure state vector
        integer(pr), allocatable                 :: h_d(:)
        integer(pr)                              :: ii, jj, kk !auxiliary variable
            
            
        allocate(h_d(d**n_particles))
        allocate(h_m(d**n_particles,d**n_particles))
        
        ! Building diagonal term 
        h_d = 0
        do ii=1,n_particles
            kk = d**(n_particles+1-ii)
            do jj=0,d**n_particles-1,kk
                h_d((jj+1):(jj+kk/2)) = h_d((jj+1):(jj+kk/2))+1
                h_d((jj+kk/2+1):(jj+kk)) = h_d((jj+kk/2+1):(jj+kk))-1
            enddo
        enddo

    
    !Printing diagonal term
        do ii = 1, size(h_d)
            write(*,*) (h_d(ii))
        end do
        

        h_m = 0
        do ii = 1, d**n_particles
                h_m(ii,ii)=h_d(ii)
        enddo
        
            
        
        deallocate(h_d)

    end subroutine
    
    
    
        subroutine ising_model_J(h_j, d, n_particles)
        !Initialize the magnetic part of the Ising model
        implicit none
        integer(pr), intent(in)                  :: d, n_particles !Dimension of the pure vector, number of particles
        complex (pr), allocatable, intent(inout) :: h_j(:,:) !Pure state vector
        complex(pr), allocatable                :: temp(:,:) !support matrix
        integer(pr)                              :: ii, jj, kk !auxiliary variable
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
                
        !initializing matrix for interaction term 
   
        allocate(h_j(d**n_particles,d**n_particles))
        allocate(temp(d**n_particles,d**n_particles))

        
        h_j = 0
            do ii=1,n_particles-1
                temp = 0
                do jj=1,n_particles
                
                    if (jj==1) then
                    
                        if (((n_particles+1-jj).eq.ii).or.((n_particles+1-jj).eq.(ii+1))) then
                        
                            temp(1:d**(jj),1:d**(jj)) = sigma_x(:,:) !sigma x
                        
                        else
                        
                            temp(1:d**(jj),1:d**(jj)) = I(:,:) !identity
                        
                        endif
                        
                    else

                        if (((n_particles+1-jj).eq.ii).or.((n_particles+1-jj).eq.(ii+1))) then
                            !perform tensor product
                            temp(1:d**(jj),1:d**(jj)) = tensor_product(sigma_x, temp( 1:d**(jj-1),1:d**(jj-1) ))
                        
                        else 
                            temp(1:d**(jj),1:d**(jj)) = tensor_product(I, temp( 1:d**(jj-1),1:d**(jj-1) ))
                            
                        endif
                
                    endif
                
                enddo
                
                h_j(:,:) = h_j(:,:) + temp(:,:)
                
            enddo
                
            deallocate(temp)

    end subroutine
    
    subroutine intializing_h_AB(h_left, h_right, d, n_particles)
        !Initialize the magnetic part of the Ising model
        implicit none
        integer(pr), intent(in)                  :: d, n_particles !Dimension of the pure vector, number of particles
        complex (pr),              intent(inout) :: h_left(:,:), h_right(:,:) !Pure state vector
        integer(pr)                              :: ii, jj, kk !auxiliary variable
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
                
        !initializing matrix for interaction term 
   
        h_left = 0._pr
        h_right = 0._pr
        
        do jj = 0, (n_particles- 1)
        
            if (jj == 0) then
            
                h_left(1:d**(jj+1), 1:d**(jj+1)) = I(:,:)
                h_right(1:d**(jj+1), 1:d**(jj+1)) = sigma_x(:,:)
            
            else
                
                if (jj == (n_particles-1)) then
                    h_left(1:d**(jj+1), 1:d**(jj+1)) = tensor_product(h_left(1:d**(jj), 1:d**(jj)), sigma_x(:,:))
                    h_right(1:d**(jj+1), 1:d**(jj+1)) = tensor_product(h_right(1:d**(jj), 1:d**(jj)), I(:,:))
                else 
                    h_left(1:d**(jj+1), 1:d**(jj+1)) = tensor_product(h_left(1:d**(jj), 1:d**(jj)), I(:,:))
                    h_right(1:d**(jj+1), 1:d**(jj+1)) = tensor_product(h_right(1:d**(jj), 1:d**(jj)), I(:,:))
                endif
            endif
        enddo

    end subroutine
    
    
    
end module
