module utilities

    use precision, pr=>dp

    contains 
    
    subroutine random_pure_state(d,pure_vector, n_particles)
        !Generates a pure **RANDOM** quantum states **NON separable** of d**n_particles
        implicit none
        integer(pr), intent(in)                  :: d, n_particles !Dimension of the pure vector, number of particles
        complex (pr), allocatable, intent(inout) :: pure_vector(:) !Pure state vector
        integer(pr)                              :: i, ii, jj, kk !auxiliary variable
        real(pr)                                 :: norm !Normalization
        real(pr), allocatable                    :: rand1(:),  rand2(:) !Random vectors 
        double precision                         :: tempp, summ
        
        
        !Check if any of the dimensions are negative
        if (d <= 0 .or. n_particles <= 0) then
            write(*,*) 'invalid dimension or number of particles'
        end if
    
        allocate(pure_vector(d**n_particles), rand1(d**n_particles), rand2(d**n_particles))
    
        if (allocated(pure_vector)) then
             print *, 'Array is allocated'
         endif
        
        !Generating random vectors 
        call random_number(rand1)
        call random_number(rand2)
        
        !Computing norm of the vector 
        norm = sqrt(sum(abs(rand1**2+rand2**2)))
        
        !Generating pure random vector
        pure_vector = rand1/norm + rand2*complex(0,1/norm)
        
        !Checking normalization
        summ = 0._pr
        
        do kk = 1, d**n_particles
            tempp = abs(pure_vector(kk))**2
            summ = summ + tempp
        enddo
        
        summ = sqrt(summ)
        
        write(*,*) 'Normalization of pure vector=', summ

    end subroutine
    
    
    subroutine random_pure_separable_state(d,pure_vector, n_particles)
        !Generates a pure **RANDOM SEPARABLE** quantum states of size d*n_particles
        implicit none
        integer(pr), intent(in)                  :: d, n_particles !Dimension of the pure vector, number of particles
        complex (pr), allocatable, intent(inout) :: pure_vector(:) !Pure state vector
        double precision                         :: tempp, summ 
        integer(pr)                              :: ii, jj, kk
        double precision                         ::rand1, rand2 !rand vectors
        
        !Check if any dimensions are negative
        if (d <= 0 .or. n_particles <= 0) then
            write(*,*) 'invalid dimension or number of particles'
        end if
        
        allocate(pure_vector(d*n_particles))
        
        !Check if array is allocated
        if (allocated(pure_vector)) then
            print *, 'Array is allocated'
        endif
        
        !Generate separable states
        do ii = 1, n_particles  
            do jj = 1, d
                call random_number(rand1)
                if (jj ==1) then
                    pure_vector(jj+(ii-1)*d) = cmplx(rand1, 0._pr)
                else
                    call random_number(rand2)
                    pure_vector(jj+(ii-1)*d) = cmplx(rand1, rand2)
                end if
            enddo
    
            !Normalization
            summ = 0._pr            
            do kk = 1, d
                tempp= abs(pure_vector(kk+(ii-1)*d))**2
                summ = summ + tempp
            enddo
            summ = sqrt(summ)
            
            !Normalizing Vector
            do kk = 1, d
                pure_vector(kk+(ii-1)*d) = pure_vector(kk+(ii-1)*d)/(summ) 
            enddo
        
            !Check normalization
            summ=0.0
            do kk= 1, d*n_particles
                tempp= abs(pure_vector(kk+(ii-1)*d))**2
                summ=summ + tempp
            end do
            summ=sqrt(summ)
            
            write(*, *) 'Normalization of separable state=', summ
            
        enddo


    end subroutine
    
    subroutine print_complex_matrix(matriz)
        !Prints the matrix in a aesthetic way in the bash
        implicit none
        complex(pr) :: matriz(:, :)
        integer(pr) :: ii, jj
        
        100 format (*('('sf6.2xspf6.2x'i)':x))
        
        do ii = 1, size(matriz, 1)
            write(*,100) (matriz(ii, jj), jj = 1, size(matriz, 2))
        end do
        
    end subroutine

    SUBROUTINE writing_file(txt_name, mat)
        !Save the matrix and the main properties in a file 
        use precision, pr=>dp
        implicit none
        character(:), allocatable :: txt_name
        complex(pr),  intent(in)  :: mat(:,:)
        integer(pr)               ::ii, jj
        complex(pr)               :: mat_trace
        integer(pr)               :: nrow, ncol
        
        open(22, file = txt_name, status = 'unknown', action = 'write')
        
        write(22, 120)
        120 format ('#' 25x, 'Matrix values ',/)
        
        100 format (*('('sf6.2xspf6.2x'i)':x))
        do ii = 1, size(mat, 1)
             write(22,100) (mat(ii, jj), jj = 1, size(mat, 2))
        end do
        nrow = size(mat, 1)
        ncol = size(mat, 2)
        
        mat_trace = trace(mat, nrow,ncol )
        write(22, *) 'Number of rows:', size(mat, 1)
        write(22, *) 'Number of cols:',size(mat, 2)       
        write(22, *) 'Trace:', mat_trace
        
        close(22)
    end SUBROUTINE
    
    
    function trace(matt, nrow, ncol)result(diagonal_sum)
        !computes the trace of the matrix
        implicit none
        complex(pr), intent(in) :: matt(:,:) !matrix
        integer(pr)             :: kk ! auxiliary variable
        integer(pr), intent(in) :: nrow, ncol !number of rows, number of colums of the matrix
        complex(pr)             :: diagonal_sum !trace
        
        diagonal_sum = (0._pr, 0._pr)
        
        !the trace is defined for square matrices only
        if (abs(nrow - ncol) < 0.5) then
            do kk = 1, nrow 
                diagonal_sum = diagonal_sum + matt(kk,kk)
            end do
        end if
    
    end function
    
    function partial_trace(rho, nnum_part, hham_dim,particle)result(rho_reduce)
        !compute the partial trace of a generic matrix and provide the reduce matrix of the rho.
        implicit none
        complex(pr), intent(in)  :: rho(:,:) !density matrix 
        integer(pr), intent(in)  :: nnum_part, hham_dim, particle !number of particle, dimension, particle is a variable that relates the left or right particle (a or b)
        complex(pr), allocatable :: rho_reduce(:,:) ! reduce matrix 
        integer(pr)              :: sz,sz_reduce, max_right_ind, max_left_ind !sizes and indeces
        integer(pr)              ::ii, jj, ll, kk, xx, yy, xxy, yyx, ind !auxiliary variable
        complex(pr)              :: temp !temporaly matrix
        
        !Dimensionality
        sz = hham_dim**nnum_part
        sz_reduce = hham_dim**(nnum_part-1)
        
        allocate(rho_reduce(sz_reduce, sz_reduce))
    
        !indices to loop
        max_right_ind = hham_dim**(particle-1)
        max_left_ind = hham_dim**(nnum_part-particle)
        
        do ii=0, max_right_ind-1
 
            do jj=0, max_left_ind-1
            
                do kk=0, max_right_ind-1
                    
                    do ll =0, max_left_ind-1
                        
                        xxy = jj*hham_dim**(particle-1)+ii+1
                        
                        yyx = ll*hham_dim**(particle-1)+kk+1
                        temp = complex(0.0,0.0)
                        
                        do ind=0,hham_dim-1
                            xx =(jj*hham_dim+ind)*hham_dim**(particle-1)+ii+1
                            yy =(ll*hham_dim+ind)*hham_dim**(particle-1)+kk+1
                            temp = temp + rho(xx,yy)
                        end do
                        
                        rho_reduce(xxy,yyx)=temp
                    
                    end do
                    
                end do
                
            end do
            
        end do
        
    end function


    subroutine check_hermitian(mat, mat_hermitian)
        !check if a matrix is hermitian by observing if the congugate is the same or not. 
        implicit none
        complex(pr), intent(in) :: mat(:,:)
        logical, intent(inout)  :: mat_hermitian
        integer(pr)             ::ii, jj
        
        do ii=1,size(mat, 1)
            do jj=ii,size(mat, 2)
                if (mat(ii,jj).ne.conjg(mat(jj,ii))) THEN
                    mat_hermitian = .FALSE.
                    write(*,*)'Matrix is not hermitian -- something is wrong'
                end if
            enddo
        enddo
        
    end subroutine
    
    
    subroutine transform_separable(psi,psi_large, dd, nn)
        !Converts a separable state d*n and convert a vector with dd**nn
        complex(pr), intent(in)                 :: psi(:)
        complex(pr), allocatable, intent(inout) :: psi_large(:)
        integer(pr), intent(in)                 :: dd, nn
        integer(pr)                             :: ii, jj, ll
        integer(pr)                             :: cc, kk
        integer(pr), allocatable                :: temp_vect(:)
        double precision                        :: summ, tempp
        double complex                          :: temp

        !call init_state(psi_large, psi%nsub, psi%ndim, debug, sepp, purr)
        call random_pure_state(dd,psi_large, nn)
        allocate(temp_vect(nn))
            
        
        do ii= 1, dd**nn
            kk=ii-1
            do jj= 1, nn
                cc=nn-jj
                temp_vect(cc+1)=int((kk-mod(kk,dd**cc))/dd**cc +0.1) +1		
                kk=kk-(temp_vect(cc+1)-1)*(dd**cc)
            end do

            temp=(1.,0.)
            do jj= 1, nn
                temp= temp*psi(temp_vect(jj)+(jj-1)*dd )
            end do
            psi_large(ii)=temp
        end do
    end subroutine 

    
    function tensor_product(mat1, mat2)result(tensor_result)
    
    implicit none
        complex(pr):: mat1(:,:), mat2(:,:)
        
        complex(pr):: tensor_result(SIZE(mat1,1)*SIZE(mat2,1),SIZE(mat1,2)*SIZE(mat2,2))
        
        integer(pr) :: i, j , sz1_1, sz1_2, sz2_1, sz2_2
        
        sz1_1 = size(mat1, 1)
        sz1_2 = size(mat1, 2)
        
        
        sz2_1 = size(mat2, 1)
        sz2_2 = size(mat2, 2)
        
        
        
        do i = 0, sz1_1 -1
            do j = 0, sz1_2 -1
            
                tensor_result(i*sz2_1+1:(i+1)*sz2_1, j*sz2_2+1:(j+1)*sz2_2) = mat1(i+1,j+1)*mat2(:,:)
            
            enddo
        
        enddo
    
    end function
    
    
    subroutine eigenvalues(matr, n, work, eig, rwork,  lwork, info)
        !Computes the eigenvalues of a matrix n x n using the LAPACK function zheev
        implicit none
        complex (8), allocatable  :: work(:) 
        real (pr), allocatable    ::  rwork(:)
        integer(pr)              :: lwork
        integer(pr), intent(in)              :: info
        integer(pr), intent(in )             :: n
        complex (8), intent(inout)  :: matr(1:n,1:n)
        real (pr), allocatable, intent(inout)   :: eig(:)
        
        ! optimal lwork
        lwork=-1
        allocate(work(1))
        allocate(rwork(max(1, 3*n-2)))
        
        if(.not.allocated(eig))allocate(eig(n))
        
        call zheev('N','U',n,matr,n,eig,work,lwork,rwork,info)
        
        lwork = int(real(work(1)))
        
        deallocate(work)
        deallocate(rwork)


        ! actual diag
        allocate(work(max(1,lwork)))
        allocate(rwork(max(1, 3*n-2)))
        
        call zheev('N','U',n,matr,n,eig,work,lwork,rwork,info)
        
        deallocate(work)
        deallocate(rwork)
        
    end subroutine
    
end module
