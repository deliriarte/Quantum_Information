module utilities
    use precision, pr=>dp

    contains 
    subroutine random_hermitian(mat, nn)
        !initializate a complex hermitian matrix of size n x n 
        implicit none
        integer(pr), intent(in)      :: nn
        complex (pr), intent(inout) :: mat(1:nn,1:nn)
        integer(pr):: ii, jj
        real(pr) :: xx, yy, r, theta
        
        if (nn <= 0) then
            write(*,*) 'invalid dimension'
        end if
        
        do ii  = 1, nn
            do jj= 1, nn
            call random_number(xx)
            call random_number(yy)
            
            r = sqrt(2*(-log(xx)))
            
            theta = 2*(2*asin(1._pr)) * yy
            
            xx = r*cos(theta)
            yy = r*sin(theta)
            
            mat(ii, jj ) = cmplx (xx, yy)
            mat(jj, ii ) = conjg (mat (ii,jj))
            
            if (ii == jj ) then 
                mat(ii, jj )= cmplx(xx, 0)
            end if
            
            end do
        end do
        
    end subroutine

    subroutine diagonalize_hermitian(mat, nn)
        implicit none
        !initializate a diagonalize matrix of size n x n 
        integer(pr), intent(in)      :: nn
        complex (pr), intent(inout) :: mat(1:nn,1:nn)
        integer(pr):: ii, jj
        real(pr) :: xx, yy, r, theta
        
        
        do ii = 1, nn
            do jj = 1, nn
                call random_number(xx)
                call random_number(yy)
                
                r = sqrt(2*(-log(xx)))
                ! π = 2*arcsin(1)
                theta = 2*(2*asin(1._pr)) * yy
                
                xx = r*cos(theta)
                yy = r*sin(theta)
                
                if (ii == jj ) then 
                    mat(ii, jj )= cmplx(xx, 0._pr)
                else 
                    mat(ii,jj) = complex(0._pr,0._pr)
                end if
            end do
        end do
        
    end subroutine

    
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
    
    
    function ratio(lambda, siz) result(r)
        !Computes the ratio 
        implicit none
        real(pr), dimension(:),intent(in) :: lambda
        integer(pr), intent(in) :: siz
        integer(pr) :: ii
        real(pr), dimension(:), allocatable :: r
        
        allocate(r(siz-1))
        
        do ii=2,siz
            r(ii-1)=min(lambda(ii), lambda(ii-1))/max(lambda(ii),lambda(ii-1))
        end do
    
    end function
    
    subroutine calculate_spacing(spac, i_space, eigenvalues, i_eigen, num_eigen, resuls)
        !Computes the local spacing 
        implicit none
        integer(pr), intent(in) :: spac(:)
        integer(pr), intent(in) :: i_space, i_eigen, num_eigen
        real(pr), intent (in) :: eigenvalues(:)
        real(pr) :: summ, delta, t_size
        real (pr), intent(out) ::  resuls
        integer(pr) :: temp_s, min_index, max_index
        
        temp_s = spac(i_space)
        
        min_index = MAX(1, i_eigen - (temp_s/2))
        max_index = min ( i_eigen + (temp_s/2), num_eigen-1)
        
        delta = eigenvalues(i_eigen+1) - eigenvalues(i_eigen)
        summ = SUM(eigenvalues((min_index +1):(max_index +1)) - eigenvalues(min_index:max_index)  )

        t_size = size(eigenvalues((min_index +1) :(max_index +1)) - eigenvalues((min_index) :(max_index))  )
        resuls = delta/(summ/t_size)

    end subroutine

    
end module
