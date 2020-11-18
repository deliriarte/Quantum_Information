program quantum_oscillator_1d
    use precision, pr=>dp
    use utilities
 
    implicit none
    real(pr) :: dx, L, fact_const
    integer(pr) :: n_points, i, j, ii, jj
    real(pr), allocatable :: v(:), x(:), eigen(:), energy(:), eigen_theo(:, :)
    complex(pr), allocatable :: hamiltonian(:,:)
    complex (8), allocatable  :: work(:) 
    real (pr), allocatable    ::  rwork(:)
    integer(pr)              :: lwork, info, n
    real (pr), parameter :: mass = 1._pr, omega = 1._pr, hbar = 1._pr
   
   OPEN(12, file="evect_exp.txt", status='old', action='WRITE')
   open(13, file = 'eigenvalues.txt', status = 'old', action = 'write')
   open(14, file = 'energy.txt', status = 'old', action = 'write')
   OPEN(17, file="evect_exp_imag.txt", status='old', action='WRITE')

   
    !The idea is that potential and wavefunctions are discretized and the second derivative in the kinetcis energy is approximated as a dinite difference. this the hamiltonian becomes a matrix whose eigenvalues can be found using LAPACK.
    
    !mesh
    open(16, file = "mesh.txt", status = 'old', action= 'read')
    read(16, *) n_points
    close(16)
    
    
    allocate(v(n_points), x(n_points), hamiltonian(n_points,n_points), eigen(n_points), energy(n_points))
    
    !Lets started by writting the potential
    ! size of simulation box
    open(17, file = "L.txt", status = 'old', action= 'read')
    read(17, *) L
    close(17)

    
    dx = 2._pr *L/dble(n_points-1)
    
    !potential
    do i = 1, n_points
        x(i) =-L + ( float(i) - 1) * dx
        v(i) = 0.5_pr * mass *  omega *omega* x(i) * x(i)
    end do
    
    
    !kinetic matrix +potential matrix
    hamiltonian = 0
    do i = 1, n_points 
        do j = 1, n_points
            if (i == j) then
                hamiltonian(i,j) = 2._pr *hbar* hbar/(2._pr*mass *dx**2) + v(i)
            elseif (i == j+1) then
                hamiltonian(i, j) = -1._pr*hbar* hbar/(2._pr*mass *dx**2)
            elseif(i == j-1) then
                hamiltonian(i, j) = -1._pr*hbar* hbar/(2._pr*mass *dx**2)
            end if
        end do
    end do
       
    !compute eigenvalues
    call eigenvalues(hamiltonian, n_points, work, eigen, rwork,  lwork, info)
    
    if(info.NE.0)then
        print*,'DIAGONALIZATION FAILED'
        stop
    end if
    
    !theoretical eigenvalues
    do i= 1, n_points
        energy(i) = (float(i-1) + 0.5_pr) *hbar *omega
    enddo
    
    !saving data
    do i=1,n_points
        write(13,*) eigen(i)
    enddo
    
    !saving data
    do i=1,n_points
        write(14,*) energy(i)
    enddo
    
    ! Writing approximated and "NORMALIZED" (scaled so that they are normalized) eigenvectors
    !	(which are now in the Hamiltonian matrix because lapack's subroutine overwrite the matrix)
    ! (writing them in columns, so are easier to plot using gnuplot)

!     
!    100 format (*('('sf6.2xspf6.2x'i)':x))
!     do ii = 1, size(hamiltonian, 1)
!         write(*,100) (hamiltonian(ii, jj), jj = 1, size(hamiltonian, 1))
!     end do

    do j = 1, n_points
        do i = 1, n_points
            write(12,*) real(hamiltonian(i,j))
            write(17, *) aimag(hamiltonian(i,j))
        enddo
    enddo
    
    CLOSE(12)
    CLOSE(13)
    CLOSE(14)
    
end program
