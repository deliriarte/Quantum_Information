program quantum_oscillator_1d_evolution
    !compile using gfortran -o ./p7 prec_mod.f90 Ex5_mod.f95 Ex7_IRIARTE_CODE.f95 -llapack  -lfftw3
    !Computes the time evolution of the quantum harmonic oscillator using the fftw library.
    use precision, pr=>dp
    use utilities
    use debugging
    use ISO_C_BINDING
 
    implicit none
    include 'fftw-3.3.8/api/fftw3.f03'
    include 'fftw-3.3.8/api/fftw3l.f03'
    
    real(pr)                 :: dx, L, fact_const, total_time,  delta_t, V_max
    real(pr)                 :: max_p, pi, mean, norm, sigma
    integer(pr)              :: time_split,n_points, i, j, ii, jj
    real(pr), allocatable    :: v(:), x(:), eigen(:), energy(:), eigen_theo(:, :), moments(:), grid_psi_square(:,:)
    complex(pr), allocatable :: hamiltonian(:,:) , psi_0(:), fft1(:), fft2(:), temp(:)
    complex(pr), allocatable :: temp2(:), grid_psi(:,:)
    complex (8), allocatable :: work(:) 
    real (pr), allocatable   ::  rwork(:)
    integer(pr)              :: lwork, info, n, dim_space, dim_time
    type(C_PTR)              :: plan, planb! to be used in "fft" subroutines
    real (pr), parameter     :: mass = 1._pr, omega = 1._pr, hbar = 1._pr
    logical                  :: debug
   
    OPEN(12, file="results.txt", status='old', action='write')
   
    OPEN(21, file="results_square.txt", status='old', action='write')
   
    OPEN(22, file="results_potential.txt", status='old', action='write')
    
    OPEN(23, file="results_mean.txt", status='old', action='write')
    
    debug = .False.
   
   !---------------------------- Computation of eigenvalues and eigenstates of --------------------------------------------
   !------------------------------------- time independent hamiltonian ----------------------------------------------------
    
    !The idea is that potential and wavefunctions are discretized and the second derivative in the kinetcis energy is approximated as a dinite difference. this the hamiltonian becomes a matrix whose eigenvalues can be found using LAPACK.
    
    !mesh
    open(16, file = "mesh.txt", status = 'old', action= 'read')
    read(16, *) n_points
    close(16)
    
    !n_points = 10
    
    
    !Lets started by writting the potential
    ! size of simulation box
    open(17, file = "L.txt", status = 'old', action= 'read')
    read(17, *) L
    close(17)

    dx = 2._pr *L/dble(n_points)
    
    dim_space = n_points + 1
    
    !Now we also need to add something more additional that is the time
    open(17, file = "total_time.txt", status = 'old', action= 'read')
    read(17, *) total_time
    close(17)
    
    !total_time = 1
    !Now we also need to add something more additional that is the time
    open(17, file = "time_split.txt", status = 'old', action= 'read')
    read(17, *) time_split
    close(17)
   ! time_split = 4
    dim_time = time_split+1
    !Let's compute the time intervarls spacing
    delta_t = total_time/(time_split)
    write(*,*) delta_t
    allocate(v(dim_space), x(dim_space), hamiltonian(dim_space,dim_space), grid_psi(dim_space,dim_time), moments(dim_space))
    
    allocate (grid_psi_square(dim_space, dim_time))
    allocate( eigen(dim_space), energy(dim_space), psi_0(dim_space))
    allocate(fft1(dim_space), fft2(dim_space), temp(dim_space), temp2(dim_space))
    
    !checks on dimension and if the matrix is square times square
    
    call check_dimension(debug, dim_space, dim_space)
    call check_square_matrix(debug, hamiltonian)
    !potential
    do i = 1, dim_space
        x(i) =-L + ( float(i) - 1) * dx
        v(i) = 0.5_pr * mass *  omega *omega* x(i) * x(i)
    end do
    
    !kinetic matrix +potential matrix
    hamiltonian = 0
    do i = 1, dim_space 
        do j = 1, dim_space
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
    call eigenvalues(hamiltonian, dim_space, work, eigen, rwork,  lwork, info)
    
    if(info.NE.0)then
        print*,'DIAGONALIZATION FAILED'
        stop
    end if
    
    !ground state normalized
    psi_0 = hamiltonian(:,1)/sqrt(dx)
    i=1
 
 !---------------------------- Computation of time evolution of ground state --------------------------------------------
     
    !now we have to propagate them 
    !Create grid with first state beign the ground state
    grid_psi(:,1) = psi_0
    
    ! Evaluating pi
	pi = 2._pr*asin(1._pr)
	
    ! Evaluating the max
	max_p = 2._pr*pi/dx
    
    ! Initializing vector of moments
    do i=1,(dim_space-1)/2
        moments(i) = max_p/(dim_space) * (float(i)-1)
    enddo

    do i=(dim_space+1)/2, dim_space
        moments(i) = max_p/(dim_space) * (float(i)-1) - max_p 
    enddo
    
    
    do i = 2, dim_time
        !exp(-i*V(x)*delta_t/2)
        temp(:) = exp(cmplx(0._pr, -0.5_pr*delta_t*0.5_pr*(x(:) - i * delta_t/total_time)**2)) * grid_psi(:,i-1)
        
        !Fourier from x to p 
        call dfftw_plan_dft_1d(plan,dim_space,temp,fft1,FFTW_FORWARD,FFTW_ESTIMATE)
        call dfftw_execute_dft(plan, temp, fft1)
        call dfftw_destroy_plan(plan)
        
        !exp(-i*kinetic*delta_t)
        temp2(:) = exp(cmplx(0._pr, -delta_t*(moments(:)**2)*0.5_pr)) * fft1(:)
        
        call dfftw_plan_dft_1d(planb, dim_space, temp2, fft2, FFTW_BACKWARD, FFTW_ESTIMATE)
        call dfftw_execute_dft(planb, temp2, fft2)
        call dfftw_destroy_plan(planb)   
        
        !issues with normalization
        fft2(:) = fft2(:) / dim_space
        
        !exp(-i*V(x)*delta_t/2)
         grid_psi(:,i) = exp(cmplx(0_pr, -0.5_pr*delta_t*0.5_pr*(x(:)-i*delta_t/total_time)**2)) * fft2(:)
        
    enddo
    
    !Square of the eigenfunctions
    do ii=1,dim_space
        grid_psi_square(ii,:) = abs(grid_psi(ii,:))**2
    enddo
    
    do i = 1, dim_space
        write(12, "(F14.8)",advance="No") x(i)
        do j = 1, dim_time -1 
            write(12, "(F14.8, F14.8)",advance="No") grid_psi(i,j)
        enddo
        write(12,"(F14.8, F14.8)",advance="Yes") grid_psi(i,dim_time)

    enddo
    
    do i = 1,dim_space
        write(21, "(F14.8)",advance="No") x(i)
        do j = 1, dim_time -1 
            write(21, "(F14.8)",advance="No") grid_psi_square(i, j)
        enddo
       write(21,"(F14.8)",advance="Yes") grid_psi_square(i,dim_time) 
    
    enddo
    
    
    !lets now compute the average
    do i = 1, dim_time
        mean = 0._pr
        sigma = 0._pr
        norm = 0._pr
        do j = 1, dim_space
            mean = mean + (-L + (float(j) - 1)*dx)*dx * abs(grid_psi(j, i))**2
            norm = norm + dx * abs(grid_psi(j,i))**2
        enddo
        
        mean = mean/norm
        sigma = sigma/norm
        
        write(23,*) (i-1)*delta_t, mean , norm 
    enddo

    close(12)
    close(21)
    close(22)
    close(23)
    
end program
