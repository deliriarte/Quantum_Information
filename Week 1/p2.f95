  !number precision

program precision
  implicit none
  
  integer(kind= 2):: int_a_sp, int_b_sp
  integer(kind = 4) :: int_a_dp, int_b_dp
  
  integer, parameter :: sp = kind(1.0), &
                        dp = kind(1.d0)
  
  real(sp),  parameter :: PI_sp  = 4 * atan (1.0_sp)
  real(dp),  parameter :: PI_dp  = 4 * atan (1.0_dp)
  real(sp) :: r_a_sp, r_b_sp
  real(dp) :: r_a_dp, r_b_dp
  

  write(*,*)'When using an integer of kind', 2, 'the value of a overflow'
  
  int_a_dp = 2000000
  int_b_dp = 1
  
  write(*,'(A,I8,A,I2)')'The sum of', int_a_dp, ' and', int_b_dp,' which are integer with double precision is', int_a_dp + int_b_dp
  
  write(*,*)'The results from part b are'
  write(*,*) 'single precision', PI_sp * (10._sp **32)
  write(*,*) 'single precision', sqrt(2._sp) * (10._sp **21)
  
  write(*,*) 'double precision', PI_dp * (10._dp **32)
  write(*,*) 'double precision',sqrt(2._dp) * (10._dp **21)

 end program
