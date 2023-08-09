PROGRAM test_shape
	use omp_lib
	implicit none 
  integer :: i 
  real :: arr(100, 100), k(100, 100)

  arr= 30 
	
	
	call omp_set_num_threads(4) ! set number of threads 
	
	!$OMP PARALLEL PRIVATE(i, k, arr )


  PRINT '("Thread number: ", i0)', omp_get_thread_num()

  k =0 

  DO i =1, 1000000
    k = matmul(arr, arr) **i + k 
  end do 


  print*, k(1, 1)
  
	
	!$OMP END PARALLEL 
  
    
END PROGRAM