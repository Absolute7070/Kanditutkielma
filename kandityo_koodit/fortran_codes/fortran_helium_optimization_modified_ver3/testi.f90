program testi    
  USE omp_lib
  real :: params(4) 
  integer:: numberOfParams 
  	

  
	call omp_set_num_threads(4) ! set number of threads 
	


  params = [0,0,1,1]
  numberOfParams= 10 

  !$OMP PARALLEL DO 
  do i =1, numberOfParams
                
    print*, params  
  end do 
  !$OMP END PARALLEL DO 



 
  
end program testi 