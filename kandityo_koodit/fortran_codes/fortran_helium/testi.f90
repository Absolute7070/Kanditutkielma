program testi 
  use omp_lib
  implicit none
  integer :: arr(2, 2), deltai(2)

  arr(:, 1) = 1 
  arr(:, 2) = 2 


  deltai(1) = 1 
  deltai(2) = 2 

  print*, arr(:, 1)+ deltai


 
  
end program testi 