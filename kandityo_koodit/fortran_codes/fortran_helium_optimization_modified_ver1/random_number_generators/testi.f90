program testi 
    use mt19937
    implicit none
    
    
    integer :: i, j , iseed 


    iseed = time()
    call sgrnd(iseed)
    print*, grnd()

    print*, iseed

    
end program testi         

