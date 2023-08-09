PROGRAM Parallel_Stored_Hello
    ! USE OMP_LIB
    
    INTEGER ::   i , j  
    real :: thread_id1

    thread_id1 = 1 


    
    do i =1, 10 
        call testi(thread_id1)
    end do 


    print*, thread_id1   



    contains 

        subroutine testi(result)
            real, intent(out):: result 


            result = 10 
        end subroutine 


    
END