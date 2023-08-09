program optimization 

    ! version 2: thermalisation only in certain interval of optimization steps 
    ! date: 20.4.



    !!!!!!!!!!!!!!!!!THINGS TO KNOW !!!!!!!!!!!!!!!!!!!!!!!!!!
    ! - Number of configuration subspaces is determined, how much you set thread to be 
    ! - when using evaluateexpressions, remember to trim() the expression (fparser)
    ! - thermalisation_samples.txt: each row is different configuration, and each three consecutive number represents 
    !   one particle's x,y,z-coordinates



    ! Known error causing problems and solutions: 
    !     1. Evaluation error due to string length: Change function string length MAX_FUN_LENGTH  in fparser, when string length exceeded




    USE mt19937
    USE omp_lib
    use FortranParser, only : EquationParser
    use FortranParser_parameters, only: rn
    

    
    implicit none

    ! integer kinds 
    integer, parameter :: ki8 = selected_int_kind(18)          ! double precision integer



    ! measuring performance 
    real(rn) :: start_time, end_time 
    real(rn):: total_start_time, total_end_time



    !!!!!!!WHEN CHANGE ATOM, CHANGE THESES!!!!!!!!!!!




    integer(ki8), parameter :: characterlengthEnergy = 50000, characterlengthProb = 300      ! check based on wc-command in bash 
    
    
    integer, parameter :: numberOfCoordinates = 6, numberOfParams = 2  ! number of coordinates and number of parameters 
    integer, parameter :: numberOfIterations = 5, numberOfConfig = 10, numberOfParticles = 2     ! number of iterations and number of configurations 


    ! Symbolic variables 
    character(len=*), dimension(numberOfCoordinates+numberOfParams), parameter :: variables = [character(len=5) :: 'x1', 'y1', 'z1', 'x2', 'y2', 'z2', 'A1', 'A2']

    ! coordinate value range 
    real(rn), dimension(2), parameter :: coordinateValueRange =[-6, 6]



    ! thermalisation related 
    integer, parameter :: numberOfParallel = 8 ! number of cores to be used 
    integer, parameter :: intervalForThermals = 10 ! interval for the next thermalisation 
    integer :: iseed 


    ! optimization related 
    real(rn), dimension(numberOfParams) :: parameterValues = [0.21336979283896812 ,    3.9996145624259256]  
    real(rn), parameter :: learning_rate=0.001, tolerance =0.00000000000000000000001
    integer, parameter :: maximumsteps = 50


    !!!!!!!END OF: WHEN CHANGE ATOM, CHANGE THESES!!!!!!!!!!!

    

    !!!!!!!! USER DEFINED FROM COMMAND-LINE !!!!!!!!!!!!!
    ! expressions for probability and local energy, user defined 
    character(len=characterlengthEnergy) :: localenergyexpr
    character(len =characterlengthProb) :: probabilityexpr  
                 

    !!!!!!!! END OF: USER DEFINED FROM COMMAND-LINE !!!!!!!!!!!!!
 


    !!!!!!!! Command line stuff !!!!!!!!!!!!

    ! dealing with command line stuff 
    integer :: i, iarg                     ! number of command line arguments                       
    character(len=80) :: arg                    ! command line argument 
    integer :: ios
    logical :: there              ! whether the thermalised samples already exists 
    !!!!!!!! END OF: Command line stuff !!!!!!!!!!!!    


    





    !!!!!!!! THERMALISATION VARIABLES !!!!!!!!!!!!        
    real(rn), allocatable :: samplesInConfigSubspace(:, :, :)
    integer :: counter 


    !!!!!!!! END OF: THERMALISATION VARIABLES !!!!!!!!!!!! 


 


    








    call omp_set_nested(.true.)     ! enable nested parallelism 
      ! set number of threads to be used 
    call omp_set_num_threads(numberOfParallel)


    call cpu_time(total_start_time)


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






    !!!!!! USER (COMMAND-LINE) INTERFACE !!!!!!!!


    !  1. filepath of probability expression 
    !  2. filepath of local energy expression


    ! check: correct number of command-line arguments 
    iarg = command_argument_count()
    if (iarg/=2) then 
        print*, "ERROR: invalid number of arguments"
        print*, "Usage: ./a.out ", & 
        "[filepath of probability expression] [filepath of local energy expr]"
        stop 
    end if 
    
        
    !1. argument: filepath probability expr
    call get_command_argument(1, arg )
    open(unit=1, file = arg, iostat= ios, status = 'old')
    if (ios/=0) then 
        print '(a, a, 5x, i0)', 'ERROR: failed in opening file named ', trim(arg), ios 
        stop 
    end if 


    ! assuming one-line expression file
    read(1, '(a)', iostat = ios) probabilityexpr 
    if (ios/=0) then 
        print '(a)', 'ERROR: reading probabilityexpr failed!'
        stop 
    end if 
    

    ! 2. argument: local energy file 
    call get_command_argument(2, arg )
    open(unit=1, file = arg, iostat= ios, status = 'old')
    if (ios/=0) then 
        print '(a, a, 5x, i0)', 'ERROR: failed in opening file named ', trim(arg), ios 
        stop 
    end if 



    ! assuming one-line expression file
    read(1, '(a)', iostat = ios) localenergyexpr
    if (ios/=0) then 
        print '(a)', 'ERROR: reading localenergyexpr failed! '
        stop 
    end if 

  
    !!!!!! END OF USER (COMMAND-LINE) INTERFACE !!!!!!!!



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!OPTIMIZATION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! initiate random number for Mersenne Twister 
    iseed = time()
    call sgrnd(iseed)


    ! ! get thermalisation_samples from txt file, if previous exists 
    ! inquire(file = 'thermalisation_samples.txt', exist= there )
	! if (there) then 
    !     allocate(samplesInConfigSubspace(numberOfConfig, 3, numberOfParticles)) 


    !     ! read in thermalisation data
    !     open(unit=1,file='thermalisation_samples.txt',iostat=ios,status='old')
    !     if (ios/=0) then
    !         print '(a,a,5x,i0)','*** Error in opening file thermalisation_samples.txt.',ios
    !         stop
    !     end if


    !     do i =1, numberOfConfig
    !         read(1, *)  samplesInConfigSubspace(i, :, :)
    !     end do 

    !     close(1)

    


        


    ! else 
    !     call cpu_time(start_time)

    !     !Generating samples and thermalisation 
    !     print*, '----------Thermalisation started------------------'
    !     call thermalisation(numberOfIterations, samplesInConfigSubspace) 
    !     print*, '----------Thermalisation ended--------------------'

    !     call cpu_time(end_time)
    !     print*, 'TIME USED IN THERMALISATION: ', end_time-start_time


    !     ! save data for later faster optimization 
    !     open(unit=1,file='thermalisation_samples.txt',iostat=ios,status='new')
    !     ! prepared for an error
    !     if (ios/=0) then
    !         print '(a,a)','*** Error in opening file thermalisation_samples.txt.'
    !         stop
    !     end if
    
    !     do i =1, numberOfConfig 
    !         write(1, *) samplesInConfigSubspace(i, :, :)
    !     end do 

    !     close(1)

	! end if



    

    ! initiate counter 
    counter = 0 

    do i = 1, maximumsteps


        ! thermalisation only after certain number of time of optimization 
        if (counter == 0) then 
            call cpu_time(start_time)

            !Generating samples and thermalisation 
            call thermalisation(numberOfIterations, samplesInConfigSubspace) 
            call cpu_time(end_time)
            print*, 'TIME USED IN THERMALISATION: ', end_time-start_time

        end if 



        print*, 'Counter: ', counter

        ! print*, 'ennen gradient'


        call cpu_time(start_time)
        ! gradient descent 
        call gradientdescentsub(learning_rate, samplesInConfigSubspace, parameterValues, tolerance, 1 )

        
        call cpu_time(end_time)
        print*, 'TIME USED IN GRADIENT DESCENT ', end_time-start_time
       
       
        ! print*, 'jälkeen gradient gradient'

        PRINT*, 'Final parameters', parameterValues
        print '(70("-"))'


        

        ! mechanism: thermalisation only after certain number of steps of optimization 
        if (counter < intervalForThermals) then 

            counter = counter + 1 
        elseif (counter >= intervalForThermals) then 
            counter = 0 
            deallocate(samplesInConfigSubspace)
        end if 



        

    end do 


  



    







    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!END OF: OPTIMIZATION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    call cpu_time(total_end_time)
    print*, 'TOTAL TIME USED: ', total_end_time-total_start_time
    



    



        




        


   




    contains 


        subroutine gradientdescentsub(learning_rate, configurationspace, params, tolerance, maximumsteps )
            real(rn), intent(in) :: learning_rate, tolerance
            integer, intent(in) :: maximumsteps 
            real(rn), dimension(numberOfConfig, 3, numberOfParticles), intent(in):: configurationspace 
            real(rn), dimension(numberOfParams), intent(inout) :: params

            
            real(rn), dimension(numberOfParams) :: gradient
            integer :: counter 

            gradient = gradient_atapoint(configurationspace, params) ! compute gradient at initial point 

            
            
            counter = 0 
            DO WHILE (ANY( abs(gradient) > tolerance) )
                print*, 'Gradient: ', gradient 
               
                print*, 'Parameters: ', params 
                


                params = params - learning_rate * gradient                ! gradient descent: update parameters 
                gradient = gradient_atapoint(configurationspace, params)   ! compute gradient at current point 




                counter = counter + 1 

                ! exit if maximumsteps (user-defined) exceeded
                if (counter >= maximumsteps) then 
                    exit 
                end if 
            END DO 

            
        end subroutine gradientdescentsub

        function gradient_atapoint(configurationspace,  params)
            implicit none
            real(rn), dimension(numberOfParams) :: gradient_atapoint 

            real(rn), dimension(numberOfConfig, 3, numberOfParticles), intent(in):: configurationspace 
            real(rn), dimension(numberOfParams), intent(in) :: params

            integer :: i 


            
            !$OMP DO PRIVATE(i)
            do i =1, numberOfParams
                
                 gradient_atapoint(i ) = fivePointDiff(configurationspace, params, i )
            end do 
            !$OMP END DO 
            
            


           
            
        end function gradient_atapoint



        ! five point difference method 

        real(rn) function fivePointDiff(configurationspace,  params, whichParameter)
        implicit none
        real(rn), dimension(numberOfConfig, 3, numberOfParticles), intent(in):: configurationspace 
        real(rn), dimension(numberOfParams), intent(in) :: params
        integer, intent(in) :: whichParameter 

        real(rn), parameter :: h = 0.0001 ! diff. size  
        real(rn) :: numerator, denominator 
        real(rn), dimension(numberOfParams) :: parameters_minus2h, parameters_plus2h, parameters_minus1h, parameters_plus1h
        real(rn) :: energy_plus1h, energy_minus1h, energy_plus2h, energy_minus2h

        ! make copy of parameters list 
        parameters_minus1h = params 
        parameters_plus1h = params 
        parameters_minus2h = params 
        parameters_plus2h = params 



        parameters_minus1h(whichParameter)  = parameters_minus1h(whichParameter) -h 
        parameters_plus1h(whichParameter)  = parameters_plus1h(whichParameter) +h 
        parameters_plus2h(whichParameter) = parameters_plus2h(whichParameter) + 2*h 
        parameters_minus2h(whichParameter) = parameters_minus2h(whichParameter) - 2*h 




        
        !$OMP SECTIONS 
        !$OMP SECTION 
        call mcenergy(configurationspace, parameters_plus1h,energy_plus1h )               ! E(alpha+1h )
        !$OMP SECTION 
        call mcenergy(configurationspace, parameters_minus1h,energy_minus1h)              ! E(alpha-1h )
        !$OMP SECTION 
        call mcenergy(configurationspace, parameters_minus2h,energy_minus2h )              ! E(alpha-2h )
        !$OMP SECTION 
        call mcenergy(configurationspace, parameters_plus2h,energy_plus2h )              ! E(alpha + 2h )
        !$OMP END SECTIONS 

        numerator = energy_minus2h- 8*energy_minus1h+ 8*energy_plus1h-energy_plus2h
        denominator = 12*h 


        fivePointDiff = numerator/denominator



        end function fivePointDiff



        ! two point difference method
        real(rn) function twoPointDiff(configurationspace,  params, whichParameter)
            implicit none
            real(rn), dimension(numberOfConfig, 3, numberOfParticles), intent(in):: configurationspace 
            real(rn), dimension(numberOfParams), intent(in) :: params
            integer, intent(in) :: whichParameter 

            real(rn), parameter :: h = 0.0001 ! diff. size  
            real(rn) :: numerator, denominator 
            real(rn), dimension(numberOfParams) :: parameters_minush, parameters_plush
            real(rn) :: energy1, energy2 

            ! make copy of parameters list 
            parameters_minush = params 
            parameters_plush = params 


            parameters_minush(whichParameter)  = parameters_minush(whichParameter) -h 
            parameters_plush(whichParameter)  = parameters_plush(whichParameter) +h 

            
            !$OMP SECTIONS 
            !$OMP SECTION 
            call mcenergy(configurationspace, parameters_plush,energy1 )               ! E(alpha+1 )
            !$OMP SECTION 
            call mcenergy(configurationspace, parameters_minush,energy2 )              ! E(alpha-1 )
            !$OMP END SECTIONS 

            numerator = energy1- energy2
            denominator = 2*h 


            twoPointDiff = numerator/denominator



        end function twoPointDiff


        subroutine mcenergy(configurationspace, params, energy)
            implicit none
            real(rn), dimension(numberOfConfig, 3, numberOfParticles), intent(in):: configurationspace 
            real(rn), dimension(numberOfParams), intent(in) :: params
            real(rn), intent(out) :: energy 
            


            integer :: i 
            real(rn), dimension(numberOfCoordinates+numberOfParams) :: flattened_configurationAndParams 
            real(rn) :: localEnergy 
            real(rn), dimension(numberOfConfig) :: localEnergyList 


            !$OMP DO PRIVATE(i,  flattened_configurationAndParams, localEnergy)
            configurationdo: DO i=1, numberOfConfig 
                

                ! get flattened configuration with parameters at the tail 
                flattened_configurationAndParams(: numberOfCoordinates)  = reshape(configurationspace(i, :, :), [numberOfCoordinates])
                flattened_configurationAndParams(numberOfCoordinates+1:) = params 
                
          
                
                call evaluateexpressions(trim(localenergyexpr), variables, flattened_configurationAndParams, localEnergy )
       


                localEnergyList(i) = localEnergy

            END DO configurationdo
            !$OMP END DO

            
            energy = (1.0_rn /size(localEnergyList)) * SUM(localEnergyList )

          




        end subroutine mcenergy





        subroutine thermalisation(iterations, samples)
            implicit none
            integer, intent(in) :: iterations 
            real(rn), allocatable, intent(out) :: samples(:, :, :)


            real(rn) :: delta = 0.5
            real(rn), allocatable :: deltai(:)
            real(rn), allocatable :: trial_configuration(:, :), current_configuration(:, : )
            real(rn) :: trial_coordinates(3)
            real(rn), allocatable :: trial_configuration_vec(:), current_configuration_vec(: )
            real(rn) :: trial_probability, old_probability 

            

            integer :: i, j, k
            real(rn) :: r, w 

            ! each table represent one configuration, each row one particle's xyz-coordinates. 
            ! different tables represent different independent sample  
            allocate(samples(numberOfConfig, 3, numberOfParticles))  

            
            
            call random_uniform_tensor(coordinateValueRange(1), coordinateValueRange(2), samples, shape(samples) )
            ! print*, samples 
            ! stop




            allocate(trial_configuration(3, numberOfParticles))  ! initialize trial configuration table 
            allocate(current_configuration(3, numberOfParticles)) ! initialize current configuration table 

            ! for flattening the configuration tables into vectors 
            allocate(trial_configuration_vec(size(trial_configuration)+numberOfParams))  
            allocate(current_configuration_vec(size(current_configuration)+numberOfParams))

            ! assign parameters' values to the end 
            trial_configuration_vec(3*numberOfParticles+1:) = parameterValues
            current_configuration_vec(3*numberOfParticles+1:) = parameterValues



            ! thermalisation 
            DO i=1, iterations

                DO j=1, numberOfConfig 
                    DO k=1,  numberOfParticles 
                        ! create step size 
                        call random_uniform_vector(-delta, delta , deltai, 3) 

                        trial_coordinates = samples(j, :, k) + deltai 

                        trial_configuration = samples(j, :, :)  ! make copy of initial samples 
                        trial_configuration(:, k )  =  trial_coordinates      ! make trial move of k-th particle

                        current_configuration = samples(j, :, : )  ! make copy of current configuration 




                        ! vectorize so that we can use them to evaluate probability expression 
                        trial_configuration_vec(: 3*numberOfParticles) = reshape(trial_configuration, shape(trial_configuration_vec(1:3*numberOfParticles)))    
                        current_configuration_vec(: 3*numberOfParticles)= reshape(current_configuration, shape(current_configuration_vec(1:3*numberOfParticles)))

                        


                        
                        
                        
                        !$OMP SECTIONS 
                        !$OMP SECTION 
                        call evaluateexpressions(trim(probabilityexpr), variables, trial_configuration_vec, trial_probability )
                        !$OMP SECTION 
                        call evaluateexpressions(trim(probabilityexpr), variables, current_configuration_vec, old_probability ) 
                        !$OMP END SECTIONS
                        
                       

                        w= trial_probability/old_probability 

                        if (w>=1) then 
                            current_configuration = trial_configuration
                        else if (w<1) then 
                            call random_number(r) ! generate uniformly distributed random number 0<= x <1 
                            if (r<=w) then 
                                current_configuration = trial_configuration 
                            end if 
                        end if 

                        samples(j, :, :) = current_configuration


                    END DO 

                END DO 

            END DO 






        end subroutine thermalisation










        subroutine evaluateexpressions(func, var, val, res ) 
            implicit none 
      
            type(EquationParser) :: eqParser
            character(len=*),   intent(in) :: func
            INTEGER,                             PARAMETER :: nfunc = 1
      
            INTEGER,                             PARAMETER :: nvar = numberOfCoordinates+ numberOfParams
            CHARACTER (LEN=*), DIMENSION(nvar),  intent(in):: var  
            REAL(rn),          DIMENSION(nvar),  intent(in) :: val  
            REAL(rn), intent(out)                                      :: res
      
            eqParser = EquationParser(func, var)
      
            res = eqParser%evaluate(val)
         
      
            
        end subroutine evaluateexpressions





        !!!!!!!!!!!!!!RANDOM  NUMBER GENERATORS !!!!!!!!!!!!!!!!!!!!


        subroutine random_uniform_tensor(lower, upper, arr, arr_shape )
            implicit none 
            real(rn),intent(in) :: lower, upper     ! lower and upper limit 
            integer, dimension(3), intent(in) :: arr_shape 
            real(rn), allocatable, intent(out) :: arr(:, :, : )

            integer :: i, j, k 
            real(rn) :: random_number 


            allocate(arr(arr_shape(1), arr_shape(2), arr_shape(3)))


            do i= 1, arr_shape(1)
                do j =1, arr_shape(2)
                    do k=1, arr_shape(3)
                        call random_uniform(lower, upper, random_number)
                        arr(i, j, k  ) = random_number
                    end do 
                end do 
            end do 
            



        end subroutine random_uniform_tensor




        subroutine random_uniform_array(lower, upper, arr, arr_shape )
            implicit none 
            real(rn),intent(in) :: lower, upper     ! lower and upper limit 
            integer, dimension(2), intent(in) :: arr_shape 
            real(rn), allocatable, intent(out) :: arr(:, :)

            integer :: i, j 
            real(rn) :: random_number 


            allocate(arr(arr_shape(1), arr_shape(2)))


            do i= 1, arr_shape(1)
                do j =1, arr_shape(2)
                        call random_uniform(lower, upper, random_number)
                        arr(i, j ) = random_number
                end do 
            end do 
            



        end subroutine random_uniform_array

        subroutine random_uniform_vector(lower, upper, vec, vec_length )
            implicit none 
            real(rn),intent(in) :: lower, upper     ! lower and upper limit 
            integer, intent(in) :: vec_length
            real(rn), allocatable, intent(out) :: vec(:)

            integer :: i, j 
            real(rn) :: random_number 


            allocate(vec(vec_length))


            do i= 1, vec_length
                call random_uniform(lower, upper, random_number)
                vec(i) = random_number
            end do 
            



        end subroutine random_uniform_vector




        ! assuming a<b, generate uniform distribution a<x<=b 
        subroutine random_uniform(a,b,x)
            implicit none
            real(rn),intent(in) :: a,b
            real(rn),intent(out) :: x
            ! integer, optional, intent(inout) :: iseed 
            real(rn) :: u


            ! if (present(iseed).eqv..false.) then 
            !     iseed = time()
            ! end if 

            call random_stduniform(u )
            x = (b-a)*u + a
        end subroutine random_uniform



        subroutine random_stduniform(u)
            implicit none
            real(rn),intent(out) :: u
            ! integer, optional, intent(inout) :: iseed 
            real(rn) :: r
            

            ! if (present(iseed).eqv..false.) then 
            !     iseed = time()
            ! end if 


            call random_number(r)

            u = 1 - r
         end subroutine random_stduniform

        subroutine random_number(r )
            implicit none
            real(rn), intent(out) :: r 
            ! integer, optional, intent(inout) :: iseed 

            ! if (present(iseed).eqv..false.) then 
            !     iseed = time()
            ! end if 
        
            ! call sgrnd(iseed)

            r = grnd()

        end subroutine random_number




        !!!!!!!!!!!!!!END OF: RANDOM  NUMBER GENERATORS !!!!!!!!!!!!!!!!!!!!





        







    
end program optimization 