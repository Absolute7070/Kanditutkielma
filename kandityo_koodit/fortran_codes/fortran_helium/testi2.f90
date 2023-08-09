PROGRAM test_shape
    use omp_lib
    use fparser 
    use parameters
    implicit none 
    ! Symbolic variables 
    character(len=*), dimension(2), parameter :: variables = [character(len=5) :: 'x1', 'y1']
     
  
    !parameters as values in the same order as in symbolic variables 
    real(rn), dimension(2) :: parameterValues = [1,2 ] 
    character(len= 20) :: func 
    real(rn) :: result 
  
    func = "exp(x1+y1 )"
  
  
    
    !$OMP PARALLEL 
      
      call evaluateexpressions(func, variables, parameterValues, result)
  
    !$OMP END PARALLEL 
    
  
  
  
  
   
  
  
  
      contains 
      subroutine evaluateexpressions(func, var, val, res ) 
        implicit none 
  
  
        character(len=*),   intent(in) :: func
        INTEGER,                             PARAMETER :: nfunc = 1
  
        INTEGER,                             PARAMETER :: nvar = 2
        CHARACTER (LEN=*), DIMENSION(nvar),  intent(in):: var  
        REAL(rn),          DIMENSION(nvar),  intent(in) :: val  
        REAL(rn), intent(out)                                      :: res
     
  
        
        CALL initf (nfunc)                      ! Initialize function parser for nfunc functions
        
        CALL parsef (1, func, var)        ! Parse and bytecompile ith function string 
        
        res = evalf (1, val)                 ! Interprete bytecode representation of ith function
        IF (EvalErrType > 0) WRITE(*,*)'*** Error: ',EvalErrMsg ()
        ! WRITE(*,*)'res=',res  
        
        
  
        
    end subroutine evaluateexpressions
  
    
      
    END PROGRAM