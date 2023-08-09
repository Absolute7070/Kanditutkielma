program testi    


  use FortranParser, only : EquationParser
  use FortranParser_parameters, only: rn
 
  implicit none

  character, dimension(2) :: variables  =  [character(len=5) :: 'x', 'y']
  real(rn) :: result 

  call evaluateexpressions("x/y", variables, [2, 1], result )

  print*, result 

  contains 

  subroutine evaluateexpressions(func, var, val, res ) 
    implicit none 

    type(EquationParser) :: eqParser
    character(len=*),   intent(in) :: func
    INTEGER,                             PARAMETER :: nfunc = 1

    INTEGER,                             PARAMETER :: nvar = 2 
    CHARACTER (LEN=*), DIMENSION(nvar),  intent(in):: var  
    REAL(rn),          DIMENSION(nvar),  intent(in) :: val  
    REAL(rn), intent(out)                                      :: res

    eqParser = EquationParser(func, var)

    res = eqParser%evaluate(val)
 

    
end subroutine evaluateexpressions


 
  
end program testi 