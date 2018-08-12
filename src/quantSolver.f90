program quantSolver 
use laplace  
implicit none 
real :: start, finish


             write(*,*) "************************************************"
             write(*,*) "************************************************"
             write(*,*) "******PARTIAL DIFFERENTIAL EQUATIONS************"
             write(*,*) "****************SOLVER**************************"
             write(*,*) "************************************************"

             call cpu_time(start)

             call start_laplace()

             call cpu_time(finish)
 
             print '("Time = ",f6.3," seconds.")',finish-start
             write(*,*) "end of the program"


               

end program quantsolver 

