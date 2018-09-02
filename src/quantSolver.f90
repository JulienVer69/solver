program quantSolver 
use heatEquation  
implicit none 
REAL :: start, finish
CHARACTER (len =100) cal_type 
INTEGER :: n,m
REAL*8  :: lambda


             write(*,*) "************************************************"
             write(*,*) "************************************************"
             write(*,*) "******PARTIAL DIFFERENTIAL EQUATIONS************"
             write(*,*) "****************SOLVER**************************"
             write(*,*) "************************************************"

             open(unit=12, file="fichier.cal", form="formatted")
             
             
             read(12,*) cal_type

              
             if ( cal_type == "HEAT_EQUATION") then 
                   read(12,*) n 
                   read(12,*) m
                   read(12,*) lambda 
             else   
                   write(*,*) "calculation label unknown"
                   CALL EXIT(1)  
             endif 

             close(12) 

             call cpu_time(start)
             call start_solver(n,m,lambda)

             call cpu_time(finish)
 
             print '("Time = ",f6.3," seconds.")',finish-start
             write(*,*) "end of the program"


               

end program quantsolver 

