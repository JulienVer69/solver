program quantSolver 
use heatEquation  
implicit none 
REAL :: start, finish
CHARACTER (len =100) cal_type
character(len=100) :: file_name
!INTEGER :: n,m
!REAL*8  :: lambda

call get_command_argument(1, file_name)

open(unit=12, file=file_name, form="formatted")

             write(*,*) "*****************************************************"
             write(*,*) "***********PARTIAL DIFFERENTIAL EQUATIONS************"
             write(*,*) "*********************SOLVER**************************"
             write(*,*) "************************************************"

             
         call cpu_time(start)
             
             read(12,*) cal_type

              
             if ( cal_type == "HEAT_EQUATION") then  
                CALL read_data_heq()
                call start_solver()



             else   
                   write(*,*) "calculation label unknown"
                   CALL EXIT(1)  
             endif 

             call cpu_time(finish)
             close(12) 

              
             print '("Time = ",f6.3," seconds.")',finish-start
             write(*,*) "end of the program"


               

end program quantsolver 

