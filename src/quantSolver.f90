program quantSolver 
!use heatEquation 
use linearAlgebra
use mpi
use para 
implicit none 
REAL :: start, finish
CHARACTER (len =100) cal_type
character(len=100) :: file_name
INTEGER :: READ_UNIT  



call mpi_init(ierr)

call get_command_argument(1, file_name)

call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

READ_UNIT = rank+100

open(unit=READ_UNIT, file=file_name, form="formatted")


if ( rank == 0 ) then 
             write(*,*) "------------------------------------------------------------------"
             write(*,*) "                       Linear Algebra Solver                      "
             write(*,*) "                           version 1.0                            "      
             write(*,*) "                                                                  "      
             write(*,*) " created and implemented by : Julien Versaci                      "
             write(*,*) " demonstration                                                    "
             write(*,*) "------------------------------------------------------------------"

             
         start= MPI_Wtime()     !cpu_time(start)
endif              
             


             read(READ_UNIT,*) cal_type



              
             if ( cal_type == "HEAT_EQUATION") then        
                CALL read_data(READ_UNIT)
                call start_solver()



             else   
                   write(*,*) "calculation label unknown"
                   CALL EXIT(1)  
             endif 

  if ( rank == 0 ) then                     
         finish= MPI_Wtime()     !cpu_time(start)

             
            write(*,*) "------------------------------------------------------------------"
            write(*,*) " Simulation time :                                                " 
            print '("Time = ",f6.6," seconds.")',finish-start
            write(*,*) "------------------------------------------------------------------"
            write(*,*) "                       END OF PROGRAM                             "
            write(*,*) "------------------------------------------------------------------"

endif

call mpi_finalize(ierr)
               

end program quantsolver 

