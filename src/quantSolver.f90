program quantSolver 
use linearAlgebra
use mpi
use para
use writeData
USE ISO_FORTRAN_ENV, ONLY : ERROR_UNIT,OUTPUT_UNIT ! access computing environment
implicit none 
REAL :: start, finish
CHARACTER (len =100) cal_type
character(len=100) :: file_name_input  
character(len=100) :: file_name_output
INTEGER :: READ_UNIT
INTEGER :: WRITE_UNIT 
character (len=100) :: argv
INTEGER*4 i, iargc, numarg


call mpi_init(ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
WRITE_UNIT = OUTPUT_UNIT 


          numarg = iargc()
i=1


        if ( numarg == 0 ) then 

        if (rank == 0 ) then 
        write(*,*) " usage : linear_solver -opt <o1> <o2> [...]" 
        endif

        call mpi_finalize(ierr)
        call exit()

endif 

    do while ( i <= numarg )
    call getarg( i, argv )

       select case (argv)

       case ('-in')
              i=i+1 
              call getarg(i,file_name_input)

        case ('-out') 
              i=i+1
              call getarg(i,file_name_output)
              WRITE_UNIT = 2589 
              open(unit=WRITE_UNIT, file=file_name_output, form="formatted")
 
                
      case default
              if (rank == 0 ) then      
          WRITE(ERROR_UNIT,*)"argument : ", argv, " unknown" 
              endif 

          call mpi_finalize(ierr)
          call exit() 

       endselect     

          
        i=i+1

    enddo 


    if ( rank == 0 ) then 

             write(*,*) "------------------------------------------------------------------"
             write(*,*) "                       Linear Algebra Solver                      "
             write(*,*) "                           version 1.0                            "      
             write(*,*) "                                                                  "      
             write(*,*) " created and implemented by : Julien Versaci                      "
             write(*,*) " demonstration                                                    "
             write(*,*) "------------------------------------------------------------------"
endif


READ_UNIT = rank+100            
open(unit=READ_UNIT, file=file_name_input, form="formatted")

          start= MPI_Wtime()     
         


             read(READ_UNIT,*) cal_type

              
             if ( cal_type == "HEAT_EQUATION") then        
                CALL read_data(READ_UNIT)
                call start_solver()



             else   
                   write(*,*) "calculation label unknown"
                   CALL EXIT(1)  
             endif 


             finish= MPI_Wtime()     
             
 if ( rank == 0 ) then         

 call write_data("testfile",WRITE_UNIT) 


             
            write(*,*) "------------------------------------------------------------------"
            write(*,*) " Simulation time :                                                " 
            print '("Time = ",f6.6," seconds.")',finish-start
            write(*,*) "------------------------------------------------------------------"
            write(*,*) "                       END OF PROGRAM                             "
            write(*,*) "------------------------------------------------------------------"

endif


call mpi_finalize(ierr)
               

end program quantsolver 

