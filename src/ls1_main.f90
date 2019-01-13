program main  
use adi2DMethod 
use mpi
use globVariables 
use writeData
USE ISO_FORTRAN_ENV, ONLY : ERROR_UNIT,OUTPUT_UNIT ! access computing environment
implicit none 
REAL :: start, finish
CHARACTER (len =100) cal_type
character(len=100) :: file_name_input  
character(len=100) :: file_name_output
character (len=100) :: argv
INTEGER*4 i, iargc, numarg


call mpi_init(ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
WRITE_UNIT = OUTPUT_UNIT 


!*****************************************************
!            READING PROGRAM'S ARGUEMENTS 
!*****************************************************

          numarg = iargc()
i=1


        if ( numarg == 0 ) then 

        if (rank == 0 ) then 
        write(ERROR_UNIT,*) " usage : linear_solver -opt <o1> <o2> [...]" 
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

!*********************************************************************    

    if ( rank == 0 ) then 

             write(WRITE_UNIT,*) "------------------------------------------------------------------"
             write(WRITE_UNIT,*) "                       Linear Algebra Solver                      "
             write(WRITE_UNIT,*) "                           version 1.0                            "      
             write(WRITE_UNIT,*) "                                                                  "      
             write(WRITE_UNIT,*) " created and implemented by : Julien Versaci                      "
             write(WRITE_UNIT,*) " demonstration                                                    "
             write(WRITE_UNIT,*) "------------------------------------------------------------------"
endif


READ_UNIT = rank+100            
open(unit=READ_UNIT, file=file_name_input, form="formatted")

          start= MPI_Wtime()     
     



              
                CALL adi_2d_read_data()
                call adi_2d_start_solver()





             finish= MPI_Wtime()     
             
 if ( rank == 0 ) then         

 call write_data("testfile") 


             
            write(WRITE_UNIT,*) "------------------------------------------------------------------"
            write(WRITE_UNIT,*) " Simulation time :                                                " 
           ! print '("Time = ",f6.6," seconds.")',finish-start
            write(WRITE_UNIT,*) "------------------------------------------------------------------"
            write(WRITE_UNIT,*) "                       END OF PROGRAM                             "
            write(WRITE_UNIT,*) "------------------------------------------------------------------"


            
  CLOSE(WRITE_UNIT) 
endif


call mpi_finalize(ierr)
               

end program  

