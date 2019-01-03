module globvariables 
       use mpi  
implicit none 


! parallel variables 
integer :: rank 
integer :: numprocs
integer :: tag
integer :: ierr 
integer :: status(MPI_STATUS_SIZE)

! ouput variables
INTEGER :: WRITE_UNIT 

end module  
