module writeData  
    use globvariables
      IMPLICIT NONE 

        CONTAINS 

 subroutine write_data(file_name_in) 
  CHARACTER (len=*), INTENT(IN) :: file_name_in                     
  REAL, dimension(:,:), allocatable :: data_write   
  INTEGER :: TS,NUMB,i, i_step,j
  INTEGER :: nc,nl 
  INTEGER :: BUFSIZE, thefile   

 

call MPI_FILE_OPEN(MPI_COMM_SELF,file_name_in, & 
             MPI_MODE_RDONLY, & 
                     MPI_INFO_NULL, thefile, ierr)

call MPI_FILE_READ(thefile,TS,1,MPI_INTEGER,& 
    MPI_STATUS_IGNORE, ierr )


call MPI_FILE_READ(thefile,NUMB,1,MPI_INTEGER,& 
    MPI_STATUS_IGNORE, ierr )

 do i_step = 1,TS

          write(WRITE_UNIT,*) " time step :", i_step 
              
do i = 0,NUMB 

call MPI_FILE_READ(thefile,nl,1,MPI_INTEGER,& 
    MPI_STATUS_IGNORE, ierr ) 

call MPI_FILE_READ(thefile,nc,1,MPI_INTEGER,& 
    MPI_STATUS_IGNORE, ierr )  

 
 BUFSIZE = nl*nc 
 allocate(data_write(nl,nc))

call MPI_FILE_READ(thefile,data_write,BUFSIZE,MPI_REAL,& 
    MPI_STATUS_IGNORE, ierr )             

     do j=1,nl
      write(WRITE_UNIT,*) data_write(j,:)
     enddo 


  deallocate(data_write) 

enddo 

              write(WRITE_UNIT,*) "------------------------------------------------------------------"
            
enddo 

  call MPI_FILE_CLOSE(thefile, ierr)
endsubroutine       


end module  
