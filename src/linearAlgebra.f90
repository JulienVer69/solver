module linearAlgebra  
    use linearSystemSolver
    use para 
      IMPLICIT NONE  


      INTEGER  :: n,m,Nt
      INTEGER :: np
      INTEGER :: niter  
      REAL, dIMENSION(:,:), ALLOCATABLE :: mesh_init  
      REAL, dIMENSION(:,:), ALLOCATABLE :: mesh_sol   
      REAL  ::  cl_north,cl_south,cl_east,cl_west
      REAL  :: alpha, lambda,rho,dif,cap
      REAL  ::  dim_x , delta_x, dim_y,  delta_y , dim_t, delta_t, time_step 
      INTEGER :: WRITE_UNIT
      CHARACTER (len=80) :: file_name_output         

        CONTAINS 

                
!*************************************************************************
! subroutine 
!  name : read data 
!*************************************************************************
                
subroutine read_data(READ_UNIT)  
         INTEGER, INTENT(IN) :: READ_UNIT  
      
      read(READ_UNIT,*) dim_x , n 
      read(READ_UnIT,*) dim_y,  m 
      read(READ_UNIT,*) dim_t, time_step 
      read(READ_UNIT,*) lambda 
      read(READ_UNIT,*) rho
      read(READ_UNIT,*) cap
      read(READ_UNIT,*) cl_north, cl_south, cl_east, cl_west



      close(READ_UNIT)  
 
      delta_x = dim_x/(1.0*n) 
      delta_y = dim_y/(1.0*m) 
      delta_t = dim_t/(1.0*time_step) 

      dif = lambda/(rho*cap) 
      alpha = (dif*delta_t)/(2.0*delta_x*delta_x) 

      if ( rank == 0 ) then 
      write(*,*) "HEAT EQUATION SOLVER"
      write(*,*) "parameters of the equation :"
      write(*,*) "thermal diffusivity :" , dif , "m^2/s"
      write(*,*) "coeff :" , alpha
      write(*,*) "------------------------------------------------------------------"
      endif  
       
      cl_north = cl_north +273
      cl_south = cl_south +273
      cl_east  = cl_east + 273
      cl_west  = cl_west + 273

    

end subroutine 

!*************************************************************************
! subroutine 
! name : start_solver  
!*************************************************************************

 subroutine start_solver()
INTEGER :: i 
INTEGER :: rest 

if ( rank == 0 ) then 
write(*,*) "         number of points to compute : ", (n-2)*(m-2) 
write(*,*) "         number of process : " , numprocs   
write(*,*) "------------------------------------------------------------------"
                   
endif 
         
         niter = (n-2)/numprocs
         rest = mod((n-2),numprocs) 

         if ( rank .lt. rest ) then 
            niter = niter +1
         endif 
      
         write(*,*) "nb d'it : ", niter   
         np = niter+2 


         allocate(mesh_init(np,m))
           
         call mesh_initialization()  
          

           call adi_method() 

   
         

         close(WRITE_UNIT) 
         deallocate(mesh_init)


end subroutine

!*************************************************************************
! subroutine 
! name : mesh_initilialization  
!*************************************************************************

subroutine mesh_initialization()

                mesh_init = 293.
                mesh_init(:,1) = cl_west                 
                mesh_init(:,m) = cl_east                 
           
           if ( rank == 0 ) then 
                mesh_init(1,:) = cl_north                        
           endif         

           if ( rank == numprocs-1 ) then                    
                mesh_init(np,:) = cl_south 
           endif 

                 

             
end subroutine


!*************************************************************************
! subroutine : 
! name : adi method  
!*************************************************************************


subroutine adi_method() 
!integer i,BUFSIZE, thefile 
!    parameter (BUFSIZE=100) 
!    character buf(BUFSIZE) 
!    integer(kind=MPI_OFFSET_KIND) disp  
REAL, dimension(:,:), allocatable ::  local_mesh
REAL, dimension(:,:), allocatable ::  matrix 
INTEGER :: i,j,i_step,mini,maxi  
REAL, dimension(:), allocatable ::  tab
INTEGER :: nc, nb,ne   

if ( rank == 0 ) then 
write(*,*) "method of resolution : Alternating direction implicit method (2-D)"
write(*,*) "------------------------------------------------------------------"
endif 

  WRITE_UNIT = rank+100   
         write(File_name_output,'(A,I3.3,A)') 'init/GenericName',rank,'.ext'
       
         open(unit=WRITE_UNIT,file = File_name_output)

allocate ( matrix(m-2,m-2) )
allocate (tab(m-2))

        do i_step = 1,time_step


allocate ( local_mesh(np,m) )


          do i =2,niter+1
               CALL adi_method_init(1,i,tab,matrix)
               CALL lapack_solver(matrix,tab,n-2,1)
               local_mesh(i,2:m-1) = tab
           enddo




 local_mesh(1,:) = mesh_init(1,:)
 local_mesh(np,:) = mesh_init(np,:)
 local_mesh(2:np-1,m) = mesh_init(2:np-1,m)
 local_mesh(2:np-1,1) = mesh_init(2:np-1,1)



deallocate(mesh_init)
allocate(mesh_init(m,np))

mesh_init(1,:) = cl_north
mesh_init(m,:) = cl_south 


do i=0,numprocs-1

        nb = i*niter +2
        ne = nb + niter-1 

        if ( i == rank ) then  
         mesh_init(nb:ne,1:np) = local_mesh(2:np-1,nb-1:ne+1)
        endif  

   if (rank /= i) then 
  call MPI_Sendrecv(local_mesh(2:np-1,nb-1:ne+1),niter*np,MPI_REAL,i,rank+101,&
  mesh_init(nb:ne,1:np),niter*np, MPI_REAL,i,i+101,MPI_COMM_WORLD,status,ierr)
   endif        
enddo 


deallocate(local_mesh)
allocate(local_mesh(m,np))

           do i =2,niter+1
               CALL adi_method_init(2,i,tab,matrix)
               CALL lapack_solver(matrix,tab,m-2,1)
               local_mesh(2:m-1,i) = tab
           enddo


 local_mesh(1,:) = mesh_init(1,:)
 local_mesh(m,:) = mesh_init(m,:)
 local_mesh(2:np-1,np) = mesh_init(2:np-1,np)
 local_mesh(2:np-1,1) = mesh_init(2:np-1,1)


deallocate(mesh_init)
allocate(mesh_init(np,m))

mesh_init(:,np) = cl_east 
mesh_init(:,1)  = cl_west 


do i=0,numprocs-1

        nb = i*niter +2
        ne = nb + niter-1 

        if ( i == rank ) then   
         mesh_init(1:np,nb:ne) = local_mesh(nb-1:ne+1,2:np-1)
          endif  

   if (rank /= i) then 
  call MPI_Sendrecv(local_mesh(nb-1:ne+1,2:np-1),niter*np,MPI_REAL,i,rank+101,&
  mesh_init(1:np,nb:ne),niter*np, MPI_REAL,i,i+101,MPI_COMM_WORLD,status,ierr)
   endif        
enddo


if (rank == 0 ) then 
 do i = 1,np 
write(*,*) mesh_init(i,:)
 enddo
 endif 

 
deallocate( local_mesh )

 if ( rank == 0 ) then 
 write(*,*) "calculation time step", i_step , "finished"
endif 

enddo 


 deallocate( matrix )
 deallocate( tab )

!write(*,*) "je suis le proc : " , rank, " je calcul de : ", begin_iter, "à ", end_iter 
        
 !       do i = 0, BUFSIZE 
 !       buf(i) = 'b' 
 !   enddo 
 !   call MPI_FILE_OPEN(MPI_COMM_WORLD, 'testfile', & 
 !                      MPI_MODE_WRONLY + MPI_MODE_CREATE, & 
 !                      MPI_INFO_NULL, thefile, ierr) 
    ! assume 4-byte integers 
 !   disp = rank * BUFSIZE * 4 
 !   call MPI_FILE_SET_VIEW(thefile, disp, MPI_INTEGER, & 
 !                          MPI_CHAR, 'native', & 
 !                          MPI_INFO_NULL, ierr) 
 !   call MPI_FILE_WRITE(thefile, buf, BUFSIZE, MPI_CHAR, & 
 !                       MPI_STATUS_IGNORE, ierr) 
 !   call MPI_FILE_CLOSE(thefile, ierr) 


! print message to screen

!if (rank == 0 ) then 
!      write(*,*) 'Hello World!'
!endif




end subroutine

!*************************************************************************
! subroutine : 
! name : adi init method  
!*************************************************************************


subroutine adi_method_init(step,line,tab,matrix) 
INTEGER :: step
REAL :: tab(:) 
REAL :: matrix(:,:) 
INTEGER, INTENT(IN)  :: line 
INTEGER :: i,j  

matrix =0. 
tab =0.

do i = 1,m-2           
   do j=1,m-2
      ! A matrix
      if ( i == j ) then 
         matrix(i,j) = 1.0 +4.0*alpha  
      end if
      if (j == i+1 .and. i < n-2 ) then  
         matrix(i,j) = -2.0*alpha 
      end if
      if (  j == i-1 .and. i>1 ) then 
         matrix(i,j) = -2.0*alpha 
      end if
    enddo
enddo


if ( step == 1) then  
! time 
     tab(:) = (1-4.0*alpha)*mesh_init(line,2:m-1) 

! top side  
     tab(:) = tab(:) + 2.0*alpha*mesh_init(line-1,2:m-1) 

! lower side  
     tab(:) = tab(:) + 2.0*alpha*mesh_init(line+1,2:m-1)

! left side 
     tab(1) = tab(1) +  2.0*alpha*mesh_init(line,1)

! right side 
     tab(m-2) = tab(m-2) +  2.0*alpha*mesh_init(line,m)
endif 

if ( step == 2) then  
! time 
     tab(:) = (1-4.0*alpha)*mesh_init(2:m-1,line) 
! left side  
     tab(:) = tab(:) + 2.0*alpha*mesh_init(2:m-1,line-1) 
! right side  
     tab(:) = tab(:) + 2.0*alpha*mesh_init(2:m-1,line+1)
! top of the grid  
     tab(1) = tab(1) +  2.0*alpha*mesh_init(1,line)
! bottom of the grid 
     tab(m-2) = tab(m-2) +  2.0*alpha*mesh_init(m,line)
endif 

end subroutine 


end module 
