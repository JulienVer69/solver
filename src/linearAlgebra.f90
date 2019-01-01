module linearAlgebra  
    use linearSystemSolver
    use para
    use writeData
      IMPLICIT NONE  


      INTEGER  :: n,m,Nt
      INTEGER :: np
      INTEGER :: niter,niter_glob  
      REAL, dIMENSION(:,:), ALLOCATABLE :: mesh_init  
      REAL, dIMENSION(:,:), ALLOCATABLE :: mesh_sol   
      REAL  ::  cl_north,cl_south,cl_east,cl_west
      REAL  :: alpha, lambda,rho,dif,cap
      REAL  ::  dim_x , delta_x, dim_y,  delta_y , dim_t, delta_t
      INTEGER :: time_step  
      INTEGER :: rest 
      
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

if ( rank == 0 ) then 
write(*,*) "         number of points to compute : ", (n-2)*(m-2) 
write(*,*) "         number of process : " , numprocs   
write(*,*) "------------------------------------------------------------------"
                   
endif 
         
         niter_glob = (n-2)/numprocs
         rest = mod((n-2),numprocs) 

         if ( rank .lt. rest ) then 
            niter = niter_glob +1
         else
            niter = niter_glob 
         endif   
      
         np = niter+2 


         allocate(mesh_init(np,m))
           
         call mesh_initialization()  
          

           call adi_method() 

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
integer BUFSIZE, thefile 
integer(kind=MPI_OFFSET_KIND) disp  
REAL, dimension(:,:), allocatable ::  local_mesh
REAL, dimension(:,:), allocatable ::  mesh_write 
REAL, dimension(:,:), allocatable ::  matrix 
INTEGER :: i,j,i_step,mini,maxi  
REAL, dimension(:), allocatable ::  tab
INTEGER :: nb,ne,mb,me  
INTEGER :: niter_send,np_send 
INTEGER :: mesh_init_cursor, local_mesh_cursor, write_cursor  
INTEGER :: offset  
INTEGER :: nl,nc,TS,NUMB 

if ( rank == 0 ) then 
write(*,*) "method of resolution : Alternating direction implicit method (2-D)"
write(*,*) "------------------------------------------------------------------"
endif 

! We create a binary file to save the results of the calculation in parallel 
! each task will write his results in a specific location of the file :

call MPI_FILE_OPEN(MPI_COMM_WORLD, 'testfile', & 
             MPI_MODE_WRONLY + MPI_MODE_CREATE, & 
                     MPI_INFO_NULL, thefile, ierr) 

disp = 0 

! meta data  -> time step 
    call MPI_FILE_WRITE_AT_ALL(thefile,disp,time_step,1, MPI_INTEGER, & 
                        MPI_STATUS_IGNORE, ierr) 

disp = disp + 4

! meta data  -> number of blocks by time edition

    call MPI_FILE_WRITE_AT_ALL(thefile,disp,numprocs-1,1, MPI_INTEGER, & 
                        MPI_STATUS_IGNORE, ierr) 


                
allocate ( matrix(m-2,m-2) )
allocate (tab(m-2))


do i_step = 1,time_step 

offset = (i_step-1)*(n*m + numprocs*2) + 2
offset = offset*4


allocate ( local_mesh(np,m) )

          do i =2,niter+1
               CALL adi_method_init(1,i,tab,matrix)
               CALL lapack_solver(matrix,tab,n-2,1)
               local_mesh(i,2:m-1) = tab
           enddo


! We set the local boundary conditions : 

 local_mesh(1,:) = mesh_init(1,:)  
 local_mesh(np,:) = mesh_init(np,:) 
 local_mesh(2:np-1,m) = mesh_init(2:np-1,m)
 local_mesh(2:np-1,1) = mesh_init(2:np-1,1)


deallocate(mesh_init)
allocate(mesh_init(m,np))
mesh_init =.0
mesh_init(1,:) = cl_north
mesh_init(m,:) = cl_south 

local_mesh_cursor = 1 
mesh_init_cursor = 2

do i=0,numprocs-1

      if ( rest /= 0 ) then 
        if ( i >= rest ) then 
           niter_send = niter_glob
        else
           niter_send = niter_glob +1 
        endif      
      else
           niter_send = niter_glob 
      endif  
 

! Déplacement horizontal sur la grille : local_mesh  
        
         mb = local_mesh_cursor 
         local_mesh_cursor = mb + niter_send  
         me = local_mesh_cursor + 1  
        
! Déplacement vertical sur la grille : mesh_init 

         nb = mesh_init_cursor    
         mesh_init_cursor = nb + niter_send  
         ne =  mesh_init_cursor -1 


        if ( i == rank ) then
                 mesh_init(nb:ne,1:np) = local_mesh(2:np-1,mb:me)
                 write_cursor = nb 
        else           
                call MPI_Sendrecv(local_mesh(2:np-1,mb:me),niter*(niter_send+2),MPI_REAL,i,rank+101,&
                mesh_init(nb:ne,1:niter+2),niter_send*(niter+2), MPI_REAL,i,i+101,MPI_COMM_WORLD,status,ierr)
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
local_mesh(2:m-1,np) = mesh_init(2:m-1,np)
local_mesh(2:m-1,1) = mesh_init(2:m-1,1)


deallocate(mesh_init)
allocate(mesh_init(np,m))

mesh_init(:,m) = cl_east 
mesh_init(:,1)  = cl_west 

local_mesh_cursor = 2 
mesh_init_cursor = 1

do i=0,numprocs-1

   if ( rest /= 0 ) then 
        if ( i >= rest ) then 
           niter_send = niter_glob
        else
           niter_send = niter_glob +1 
        endif      
      else
           niter_send = niter_glob 
      endif  

! Déplacement horizontal sur la grille : mesh_init  
        
         mb = local_mesh_cursor 
         local_mesh_cursor = mb + niter_send  
         me = local_mesh_cursor - 1  
        
! Déplacement vertical sur la grille : local_mesh  

         nb = mesh_init_cursor    
         mesh_init_cursor = nb + niter_send  
         ne =  mesh_init_cursor +1 

        if ( i == rank ) then   
           mesh_init(1:np,mb:me) = local_mesh(nb:ne,2:np-1)
        else  
              call MPI_Sendrecv(local_mesh(nb:ne,2:np-1),niter*(niter_send+2),MPI_REAL,i,rank+101,&
              mesh_init(1:np,mb:me),niter_send*(niter+2), MPI_REAL,i,i+101,MPI_COMM_WORLD,status,ierr)
        endif 
enddo  

deallocate( local_mesh )

! now we save the results in the binary file 
 
     BUFSIZE = niter
     nb = 2
     ne = np-1

    if ( rank == 0 ) then 
       BUFSIZE = BUFSIZE + 1
       nb = 1
       write_cursor = 1
    endif 
   if ( rank == numprocs-1) then    
        BUFSIZE = BUFSIZE + 1 
        ne = np
   endif 
 
      disp = (write_cursor-1)*m*4 + offset + rank*8    
   
! meta data
! number of lines in block
    call MPI_FILE_SET_VIEW(thefile,disp, MPI_INTEGER, & 
                           MPI_INTEGER, 'native', & 
                           MPI_INFO_NULL, ierr) 
    call MPI_FILE_WRITE(thefile,BUFSIZE,1, MPI_INTEGER, & 
                        MPI_STATUS_IGNORE, ierr) 

                disp = disp + 4

! number of columns in block 

    call MPI_FILE_SET_VIEW(thefile,disp, MPI_INTEGER, & 
                           MPI_INTEGER, 'native', & 
                           MPI_INFO_NULL, ierr) 
    call MPI_FILE_WRITE(thefile,m,1, MPI_INTEGER, & 
                        MPI_STATUS_IGNORE, ierr) 

                disp = disp + 4


! data
       BUFSIZE = BUFSIZE*m 

           
     call MPI_FILE_SET_VIEW(thefile, disp, MPI_REAL, & 
                           MPI_REAL, 'native', & 
                           MPI_INFO_NULL, ierr) 
    call MPI_FILE_WRITE(thefile,mesh_init(nb:ne,:),BUFSIZE, MPI_REAL, & 
                        MPI_STATUS_IGNORE, ierr) 


 if ( rank == 0 ) then 
 write(*,*) "calculation time step", i_step , "finished"
endif 

enddo 

 deallocate( matrix )
 deallocate( tab )

 call MPI_FILE_CLOSE(thefile, ierr)
 
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
