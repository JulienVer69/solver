module adi2DMethod  
    use linearSystemSolver
    use globVariables 
      IMPLICIT NONE  


      INTEGER  :: n,m,Nt
      INTEGER :: np,mp 
      INTEGER :: niter_x,niter_x_glob  
      INTEGER :: niter_y,niter_y_glob  
      REAL, dIMENSION(:,:), ALLOCATABLE :: mesh_init  
      REAL, dIMENSION(:,:), ALLOCATABLE :: mesh_sol   
      REAL  ::  cl_north,cl_south,cl_east,cl_west
      REAL  :: alpha, lambda,rho,dif,cap
      REAL :: init_value
      REAL  ::  dim_x , delta_x, dim_y,  delta_y , dim_t, delta_t
      INTEGER :: time_step  
      INTEGER :: rest_x,rest_y  
      
      CONTAINS 

                
!*************************************************************************
! subroutine 
!  name : read data 
!*************************************************************************
                
subroutine adi_2d_read_data()  
 INTEGER :: stat 
 character(len=100) :: line
      
      dim_x = .0
      dim_y = .0 
      dim_t = .0 
      n = 0 
      m = 0
      rho = .0
      cap = .0 
      lambda = .0 
      cl_east = .0
      cl_west = .0 
      cl_south = .0 
      cl_north = .0 
      
          

      
  DO WHILE (.TRUE.) 
              READ(READ_UNIT, *, IOSTAT=stat) line 
              IF(IS_IOSTAT_END(stat)) EXIT 

           SELECT CASE (line) 
          
           
               case("GEOMETRY_RECTANGLE") 
                  
                      READ(READ_UNIT,*) dim_x,dim_y

               case("CARTESIAN_GRID") 
                  
                      READ(READ_UNIT,*) n,m

               case("TIME_DEPENDANT") 
                 
                      READ(READ_UNIT,*) dim_t, time_step  

                        


               case("BOUNDARY_CONDITION") 
                 
                      READ(READ_UNIT,*) cl_east,cl_west,cl_north, cl_south
                              

               case ("INIT_ALL_GRID") 

                      READ(READ_UNIT,*) init_value 

                case ("PHYSICS_HEAT") 

                      READ(READ_UNIT,*) lambda,rho,cap

                        
                case default  

                      write(WRITE_UNIT,*) "unknown data set element" 

                end select  
   
      enddo 

      close(READ_UNIT)  
! control

  if ( lambda == .0 .or. rho == .0 .or. cap == .0 ) then 
                                write(*,*) "bad constants"
                                call mpi_finalize(ierr)
                                call exit   
  endif


 if ( dim_t == .0 .or. time_step == 0 ) then 

                                write(*,*) "bad time discretization"
                                call mpi_finalize(ierr)
                                call exit   
      endif



 if ( n == 0 .or. m == 0 ) then 
                              write(*,*) "bad domain discretization"
                              call mpi_finalize(ierr)
                              call exit   
 endif

  if ( dim_x == .0 .or. dim_y == .0 ) then 
                              write(*,*) "bad domain definition"
                              call mpi_finalize(ierr)
                              call exit   
   endif




 
      delta_x = dim_x/(1.0*n) 
      delta_y = dim_y/(1.0*m) 
      delta_t = dim_t/(1.0*time_step) 
      dif = lambda/(rho*cap) 

      if ( rank == 0 ) then 
      write(WRITE_UNIT,*) "HEAT EQUATION SOLVER"
      write(WRITE_UNIT,*) "parameters of the equation :"
      write(WRITE_UNIT,*) "thermal diffusivity :" , dif , "m^2/s"
      write(WRITE_UNIT,*) "------------------------------------------------------------------"
      endif  
       
end subroutine 

!*************************************************************************
! subroutine 
! name : start_solver  
!*************************************************************************

 subroutine adi_2d_start_solver()

if ( rank == 0 ) then 
write(WRITE_UNIT,*) "         number of points to compute : ", (n-2)*(m-2) 
write(WRITE_UNIT,*) "         number of process : " , numprocs   
write(WRITE_UNIT,*) "------------------------------------------------------------------"
                   
endif 
  
! x direction 
         niter_x_glob = (m-2)/numprocs
         rest_x = mod((m-2),numprocs) 

         if ( rank .lt. rest_x ) then 
            niter_x = niter_x_glob +1
         else
            niter_x = niter_x_glob 
         endif   
      
         mp = niter_x+2


! y direction 


         niter_y_glob = (n-2)/numprocs
         rest_y = mod((n-2),numprocs) 

         if ( rank .lt. rest_y ) then 
            niter_y = niter_y_glob +1
         else
            niter_y = niter_y_glob 
         endif   
      
         np = niter_y+2 

         call adi_2d_simulation() 

end subroutine

!*************************************************************************
! subroutine 
! name : mesh_initilialization  
!*************************************************************************

subroutine adi_2d_mesh_initialization()

                mesh_init = init_value 
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


subroutine adi_2d_simulation() 
integer BUFSIZE, thefile 
integer(kind=MPI_OFFSET_KIND) disp  
REAL, dimension(:,:), allocatable ::  local_mesh
REAL, dimension(:,:), allocatable ::  mesh_write 
REAL, dimension(:,:), allocatable ::  matrix 
INTEGER :: i,j,i_step,mini,maxi  
REAL, dimension(:), allocatable ::  tab
INTEGER :: nb,ne,mb,me  
INTEGER :: niter_send,np_send, niter_receiv  
INTEGER :: mesh_init_cursor, local_mesh_cursor, write_cursor  
INTEGER :: offset  
INTEGER :: nl,nc,TS,NUMB 


if ( rank == 0 ) then 
write(WRITE_UNIT,*) "method of resolution : Alternating direction implicit method (2-D)"
write(WRITE_UNIT,*) "------------------------------------------------------------------"
endif 



allocate(mesh_init(np,m))

call adi_2d_mesh_initialization() 

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

 !meta data  -> number of blocks by time edition

    call MPI_FILE_WRITE_AT_ALL(thefile,disp,numprocs-1,1, MPI_INTEGER, & 
                        MPI_STATUS_IGNORE, ierr) 


                
do i_step = 1,time_step 

offset = (i_step-1)*(n*m + numprocs*2) + 2
offset = offset*4


! x direction 

allocate ( matrix(m-2,m-2) )
allocate (tab(m-2))
allocate ( local_mesh(np,m) )

local_mesh(1,:) = mesh_init(1,:)  
local_mesh(np,:) = mesh_init(np,:) 
local_mesh(2:np-1,m) = mesh_init(2:np-1,m)
local_mesh(2:np-1,1) = mesh_init(2:np-1,1)

           do i =2,niter_y+1
               CALL adi_2d_arrays_init(1,i,tab,matrix)
               CALL lapack_solver(matrix,tab,m-2,1)
               local_mesh(i,2:m-1) = tab
           enddo

deallocate(matrix)
deallocate(tab)


! y direction calculation initialization 
! We set the local boundary conditions : 

 

deallocate(mesh_init)
allocate(mesh_init(n,mp))
mesh_init(1,:) = cl_north
mesh_init(n,:) = cl_south 

local_mesh_cursor = 1 
mesh_init_cursor = 2

do i=0,numprocs-1

      if ( rest_x /= 0 ) then 
        if ( i >= rest_x ) then 
           niter_send = niter_x_glob
        else
           niter_send = niter_x_glob +1 
        endif      
      else
           niter_send = niter_x_glob 
      endif 

         if ( rest_y /= 0 ) then 
        if ( i >= rest_y ) then 
           niter_receiv = niter_y_glob
        else
           niter_receiv = niter_y_glob +1 
        endif      
      else
           niter_receiv = niter_y_glob 
      endif 


! Déplacement horizontal sur la grille : local_mesh  
        
         mb = local_mesh_cursor 
         local_mesh_cursor = local_mesh_cursor + niter_send  
         me = local_mesh_cursor + 1  
        
! Déplacement vertical sur la grille : mesh_init 

         nb = mesh_init_cursor    
         mesh_init_cursor = mesh_init_cursor + niter_receiv  
         ne =  mesh_init_cursor -1 



        if ( i == rank ) then
                 mesh_init(nb:ne,1:mp) = local_mesh(2:np-1,mb:me)
                 write_cursor = nb
        else           
                call MPI_Sendrecv(local_mesh(2:np-1,mb:me),niter_y*(niter_send+2),MPI_REAL,i,rank+101,&
                mesh_init(nb:ne,1:mp),niter_receiv*mp, MPI_REAL,i,i+101,MPI_COMM_WORLD,status,ierr)
        endif


         
enddo

deallocate(local_mesh)



allocate ( matrix(n-2,n-2) )


allocate (tab(n-2))
allocate(local_mesh(n,mp))

local_mesh(1,:) = mesh_init(1,:)
local_mesh(n,:) = mesh_init(n,:)
local_mesh(2:n-1,mp) = mesh_init(2:n-1,mp)
local_mesh(2:n-1,1) = mesh_init(2:n-1,1)
       
            do i =2,niter_x+1
               CALL adi_2d_arrays_init(2,i,tab,matrix)
               CALL lapack_solver(matrix,tab,n-2,1)
               local_mesh(2:n-1,i) = tab
            enddo

            

deallocate( matrix )
deallocate( tab )

deallocate(mesh_init)
allocate(mesh_init(np,m))
mesh_init(:,m) = cl_east 
mesh_init(:,1)  = cl_west 

local_mesh_cursor = 2 
mesh_init_cursor = 1

do i=0,numprocs-1

   if ( rest_y /= 0 ) then 
        if ( i >= rest_y ) then 
           niter_send = niter_y_glob
        else
           niter_send = niter_y_glob +1 
        endif      
      else
           niter_send = niter_y_glob 
      endif  

          if ( rest_x /= 0 ) then 
        if ( i >= rest_x ) then 
           niter_receiv = niter_x_glob
        else
           niter_receiv = niter_x_glob +1 
        endif      
      else
           niter_receiv = niter_x_glob 
      endif 

! Déplacement horizontal sur la grille : mesh_init  
        
         mb = local_mesh_cursor 
         local_mesh_cursor = local_mesh_cursor + niter_receiv  
         me = local_mesh_cursor - 1  
        
! Déplacement vertical sur la grille : local_mesh  

         nb = mesh_init_cursor    
         mesh_init_cursor = mesh_init_cursor + niter_send   
         ne =  mesh_init_cursor +1 

        if ( i == rank ) then   
           mesh_init(1:np,mb:me) = local_mesh(nb:ne,2:mp-1)
        else  
              call MPI_Sendrecv(local_mesh(nb:ne,2:mp-1),niter_x*(niter_send+2),MPI_REAL,i,rank+101,&
              mesh_init(1:np,mb:me),niter_receiv*np, MPI_REAL,i,i+101,MPI_COMM_WORLD,status,ierr)
        endif

enddo  

deallocate( local_mesh )


! now we save the results in the binary file 
 
     BUFSIZE = niter_y
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
 write(WRITE_UNIT,*) "calculation time step", i_step , "finished"
endif 

enddo 

 deallocate(mesh_init)
 

 call MPI_FILE_CLOSE(thefile, ierr)
 

end subroutine

!*************************************************************************
! subroutine : 
! name : adi init method  
!*************************************************************************


subroutine adi_2d_arrays_init(step,line,tab,matrix) 
INTEGER :: step
REAL :: tab(:) 
REAL :: matrix(:,:) 
INTEGER, INTENT(IN)  :: line 
INTEGER :: i,j   

matrix =0. 
tab =0.

if ( step == 1) then

 
        alpha = (dif*delta_t)/(2.0*delta_x*delta_x) 

do i = 1,m-2           
   do j=1,m-2
      
   if ( i == j ) then 
         matrix(i,j) = 1.0 +2.0*alpha  
      end if
      if (j == i+1 .and. i < m-2 ) then  
         matrix(i,j) = -1.0*alpha 
      end if
      if (  j == i-1 .and. i>1 ) then 
         matrix(i,j) = -1.0*alpha 
      end if
    enddo
enddo


! time 
     tab(:) = (1-2.0*alpha)*mesh_init(line,2:m-1) 

! top side  
     tab(:) = tab(:) + 1.0*alpha*mesh_init(line-1,2:m-1) 

! lower side  
     tab(:) = tab(:) + 1.0*alpha*mesh_init(line+1,2:m-1)

! left side 
     tab(1) = tab(1) +  1.0*alpha*mesh_init(line,1)

! right side 
     tab(m-2) = tab(m-2) +  1.0*alpha*mesh_init(line,m)
endif 

if ( step == 2) then  

        alpha = (dif*delta_t)/(2.0*delta_y*delta_y) 

 do i = 1,n-2           
   do j=1,n-2
      
    if ( i == j ) then 
         matrix(i,j) = 1.0 +2.0*alpha  
      end if
      if (j == i+1 .and. i < n-2 ) then  
         matrix(i,j) = -1.0*alpha 
      end if
      if (  j == i-1 .and. i>1 ) then 
         matrix(i,j) = -1.0*alpha 
      end if
    enddo
enddo
                        
! time 
     tab(:) = (1-2.0*alpha)*mesh_init(2:n-1,line) 
! left side  
     tab(:) = tab(:) + 1.0*alpha*mesh_init(2:n-1,line-1) 
! right side  
     tab(:) = tab(:) + 1.0*alpha*mesh_init(2:n-1,line+1)
! top of the grid  
     tab(1) = tab(1) +  1.0*alpha*mesh_init(1,line)
! bottom of the grid 
     tab(n-2) = tab(n-2) +  1.0*alpha*mesh_init(n,line)
endif 

end subroutine 

end module 
