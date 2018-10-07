module heatEquation 
    use linearSystemSolver
      implicit none 

         
         REAL  ::  dim_x , delta_x,  dim_y,  delta_y , dim_t, delta_t 
         INTEGER  :: n,m,Nt,time_step
         REAL  :: alpha, lambda,rho,dif,cap 
         REAL, dimension(:,:), allocatable :: matrix
         REAL, dimension(:), allocatable  :: cl_tab
         REAL, dimension(:), allocatable :: x
         REAL  ::  cl_north,cl_south,cl_east,cl_west
 
contains


subroutine read_data_heq( ) 
      read(12,*) dim_x , n 
      read(12,*) dim_y,  m 
      read(12,*) dim_t, time_step 
      read(12,*) lambda 
      read(12,*) rho
      read(12,*) cap
      read(12,*) cl_north, cl_south, cl_east, cl_west
 
      n = n-2
      m = m-2 

      !n = int(dim_x/delta_x) -2 
      !m = int(dim_y/delta_y) -2
      !time_step = int(dim_t/delta_t) 
   
       delta_x = dim_x/(1.0*n) 
       delta_y = dim_y/(1.0*m) 
       delta_t = dim_t/(1.0*time_step) 

      dif = lambda/(rho*cap) 
      alpha = (dif*delta_t)/(2.0*delta_x*delta_x)  
  
 write(*,*) "thermal diffusivity :" , dif , "m^2/s"


 write(*,*) "coeff :" , alpha 

 


end subroutine 

        
  subroutine start_solver()
    
             write(*,*) "***********************************************************"
             write(*,*) "*******************2D HEAT EQUATION************************"
             write(*,*) "***********************************************************"
          
              write(*,*) "taille :",n,m 

              Nt = n*m  

            !call adi_opt_method()

            call adi_method()

            call clear_memory()

     end subroutine 


!*************************************************************************
! We use the Crank-Nicolson scheme for two spatial heat equation. More 
! precisely we use the adi method because of the numerical instabilities
! of the previous method. 
!*************************************************************************

subroutine adi_method() 
INTEGER :: i,j,k,p,i_step 
CHARACTER (len=80) :: file_name_output  
REAL*8 :: temp 



        allocate ( matrix(Nt,Nt) )
        allocate ( cl_tab(Nt) )
        allocate ( x(Nt) )

        
        x = 293.
        cl_tab =0.

        cl_north = cl_north +273
        cl_south = cl_south +273
        cl_east  = cl_east + 273
        cl_west  = cl_west + 273


        do i_step = 0,time_step 

         write(File_name_output,'(A,I3.3,A)') 'data2/GenericName',i_step,'.ext'
         open(unit=13,file = File_name_output)  

         if (i_step == 0 ) then 
                 CALL  write_mesh(13)

                 write(*,*) "calculation time step", i_step , "finished" 
         else  



        call adi_init(1)
        !call gaussJordanMethod(matrix,cl_tab,Nt)
        call lapack_solver(matrix,cl_tab,Nt,1) 
        
        x = cl_tab
        

        call adi_init(2) 
       ! call gaussJordanMethod(matrix,cl_tab,Nt)

        call lapack_solver(matrix,cl_tab,Nt,1)



        do i=1,n 
          do j=1,m
           
              k = (i-1)*m + j 
              p = (j-1)*n + i

              x(k) = cl_tab(p) 
         enddo 
        enddo       
 
      call write_mesh(13) 
      
           write(*,*) "calculation time step", i_step , "finished"
            
    endif

 enddo 

end subroutine


!*************************************************************************
! adi_init() 
!this subroutine intialize all the matrices that we need to compute the
!solution. First we solve the equation for time 
!step t=t+1/2 and with the result of this operation we compute the solution 
!**************************************************************************

subroutine adi_init(step)           
       INTEGER :: i,j,k,p
       INTEGER,INTENT(IN) :: step 
     !  write(*,*) "***********initialization step***************" 
     !  write(*,*) "****************ADI METHOD*******************" 

     if ( step /= 1 .and. step /=2 ) then 

             write(0,*) " error :" 
             write(0,*) "step not defined"
             write(0,*) "in adi_init subroutine : heatEquation.f90 "
             CALL exit(0) 
     endif 

        matrix =0.
        cl_tab =0.
       


    do i=1,Nt
      do j=1,Nt

         if ( j == i ) then 
             matrix(i,j) = 1.0 +4.0*alpha  
          end if
         if (j == i+1 .and. i < Nt .and. modulo(i,m) /= 0 ) then  
         matrix(i,j) = -2.0*alpha 
          end if
          if (  j == i-1 .and. i>1 .and. modulo(j,m) /=0) then 
          matrix(i,j) = -2.0*alpha 
          end if
     enddo
    enddo
 


    do i =1,n
       do j=1,m 
      

          k = (i-1)*m +j

     if ( step == 1) then 

          cl_tab(k) = (1 - 4.0*alpha)*x(k) 

! top side of the point  
          if ( i == 1) then       
          cl_tab(k) = cl_tab(k) + 2.0*alpha*cl_north 
          else
          cl_tab(k) = cl_tab(k) + 2.0*alpha*x(k-m)
          endif  
! lower side of the point 
          if ( i == n ) then 
          cl_tab(k) = cl_tab(k) + 2.0*alpha*cl_south  
          else
          cl_tab(k) = cl_tab(k) + 2.0*alpha*x(k+m)
          endif  
! left side of the point
          if ( j == 1 ) then 
          cl_tab(k) = cl_tab(k) +  2.0*alpha*cl_west
          endif 
! right side of the point            
          if ( j == m ) then  
          cl_tab(k) = cl_tab(k)  + 2.0*alpha*cl_east  
          endif 
  endif 

        if ( step == 2 ) then 

           p = (j-1)*n + i

           cl_tab(p) = (1 - 4.0*alpha)*x(k)

! top side of the point 
          if ( i == 1) then             
          cl_tab(p) = cl_tab(p) + 2.0*alpha*cl_north  
          endif  
! lower side of the point  
          if ( i == n ) then 
          cl_tab(p) = cl_tab(p)  + 2.0*alpha*cl_south 
          endif 
! left side of the point 
          if ( j == 1 ) then 
          cl_tab(p) = cl_tab(p) + 2.0*alpha*cl_west
          else 
          cl_tab(p) = cl_tab(p) + 2.0*alpha*x(k-1) 
          endif 
! right side of the point            
          if ( j == m ) then  
          cl_tab(p) = cl_tab(p) + 2.0*alpha*cl_east    
          else
          cl_tab(p) = cl_tab(p) + 2.0*alpha*x(k+1)     
          endif 
  
     endif 
   
           
         enddo    
           enddo

 
           
end subroutine


subroutine adi_opt_method() 
INTEGER :: i,j,i_step,mini,maxi 
CHARACTER (len=80) :: file_name_output  
REAL, dimension(:), allocatable ::  tab
INTEGER :: k,p 


allocate ( matrix(n,m) )
allocate ( cl_tab(Nt) )
allocate ( x(Nt) )
allocate (tab(m))

x = 293.

cl_north = cl_north +273
cl_south = cl_south +273
cl_east  = cl_east + 273
cl_west  = cl_west + 273
 
write(File_name_output,'(A,I3.3,A)') 'data/GenericName',0,'.ext'
open(unit=13,file = File_name_output) 
CALL  write_mesh(13)

        do i_step = 1,time_step
           write(File_name_output,'(A,I3.3,A)') 'data/GenericName',i_step,'.ext'

           
           open(unit=13,file = File_name_output)
           do i=1,n
                mini = (i-1)*m +1 
                maxi = i*m 
                CALL adi_opt_init(1,i,tab)
                CALL lapack_solver(matrix,tab,n,1)
                cl_tab(mini:maxi) = tab
           enddo
 
 x = cl_tab
 cl_tab=0.
           do i=1,n
                mini = (i-1)*m +1 
                maxi = i*m   
               CALL adi_opt_init(2,i,tab)
               CALL lapack_solver(matrix,tab,n,1)
               cl_tab(mini:maxi) = tab
           enddo

           
     do i=1,n 
          do j=1,m
           
              k = (i-1)*m + j 
              p = (j-1)*n + i

              x(k) = cl_tab(p) 
         enddo 
        enddo       



            CALL write_mesh(13) 
          enddo  

deallocate(tab)

end subroutine

subroutine adi_opt_init(step,line,tab)
INTEGER :: i,j,k,p, tab_ind
INTEGER :: step
INTEGER :: line
REAL :: tab(:) 
!INTEGER :: mini,maxi 

matrix =0. 
tab =0.

do i = 1,n           
   do j=1,m
      ! A matrix
      if ( i == j ) then 
         matrix(i,j) = 1.0 +4.0*alpha  
      end if
      if (j == i+1 .and. i < n ) then  
         matrix(i,j) = -2.0*alpha 
      end if
      if (  j == i-1 .and. i>1 ) then 
         matrix(i,j) = -2.0*alpha 
      end if
    enddo
enddo


     if ( step == 1) then

do j=1,m
   i=line
   
   k = (i-1)*m +j


          tab(j) = (1 - 4.0*alpha)*x(k)
! top side of the point  
          if ( i == 1) then       
            tab(j) = tab(j) + 2.0*alpha*cl_north 
          else
            tab(j) = tab(j) + 2.0*alpha*x(k-m)
          endif 


! lower side of the point 
          if ( i == n ) then 
            tab(j) =  tab(j) + 2.0*alpha*cl_south  
          else
            tab(j) =  tab(j) + 2.0*alpha*x(k+m)
          endif  

          
! left side of the point
          if ( j == 1 ) then 
            tab(j) =  tab(j) +  2.0*alpha*cl_west
          endif 


! right side of the point            
          if ( j == m ) then  
            tab(j) = tab(j)  + 2.0*alpha*cl_east  
          endif


   enddo
endif 

        if ( step == 2 ) then 

do tab_ind =1,n
       j = line  

           p = (tab_ind-1)*n + j

           tab(tab_ind) = (1 - 4.0*alpha)*x(p)

! top side of the point 
          if ( tab_ind == 1) then             
           tab(tab_ind) =  tab(tab_ind) + 2.0*alpha*cl_north  
          endif  
! lower side of the point  
          if ( tab_ind == n ) then 
           tab(tab_ind) =  tab(tab_ind)  + 2.0*alpha*cl_south 
          endif 
! left side of the point 
          if ( j == 1 ) then 
            tab(tab_ind) =  tab(tab_ind) + 2.0*alpha*cl_west
          else 
            tab(tab_ind) =  tab(tab_ind) + 2.0*alpha*x(p-1) 
          endif 
! right side of the point            
          if ( j == m ) then  
            tab(tab_ind) =  tab(tab_ind) + 2.0*alpha*cl_east    
          else
            tab(tab_ind) =  tab(tab_ind) + 2.0*alpha*x(p+1)     
          endif 
  
enddo
     endif 




end subroutine 


                                    

       subroutine write_mesh(out_unit) 
       INTEGER :: i,j,k
       INTEGER, INTENT(IN) :: out_unit  

           do i=0,n+1 
             do j=0,m+1

                    k = (i-1)*m + j 

                  if (i == 0 )  then  
                  write(out_unit,*)  i,j, cl_north 

                  else if ( i == n+1 ) then   
                         
                  write(out_unit,*) i,j, cl_south 

                 else if ( j == 0 ) then 

                  write(out_unit,*) i,j, cl_west 

                 else if ( j == m+1 ) then 

                  write(out_unit,*) i,j, cl_east 

                 else 

                  write(out_unit,*)  i,j, x(k)

               endif

             enddo
           enddo   
                    


       end subroutine         
 

    subroutine clear_memory() 
     deallocate(x)
     deallocate(cl_tab) 
     deallocate(matrix)

  !   deallocate( cl_south ) 
  !   deallocate( cl_north ) 
  !   deallocate( cl_east ) 
  !   deallocate( cl_west )
    end subroutine  


end module  
