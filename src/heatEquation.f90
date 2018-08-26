module heatEquation 
      implicit none 

         integer   :: n,m,Nt
         REAL*8, dimension(:,:), allocatable :: matrix
         REAL*8, dimension(:), allocatable  :: cl_tab
         REAL*8, dimension(:), allocatable :: x
         REAL*8, dimension(:), allocatable :: cl_north
         REAL*8, dimension(:), allocatable :: cl_south
         REAL*8, dimension(:), allocatable :: cl_east
         REAL*8, dimension(:), allocatable :: cl_west
         REAL*8  :: lambda  

! External Subroutines ..
      EXTERNAL         SGESV
contains 


        
  subroutine start_solver()
       INTEGER :: i
           
    
             write(*,*) "****************************************"
             write(*,*) "****************************************"
             write(*,*) "**********2D HEAT EQUATION**************"
             write(*,*) "****************************************"
             write(*,*) "****************************************"
     
             n=3
             m=3
             Nt = n*m 

             call adi_init()

             call clear_memory()

     end subroutine 
  


   subroutine adi_init()
       integer :: i,j,k
       write(*,*) "***********initialization step***************" 
       write(*,*) "****************ADI METHOD*******************" 

        allocate ( matrix(Nt,Nt) )
        allocate ( cl_north(m) )
        allocate ( cl_south(m) )
        allocate ( cl_east(n) )
        allocate ( cl_west(n) )
        allocate ( cl_tab(Nt) )
        allocate ( x(Nt) )
         
        x = 0.
        matrix =0.
        cl_tab =0.
        cl_east =0. 
        cl_west =0. 
        cl_south =0. 
        cl_north =0. 

        do i=1,m 
        cl_north(i) = 20 
        cl_south(i) = 12
        enddo

        do i=1,n
        cl_east(i) = 15
        cl_west(i) = 14
        enddo 


!  Fist we resolve the second x derivative in the explicitly way
!  We store the triadiagonal matrix for this case 
!  ( (1 + 4)*lambda      -2*lambda          0     ---------   0    -----   0  )    
!  ( -2*lambda         (1 + 4)*lambda   -2*lambda ---------   0    .....   0  )
!  (    -                                   --                                )
!  (    -                                             --                      )
!  (    -                                                    --               )
!  (    0 -------------------                     -2*lambda    (1 + 4)*lambda )    


        k=-1
       do i=1,Nt
          do j=1,Nt
           k=k+1
           if (k>0) then 
                  matrix(k,k) = 2.0*(1.0+lambda)
                  matrix(k,k-1) = -lambda 
           endif 
           
           if (k < Nt) then 
                  matrix(k,k+1) = -lambda
           endif
         enddo
        enddo  


! And the right-hand side of the equation : 

            
     end subroutine 



    subroutine clear_memory() 
     deallocate(x)
     deallocate(cl_tab) 
     deallocate(matrix) 
     deallocate( cl_south ) 
     deallocate( cl_north ) 
     deallocate( cl_east ) 
     deallocate( cl_west )
    end subroutine  


end module  
