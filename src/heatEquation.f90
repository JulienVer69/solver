module heatEquation 
    use linearSystemSolver
      implicit none 

         INTEGER  :: Nt
         REAL*8, dimension(:,:), allocatable :: matrix
         REAL*8, dimension(:), allocatable  :: cl_tab
         REAL*8, dimension(:), allocatable :: x
         REAL*8, dimension(:), allocatable :: cl_north
         REAL*8, dimension(:), allocatable :: cl_south
         REAL*8, dimension(:), allocatable :: cl_east
         REAL*8, dimension(:), allocatable :: cl_west
 
contains 

        
  subroutine start_solver(n,m,lambda)
       INTEGER, INTENT(IN) :: n 
       INTEGER, INTENT(IN) :: m 
       REAL*8, INTENT(IN) :: lambda
    
             write(*,*) "****************************************"
             write(*,*) "****************************************"
             write(*,*) "**********2D HEAT EQUATION**************"
             write(*,*) "****************************************"
             write(*,*) "****************************************"
     
             Nt = n*m


            call adi_method(n,m,lambda)

            call clear_memory()

     end subroutine 


!*************************************************************************
! We use the Nicholson scheme for two spatial heat equation. More 
! preciesly we use the adi method because of the numercal instabilities
! of the previous method. 
!*************************************************************************

subroutine adi_method(n,m,lambda) 
INTEGER :: i
INTEGER, INTENT(IN) ::n,m
REAL*8, INTENT(IN) :: lambda

        allocate ( matrix(Nt,Nt) )
        allocate ( cl_north(m) )
        allocate ( cl_south(m) )
        allocate ( cl_east(n) )
        allocate ( cl_west(n) )
        allocate ( cl_tab(Nt) )
        allocate ( x(Nt) )

        
        x = 0.
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
          
        CALL adi_init(n,m,lambda,1) 
        CALL lapack_solver(matrix,cl_tab,Nt,1)
        CALL adi_init(n,m,lambda,2)
        CALL lapack_solver(matrix,cl_tab,Nt,1)

        x = cl_tab

        write(*,*) x(:)

end subroutine


!*************************************************************************
! adi_init() 
!this subroutine intialize all the matrices that we need to compute the
!solution. First we solve the equation for time 
!step t=t+1/2 and with the result of this operation we compute the solution 
!**************************************************************************

subroutine adi_init(n,m,lambda,time_step)           
       INTEGER :: i,j,k
       INTEGER, INTENT(IN) ::n,m
       REAL*8, INTENT(IN) :: lambda
       INTEGER,INTENT(IN) :: time_step 
       write(*,*) "***********initialization step***************" 
       write(*,*) "****************ADI METHOD*******************" 

                
        matrix =0.
       

     


!  Fist we resolve the second x derivative in the explicitly way
!  We store the triadiagonal matrix for this case 
!  ( 1 + 4*lambda      -2*lambda          0     ---------   0    -----   0    )    
!  ( -2*lambda         1 + 4*lambda   -2*lambda ---------   0    .....   0    )
!  (    -                                   --                                )
!  (    -                                             --                      )
!  (    -                                                    --               )
!  (    0 -------------------                     -2*lambda      1 + 4*lambda )    
!
! multiple of lambda -> lambda = 2*lambda 

    do i=1,Nt
      do j=1,Nt

        ! matrix A 
         if ( j == i ) then 
             matrix(i,j) = 1.0 +4.0*lambda  
          end if
         if (j == i+1 .and. i < Nt ) then  
         matrix(i,j) = -2.0*lambda 
          end if
          if (  j == i-1 .and. i>1) then 
          matrix(i,j) = -2.0*lambda 
          end if
     enddo
    enddo

         
! And the right-hand side of the equation :
        do i=1,n 
          do j=1,m 
            k = (i-1)*m +j 
! time part  
            cl_tab(k) = cl_tab(k) + (1.0 - 2.0*lambda)*x(k) 

           if (time_step == 1) then             
! right side of the grid 
                if (j == m) then
                cl_tab(k) = cl_tab(k) + 2.0*lambda*cl_east(i) + x(k-1)
! left side of the grid 
                else if (j == 1) then
           cl_tab(k) = cl_tab(k) + 2.0*lambda*cl_west(i) + x(k+1) 
! center of the gride 
                else 
           cl_tab(k) = cl_tab(k) + 2.0*lambda*(x(k+1) + x(k-1)) 
                endif  
           endif
        if (time_step == 2) then 
! upper side of the grid 
                if (i == 1) then
                cl_tab(k) = cl_tab(k) + 2.0*lambda*cl_north(j) + x(k+m)
! lower side of the grid 
                else if (i == n) then
           cl_tab(k) = cl_tab(k) + 2.0*lambda*cl_south(j) + x(k-m) 
! center of the gride 
                else 
           cl_tab(k) = cl_tab(k) + 2.0*lambda*(x(k+m) + x(k-m)) 
                endif  
           endif
          enddo
        enddo
           
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
