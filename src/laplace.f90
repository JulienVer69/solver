module laplace
   use gaussJordan  
      implicit none 

         integer   :: n,m,Nt
         real, dimension(:,:), allocatable :: matrix
         real, dimension(:), allocatable  :: cl_tab
         real, dimension(:), allocatable :: x
         real, dimension(:), allocatable :: cl_north
         real, dimension(:), allocatable :: cl_south
         real, dimension(:), allocatable :: cl_east
         real, dimension(:), allocatable :: cl_west


! External Subroutines ..
      EXTERNAL         SGESV
contains 


        
  subroutine start_laplace() 
             write(*,*) "*********************************"
             write(*,*) "*********************************"
             write(*,*) "******LAPLACE EQUATION***********"
             write(*,*) "*********************************"
             write(*,*) "*********************************"
     
             n=100
             m=100
             Nt = n*m 
                   
             call laplace_init()

             
     !       x=gaussJordanMethod(n,m,matrix,cl_tab) 
            call lapack_solver()  

            x = - cl_tab 
            
            call writeData()
            call clear_memory()

     end subroutine 
  


   subroutine laplace_init()
       integer :: i,j,k
       write(*,*) "...initialization step..." 


        allocate ( matrix(Nt,Nt) )
        allocate ( cl_north(m) )
        allocate ( cl_south(m) )
        allocate ( cl_east(n) )
        allocate ( cl_west(n) )
        allocate ( cl_tab(Nt) )
        allocate ( x(Nt) )

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


       do i=1,Nt
      do j=1,Nt

        ! matrix A 
         if ( j == i ) then 
             matrix(i,j) = -4
          end if
          if ( ( j == i+m .and. i < Nt-m+1  ) .or. ( j == i-m .and. i > m ) ) then  
             matrix(i,j) = 1
          end if
          if ( ( j == i+1 .and. modulo(i,m) /= 0  ) .or. ( j == i-1 .and. modulo(j,m) /= 0   ) ) then 
             matrix(i,j) = 1
          end if

        ! matrix B 
        if (i<=n .and. j<=m) then 

         k = (i-1)*m +j 
        
        if ( i == 1) then 
           cl_tab(k)  = cl_tab(k) +cl_north(j)
        end if 
        if ( i == n) then      
           cl_tab(k) = cl_tab(k) +  cl_south(j)
        end if 
         if (j == m) then
           cl_tab(k) = cl_tab(k) + cl_east(i)
         end if
        if (j == 1) then
           cl_tab(k) = cl_tab(k) + cl_west(i)
        end if  
       end if 
     
      enddo 
        enddo
            
     end subroutine 


     subroutine writeData() 
     integer :: i,j 


     do i=0,n+1 
          do j=0,m+1
           
           if (i == 0 .and. j < m+1 .and. j>0) then 
                write(*,*) i,j,cl_north(j)
           end if
           if (i == n+1 .and. j < m+1 .and. j>0 ) then 
                write(*,*) i,j,cl_south(j)
           end if    
           if (j == m+1 .and. i < n+1 .and. i>0 ) then 
                write(*,*) i,j,cl_east(i)
           end if
           if (j == 0 .and. i < n+1 .and. i>0 ) then 
                write(*,*) i,j,cl_west(i)
           end if    
           if ( i == 0 .and. j == 0 ) then 
                write(*,*) i,j,(cl_north(1)+cl_west(1))/2.0
           end if
           if ( i == 0 .and. j == m+1 ) then 
                write(*,*) i,j,(cl_north(m)+cl_east(1))/2.0
           end if 
           if ( i == n+1 .and. j == 0 ) then 
                write(*,*) i,j,(cl_south(1)+cl_west(n))/2.0
            end if 
            if ( i == n+1 .and. j == m+1 ) then 
                write(*,*) i,j,(cl_south(m)+cl_east(n))/2.0
            end if 
             if ( ( i >0 .and. i<n+1) .and. ( j>0 .and. j<m+1) ) then  
               write(*,*) i,j,x(i)
           end if 
       end do 
     end do  

     
          end subroutine


    subroutine lapack_solver 
      ! Parameters ..
      integer ::     NRHS, LDA , LDB 

      ! local scalar
      INTEGER          INFO 

      ! local array 
      INTEGER          IPIV(Nt)

      

     NRHS = 1
     LDA = Nt
     LDB = Nt

      !!  Executable Statements ..
      WRITE(*,*)'DGESV Example Program Results'
!
!     Solve the equations A*X = B.
!
      CALL SGESV( Nt, NRHS, matrix, LDA, IPIV, cl_tab, LDB, INFO )

  
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
