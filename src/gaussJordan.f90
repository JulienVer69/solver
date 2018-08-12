module gaussJordan 
implicit none 
 
 contains

     function gaussJordanMethod(n,m,matrix,cl_tab) 
          integer :: i,j
          integer :: k,r,p
          real :: max_,tmp,tmp2,pivot
          integer, intent(in) :: n,m
          integer :: Nt 
          real :: matrix(:,:) 
          real :: cl_tab(:) 
          real, dimension(:), allocatable ::  gaussJordanMethod 
          

          write(*,*) "start matrix inversion"   


        Nt=n*m 
 
         gaussJordanMethod = -cl_tab 
          r=0 ! Initialization of the pivot row 
       
        do j=1,Nt ! loop of all column   
           max_ = 0.0

           !step one, find the largest value of the row j  

           do i=r+1,Nt ! loop of the row j to find the max 
              
               if ( abs(matrix(i,j)) > abs(max_) ) then 
               max_ = matrix(i,j) 
               k=i
               end if
           end do 
         
           ! step two, pivoting 1
           ! find a row with the largest pivoting element      

           if ( max_ /= 0) then 
           r=r+1 
           matrix(k,:) = matrix(k,:)/max_ 
            gaussJordanMethod(k) =  gaussJordanMethod(k)/max_

            do p=1,Nt  
            tmp = matrix(k,p)
            matrix(k,p) = matrix(r,p)
            matrix(r,p) = tmp
            end do 
          
            tmp2 =  gaussJordanMethod(k)
             gaussJordanMethod(k) = gaussJordanMethod(r)
             gaussJordanMethod(r) =tmp2
            
           do i=1,Nt 
             if (i /= r) then
             pivot = matrix(i,j) 
             matrix(i,:) = matrix(i,:) - pivot*matrix(r,:)
              gaussJordanMethod(i) =  gaussJordanMethod(i) -  pivot* gaussJordanMethod(r)

            end if  
           end do 
         end if 
       end do 

end function 
end module  
