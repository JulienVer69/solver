module sparseStorage 


    

    contains 

    subroutine storage_csr3(matrix,Nt,a,ia,ja) 
    INTEGER  :: Nt
    INTEGER :: i,j
    INTEGER :: k=1
    REAL*8 :: a(:)
    INTEGER :: ia(:)
    INTEGER :: ja(:)
    REAL*8 :: matrix(:,:) 
    logical :: ind  
    
    do i=1,Nt
      ind = .TRUE. 
        do j=1,Nt
           if ( matrix(i,j) /= 0. ) then 
           a(k) = matrix(i,j) 
           ja(k) = j 
                if (ind) then 
                ia(i) = j
                ind = .FALSE.  
                endif 
           k=k+1        
           endif
        enddo
      enddo 
                ia(Nt+1) = k-1
        end subroutine 

end module 
