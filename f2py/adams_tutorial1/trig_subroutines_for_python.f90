!module trig_subroutines_md
!
!contains
 
   !A subroutine to compute the area of a circle, 
   !from given radius      
   subroutine compute_area(area, radius)
   	real, intent(in)   :: radius
   	real, intent(out)  :: area
   	real, parameter    :: pi = 4.*atan(1.)
           area = pi*radius**2
   end subroutine compute_area
   
   !A subroutine to compute the total area of many circles, 
   !from given array of radii
   !(calls previous subroutine)
   subroutine compute_total_area(totarea, radii)
   	real, dimension(:), intent(in)    :: radii
   	real, intent(out)                 :: totarea
        real                              :: area
        totarea = 0.0
        do i = 1, size(radii)
                  call compute_area(area, radii(i))
                  totarea = totarea + area
        end do                  
   end subroutine compute_total_area 

!end module trig_subroutines_md
