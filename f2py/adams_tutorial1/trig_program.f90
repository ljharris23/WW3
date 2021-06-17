program some_trig  

     use trig_subroutines_md
     
     implicit none
     
     !-------------------------------------------------------------!
       real                :: area, areatot
       real                :: radius = 10
       real, dimension(3)  :: radii
     !-------------------------------------------------------------!
     
     call compute_area(area, radius)
     
     print * , "A circle of radius", radius, "has area", area
     
     radii(1)=1.0
     radii(2)=2.0
     radii(3)=3.0
     
     call compute_total_area(areatot, radii)
     
     print * , "The sum of the areas of circles with radii", radii, "is", areatot
     
end program some_trig
   
   
   
   
   
