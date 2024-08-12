# greenland_flowlines
Makes flowlines from user input cross section and subset greenland velocity data from fastice. See example.png.
https://github.com/beynoldsj/greenland_flowlines/blob/main/example.png 

This script takes in subset velocity from the fastice MEaSUREs subsetting tool (https://github.com/fastice/GrIMPNotebooks) and creates flowlines the pass through points on a user-defined cross section. The flowlines are output to a shape file in the same format as the shape files from Felikson et al. (2020).

To find the margins of the glacier, the code can:
1) Find where velocity is a user-defined fraction of the max velocity.
2) Find where velocity is a user-defined fraction of the difference between the max and the minimum velocity (for each side of the glacier).
3) Find the maximum of the derivative of the velocity magnitude across the cross section, which will tend to land right in the center of the shear margins.
