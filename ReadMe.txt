This file explains how to use the SFBM.m code. as well as supporting 
accessories. SFBM means Shear Force & Bending Moment. 

This program calculates the shear force and bending moment profiles, draw 
the free body, shear force and bending moment diagrams of the problem.

Under the free body diagram, the equations of each section is clearly 
written with Latex

To use this program, you call the function placing the arguments in cells
with keywords at the beginning of each cell except for the first 2 arguments.

First Argument
The first argument is the name of the problem as a string e.g.: 'PROB 1'.

Second Argument
-Simply supported beam
The second argument is a row vector containing length of the beam and 
location of the supports, for example, if the length of the beam is 20m and 
has 2 supports, one at 3m and the other at 17m, the second argument will 
thus be: [20, 3, 17]

-Cantilever
If the problem is a cantilever problem, then you have only one clamped 
support, at the beginning or end of the beam. In such a case, the number is
second argument contains 2 elements instead of three. For instance, fir a 
cantilever of length 20m, supported at the beginning, the second argument 
would be [20,0], and if supported at the end, we have [20,20].

-Beam on the floor
Its possible to have a problem in which the body is lying on the floor 
without any point support. In such scenario, the second argument will just 
be the length of the beam 

Third argument and on
From the third argument and onward, we use cells. The first element of the 
cell contains a keyword describing what type of load is inside the argument.
The second element is the magnitude of the load while, the third element of
 a cell argument is its location.

Keywords: Point Load       = 'CF'
          Moment           = 'M'
          Distributed Load = 'DF'

To add a downward point load of magnitude 5N at location 4m, the argument 
would be {'CF',-5,4}. Note the negative sign. If the force is acting upward
the argument would be {'CF',5,4};

To add a clockwise moment of magnitude 10N-m at location 14m, the argument 
would be {'M',-10,14}. Note the negative sign. If the moment is anticlockwise
the argument would be {'M',10,14};

To add distributed load we need to describe all of them with the minimum 
number of point required to describe the profile with the highest 
complexity. For example, a linear profile can be described as {'DF',[5,5],[2,10]}
meaning uniform force per unit length of 5N/m from point 2m to 10m. If the 
values of the profile were given at 3 points, the code will automatically 
assume it to be quadratic. If profile is uniform, the coefficient of the 
second and first degrees would be zero.Hence describing the constant 5N/m 
from 2m to 10m as {'DF',5,[2,10]}, {'DF',[5,5],[2,10]}, {'DF',[5,5,5],[2,8,10]}
will make no difference. But in case where the values in the force vector 
are different, SFBM will generate a polynomial fit for the forces as a 
function of position. For instance {'DF',[1,5,5],[2,8,10]} will generate a 
quadratic function, while {'DF',[1,4,5],[2,8,10]} will generate a linear 
expression and {'DF',[5,5,5],[2,8,10]}will generate a degree zero expression

There is no limit to the number degree of polynomial that can be used. 

its is important that all concentrated loads and torques are listed in the 
order of locations 

For example:
SFBM('Prob 200',[20,5,20],{'CF',-2,0},{'M',10,8},{'DF',[5,8,6],[1,3,5]},{'M',-10,12},{'DF',-4,[14,17]})

Name :            Prob 200
Length:           20m
Supports:         5m and 20m
Point Load:       2N at point 0
Distributed Load: Constant 5N/m from 1m to 3m and Constant -4N/m from 14m to 17m
Moment:           ACW 10Nm at point 8m and CW 10Nm at point 12


Solution:
The Freebody (with Legend), Shear Force and Bending moment diagrams are 
generated and saved in picture format titles Prob 200.png.
