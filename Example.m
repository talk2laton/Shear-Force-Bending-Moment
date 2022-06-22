%% Shear Force & Bending Moment Examples

%     This program calculates the shear force and bending moment profiles, draw 
%     the free body, shear force and bending moment diagrams of the problem.
% 
%     Under the free body diagram, the equations of each section is clearly 
%     written with Latex
% 
%% How to call the function
%     To use this program, you call the function placing the arguments in cells
%     with keywords at the beginning of each cell except for the first 2 arguments.
% 
%%     First Argument
%     The first argument is the name of the problem as a string e.g.: 'PROB 1'.
% 
%%     Second Argument
%     -Simply supported beam
%     The second argument is a row vector containing length of the beam and 
%     location of the supports, for example, if the length of the beam is 20m and 
%     has 2 supports, one at 3m and the other at 17m, the second argument will 
%     thus be: [20, 3, 17]
% 
%%     -Cantilever
%     If the problem is a cantilever problem, then you have only one clamped 
%     support, at the beginning or end of the beam. In such a case, the number is
%     second argument contains 2 elements instead of three. For instance, fir a 
%     cantilever of length 20m, supported at the beginning, the second argument 
%     would be [20,0], and if supported at the end, we have [20,20].

%%     -Beam on the floor
%     Its possible to have a problem in which the body is lying on the floor 
%     without any point support. In such scenario, the second argument will just 
%     be the length of the beam 
% 
%%     Third argument and on
%     From the third argument and onward, we use cells. The first element of the 
%     cell contains a keyword describing what type of load is inside the argument.
%     The second element is the magnitude of the load while, the third element of
%      a cell argument is its location.
% 
%     Keywords: Point Load       = 'CF'
%               Moment           = 'M'
%               Distributed Load = 'DF'
% 
%     To add a downward point load of magnitude 5N at location 4m, the argument 
%     would be {'CF',-5,4}. Note the negative sign. If the force is acting upward
%     the argument would be {'CF',5,4};

%%     Examples
%% Moment(Torque)
%     To add a clockwise moment of magnitude 10N-m at location 14m, the argument 
%     would be {'M',-10,14}. Note the negative sign. If the moment is anticlockwise
%     the argument would be {'M',10,14};
%% Concentrated Load(Torque)
%     To add a downward force of magnitude 10N at location 14m, the argument 
%     would be {'CF',-10,14}. Note the negative sign. If the moment is
%     upward the argument would be {'CF',10,14};
%% Distributed Force
%     To add distributed load we need to describe all of them with the minimum 
%     number of point required to describe the profile with the highest 
%     complexity. For example,  {'DF',[5,5],[2,10]}, or {'DF',[1,4,5],[2,8,10]} 
%     There is no limit to the number degree of polynomial that can be used. 

%% Note
%     its is important that all concentrated loads and torques are listed in the 
%     order of locations 

%Problem Name
Name = 'Example';
% Length and Supports
LengthSupport = [20,5,20]; % length  = 20m, supports at 5m and 20m;
% Concetrated Loads
F1 = {'CF',-2,0}; F2 = {'CF',12,7};  % 2N downward at point 0
% Torques
T1 = {'M',10,8}; T2 = {'M',-10,12}; % ACW 10Nm at point 8m and CW 10Nm at point 12
% Distributed Loads
D1 = {'DF',5,[1,3]}; D2 = {'DF',-4,[14,17]}; % Constant 5N/m upwards from 1m to 3m and Constant 4N/m downwards from 14m to 17m
% Call the function
SFBM(Name,LengthSupport,F1,F2,T1,D1,T2,D2);