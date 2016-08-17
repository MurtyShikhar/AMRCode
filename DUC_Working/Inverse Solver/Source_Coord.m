function [x,y] = Source_Coord(i, diameter, M);

%Calculates x-y coordinate of the ith incident source on the measurement
%boundary
theta_step=(2*pi)/M;
%%%%%%-----Finding out the polar coordinates of the source----%%%%%% 
r=diameter/2;
theta=(i-1)*theta_step;

%receiver coordinates start at -pi and move anticlockwise
%%%%%%-----Transforming the polar coordinate into Cartesian----%%%%%%
x=r*cos(theta);
y=r*sin(theta);