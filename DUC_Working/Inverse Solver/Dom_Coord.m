%%%%%%%%%------As of now, this function divides the square object domain into N^2
%%%%%%%%%------squares. The inputs are N and the dimension of the square object
%%%%%%%%%------domain (in meters)---------%%%%%%%%%%


function [coord_matrix] = Dom_Coord(side,N)

%Given i from 1 to N^2, gives x-y coordinates of the center of the
%corresponding ith square cell

grid_dim_x=(side/N); %dimension of the cell along x
grid_dim_y=(side/N); %dimension of the cell along y
coord_matrix = zeros(N^2, 2);

for count_1=1:N %row starting from the bottom. so corresponding coordinate in coord_matrix is (count_1-1)*N + count_2
	coordinate = (count_1-1)*N + 1;   
    coord_matrix(coordinate,1)=(grid_dim_x/2) - side/2;
    coord_matrix(coordinate,2)=grid_dim_y/2+(count_1-1)*grid_dim_y - side/2;
    for count_2=2:N
    	curr_coord = coordinate + count_2 -1;
        coord_matrix(curr_coord, 1)= coord_matrix(coordinate,1) + (count_2-1)*grid_dim_x;      
        coord_matrix(curr_coord,2)= coord_matrix(coordinate,2);
    end
end

%%%%----calculating the position of the cell----%%%%



