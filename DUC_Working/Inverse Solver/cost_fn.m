function[cost] = cost_fn(w,chi)

global eta_S; 
global eta_D; 
global aux_GD; % N^2 x N^2
global aux_GS;
global f;
global M;
global u_inc;

total_data_error = data_error(w)
total_obj_error = object_error(w,chi)
cost = eta_S*total_data_error + eta_D*total_obj_error;