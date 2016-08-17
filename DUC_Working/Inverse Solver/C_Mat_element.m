function [ out ] = C_Mat_element(vec,mat) 	% it multiplies each element of the vector 'vec' with each row of the matrix 'mat', size of vec= mX1 and size of mat = mXn
out=mat; 			% size initialization
for i=1:size (mat,2)
    out(:,i)=vec.*mat(:,i);
end
end

