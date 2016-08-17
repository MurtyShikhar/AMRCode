function[num] = dot_2(mat1, mat2)
aux = conj(mat2');
if (iscolumn(mat1))
    num = aux*mat1;
else
    num = mat1*aux;
end