function [status, thetas, lambda, meas, tr_x, tr_y, bl_x, bl_y, fields ] = read_synthetic_from_fwd(configfile)
%Wrapper to parse config file: determine incident field angles and fields
%
%[s,th]=system('grep theta config_fem.txt |sed "s/theta/ /" |sed "s/=/ /"');
%concatanate strings to create the linux command required to parse file
%status should be 0 if all is okay with the file read
%
%**This file will not work in Windows. Need grep and sed at command line**

%Returns an array of incident field angles (in degrees) 'thetas'
%Set of incident angles
comm1 = ['grep theta ' configfile '|sed "s/theta/ /" |sed "s/=/ /"'];
[status1,th] = system(comm1);
thetas = str2num(th);

%lambda
comm2 = ['grep lambda ' configfile '|sed "s/lambda/ /" |sed "s/=/ /"'];
[status3,th] = system(comm2);
lambda = str2num(th);

%Returns a cmp matrix of fields, F_ij, where i is tx and j is rx.
%Measured fields in complex format
comm3 = ['grep outputfile ' configfile '|sed "s/outputfile/ /" |sed "s/=/ /"'];
[status2,df] = system(comm3);
field_raw = csvread(strtrim(df),1,2);
l = length(thetas);
fields = zeros(l,l);
for i = 1:l
    for j=1:l
        fields(i,j) = complex(field_raw(j+(l-1)*i,1),field_raw(j+(l-1)*i,2));
    end
end
fields = (fields)';

csvwrite('fields.dat', fields);
csvwrite('thetas.dat', thetas);
file = 'extents.csv';
%file = '/home/uday/Dropbox/research-uday-csi/Inverse Solver/extents.csv';

if exist(file, 'file')
    ex = csvread(file,1,0);
else
    error('could not find extents');
end

csvwrite('extents.dat', ex);

meas = ex(1);
tr_x = ex(2); tr_y = ex(3);
bl_x = ex(4); bl_y = ex(5);

status = status1+status2+status3;
end
