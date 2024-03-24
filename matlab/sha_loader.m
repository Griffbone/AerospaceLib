%--------------------------------------------------------------------------
% This function loads a set of spherical harmonics
% Input: Name of spherical harmonics text file 
% Ouput: CS (matrix of C and S), C, S, n_max, m_max
%--------------------------------------------------------------------------

function [CS, n_max, m_max] = sha_loader(filename, n_max);

%--------------------------------------------------------------------------
% Select source file with set of spherical harmonics
fid = fopen(filename,'rt');

line = fgets(fid);

array = textscan(fid,'%6f,%6f,%f,%f,%f,%f', -1, 'headerlines', 0);
n=array{1,1};
m=array{1,2};
C_row=array{1,3};
S_row=array{1,4};

for z=1:length(n)
    C_norm(n(z)+1,m(z)+1) = C_row(z);
    S_norm(n(z)+1,m(z)+1) = S_row(z);
end

m_max = n_max;

%--------------------------------------------------------------------------
% Denormalizing C and S for further use

for n = 0:n_max
    for m = 0:n
        if m == 0
            k = 1;
        else
            k = 2;
        end
        C(n+1,m+1) = C_norm(n+1,m+1) / sqrt(factorial(n+m)/(k*(2*n+1)*factorial(n-m)));
        S(n+1,m+1) = S_norm(n+1,m+1) / sqrt(factorial(n+m)/(k*(2*n+1)*factorial(n-m)));
    end
end

C(1,1) = 1;
           
CS = C + [S(:,2:end) S(:,1)]';  % Combining C and S in 1 matrix
