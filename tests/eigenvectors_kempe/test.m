
% testing the algorithm of kempe for computing eigenvectors
clear
clc
close all

% number of nodes
n = 7;

% create a custom topology
cm = [ 0 1 0 1 0 0 0;
       1 0 1 1 0 0 0;
       0 1 0 0 0 0 0;
       1 1 0 0 1 1 1;
       0 0 0 1 0 1 1;
       0 0 0 1 1 0 1;
       0 0 0 1 1 1 0];
   
% number of eigenvectors
k = 3;
% line 1
Q = rand(n, k);
%Q = ones(n,k);

% "weigthed" cm
a = cm + eye(n);
% number of loops
maxloops = 100;

for i=1:maxloops,

    % line 3
    V = a * Q;
    
    % line 4-5
    K = zeros(k,k);
    for j = 1:n,
        K = K + (V(j,:)'*V(j,:));
    end
    
    K
    
    % line 6
    R = chol(K);
    
    R
    %return
    
    % line 7
    Q = V * inv(R);
end

% comparison
[tv td] = eig(cm);
Q 
tv
diag(td)'
