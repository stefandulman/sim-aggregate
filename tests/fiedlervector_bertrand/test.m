% testing the algorithm of kempe for computing eigenvectors
clear
clc
close all

maxruns = 50;
n = 7;

% create a custom topology
cm = [ 0 1 0 1 0 0 0;
       1 0 1 1 0 0 0;
       0 1 0 0 0 0 0;
       1 1 0 0 1 1 1;
       0 0 0 1 0 1 1;
       0 0 0 1 1 0 1;
       0 0 0 1 1 1 0];
   
% create laplacian
lap = diag(sum(cm,1)) - cm;

% setup alfa
alfa = max(sum(cm));

% compute m
m = eye(n) - lap / alfa;

x = rand(n,1);
for i=1:maxruns,
  v = m*x; % compute v
  v = v / norm(v); % normalize
  x = v - mean(v);
end

x

[vtmp dtmp] = eig(lap);
vtmp
diag(dtmp)'
