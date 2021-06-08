n = 10; % the size of random nxn Hermitian matrix, need to be an integer;
times = 70; % this is max iteration.
err=10^(-5);
%% the following is to generate random matrix A and vector b then we can solve Ax=b;

diff = zeros(n,times);

A = zeros(n, n); % preallocate matrix for computational efficiency;
q = zeros(n, n);
b = zeros(n,1);
x = zeros(n,1);
%generate upper triangular portion of random nxn Hermitian matrix
for i = 1:n 
%     for j = 1:i 
%         A(i,j) = 10*rand; %complex(rand, rand); % create a random complex matrix
%     end
    b(i,1) = 10*rand; %complex(rand, rand);
    x(i,1) = 10*rand; %complex(rand, rand);
end
% A = ctranspose(A) + A; % thus created a symmetric n by n matrix A 
% generate a SPD A
A = randn(n);
A = A'*A;
A = A + 0.01*eye(n);
detA=det(A)

[x,r]=SteepestDescent(A,b,x,times,err)

r(size(r))
figure(1)
hold
plot(r)
[xc,rc]=ConjugateGradient(A,b,x,n,err)
figure(2)
hold
plot(rc)

realx= A \ b
xc



