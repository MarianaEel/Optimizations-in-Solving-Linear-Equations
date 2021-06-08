% Solving Linear Equation problems using Lanczos method by zhangcb@bu.edu
% A EK501 Project
% Feel free to copy this code to any use, I hope that this may help

clc, clear all

n = 10; % the size of random nxn Hermitian matrix, need to be an integer;
times = 1; % this is a over all loop to compute mean of errors.

%% the following is to generate random matrix A and vector b then we can solve Ax=b;

diff = zeros(n,times);
for t=1:times
A = zeros(n, n); % preallocate matrix for computational efficiency;
q = zeros(n, n);
b = zeros(n,1);

%generate upper triangular portion of random nxn Hermitian matrix
for i = 1:n 
    for j = 1:i 
        A(i,j) = complex(rand-0.5, rand-0.5); % create a random complex matrix
    end
    b(i,1) = complex(rand-0.5, rand-0.5);
end
A = ctranspose(A) + A; % thus created a symmetric n by n matrix A 
det(A)

q(:,1) = b/norm(b); % choose the normed b as the 1st basis of Krylov Subsbace
alpha = zeros(n,1); % preallocate matrix for computational efficiency;
beta = zeros(n,1); 
v = zeros(n,1); 
beta(1) = 0;

err=10^(-10); % a small number represent 0 when judging if beta=0

tic
j=1;
v(:,j)=A*q(:,j);
alpha(j) = q(:,j)'*A*q(:,j) ;
v(:,j) = v(:,j) - alpha(j)*q(:,j) ;% here v_j = r_j = beta_j+2 q_j+1
beta(j+1) = norm(v(:,j));
q(:,j+1) = v(:,j)/beta(j+1);

for j = 2:n
    j;
    
    alpha(j) = q(:,j)'*A*q(:,j); 
    v(:,j) = A*q(:,j) - beta(j)*q(:,j-1) - alpha(j)*q(:,j); % here v_j = r_j = beta_j+2 q_j+1
    beta(j+1) = norm(v(:,j));

%     v(:,j) = A*q(:,j) - beta(j)*q(:,j-1) % v_j = alpha_j q_j + beta_j+2 q_j+1
%     alpha(j) = q(:,j)'*A*q(:,j) 
% %     alpha(j) = v(:,j)'*q(:,j)
%     v(:,j) = v(:,j) - alpha(j)*q(:,j) % here v_j = r_j = beta_j+2 q_j+1
%     beta(j+1) = norm(v(:,j))
    

    
    if beta(j+1)< err % break when beta j=0
        break
    end
    q(:,j+1) = v(:,j)/beta(j+1);
        
end
% this is to construct the tri-diagonal matrix T
T = diag(alpha);
T = T+ diag(beta(2:length(beta)-1),1);
T = T+ diag(beta(2:length(beta)-1),-1);
T

be = norm(b);
e1 = zeros(length(T),1);
e1(1,1) = 1;

% take out the one-more q
[r,c] = size(q);
if c<=n
    q=zeros(n,n)+q
else
    q = q(:,1:c-1)    
end

x = q* pinv(T)*be*e1
toc

x_real=A\b
diff(:,t)=x_real-x
end
mdiff=mean(diff,"all")



