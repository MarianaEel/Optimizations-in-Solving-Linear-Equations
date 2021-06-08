function [x,r]=ConjugateGradient(A,b,x,times,err)

n = size(A,1);
r = zeros(n,times);
r(:,1) = b - A*x;
d = r(:,1);
t = 1;
delta = r(:,1)'*r(:,1);
while t<times 
    if delta<err
        break
    end
    alpha = delta/(d'*A*d);
    x = x + alpha*d;
%     r(:,t+1) = r(:,t) - alpha*A*d;
    r(:,t+1) = b-A*x;
    beta = r(:,t+1)'*r(:,t+1)/delta;
    delta = r(:,t+1)'*r(:,t+1);
    d = r(:,t+1) + beta*d;
    t = t+1;
end
r=abs(r);
r=mean(r,1);

