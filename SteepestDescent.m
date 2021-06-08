function [x,r] = SteepestDescent(A,b,x,times,err)
n = size(A,2);
r=zeros(n,times);
r(:,1)=b-A*x;
% delta=r'*r
delta1=r'*r;

for t=1:times
    t;
    if abs(r)<err
        break
    end
    
    alpha=(r(:,t)'*r(:,t))/(r(:,t)'*A*r(:,t));
    x = x + alpha * r(:,t);
    r(:,t+1)=r(:,t)-alpha*A*r(:,t);
%     delta=r'*r
end

r=abs(r);
r=mean(r,1);

