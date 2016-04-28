function x=covection_diffusion2(n,u)
disp('upwind differencing scheme');
a=[];
% u=0.1;
q=1/n;
d=1;%density
c=0.1;%Diffusion coefficient
D=c/q;
F=d*u;
fia=1;
fib=0;
a=zeros(n,n);
b=zeros(n,1);
aW=zeros(n,1);
% aE=D-F/2;
aE=zeros(n,1);
% Su1=(2*D+F)*fia;
% Su2=(2*D-F)*fib;
Su=zeros(n,1);
Sp=zeros(n,1);
% Sp1=-(2*D+F);
% Sp2=-(2*D-F);
% ap=aW+aE-Sp;
ap=zeros(n,1);
for r=1:n
    if r~=1;
        aW(r,1)=D+F;
    end
    if r~=n
        aE(r,1)=D;
    end
    if r==1
        Su(r,1)=(2*D+F)*fia;
        Sp(r,1)=-(2*D+F);
    end
    if r==n
        Su(r,1)=2*D*fib;
        Sp(r,1)=-2*D;
    end
end
ap=aW+aE-Sp;
for r=1:n
    for c=1:n
        if r==c
            a(r,c)=ap(r,1);
        end
    end
end
for r=1:n-1
    a(r,r+1)=-aE(r,1);
    a(r+1,r)=-aW(r+1,1);
end
b=Su;
x=inv(a)*b;
