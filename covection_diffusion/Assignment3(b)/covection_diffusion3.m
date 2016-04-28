function x=covection_diffusion3(n,u)
disp('QUICK scheme');
a=[];
% u=0.2;
q=1/n;
d=1;%density
c=0.1;%Diffusion coefficient
D=c/q;
F=d*u;
fia=1;
fib=0;
% Da=D;
% Db=Da;
% De=D;
% Dw=De;
a=zeros(n,n);
b=zeros(n,1);
aW=ones(n,1)*(D+7/8*F);
aWW=ones(n,1)*(-1/8*F);
% aE=D-F/2;
aE=ones(n,1)*(D-3/8*F);
% Su1=(2*D+F)*fia;
% Su2=(2*D-F)*fib;
Su=zeros(n,1);
Sp=zeros(n,1);
% Sp1=-(2*D+F);
% Sp2=-(2*D-F);
% ap=aW+aE-Sp;
ap=zeros(n,1);
for r=1:n
    if r==1
        aW(r,1)=0;
        aWW(r,1)=0;
        aE(r,1)=D+1/3*D-3/8*F;
        Su(r,1)=(8/3*D+2/8*F+F)*fia;
        Sp(r,1)=-(8/3*D+2/8*F+F);
    end
     if r==2
        aWW(r,1)=0;
        aW(r,1)=D+7/8*F+1/8*F;
        Su(r,1)=-1/4*F*fia;
        Sp(r,1)=1/4*F;
     end
    if r==n
        aW(r,1)=D+1/3*D+6/8*F;
        aE(r,1)=0;
        Su(r,1)=(8/3*D-F)*fib;
        Sp(r,1)=-(8/3*D-F);
    end
end
ap=aW+aE+aWW-Sp;
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
for r=1:n-2
    a(r+2,r)=-aWW(r+2,1);
end
b=Su;
% aW
% aE
% aWW
% Su
% Sp
% ap
x=inv(a)*b;