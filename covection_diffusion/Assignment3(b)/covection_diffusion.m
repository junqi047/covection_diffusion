function x=covection_diffusion(n,u)
disp('central dierencing scheme');
% axis([0 2 0 1]);
% xplot=[0];
% yplot=[1];
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
        aW(r,1)=D+F/2;
    end
    if r~=n
        aE(r,1)=D-F/2;
    end
    if r==1
        Su(r,1)=(2*D+F)*fia;
        Sp(r,1)=-(2*D+F);
    end
    if r==n
        Su(r,1)=(2*D-F)*fib;
        Sp(r,1)=-(2*D-F);
    end
end
% aW
% aE
% Sp
% Su
ap=aW+aE-Sp;
for r=1:n
    for s=1:n
        if r==s
            a(r,s)=ap(r,1);
        end
    end
end
for r=1:n-1
    a(r,r+1)=-aE(r,1);
    a(r+1,r)=-aW(r+1,1);
end
b=Su;
% b
% a
x=inv(a)*b;
% for k=1:n
% %     
% %     (2.7183-exp(1/n/2+(k-1)*1/n))/1.7183
% %     d
% %     u
% %     c
% %     exp(d*u*1/c)
% %        (2.7183-exp(1/n/2+(k-1)*1/n))/1.7183
% %     1-((exp(d*u*(1/n/2+(k-1)*1/n)/c)-1)/(exp(d*u*1/c)-1))
%    fi(k,1)=((exp(d*u*1/c)-exp(d*u*(1/n/2+(k-1)*1/n)/c))/(exp(d*u*1/c)-1));
% end
% % x
% % fi
% Difference=-(x-fi);
% for k=1:n
%     error(k,1)=Difference(k,1)/fi(k,1)*100;
% end
% hold off
% for k=1:n
%     %     1/n/2+(k-1)*1/n
%     %     fi(k,1)
%     plot(1/n/2,fi(1,1),'sr','MarkerSize',10);
%     plot((1/n/2+(k-1)*1/n), fi(k,1),'sr','MarkerSize',10);
%     xplot=[xplot 1/n/2+(k-1)*1/n];
%     yplot=[yplot x(k,1)];
%     plot(xplot,yplot,'b','LineWidth',2);
%     hold on
% end
% xplot=[xplot 1];
% yplot=[yplot 0];
% plot(xplot,yplot,'b','LineWidth',2);
% % fi