function plot_NE(n,x,u)
xplot=[0];
yplot=[1];
d=1;%density
c=0.1;%Diffusion coefficient
for k=1:n
%     
%     (2.7183-exp(1/n/2+(k-1)*1/n))/1.7183
%     d
%     u
%     c
%     exp(d*u*1/c)
%        (2.7183-exp(1/n/2+(k-1)*1/n))/1.7183
%     1-((exp(d*u*(1/n/2+(k-1)*1/n)/c)-1)/(exp(d*u*1/c)-1))
   fi(k,1)=((exp(d*u*1/c)-exp(d*u*(1/n/2+(k-1)*1/n)/c))/(exp(d*u*1/c)-1));
end
% x
% fi
Difference=-(x-fi);
for k=1:n
    error(k,1)=Difference(k,1)/fi(k,1)*100;
end
hold off
for k=1:n
    %     1/n/2+(k-1)*1/n
    %     fi(k,1)
    plot(1/n/2,fi(1,1),'sr','MarkerSize',10);
    plot((1/n/2+(k-1)*1/n), fi(k,1),'sr','MarkerSize',10);
    xplot=[xplot 1/n/2+(k-1)*1/n];
    yplot=[yplot x(k,1)];
    plot(xplot,yplot,'b','LineWidth',2);
    hold on
end
xplot=[xplot 1];
yplot=[yplot 0];
plot(xplot,yplot,'b','LineWidth',2);
% fi
end