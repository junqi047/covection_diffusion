
clc
disp('Give the input to solve the set of equations AX=B')
m=length(a);
%z is a two dimensional array in which row corresponds to values of X in a
%specific iteration and the column corresponds to values of specific
%element of X in different iterations
c=0;%random assignment
e=1;%'e' represents the maximum error
d=0;%random assignment
for u=1:m
    x(u)=b(u,1)/a(u,u);
    z(1,u)=0;%initializing the values for matrix X(x1;x2;...xm)
end
l=2;%'l' represents the iteration no.
%loop for finding the convergence factor (C.F)
for r = 1:m
    for s = 1:m
        if r~=s
           p(r)=abs(a(r,s)/a(r,r))+d;%p(r) is the C.F for equation no. r
           d=p(r);
        end
    end
    d=0;
end
if min(p)>=1 %at least one equation must satisfy the condition p<1
   fprintf('Roots will not converge for this set of equations')
else
    while(e>=1e-4)
        j1=1;%while calculating elements in first column we consider only the old values of X
        for i1=2:m
            q(j1)=(a(j1,i1)/a(j1,j1))*z(l-1,i1)+c;
            c=q(j1);
        end
        c=0;
        z(l,j1)=x(j1)-q(j1);%elements of z in the iteration no. l
        x(j1)=z(l,j1);
        for u=1:m
            x(u)=b(u,1)/a(u,u);
            z(1,u)=0;
        end
        for j1=2:m-1%for intermediate columns between 1 and m, we use the updated values of X 
            for j1=1:i1-1
                q(j1)=(a(i1,j1)/a(i1,i1))*z(l,i1)+c;
                c=q(j1);
            end
            for i1=j1+1:m
                q(j1)=(a(j1,i1)/a(j1,j1))*z(l-1,i1)+c;
                c=q(j1);
            end
            c=0;
            z(l,j1)=x(j1)-q(j1);
            x(j1)=z(l,j1);
            for u=1:m
                x(u)=b(u,1)/a(u,u);
                z(1,u)=0;
            end
        end
        j1=m;%for the last column, we use only the updated values of X
        for i1=1:m-1
            q(j1)=(a(j1,i1)/a(j1,j1))*z(l,i1)+c;
            c=q(j1);
        end
        c=0;
        z(l,j1)=x(j1)-q(j1);
        for v=1:m
            t=abs(z(l,v)-z(l-1,v));%calculates the error 
        end 
        e=max(t);%evaluates the maximum error out of errors of all elements of X 
        l=l+1;%iteration no. gets updated
        for i=1:m 
            X(1,i)=z(l-1,i);%the final solution X 
        end
    end
    %loop to show iteration number along with the values of z 
    for i=1:l-1    
        for j=1:m        
            w(i,j+1)=z(i,j);    
        end
        w(i,1)=i;
    end
    disp('   It. no.      x1        x2       x3        x4 ') 
    disp(w) 
    disp('The final solution is ') 
    disp(X) 
    fprintf('The total number of iterations is %d',l-1)
end
