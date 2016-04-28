function main_box
    function x=Gauss_seidel(a,b)
        er=0.75;
        iterationNo=0;
        L=length(a);
        n = length(b);
        x = zeros(n,1);
        e = ones(n,1);
        d=0;
        % Check if the matrix A is diagonally dominant
        for i = 1:n
            j = 1:n;
            j(i) = [];
            b = abs(a(i,j));
            Check(i) = abs(a(i,i)) - sum(b); % Is the diagonal value greater than the remaining row values combined?
            if Check(i) < 0
                fprintf('The matrix is not strictly diagonally dominant at row %2i\n\n',i)
            end
        end
        
        
        %     while max(e) > 0.001
        %         iterationNo=iterationNo+1;
        %         z=x;
        %         for i=1:n
        %             j=1:n;
        % %                j(i) = [] ;
        %             temp=x;
        %                 temp(i)=[];
        % %                 if i>j
        % %                     sum1=sum1+(-a(i,j)/a(i,i));
        % %                 elseif i<j
        % %                     sum2=sum2+(-a(i,j)/a(i,i));
        % %                 else
        % %                     sum2=0;
        % %                     sum1=sum2;
        % %                 end
        % %                 x(i)=sum1*temp+b(i)/a(i,i)+sum2*temp;
        %                   x(i) = (b(i) - sum(a(i,j) * temp)) / a(i,i);
        %
        %             end
        %         end
        %         Xsolution(:,iterationNo) = x
        %         e = sqrt((x - z).^2);
        %     end
        
        while max(e) > 0.001
            iterationNo = iterationNo + 1;
            z = x;  % save current values to calculate error later
            for i = 1:L
                j = 1:L; % define an array of the coefficients' elements
                j(i) = [];  % eliminate the unknow's coefficient from the remaining coefficients
                temp = x;  % copy the unknows to a new variable
                temp(i) = [];  % eliminate the unknown under question from the set of values
                %                 x(i) = (b(i) - sum(a(i,j) * temp) )/ a(i,i);
                x(i,1) = (1-er)*temp+er*(b(i) - sum(a(i,j) )* temp) / a(i,i);
                %                  X(i) = (C(i) - sum(A(i,j) * Xtemp)) / A(i,i);
            end
            S(:,iterationNo) = x;
            e = sqrt((x - z).^2);
        end
        x = [1:iterationNo;S];
        
    end
    function plot_NE(n,x,u)
        %         analysisS=[1];
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
            xplot=[xplot 1/n/2+(k-1)*1/n];
            yplot=[yplot x(k,1)];
            %     1/n/2+(k-1)*1/n
            %     fi(k,1)
            plot(1/n/2,fi(1,1),'sr','MarkerSize',10);
            %             analysisS=[analysisS fi(k,1)];
            plot(1/n/2+(k-1)*1/n, fi(k,1),'sr','MarkerSize',10);
            hold on
            plot(xplot,yplot,'b','LineWidth',2);
        end
        xplot=[xplot 1];
        yplot=[yplot 0];
        plot(xplot,yplot,'b','LineWidth',2);
        % fi
        %         axis([0 1.1 0 1.1])
    end
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
        x=Gauss_seidel(a,b);
    end
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
        x=Gauss_seidel(a,b);
    end
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
        x=Gauss_seidel(a,b);
    end
clc
clear
number=0;
while true
    %     number=0;
    clc
    if number==4
        break
    end
    disp('Which scheme you want to use?');
    disp('1, central differencing scheme.');
    disp('2, upwind scheme.');
    disp('3, QUICK scheme.');
    disp('4, Quit.');
    number=input('Please enter the number: ');
    n=input('Enter the node: ');
    u=input('Enter the velocity(m/s): ');
    switch number
        case 1
            x=covection_diffusion(n,u);
        case 2
            x=covection_diffusion2(n,u);
        case 3
            x=covection_diffusion3(n,u);
        otherwise
            disp('Error');
    end
    if number==1||number==2 ||number==3
        plot_NE(n,x,u);
    end
end
end