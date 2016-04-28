function x=Gauss_seidel(a,b)
er=0.75;
iterationNo=0;
L=length(a);
n = length(b);
x = zeros(n,1);
e = ones(n,1);
sum1=0;
sum2=0;
% Check if the matrix A is diagonally dominant
for i = 1:n
    j = 1:n;
    j(i) = [];
    Check(i) = abs(a(i,i)) - sum(abs(a(i,j))); % Is the diagonal value greater than the remaining row values combined?
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
    
%     while max(e) > 0.001
%         iterationNo = iterationNo + 1;
%         z = x;  % save current values to calculate error later
%         for i = 1:L
%             j = 1:L; % define an array of the coefficients' elements
%             j(i) = [];  % eliminate the unknow's coefficient from the remaining coefficients
%             temp = x;  % copy the unknows to a new variable
%             temp(i) = [];  % eliminate the unknown under question from the set of values
%                                x(i) = (b(i) - sum(a(i,j) * temp) )/ a(i,i);
%                 %x(i) = (b(i) - sum(a(i,j) )* temp) / a(i,i);
%                 %                  X(i) = (C(i) - sum(A(i,j) * Xtemp)) / A(i,i);
%         end
%         S(:,iterationNo) = x;
%             e = sqrt((x - z).^2);
%     end
while max(e) > 0.0001
    iterationNo = iterationNo + 1;
    z = x;  % save current values to calculate error later
    for i = 1:n
       for j = 1:i-1; % define an array of the coefficients' elements
%         j(i) = [];  % eliminate the unknow's coefficient from the remaining coefficients
        sum1=sum1-a(i,j)/a(i,i);
        T1(i,j) =sum1;
       end
    end
    for i = 1:n
        for j = i+1:n; % define an array of the coefficients' elements
%         j(i) = [];  % eliminate the unknow's coefficient from the remaining coefficients
        %temp = x;  % copy the unknows to a new variable
       % temp(i) = [];  % eliminate the unknown under question from the set of values
       sum2=sum2-a(i,j)/ a(i,i);
       T2(i,j) =  sum2;
        end
    end
    for i = 1:n
        j = 1:n; % define an array of the coefficients' elements
%         j(i) = [];  % eliminate the unknow's coefficient from the remaining coefficients
        %temp = x;  % copy the unknows to a new variable
        %temp(i) = [];  % eliminate the unknown under question from the set of values
        c(i,j) = b(i) / a(i,i);
     end
%     x = (1-er)*temp+er*(c(i)-T1*temp-T2*temp) ;
    temp=x;
    T1
    T2
    c
    x = (1-er)*temp+er*(c-T1*x-T2*temp) ;
    s(:,iterationNo) = x;
    e = sqrt((z - x).^2);
end

