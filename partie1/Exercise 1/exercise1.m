%% point 1)
clear all
close all
path="Z:\archivedwl-528\Second report\partie1\Exercise 1\";
%--------------------------- 1)
A1=[1,1,1,1;2,3,0,-1;-3,4,1,2;1,2,-1,1]; %1)
[L1,U1,P1]=lu_decomposition(A1);
b=[13;-1;10;1];
n=length(A1);
x=zeros(1,n);
y=zeros(1,n);

%solving the system:-------------------------------------------
b_=P1*b;
%------------forward substitution
y(1)=1/L1(1,1)*b_(1);
for i=2:n
    y(i)=1/L1(i,i)*(b_(i)-sum(L1(i,1:i-1).*y(1:i-1)));
end
%------------backward substitution
x(n)=1/U1(n,n)*y(n);
for i=1:n-1
    i=n-i;
    x(i)=1/U1(i,i)*(y(i)-sum(U1(i,i+1:n).*x(i+1:n)));
end
%solution:
x=x'
%verification:
xmat=A1\b
%----------------------------------------------------------


%%
clear all
close all
%------------------------ point 2)
range=150;
for m=1:range
    for j= 1:10
        A2=magic(m);
        [L2mat,U2mat,P2mat]=lu(A2);
        [L2,U2]=lu_withoutpiv(P2mat*A2);
        tempL(j)=abs(max(max(L2mat-L2)));
        tempU(j)=abs(max(max(U2mat-U2)));
    end
    diffL(m)=mean(tempL);
    diffU(m)=mean(tempU);
end
%take the maximum difference over all the indices of U between U from
%"lu_withoutpiv" and matlab's "lu" and take the mean of this over 10 different
%matrices of same size

figure
hold on
plot(1:range,diffL,'*--b')
xlabel('matrix size $\mathrm{[\,]}$','Interpreter','latex','FontSize',18);
ylabel('max component difference $\mathrm{[\,]}$','Interpreter','latex','FontSize',18);
set(gca,'FontSize',14)
saveas(gcf,'comparison_lu','png');
%----------------------------


%%
clear all
close all
%----------------------3)
A3=[0,-7,0;-3,2,6;5,-1,5]; %3)
[L3without,U3without]=lu_withoutpiv(A3)
[L3mat,U3mat,P3mat]=lu(A3)
[L3with,U3with]=lu_decomposition(P3mat*A3)

