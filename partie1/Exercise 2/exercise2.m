clear all
close all
alpha=[0;6];
beta=[9;0];
gamma=[3;5];
delta=[3;4];
epsilon=[10;12];
zeta=[7;0];
A=[alpha/12, beta/12, [0;0], [0;0], [0;0], [0;0]; [0;0], [0;0], gamma/14, [0;0], epsilon/14, [0;0]; alpha/32, [0;0], [0;0], delta/32, [0;0], zeta/32];
b=[6;2;4;6;25/8;15/4];
% A2=[1,1,0,0,0,0;0,0,1,0,1,0;1,0,0,1,0,1];
% b2=[12;14;32];

[L1,U1,P1]=lu_decomposition(A);
n=length(A);
b_=P1*b;
%------------forward
y(1)=1/L1(1,1)*b_(1);
for i=2:n
    y(i)=1/L1(i,i)*(b_(i)-sum(L1(i,1:i-1).*y(1:i-1)));
end
%-----------------

%------------backward
x(n)=1/U1(n,n)*y(n);
for i=1:n-1
    i=n-i;
    x(i)=1/U1(i,i)*(y(i)-sum(U1(i,i+1:n).*x(i+1:n)));
end
x


%verification:
%x=A\b %alpha, beta, gamma,  delta, epsilon, zeta
earthmass=5.972e24;
x*earthmass

% barycenter-----------------
lightminute=17987547.5; %in kilometer
%R in kilometer
R=lightminute/sum(x)*[(x(1)*alpha(1)+x(2)*beta(1)+x(3)*gamma(1)+x(4)*delta(1)+x(5)*epsilon(1)+x(6)*zeta(1));(x(1)*alpha(2)+x(2)*beta(2)+x(3)*gamma(2)+x(4)*delta(2)+x(5)*epsilon(2)+x(6)*zeta(2))];



