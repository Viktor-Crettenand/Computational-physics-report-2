function H = ham(c)
%------------------constants----------------------------------------------
N=100;
r=1e-9;
a=5*r;
eV=1.60217656535e-19;
hbar=1.05457172647e-34;
mel=9.10938291e-31;
V0=-1.5*eV;
dx=a/(N+1);
%-------------------------------------------------------------------------

%-------------------------------------V(x)--------------------------------
V=zeros(N^2);
x=@(i) mod(i,N)*dx-dx*N/2; % gives x coordinate with vector index i (assumes (0,0) coordinate is in the midle)
y=@(i) ceil(i/N)*dx-dx*N/2; %y coordinate from i
for i=1:N^2 
    if (y(i)/c)^2+x(i)^2<r^2
        V(i,i)=V0;
    end
end
%-----------------test
% test=(diag(V)<0);
% for i=1:N
%     test_(i,1:N)=test((i-1)*N+1:i*N);
% end
% pcolor(test_)
%---------------------
%-------------------------------------------------------------------------

% ------------------------delta-------------------------------------------
delta=zeros(N^2);
for i=1:N^2
    delta(i,i)=-4;
    if i+1<N^2
        delta(i,i+1)=1;
    end
    if i-1>0
        delta(i,i-1)=1;
    end
    if i+N<N^2
        delta(i,i+N)=1;
    end
    if i-N>0
        delta(i,i-N)=1;
    end
end
delta=delta/dx^2;
%--------------------------------------------------------------------------

%------------construction of H (hamiltonian)-------------------------------
H=sparse(-hbar^2/2/mel*delta + V);
%--------------------------------------------------------------------------

end

