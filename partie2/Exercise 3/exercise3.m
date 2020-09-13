close all
clear all
path='Z:\archivedwl-528\Second report\partie2\Exercise 3\';
a=0.4; %nm
gamma= -2.7; %eV
for N=2:7
f=figure('visible','off')
xlim([-2.5,3.5])
ylim([-3.5,2.5])
hold on
for m=1:2*N
    if mod(m-1,4)<2
        pos{N,m}=[mod(m-1,2)*(a/sqrt(3)),-fix((m-1)/4)*a];
    else
        pos{N,m}=[-a/2/sqrt(3)+mod(m-1,2)*2/sqrt(3)*a, -a/2-a*fix((m-1)/4)];
    end
    plot(pos{N,m}(1),pos{N,m}(2),'*')
    plot(pos{N,m}(1)+3/sqrt(3)*a,pos{N,m}(2),'*')
    plot(pos{N,m}(1)-3/sqrt(3)*a,pos{N,m}(2),'*')
end 
saveas(f,strcat(path,'atompositionN=',num2str(N)),'png');
end

dist(pos{7,1},pos{7,2});
nearest=dist(pos{7,1},pos{7,3});

for N=2:7
for j=1:2*N
for i=1:2*N
    if (dist(pos{N,i},pos{N,j})<nearest+0.01 && dist(pos{N,i},pos{N,j})>0.01)
        H{N}{i,j}= @(k) gamma;
    elseif dist((pos{N,i}+[3/sqrt(3)*a,0]),pos{N,j})<nearest+0.01
        H{N}{i,j}= @(k) gamma*exp(-1i*k*a);  
    elseif dist((pos{N,i}-[3/sqrt(3)*a,0]),pos{N,j})<nearest+0.01
        H{N}{i,j}= @(k) gamma*exp(+1i*k*a);
    else
        H{N}{i,j}=@(k) 0;
    end
end
end
end
%-------------------------------- cellfunction discovery-------
% k=1;
% resultsCell =cellfun(@(x) x(k),H{2}) %'UniformOutput',true)        @(x) x(k) is evaluated at x being the function in each cell
%--------------------------------------------------------------

numbpoint=50;
ks=linspace(-pi/a,pi/a,numbpoint);
for N=2:7
    f=figure('visible','off')
for k=1:numbpoint
    E{N}(k,:)=eig(cellfun(@(x) x(ks(k)),H{N}));
end
for j=1:size(E{N},2)
plot(ks,E{N}(:,j))
hold on
end
xlabel('$k\mathrm{[nm^{-1}]}$','Interpreter','latex','FontSize',18);
ylabel('$E\mathrm{[eV]}$','Interpreter','latex','FontSize',18);
set(gca,'FontSize',14)
hold off
saveas(f,strcat(path,'bandstructN=',num2str(N)),'png');
Eg(N)=min(min(E{N}(:,N+1:2*N)))-max(max(E{N}(:,1:N)));
end

f=figure('visible','off')
plot(2:N,Eg(2:N))
%title('Imaginary part of FID as a function of time','FontSize',10)
xlabel('$N\mathrm{[\,]}$','Interpreter','latex','FontSize',18);
ylabel('$E_g\mathrm{[eV]}$','Interpreter','latex','FontSize',18);
set(gca,'FontSize',14)
saveas(f,strcat(path,'energygap'),'png');

