clear all
close all
path='Z:\archivedwl-528\Second report\partie2\Exercise 4\';
N=[1:9,linspace(10,50,5)];
tclas=zeros(1,14);
tcycl=zeros(1,14);
numofrot=zeros(1,14);
for i=1:14  
A=rmg(rand(1,N(i)));

tic
eig_jclas(A);
tclas(i)=toc;

tic
[~,numofrot(i)]=eig_j(A);
tcycl(i)=toc;

end
figure('visible','on');
hold on
plot(N,tclas,'g')
plot(N,tcycl,'r')
xlabel('size of matrix $\mathrm{[\,]}$','Interpreter','latex','FontSize',18);
ylabel('$t\mathrm{[s]}$','Interpreter','latex','FontSize',18);
legend('classical Jacobi','cyclic Joacobi','Location','northwest')
set(gca,'FontSize',14)


hold off
saveas(gcf,strcat(path,'comparecyclandclass'),'png')

figure('visible','on');
hold on
plot(N,numofrot,'r')
xlabel('size of matrix $\mathrm{[\,]}$','Interpreter','latex','FontSize',18);
ylabel('number of rotation $\mathrm{[\,]}$','Interpreter','latex','FontSize',18);
set(gca,'FontSize',14)

hold off
saveas(gcf,strcat(path,'numofrot'),'png')
