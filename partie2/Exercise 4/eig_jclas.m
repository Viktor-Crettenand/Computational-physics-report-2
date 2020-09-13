function val = eig_jclas(A)
off=10;
tol=1e-5;
A_=(ones(size(A))-eye(size(A))).*A;
off=sqrt(sum(A_(:).^2));

while off>tol
    %------------------
A_=(ones(size(A))-eye(size(A))).*A;
[line,row]=find(abs(A_)==max(max(abs(A_(:)))),1);
maxA_=max(abs(A_(:)));
tau=(A(row,row)-A(line,line))/(2*A(line,row));
if tau>0 || tau==0
    t=-tau+sqrt(1+tau^2);
else
    t=-tau-sqrt(1+tau^2);
end
c=1/sqrt(1+t^2);
s=t*c;
    %-----------------
% G=eye(size(A));
% G(line,line)=c;
% G(row,row)=c;
% G(row,line)=-s;
% G(line,row)=s;
Aor=A; %original A
%A__=G'*A*G;
%-------------- performing multiplication A=G'*A*G
for i=1:length(A)
red1(i)=c*A(line,i)+-s*A(row,i); %line
red2(i)=s*A(line,i)+c*A(row,i); %row
end
for j=1:length(A)
    if j==line
    A(line,j)=red1(line)*c+red1(row)*-s;
    A(row,j)=red2(line)*c+red2(row)*-s;
    elseif j==row
    A(line,j)=red1(line)*s+red1(row)*c;
    A(row,j)=red2(line)*s+red2(row)*c;
    else
    A(line,j)=red1(j);
    A(row,j)=red2(j);
    end
end

for i=1:length(A)
    if i~=line && i~=row
    A(i,line)=Aor(i,line)*c+Aor(i,row)*-s;
    A(i,row)=Aor(i,line)*s+Aor(i,row)*c;
    end
end
%assert(isequal(A,A__));




%------------
A_=(ones(size(A))-eye(size(A))).*A;
off=sqrt(sum(A_(:).^2));
end

val = diag(A);
end