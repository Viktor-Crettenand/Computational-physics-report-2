function [L,U,P]=lu_decomposition(A)
n=length(A);
L=eye(n);
P=eye(n);
U=A;
for k=1:n-1
    %---------------------------------------exchanges the kth line with the indexth line
    [~,index]=max(abs(U(k:n,k)));
    index=index+k-1;
    tempL=L(index,1:k-1);
    U([index k],:)=U([k index],:);
    P([index k],:)=P([k index],:);
    if k>1
    L(index,1:k-1)=L(k,1:k-1);
    L(k,1:k-1)=tempL;
    end
    %------------------------------------------    
    for i=k+1:n
        L(i,k)=L(i,k)+U(i,k)/U(k,k);
        for j=k:n
            U(i,j)=U(i,j)-L(i,k)*U(k,j);
        end
    end
end
end











% function [L,U]=lu_decomposition(A)
% n=length(A);
% L=eye(n);
% U=A;
% for k=1:n-1
%     for i=k+1:n
%         L(i,k)=U(i,k)/U(k,k);
%         for j=k:n
%             U(i,j)=U(i,j)-L(i,k)*U(k,j);
%         end
%     end
% end
% end

% function [L,U]=lu_decomposition(A)
% n=length(A);
% L=eye(n);
% U=zeros(n);
% Atemp=cell(1,n-1);
% Atemp{1}=A;
% for k=1:n-1
%     for i=k+1:n
%         L(i,k)=Atemp{k}(i,k)/Atemp{k}(k,k);
%         for j=k:n
%             Atemp{k+1}(i,j)=Atemp{k}(i,j)-L(i,k)*Atemp{k}(k,j);
%         end
%     end
% end
% 
% for i=1:n
%     for j=i:n
%         U(i,j)=Atemp{i}(i,j);
%     end
% end
% end