function[vec,val] = eig_rq(A,target)
vec=1/sqrt(length(A)).*ones(length(A),1);
val=vec'*A*vec;
err=1;
while err>1e-11
    tempvec=vec;
    tempval=val;
    vec=(A-eye(length(A))*tempval-eye(length(A))*target)\tempvec;
    vec=vec./norm(vec);
    val=vec'*A*vec;
    err=abs(tempval-val);
end

end


% function[vec,val] = eig_rq(A,target)
% vec=1/sqrt(length(A)).*ones(length(A),1);
% val=vec'*A*vec;
% vec=(A-eye(length(A))*val-eye(length(A))*target)\vec;
% tempval=val;
% val=vec'*A*vec; %first eigenvalue
% 
% 
% while abs(val-tempval)>1e-11
%     tempvec=vec;
%     tempval=val;
%     vec=(A-eye(length(A))*tempval-eye(length(A))*target)\tempvec;
%     vec=vec./norm(vec);
%     val=vec'*A*vec;
% end
% end