function[vec,val] = eig_power(A)
vec=1/sqrt(length(A)).*ones(length(A),1);
val=vec'*A*vec;
err=1;
while err>1e-11
    tempvec=vec;
    tempval=val;
    vec=A*tempvec;
    vec=vec./norm(vec);
    val=vec'*A*vec;
    err=abs(val-tempval);
end

end