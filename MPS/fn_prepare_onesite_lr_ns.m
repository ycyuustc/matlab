function [UA,UB,A_out,B_out,AB_out]=fn_prepare_onesite_lr_ns(A,B,AB_lr)

[D1,D2,d]=size(A);
A = permute(A,[1,3,2]);
A=reshape(A,[D1*d,D2]);

[UA,S,V]=svd2(A);
DA=size(S,1);

UA=reshape(UA,[D1,d,DA]);
UA=permute(UA,[1,3,2]);
A_out=S*V;

[D1,D2,d]=size(B);
B = permute(B,[1,3,2]);
B=reshape(B,[D1*d,D2]);

[UB,S,V]=svd2(B);
DB=size(S,1);

UB=reshape(UB,[D1,d,DB]);
UB=permute(UB,[1,3,2]);
B_out=S*V;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tensor = fn_contract(UA,3,3,conj(UB),3,3);
tensor = fn_contract(tensor,4,[1,3],AB_lr,2,[1,2]);

[u,s,v] = svd2(tensor);
v = v';

UA = fn_contract(UA,3,2,transpose(u'),2,1);
UA = permute(UA,[1,3,2]);
UB = fn_contract(UB,3,2,transpose(v'),2,1);
UB = permute(UB,[1,3,2]);

A_out = transpose(u) * A_out;
B_out = transpose(v) * B_out;

AB_out = s;

end

