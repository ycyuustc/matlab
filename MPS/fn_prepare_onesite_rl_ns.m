function [VA,VB,A_out,B_out,AB_out]=fn_prepare_onesite_rl_ns(A,B,AB_rl)

[D1,D2,d]=size(A);
A=reshape(A,[D1,d*D2]);

[U,S,VA]=svd2(A);
DA=size(S,1);

VA=reshape(VA,[DA,D2,d]);
A_out=U*S;

[D1,D2,d]=size(B);
B=reshape(B,[D1,d*D2]);

[U,S,VB]=svd2(B);
DB=size(S,1);

VB=reshape(VB,[DB,D2,d]);
B_out=U*S;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tensor = fn_contract(VA,3,3,conj(VB),3,3);
tensor = fn_contract(tensor,4,[2,4],AB_rl,2,[1,2]);

[u,s,v] = svd2(tensor);
v = v';

VA = fn_contract(u',2,2,VA,3,1);
VB = fn_contract(v',2,2,VB,3,1);

A_out = A_out * u;
B_out = B_out * v;

AB_out = s;

end

