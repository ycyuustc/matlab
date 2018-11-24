function [B,U,DB]=fn_prepare_onesite_lr_his(A)
%UNTITLED5 此处显示有关此函数的摘要
%   此处显示详细说明

[D1,D2,d]=size(A);

A=permute(A,[3,1,2]);
A=reshape(A,[D1*d,D2]);

[B,S,U]=svd2(A);
DB=size(S,1);

B=reshape(B,[d,D1,DB]);
B=permute(B,[2,3,1]);

U=S*U;

end

