function [B,U,DB]=fn_prepare_onesite_lr(A)
%UNTITLED5 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

[D1,D2,d]=size(A);

A=permute(A,[1,3,2]);
A=reshape(A,[D1*d,D2]);

[B,S,U]=svd2(A);
DB=size(S,1);

B=reshape(B,[D1,d,DB]);
B=permute(B,[1,3,2]);

U=S*U;

end

