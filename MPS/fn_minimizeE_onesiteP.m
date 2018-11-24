function [A,E]=fn_minimizeE_onesiteP(hsetj,Hleft,Hright,P)

DAl=size(Hleft{1},1);
DAr=size(Hright{1},1);
d=size(hsetj{1},1);
M=size(hsetj,1);

Heff=0;

for m=1:M
    
    Heffm=fn_contract(Hleft{m},3,2,Hright{m},3,2);
    Heffm=fn_contract(Heffm,5,5,hsetj{m},3,3);
    Heffm=permute(Heffm,[1,3,5,2,4,6]);
    Heffm=reshape(Heffm,[DAl*DAr*d,DAl*DAr*d]);
    Heff=Heff+Heffm;
    
end

Heff=P'*Heff*P;

% Heff=(Heff+Heff')/2; % �Ҽӵģ�Ϊ�˷�ֹ����ֵ������ʱ��ı���
options.disp=0;
[A,E]=eigs(Heff,1,'sr',options);  %  ��sa������sr 

A=P*A;
A=reshape(A,[DAl,DAr,d]);

end