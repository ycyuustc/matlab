function [A,E]=fn_minimizeE_onesite(hsetj,Hleft,Hright)

DAl=size(Hleft{1},1);
DAr=size(Hright{1},1);
d=size(hsetj{1},3);
% d = 2;
M=size(hsetj,1);

Heff=0;

for m=1:M
    
    L = Hleft{m}; R = Hright{m}; X = hsetj{m};
%     Heffm2 = ncon({L,R,X},{[-1,1,-2],[-3,2,-4],[1,2,-5,-6]}); 
%     Heffm2=permute(Heffm2,[1,3,5,2,4,6]);
%     Heffm2=reshape(Heffm2,[DAl*DAr*d,DAl*DAr*d]);
    
    Heffm2 = fn_contract(L,3,2,X,4,1);
    Heffm2 = fn_contract(Heffm2,5,3,R,3,2);
    Heffm2 = permute(Heffm2,[1,5,3,2,6,4]);
    Heffm2=reshape(Heffm2,[DAl*DAr*d,DAl*DAr*d]);
    
%     Heffm=fn_contract(Hleft{m},3,2,Hright{m},3,2); 
%     Heffm=fn_contract(Heffm,5,5,hsetj{m},3,3);
%     Heffm=permute(Heffm,[1,3,5,2,4,6]);
%     Heffm=reshape(Heffm,[DAl*DAr*d,DAl*DAr*d]);
    
%     Heff=Heff+Heffm;
    Heff = Heff + Heffm2;
%     disp(max(max(abs(Heffm - Heffm2))));
    
end

% Heff=(Heff+Heff')/2; % 我加的，为了防止特征值趋于零时候的报错
options.disp=0;
[A,E]=eigs(Heff,1,'sr',options);  %  用sa代替了sr 

A=reshape(A,[DAl,DAr,d]);

end