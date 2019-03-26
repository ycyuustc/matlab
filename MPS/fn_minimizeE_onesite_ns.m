function [A1,A2,E]=fn_minimizeE_onesite_ns(hsetj,Hleft,Hright,dL,dR)

dL0 = dL;
dR0 = dR;
dL = diag(sqrt(1./diag(dL)));
dR = diag(sqrt(1./diag(dR)));
DAl=size(Hleft{1},1);
DAr=size(Hright{1},1);
% d=size(hsetj{1},1);
d = size(hsetj{1},3);
M=size(hsetj,1);

Heff=0;

for m=1:M
    
    L = Hleft{m}; R = Hright{m}; X = hsetj{m};

    L = fn_contract(L,3,1,dL,2,1);
    L = permute(L,[3,1,2]);
    L = fn_contract(L,3,3,dL,2,1);
    
    R = fn_contract(R,3,1,dR,2,1);
    R = permute(R,[3,1,2]);
    R = fn_contract(R,3,3,dR,2,1);
    
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
[A1,E]=eigs(Heff,1,'lm',options);  %  用sa代替了sr 
[A2,~]=eigs(Heff',1,'lm',options);

% A1 = A1./sqrt(sum(A1.*conj(A1)));
% A2 = A2./sqrt(sum(A2.*conj(A2)));

% norm_factor1 = sum(A1.*conj(A2));
% disp(norm_factor1);

A1=reshape(A1,[DAl,DAr,d]);
A2=reshape(A2,[DAl,DAr,d]);

A1 = ncon({dL,A1,dR},{[-1,1],[1,2,-3],[2,-2]});
A2 = ncon({dL,A2,dR},{[-1,1],[1,2,-3],[2,-2]});

% norm_factor2 = ncon({A1,conj(A2),dL0,dR0},{[1,4,3],[2,5,3],[1,2],[4,5]});
% % E = E/abs(norm_factor2);
% disp(norm_factor2);
end