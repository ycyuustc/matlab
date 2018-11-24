function [E,mps]=fn_minimizeEP(hset,D,precision,mpsB)

[M,N]=size(hset);
d=size(hset{1,1},1);
mps=fn_createrandommps(N,D,d);
mps=fn_prepare_mps(mps);

Hstorage=fn_init_Hstorage(mps,hset);
Cstorage=fn_initCstorage(mps,[],mpsB,N);

Evalues=zeros(1,2*N-2);

while 1
    
   ind_E=1;
    
   for j=1:(N-1)
       
       B=mpsB{j};
       Cleft=Cstorage{j};
       Cright=Cstorage{j+1};
       P=fn_calcprojector_onesite(B,Cleft,Cright);
       
       Hleft=Hstorage(:,j);
       Hright=Hstorage(:,j+1);
       hsetj=hset(:,j);
       
       [A,E]=fn_minimizeE_onesiteP(hsetj,Hleft,Hright,P);
       [A,U]=fn_prepare_onesite_lr(A);
       
       mps{j}=A;
       Evalues(ind_E)=E;
       
       disp(E);
       
       ind_E=ind_E+1;
       
       for m=1:M
           h=reshape(hset{m,j},[1,1,d,d]);
           Hstorage{m,j+1}=fn_updateCleft(Hleft{m},A,h,A);
       end
       Cstorage{j+1}=fn_updateCleft(Cleft,A,[],B);
   end
   
   for j=N:(-1):2
       
       B=mpsB{j};
       Cleft=Cstorage{j};
       Cright=Cstorage{j+1};
       P=fn_calcprojector_onesite(B,Cleft,Cright);
       
       Hleft=Hstorage(:,j);
       Hright=Hstorage(:,j+1);
       hsetj=hset(:,j);
       
       [A,E]=fn_minimizeE_onesiteP(hsetj,Hleft,Hright,P);
       [A,U]=fn_prepare_onesite_rl(A);
       
       mps{j}=A;
       Evalues(ind_E)=E;
       ind_E=ind_E+1;
       
       for m=1:M
           h=reshape(hset{m,j},[1,1,d,d]);
           Hstorage{m,j}=fn_updateCright(Hright{m},A,h,A);
       end
       Cstorage{j}=fn_updateCright(Cright,A,[],B);
   end
   
   disp(reshape(Evalues,[1,2*N-2]));
%    pause(0.2);
   
   if std(Evalues)/abs(mean(Evalues))<precision
       mps{1}=fn_contract(mps{1},3,2,U,2,1);
       mps{1}=permute(mps{1},[1,3,2]);
       break;
   end
       
    
end

end