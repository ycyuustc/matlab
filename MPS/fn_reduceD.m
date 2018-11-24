function [mpsB,K]=fn_reduceD(mpsA,mpoX,DB,precision)

N=length(mpsA);
d=size(mpsA{1},3);
mpsB=fn_createrandommps(N,DB,d);
mpsB=fn_prepare_mps(mpsB);

Cstorage=fn_initCstorage(mpsB,mpoX,mpsA,N);
Kvalues=zeros(1,2*N-2);

while 1
   ind_K=1; 
   for j=1:(N-1)
       Cleft=Cstorage{j};
       Cright=Cstorage{j+1};
       A=mpsA{j};
       X=mpoX{j};
       [B,K]=fn_reduceD2_onesite(A,X,Cleft,Cright);
       [B,U]=fn_prepare_onesite_lr(B);
       mpsB{j}=B;
       Kvalues(ind_K)=K;
       ind_K=ind_K+1;
       Cstorage{j+1}=fn_updateCleft(Cleft,B,X,A);
   end
   
   for j=N:(-1):2
       
       Cleft=Cstorage{j};
       Cright=Cstorage{j+1};
       A=mpsA{j};
       X=mpoX{j};
       [B,K]=fn_reduceD2_onesite(A,X,Cleft,Cright);
       [B,U]=fn_prepare_onesite_rl(B);
       mpsB{j}=B;
       Kvalues(ind_K)=K;
       ind_K=ind_K+1;
       Cstorage{j}=fn_updateCright(Cright,B,X,A);
       
   end
   
   if std(Kvalues)/abs(mean(Kvalues))<precision
       mpsB{1}=fn_contract(mpsB{1},3,2,U,2,1);
       mpsB{1}=permute(mpsB{1},[1,3,2]);
       break;
   end
       
        
   
    
end

end