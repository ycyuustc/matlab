function [E,mps]=fn_minimizeE(hset,D,precision)

[M,N]=size(hset);
d=size(hset{1,1},2);

mps=fn_createrandommps(N,D,d);
mps=fn_prepare_mps(mps);

Hstorage=fn_init_Hstorage(mps,hset);

Evalues=zeros(1,2*N-2);

while 1
    
   ind_E=1;
    
   for j=1:(N-1)
       
       Hleft=Hstorage(:,j);
       Hright=Hstorage(:,j+1);
       hsetj=hset(:,j);
       
       [A,E]=fn_minimizeE_onesite(hsetj,Hleft,Hright);
       [A,U]=fn_prepare_onesite_lr(A);
       
       mps{j}=A;
       Evalues(ind_E)=E;
%        
%        disp(E);
%        
       ind_E=ind_E+1;
       
       for m=1:M
%            h=reshape(hset{m,j},[1,1,d,d]);
           h = hset{m,j};
           Hstorage{m,j+1}=fn_updateCleft(Hleft{m},A,h,A);
       end
   end
   
   for j=N:(-1):2
       
       Hleft=Hstorage(:,j);
       Hright=Hstorage(:,j+1);
       hsetj=hset(:,j);
       
       [A,E]=fn_minimizeE_onesite(hsetj,Hleft,Hright);
       [A,U]=fn_prepare_onesite_rl(A);
       
       mps{j}=A;
       Evalues(ind_E)=E;
       ind_E=ind_E+1;
       
       for m=1:M
%            h=reshape(hset{m,j},[1,1,d,d]);
           h = hset{m,j};
           Hstorage{m,j}=fn_updateCright(Hright{m},A,h,A);
       end
   end
   
   disp(reshape(Evalues,[2*N-2,1]));
%    pause(0.2);
   
   if std(Evalues)/abs(mean(Evalues))<precision
       mps{1}=fn_contract(mps{1},3,2,U,2,1);
       mps{1}=permute(mps{1},[1,3,2]);
       break;
   end
       
    
end

end