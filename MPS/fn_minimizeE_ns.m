function [E,mpsA,mpsB]=fn_minimizeE_ns(hset,D,precision,mps0)

[M,N]=size(hset);
% d=size(hset{1,1},1);
% d=size(hset{1,1},3);
% mpsA=fn_createrandommps(N,D,d);
% mpsB=fn_createrandommps(N,D,d);
mpsA = mps0;
mpsB = mps0;
[mpsA,mpsB,cs_R]=fn_prepare_mps_ns(mpsA,mpsB);
cs_L = cs_R;
cs_L{1} = 1;

Hstorage=fn_init_Hstorage_ns(mpsA,mpsB,hset);

Evalues=zeros(1,2*N-2);

while 1
    
   ind_E=1;
    
   for j=1:(N-1)
       
       Hleft=Hstorage(:,j);
       Hright=Hstorage(:,j+1);
       hsetj=hset(:,j);
       
       [A,B,E]=fn_minimizeE_onesite_ns(hsetj,Hleft,Hright,cs_L{j},cs_R{j});
       [A,B,A_out,B_out,AB_out]=fn_prepare_onesite_lr_ns(A,B,cs_L{j});
       cs_L{j+1} = AB_out;
       
       mpsA{j}=A;
       mpsB{j}=B;
       Evalues(ind_E)=E;
       
%        disp(E);
       
       ind_E=ind_E+1;
       
       for m=1:M
%            h=reshape(hset{m,j},[1,1,d,d]);
           h = hset{m,j};
           Hstorage{m,j+1}=fn_updateCleft(Hleft{m},B,h,A);
       end
   end
   
   for j=N:(-1):2
       
       Hleft=Hstorage(:,j);
       Hright=Hstorage(:,j+1);
       hsetj=hset(:,j);
       
       [A,B,E]=fn_minimizeE_onesite_ns(hsetj,Hleft,Hright,cs_L{j},cs_R{j});
       [A,B,A_out,B_out,AB_out]=fn_prepare_onesite_rl_ns(A,B,cs_R{j});
       cs_R{j-1} = AB_out;
       
       mpsA{j}=A;
       mpsB{j}=B;
       
       Evalues(ind_E)=E;
       ind_E=ind_E+1;
       
       for m=1:M
%            h=reshape(hset{m,j},[1,1,d,d]);
           h = hset{m,j};
           Hstorage{m,j}=fn_updateCright(Hright{m},B,h,A);
       end
   end
   
   disp(reshape(Evalues,[2*N-2,1]));
%    figure(9);
%    plot(real(Evalues),'-r*');
   
   if std(Evalues)/abs(mean(Evalues))<precision
       mpsA{1}=fn_contract(mpsA{1},3,2,A_out,2,1);
       mpsA{1}=permute(mpsA{1},[1,3,2]);
       mpsB{1}=fn_contract(mpsB{1},3,2,B_out,2,1);
       mpsB{1}=permute(mpsB{1},[1,3,2]);
       break;
   end
       
    
end

end