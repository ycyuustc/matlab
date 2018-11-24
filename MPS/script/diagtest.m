function [Hami,base,energy,wf]=diagtest(N)

%N=4;

Num=2^N;

% v_zero=zeros(1,N);

base=zeros(Num,N);


for n=0:Num-1
   
    strn=dec2base(n,2);   
    base(n+1,N-length(strn)+1:N)=str2num(strn(:))';
         
end

Hami=zeros(Num,Num);

for i=1:Num
    for j=1:Num
        
        v1=base(i,:);
        v2=base(j,:);
        
        h=0;
        for m=1:N-1
            h=h+op_s(m,m+1,v1,v2);
        end
        
        Hami(i,j)=h;
        
    end
end

[v,d]=eig(Hami);
energy=diag(d);

wf=cell(1,Num);

for i=1:Num
   
    wf{i}=[v(:,i),base];
    
end

end



function [r]=op_s(i,j,v1,v2)

H2=[1/2,0,0,0;0,-1/2,1,0;0,1,-1/2,0;0,0,0,1/2];   
    
v_int_1=v1([i,j]);
v_int_2=v2([i,j]);

v1([i,j])=[];
v2([i,j])=[];

ind_1=sum(v_int_1.*[2,1])+1;
ind_2=sum(v_int_2.*[2,1])+1;

if max(abs(v1-v2))<1e-6 
    r=H2(ind_1,ind_2);
%     r=0;
%     if v_int_1==[1,1]&&v_int_2==[1,1]
%         r=1/2;
%     end
%     if v_int_1==[0,0]&&v_int_2==[0,0]
%         r=1/2;
%     end
%     if v_int_1==[1,0]&&v_int_2==[1,0]
%         r=-1/2;
%     end
%     if v_int_1==[0,1]&&v_int_2==[0,1]
%         r=-1/2;
%     end
%     if v_int_1==[1,0]&&v_int_2==[0,1]
%         r=1;
%     end 
%     if v_int_1==[0,1]&&v_int_2==[1,0]
%         r=1;
%     end   

else
    r=0;
end



end





