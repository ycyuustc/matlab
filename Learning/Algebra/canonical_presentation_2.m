function [result,m_pers_order,m3_g]=canonical_presentation_2(Num)

%Num=3;

m_per=perms(1:1:Num);
tm=m_per;


for k=1:Num-1
   
    tv=sign(-m_per(:,(k+1):Num)+m_per(:,k)*ones(1,Num-k));
    tv=1/2*tv+1/2;
    
    tm(:,k)=sum(tv,2);
    
end

tm(:,Num)=zeros(factorial(Num),1);
v_factorial=factorial((Num-(1:1:Num)));
m_factorial=ones(factorial(Num),1)*v_factorial;

v_sort=sum(m_factorial.*tm,2);
v_sort=v_sort+1;


[~,vb]=sort(v_sort);
m_pers=[m_per,v_sort];
m_pers_order=m_pers(vb,:);

Nump=factorial(Num);

m3_g=zeros(Nump,Nump,Nump);

m_matrix=zeros(Nump,Nump);

for iL=1:Nump
    for iR=1:Nump
        
        vL=m_pers_order(iL,1:1:Num);
        vR=m_pers_order(iR,1:1:Num);
        v_prod=vR(vL);
        m_matrix(iL,iR)=fn_order(v_prod);
        
    end
end

m_id=eye(Nump);

for ig=1:Nump
   
    tv=m_matrix(ig,:);
    m3_g(:,:,ig)=m_id(:,tv);
           
end

result=m_matrix;

end


function result=fn_order(v)

len=length(v);
v_order=ones(1,len);

for k=1:len-1
   
    tv=sign(-v(:,(k+1):len)+v(:,k));
    tv=1/2*tv+1/2;
    
    v_order(k)=sum(tv);
    
end

v_order(len)=0;

v_factorial=factorial((len-(1:1:len)));

result=sum(v_factorial.*v_order)+1;
    
end

