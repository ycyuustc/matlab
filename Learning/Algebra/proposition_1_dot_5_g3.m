% 这段script用来用g3群的正则表示，验证propostion 1.5
% of GTM 129. 

[result,m_pers,m3]=canonical_presentation_2(3);

mg1=m3(:,:,1);
mg2=m3(:,:,2);
mg3=m3(:,:,3);
mg4=m3(:,:,4);
mg5=m3(:,:,5);
mg6=m3(:,:,6);

% step 1: definition of the innner product
fn_H0=@(v,w) sum(v.*w); %先定义内积

% step 2: construct H

% step3: find a subspace
% w1 and w2 are two vector in the subspace
[v,d]=eig(mg4);
w1=v(:,1);
w2=v(:,2);

tm=eye(6);
tm(:,1)=w1;
tm(:,2)=w2;

[mQ,mR]=qr(tm);

tv=rand()*mQ(:,3)+rand()*mQ(:,4);
% 我来验证这个垂直于原W空间的向量，在群作用之后依然是
% 垂直于原空间的。

tv2=mg6*tv;
disp(['tv2*w1=',num2str(sum(tv2.*w1))]);

temp=0;
for i=1:6
   
    tm=m3(:,:,i);
    
    temp=temp+sum((tm*tv2).*(tm*w1));
    
end

disp(['new tv2*w1=',num2str(temp)]);



