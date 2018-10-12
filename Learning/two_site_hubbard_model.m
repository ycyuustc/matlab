m2T=zeros(4,4);
m4U=zeros(4,4,4,4);
% m4M=zeros(4,4,4,4);

m2T(1,2)=1;
m2T(2,1)=1;
m2T(3,4)=1;
m2T(4,3)=1;

syms x;
syms y;
m2M=[cos(x),sin(x),0,0;-sin(x),cos(x),0,0;...
    0,0,cos(y),sin(y);0,0,-sin(y),cos(y)];
m2Mt=m2M.';

m4M=fn_contract(m2M,3,3,reshape(m2M,1,4,4),3,1);
m4M=permute(m4M,[1,3,2,4]);
m4Mt=fn_contract(m2Mt,3,3,reshape(m2Mt,1,4,4),3,1);
m4Mt=permute(m4Mt,[1,3,2,4]);

m4U(1,3,1,3)=1;
m4U(1,3,3,1)=-1;
m4U(3,1,1,3)=-1;
m4U(3,1,3,1)=1;
m4U(2,4,2,4)=1;
m4U(2,4,4,2)=-1;
m4U(4,2,2,4)=-1;
m4U(4,2,4,2)=1;

m4U=1/4*m4U;


t=fn_contract(m2M,2,2,m2T,2,1);
m2T_r=fn_contract(t,2,2,m2Mt,2,1);
m2T_r=simplify(m2T_r);
 
t=fn_contract(m4M,4,[3,4],m4U,4,[1,2]);
t=fn_contract(t,4,[3,4],m4Mt,4,[1,2]);
m4U_r=simplify(t);

syms t;
syms U;
tm=m4U_r(:,1,:,1)+m4U_r(:,3,:,3);
tm2=reshape(tm,4,4);
tm2=simplify(tm2);
result=-t*m2T_r+U*4*tm2;
disp(result);

energy=-t*(m2T_r(1,1)+m2T_r(3,3))...
    +U*(m4U_r(1,3,1,3)+m4U_r(3,1,3,1)...
    -m4U_r(1,3,3,1)-m4U_r(3,1,1,3));

energy=simplify(energy);

% t=t(:,1,:,1);
% t=reshape(t,2,2);
% t=simplify(t);
% 
% disp(simplify(m2T_r(2,1)));
% disp(simplify(m4U_r(2,1,1,1)));
% 
% syms lambda;
% result=m2T_r+lambda*t;
% disp(result);