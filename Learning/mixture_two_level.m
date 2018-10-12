m2T=zeros(2,2);
m4U=zeros(2,2,2,2);

m2T(1,1)=1;
m2T(2,2)=2;

syms x;
m2M=[cos(x),sin(x);-sin(x),cos(x)];
m2Mt=m2M.';

m4M=fn_contract(m2M,3,3,reshape(m2M,1,2,2),3,1);
m4M=permute(m4M,[1,3,2,4]);
m4Mt=fn_contract(m2Mt,3,3,reshape(m2Mt,1,2,2),3,1);
m4Mt=permute(m4Mt,[1,3,2,4]);

m4U(1,1,1,1)=1;
m4U(2,2,2,2)=1;
m4U(1,2,1,2)=1;
m4U(1,2,2,1)=-1;
m4U(2,1,1,2)=-1;
m4U(2,1,2,1)=1;

t=fn_contract(m2M,2,2,m2T,2,1);
m2T_r=fn_contract(t,2,2,m2Mt,2,1);

t=fn_contract(m4M,4,[3,4],m4U,4,[1,2]);
t=fn_contract(t,4,[3,4],m4Mt,4,[1,2]);
m4U_r=t;
t=t(:,1,:,1);
t=reshape(t,2,2);
t=simplify(t);

disp(simplify(m2T_r(2,1)));
disp(simplify(m4U_r(2,1,1,1)));

syms lambda;
result=m2T_r+lambda*t;
disp(result);