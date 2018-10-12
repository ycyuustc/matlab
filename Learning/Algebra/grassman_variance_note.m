%% The function fn_c represents the function fn_contact
fn_f=@(j,m) sqrt((j+m).*(j-m+1));
fn_c=@(x1,x2,x3,x4,x5,x6) fn_contract(x1,x2,x3,x4,x5,x6);

%%
J=1/2;

Num_J=2*J+1;

v_J=(Num_J-1)/2:-1:-(Num_J-1)/2;

% z operator
mJ_z=diag(v_J);

% plus minus x and y operator
v_m=J:-1:-J+1;
mJ_p=diag(fn_f(J,v_m),1);
mJ_m=diag(fn_f(J,v_m),-1);
mJ_x=(mJ_p+mJ_m)/2;
mJ_y=(mJ_p-mJ_m)/(2j);

% fprintf('%s\n','mJ_p=');
disp('mJ_p=');disp(mJ_p);

% identity operator
mJ_id=eye(Num_J);

%  delta matrix
m_dtx=mJ_x*2;
m_dty=mJ_y*2;
m_dtz=mJ_z*2;
m_dti=mJ_id;

%%
% next we try to generate Six gamma matrix
Num=Num_J; % dimension of the spin vector space

T_gamma1=fn_c(m_dtx,3,3,reshape(m_dti,[1,Num,Num]),3,1);
f1=T_gamma1;
T_gamma1=fn_c(T_gamma1,5,5,reshape(m_dti,[1,Num,Num]),3,1);

T_gamma2=fn_c(m_dtz,3,3,reshape(m_dtx,[1,Num,Num]),3,1);
f2=T_gamma2;
T_gamma2=fn_c(T_gamma2,5,5,reshape(m_dti,[1,Num,Num]),3,1);

T_gamma3=fn_c(m_dtz,3,3,reshape(m_dtz,[1,Num,Num]),3,1);
T_gamma3=fn_c(T_gamma3,5,5,reshape(m_dtx,[1,Num,Num]),3,1);

T_gamma4=fn_c(m_dty,3,3,reshape(m_dti,[1,Num,Num]),3,1);
f3=T_gamma4;
T_gamma4=fn_c(T_gamma4,5,5,reshape(m_dti,[1,Num,Num]),3,1);

T_gamma5=fn_c(m_dtz,3,3,reshape(m_dty,[1,Num,Num]),3,1);
f4=T_gamma5;
T_gamma5=fn_c(T_gamma5,5,5,reshape(m_dti,[1,Num,Num]),3,1);

T_gamma6=fn_c(m_dtz,3,3,reshape(m_dtz,[1,Num,Num]),3,1);
T_gamma6=fn_c(T_gamma6,5,5,reshape(m_dty,[1,Num,Num]),3,1);                              

T_gamma1=permute(T_gamma1,[1,3,5,2,4,6]);
m_gamma1=reshape(T_gamma1,Num^3,Num^3);
T_gamma2=permute(T_gamma2,[1,3,5,2,4,6]);
m_gamma2=reshape(T_gamma2,Num^3,Num^3);
T_gamma3=permute(T_gamma3,[1,3,5,2,4,6]);
m_gamma3=reshape(T_gamma3,Num^3,Num^3);
T_gamma4=permute(T_gamma4,[1,3,5,2,4,6]);
m_gamma4=reshape(T_gamma4,Num^3,Num^3);
T_gamma5=permute(T_gamma5,[1,3,5,2,4,6]);
m_gamma5=reshape(T_gamma5,Num^3,Num^3);
T_gamma6=permute(T_gamma6,[1,3,5,2,4,6]);
m_gamma6=reshape(T_gamma6,Num^3,Num^3);

f1=permute(f1,[1,3,2,4]);
f2=permute(f2,[1,3,2,4]);
f3=permute(f3,[1,3,2,4]);
f4=permute(f4,[1,3,2,4]);

f1=reshape(f1,Num^2,Num^2);
f2=reshape(f2,Num^2,Num^2);
f3=reshape(f3,Num^2,Num^2);
f4=reshape(f4,Num^2,Num^2);

x1=1/2*(f1+1j*f3);
p1=1/2*(f1-1j*f3);
x2=1/2*(f2+1j*f4);
p2=1/2*(f2-1j*f4);

% create vacuum
v_v=[1;0];
T_v=fn_c(v_v,2,2,reshape(v_v,[1,Num]),2,1);
T_v=fn_c(T_v,3,3,reshape(v_v,[1,Num]),2,1);
T_op_V=fn_c(T_v,4,4,reshape(T_v,[1,Num,Num,Num]),Num^3,1);

v_vacuum=reshape(T_v,[Num^3,1]);
m_vacuum=reshape(T_op_V,[Num^3,Num^3]);

g1=m_gamma1;
g2=m_gamma2;
g3=m_gamma3;
g4=m_gamma4;
g5=m_gamma5;
g6=m_gamma6;

a1=1/2*(g1+1j*g4);
ad1=1/2*(g1-1j*g4);

a2=1/2*(g2+1j*g5);
ad2=1/2*(g2-1j*g5);

a3=1/2*(g3+1j*g6);
ad3=1/2*(g3-1j*g6);