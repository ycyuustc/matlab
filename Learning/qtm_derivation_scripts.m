b = @(lambda) 1j./(lambda+1j);
c = @(lambda) lambda./(lambda+1j);

m4id=zeros(2,2,2,2);
m4id(1,1,1,1)=1;
m4id(1,2,1,2)=1;
m4id(2,1,2,1)=1;
m4id(2,2,2,2)=1;

m4p=zeros(2,2,2,2);
m4p(1,1,1,1)=1;
m4p(1,2,2,1)=1;
m4p(2,1,1,2)=1;
m4p(2,2,2,2)=1;
m4p2 = permute(m4p,[1,4,3,2]);

fn_R = @(lambda) b(lambda)*m4p+c(lambda)*m4id;
fn_RR = @(lambda) b(lambda)*m4p2+c(lambda)*m4id;
fn_Rh = @(lambda) b(lambda)*m4id+c(lambda)*m4p;

global Num;
Num = 8;
global v_tau;
v_tau = rand(1,Num);
lambda = 0.5;
global M;
M = 3;

%%
[~,op1] = fn_generate_QTM_mps(rand(),v_tau,Num);
[~,op2] = fn_generate_QTM_mps(rand(),v_tau,Num);

T1 = op1(:,:,1,1)+op1(:,:,2,2);
T2 = op2(:,:,1,1)+op2(:,:,2,2);
disp("==== check the communitation rule of t(lambda1) and t(lambda2)====");
disp(max(max(abs(T1*T2-T2*T1))));

%%
[~,op1] = fn_generate_QTM_mps(lambda,rand(1,Num),Num);
[mps2,op2] = fn_generate_QTM_mps(lambda,rand(1,Num),Num);

T1 = op1(:,:,1,1)+op1(:,:,2,2);
T2 = op2(:,:,1,1)+op2(:,:,2,2);

disp("====== the tau communiation rule does not hold ======");
disp(max(max(abs(T1*T2-T2*T1))));

%%
x1 = rand();
x2 = rand();
y = rand();

Rh = fn_Rh(x2-x1);
R_1 = fn_R(y-x1);
R_2 = fn_R(y-x2);
RR_1 = fn_RR(y+x1);
RR_2 = fn_RR(y+x2);

LHS = ncon({Rh,RR_1,RR_2},{[-1,-2,1,2],[1,-5,-3,3],[2,3,-4,-6]});
RHS = ncon({Rh,RR_1,RR_2},{[1,2,-3,-4],[-2,3,2,-6],[-1,-5,1,3]});
LHS = reshape(LHS,[1,numel(LHS)]);
RHS = reshape(RHS,[1,numel(RHS)]);
disp("======= check the YBE of rotated R matrix =======");
disp(max(max(abs(LHS-RHS))));

%%
LHS = ncon({Rh,R_1,R_2},{[-1,-2,1,2],[1,-5,-3,3],[2,3,-4,-6]});
RHS = ncon({Rh,R_1,R_2},{[1,2,-3,-4],[-2,3,2,-6],[-1,-5,1,3]});
LHS = reshape(LHS,[1,numel(LHS)]);
RHS = reshape(RHS,[1,numel(RHS)]);
disp("======= check the YBE of original R matrix =======");
disp(max(max(abs(LHS-RHS))));

%%
lambda_1 = rand();
lambda_2 = rand();
v_tau = rand(1,Num);

[~,T1] = fn_generate_QTM_mps(lambda_1,v_tau,Num);
[~,T2] = fn_generate_QTM_mps(lambda_2,v_tau,Num);
Rh = fn_Rh(lambda_1 - lambda_2);

LHS = ncon({Rh,T1,T2},{[-1,-2,1,2],[-5,3,1,-3],[3,-6,2,-4]});
RHS = ncon({Rh,T1,T2},{[1,2,-3,-4],[3,-6,-2,2],[-5,3,-1,1]});
LHS = reshape(LHS,[1,numel(LHS)]);
RHS = reshape(RHS,[1,numel(RHS)]);
disp("======= check the RTT relation =======");
disp(max(max(abs(LHS-RHS))));

%%
lambda = rand();
mu = rand();

LHS = A(lambda)*B(mu);
RHS = 1/c(mu-lambda)*B(mu)*A(lambda)...
    - b(mu-lambda)/c(mu-lambda)*B(lambda)*A(mu);
LHS = reshape(LHS,[1,numel(LHS)]);
RHS = reshape(RHS,[1,numel(RHS)]);
disp("======= check the communication rule for A B operator =======");
disp(max(max(abs(LHS-RHS))));

%%
LHS = D(lambda)*B(mu);
RHS = 1/c(-mu+lambda)*B(mu)*D(lambda)...
    - b(-mu+lambda)/c(-mu+lambda)*B(lambda)*D(mu);
LHS = reshape(LHS,[1,numel(LHS)]);
RHS = reshape(RHS,[1,numel(RHS)]);
disp("======= check the communication rule for D B operator =======");
disp(max(max(abs(LHS-RHS))));

%%
lambda = rand();
vacuum = fn_generate_vacuum(Num);
mA = A(lambda);
mD = D(lambda);
mv = reshape(vacuum,[numel(vacuum),1]);

EigA = 1;
for i=2:2:Num
    EigA = EigA*c(v_tau(i)-lambda);
end
EigD = 1;
for i=1:2:(Num-1)
    EigD = EigD*c(v_tau(i)+lambda);
end

LHS = mA*mv; RHS = EigA*mv;
disp("======= check the eigen function for A operator =======");
disp(max(max(abs(LHS-RHS))));
LHS = mD*mv; RHS = EigD*mv;
disp("======= check the eigen function for D operator =======");
disp(max(max(abs(LHS-RHS))));

%%
lambda = rand();
v_mu = rand(1,M);
tv = [lambda,v_mu];
cA = @(n) A(tv(n+1));
cB = @(n) B(tv(n+1));
cD = @(n) D(tv(n+1));
cmu = @(n) tv(n+1);

cf = @(m,n) 1/c(cmu(n)-cmu(m));
cg = @(m,n) -b(cmu(n)-cmu(m))/c(cmu(n)-cmu(m));

m=0;
n=3;
LHS = cA(m)*cB(n);
RHS = cf(m,n)*cB(n)*cA(m) + cg(m,n)*cB(m)*cA(n);
LHS = reshape(LHS,[1,numel(LHS)]);
RHS = reshape(RHS,[1,numel(RHS)]);
disp("======= check the communication rule for A B operator =======");
disp(max(max(abs(LHS-RHS))));

LHS = cD(m)*cB(n);
RHS = cf(n,m)*cB(n)*cD(m) + cg(n,m)*cB(m)*cD(n);
LHS = reshape(LHS,[1,numel(LHS)]);
RHS = reshape(RHS,[1,numel(RHS)]);
disp("======= check the communication rule for D B operator =======");
disp(max(max(abs(LHS-RHS))));

%%
LHS = cA(0);
for i=1:M
   LHS = LHS*cB(i); 
end
RHS = 0;
for l=0:M
   temp = cA(l);
   for k=0:M
      if k~=l
          temp = cB(k)*temp;
      end
   end
   RHS = RHS + fn_K(M,l,lambda,v_mu)*temp;
end
LHS = reshape(LHS,[1,numel(LHS)]);
RHS = reshape(RHS,[1,numel(RHS)]);
disp("======= check the IMPORTANT iteration relation for A B =======");
disp(max(max(abs(LHS-RHS))));

LHS = cD(0);
for i=1:M
   LHS = LHS*cB(i); 
end
RHS = 0;
for l=0:M
   temp = cD(l);
   for k=0:M
      if k~=l
          temp = cB(k)*temp;
      end
   end
   RHS = RHS + fn_K2(M,l,lambda,v_mu)*temp;
end
LHS = reshape(LHS,[1,numel(LHS)]);
RHS = reshape(RHS,[1,numel(RHS)]);
disp("======= check the IMPORTANT iteration relation for D B =======");
disp(max(max(abs(LHS-RHS))));

%% 
% prepare tau(lambda), v_tau
beta = 0.1;
M = 12;
Num = M/2;
v_tau = -1j*beta/M*ones(1,M);
epsilon = beta/M;

lambda = 0;

[mps1,op1] = fn_generate_QTM_mps(lambda,v_tau,M);
op_lambda = op1(:,:,1,1)+op1(:,:,2,2);
% from the script "qtm_solvver" or 
%"solving_QTM_BA_equation_Num_6_finite_temperature" 
% we obtain the BA roots. named the solution of the BA equation as 
% vector "v_ba6". 
v_mu = 1j*v_ba6;
B_op_cell = cell(1,Num);
for i=1:Num
    [mps1,op1] = fn_generate_QTM_mps(v_mu(i),v_tau,M);
    op_B = op1(:,:,1,2)+op1(:,:,1,2);
    B_op_cell{i} = op_B;
end

vacuum = fn_generate_vacuum(M);
vvacuum = reshape(vacuum,[numel(vacuum),1]);

a_lambda = ((lambda+1j*epsilon)/(lambda+1j*epsilon-1j))^(Num);
d_lambda = ((lambda-1j*epsilon)/(lambda-1j*epsilon+1j))^(Num);

tv = lambda-v_mu;

tv_a = prod(c(tv));
tv_d = prod(c(-tv));

max_eigen = a_lambda/tv_a + d_lambda/tv_d;
LHS = max_eigen;
RHS = eigs(op_lambda,1);
disp("======= check the spectrum =======");
disp(max(max(abs(LHS-RHS))));

temp = vvacuum;
for i=1:Num
    temp = B_op_cell{i}*temp;
    temp = temp/norm(temp);
end

LHS = op_lambda*temp;
RHS = max_eigen*temp;

disp("======= check the eigen function =======");
disp(max(max(abs(LHS-RHS))));

%===== END of MAIN SCRIPTS  ============================================
%%  The functions 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function res = A(lambda)
global v_tau Num;
[~,tm] = fn_generate_QTM_mps(lambda,v_tau,Num);
res = tm(:,:,1,1);
end

function res = B(lambda)
global v_tau Num;
[~,tm] = fn_generate_QTM_mps(lambda,v_tau,Num);
res = tm(:,:,1,2);
end

% function res = C(lambda)
% global v_tau Num;
% [~,tm] = fn_generate_QTM_mps(lambda,v_tau,Num);
% res = tm(:,:,2,1);
% end

function res = D(lambda)
global v_tau Num;
[~,tm] = fn_generate_QTM_mps(lambda,v_tau,Num);
res = tm(:,:,2,2);
end

function res = fn_generate_vacuum(Num)
v1 = [1;0];
v2 = [0;1];
c_tensor = cell(1,Num);
c_index = cell(1,Num);
for i=1:Num
   if mod(i,2)==1
      c_tensor{i} = v1; 
   else
      c_tensor{i} = v2; 
   end
   c_index{i} = -i;
end
res = ncon(c_tensor,c_index);
end

function res = fn_K(M,l,lambda,v_mu)
tv = [lambda,v_mu];
cmu = @(n) tv(n+1);

b = @(lambda) 1j/(lambda+1j);
c = @(lambda) lambda/(lambda+1j);
cf = @(m,n) 1/c(cmu(n)-cmu(m));
cg = @(m,n) -b(cmu(n)-cmu(m))/c(cmu(n)-cmu(m));

res = 1;
for i=1:M
    if i==l
        res = res*cg(0,l);
    else
        res = res*cf(l,i);
    end
end

end

function res = fn_K2(M,l,lambda,v_mu)
tv = [lambda,v_mu];
cmu = @(n) tv(n+1);

b = @(lambda) 1j/(lambda+1j);
c = @(lambda) lambda/(lambda+1j);
cf = @(m,n) 1/c(cmu(n)-cmu(m));
cg = @(m,n) -b(cmu(n)-cmu(m))/c(cmu(n)-cmu(m));

res = 1;
for i=1:M
    if i==l
        res = res*cg(l,0);
    else
        res = res*cf(i,l);
    end
end

end
