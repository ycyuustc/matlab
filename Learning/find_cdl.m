m2id=eye(2);
m2z=[1,0;0,-1];

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

lambda = 0.5;

p1 = 1/(1-lambda);
p2 = lambda/(lambda-1);

T1 = p1*m4p + p2*m4id;
T2 = p1*m4p2 + p2*m4id;
sT1 = T1;
sT2 = T2;

tm = mcon({T1,T1,T2,T2},...
    {[-1,-2,1,4],[1,-3,-4,2],[3,2,-5,-6],[-8,4,3,-7]});


% T1 = rand(2,2,2,2);
% T2 = rand(3,5,3,2);


[A1,D1] = fn_svd_ad(T1,6);
[B1,C1] = fn_svd_bc(T1,6);
[A2,D2] = fn_svd_ad(T2,6);
[B2,C2] = fn_svd_bc(T2,6);

%%%%%%%%% check   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t1ad = mcon({A1,D1},{[-1,-2,1],[-3,-4,1]});
t1bc = mcon({B1,C1},{[-2,-3,1],[-4,-1,1]});
t2ad = mcon({A2,D2},{[-1,-2,1],[-3,-4,1]});
t2bc = mcon({B2,C2},{[-2,-3,1],[-4,-1,1]});
tv = t1ad - T1;
tv = reshape(tv,[1,numel(tv)]);
disp(max(abs(tv)));
tv = t1bc - T1;
tv = reshape(tv,[1,numel(tv)]);
disp(max(abs(tv)));
tv = t2ad - T2;
tv = reshape(tv,[1,numel(tv)]);
disp(max(abs(tv)));
tv = t2bc - T2;
tv = reshape(tv,[1,numel(tv)]);
disp(max(abs(tv)));
%%%%%%%  end check  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T1 = fn_comb_abcd(A1,B1,C2,D2);
T2 = fn_comb_abcd(A2,B2,C1,D1);

[A1,D1] = fn_svd_ad(T1,16);
[B2,C2] = fn_svd_bc(T2,16);

T = fn_comb_abcd(A1,B2,C2,D1);

Z = ncon({T},{[1,2,1,2]});
disp(['Z=',num2str(Z)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55


Num_truncate = 32;
Num_iteration = 1;
T_cell = cell(1,Num_iteration);
v_normfactor = ones(1,Num_iteration);

normfactor = sqrt(ncon({T,conj(T)},{[1,2,3,4],[1,2,3,4]}));
T = T./normfactor;

T_cell{1} = T;
v_normfactor(1) = normfactor;

Z = ncon({T},{[1,2,1,2]});
tv = v_normfactor(1:1);
tv_p = (1-1):-1:0;
tv_p = power(2,tv_p);
t_sum = sum(tv_p.*log(tv));
t_factor = exp(t_sum);
disp(['Z=',num2str(Z*t_factor)]);

for i=2:Num_iteration
    
    disp(['being iteration-',num2str(i)]);
    
    [A,D] = fn_svd_ad(T,Num_truncate);
    [B,C] = fn_svd_bc(T,Num_truncate);
    
    T = fn_comb_abcd(A,B,C,D);
    normfactor = sqrt(ncon({T,conj(T)},{[1,2,3,4],[1,2,3,4]}));
    disp(normfactor);
    T = T./normfactor;
    T_cell{i} = T;
    v_normfactor(i) = normfactor;
    
    Z = ncon({T},{[1,2,1,2]});
    tv = v_normfactor(1:i);
    tv_p = (i-1):-1:0;
    tv_p = power(2,tv_p);
    t_sum = sum(tv_p.*log(tv));
    t_factor = exp(t_sum);
    disp(['Z=',num2str(Z*t_factor)]);

end

tm = mcon({T,T,T,T},{[-1,-2,1,4],[1,-3,-4,2],...
    [3,2,-5,-6],[-8,4,3,-7]});

% t1 = permute(T,[4,3,1,2]);
% t2 = permute(T,[1,4,2,3]);
% t3 = permute(T,[2,1,3,4]);
% t4 = permute(T,[3,2,4,1]);

t1 = rand(8,8,8,8);
t2 = rand(8,8,8,8);
t3 = rand(8,8,8,8);
t4 = rand(8,8,8,8);

%%%% original tensor %%%%%%%%%%
tm2 = mcon({t1,t2,t3,t4},{[4,1,-1,-2],[1,2,-3,-4],...
    [2,3,-5,-6],[3,4,-7,-8]});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% v_size = size(t1);
% shape1 = [v_size(1),v_size(2),v_size(3)*v_size(4)];
% shape2 = v_size;
% R0 = eye(size(t1,1),size(t1,2));
% while 1
%     
%     t4_reshape = reshape(t4,shape1);
%     [R0,t4] = fn_lq(t4_reshape,R0);
%     t4 = reshape(t4,shape2);
%     
%     t3_reshape = reshape(t3,shape1);
%     [R0,t3] = fn_lq(t3_reshape,R0);
%     t3 = reshape(t3,shape2);
%     
%     t2_reshape = reshape(t2,shape1);
%     [R0,t2] = fn_lq(t2_reshape,R0);
%     t2 = reshape(t2,shape2);
%     
%     t1_reshape = reshape(t1,shape1);
%     [R0,t1] = fn_lq(t1_reshape,R0);
%     t1 = reshape(t1,shape2);
%     
%     tm3 = mcon({t1,t2,t3,t4,R0},{[5,1,-1,-2],[1,2,-3,-4],...
%     [2,3,-5,-6],[3,4,-7,-8],[4,5]});
% 
%     disp(max(reshape(tm3-tm2,[1,numel(tm2)])));
%     
% end
% 


[L0,normfactor_L0] = fn_filter_L(t10,t20,t30,t40);

[R0,normfactor_R0] = fn_filter_R(t10,t20,t30,t40);

disp(L0);
disp(R0);
disp(normfactor_L0)
disp(normfactor_R0);



%%%%%%%%%% END MAIN SCRIPTS %%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%  working functions %%%%%%%%%%%%%%%%%%
function [R0,store_normfactor] = fn_filter_R(t10,t20,t30,t40)

R0 = eye(size(t10,1),size(t10,1));
store_normfactor = 1;

v_size = size(t10);
shape1 = [v_size(1),v_size(2),v_size(3)*v_size(4)];
t10_reshape = reshape(t10,shape1);
t20_reshape = reshape(t20,shape1);
t30_reshape = reshape(t30,shape1);
t40_reshape = reshape(t40,shape1);
iteration_round = 0;
while 1

R0_pre = R0;

[R0,~] = fn_lq(t40_reshape,R0);
disp('corss 4');
disp(R0);

[R0,~] = fn_lq(t30_reshape,R0);
disp('corss 3');
disp(R0);

[R0,~] = fn_lq(t20_reshape,R0);
disp('corss 2');
disp(R0);

[R0,~] = fn_lq(t10_reshape,R0);
disp('corss 1');
disp(R0);

normfactor = max(abs(reshape(R0,[1,numel(R0)])));
R0 = R0./normfactor;
store_normfactor = store_normfactor * normfactor;

disp(['store_normfactor=',num2str(store_normfactor)]);
iteration_round = iteration_round + 1;
disp(['finished round ',num2str(iteration_round)]);
disp(num2str(abs(R0-R0_pre)));

if max(abs(reshape(R0 - R0_pre,[1,numel(R0)])))<1e-10
    break;
end

end

end

function [L0,store_normfactor] = fn_filter_L(t10,t20,t30,t40)

L0 = eye(size(t10,1),size(t10,1));
store_normfactor = 1;

v_size = size(t10);
shape1 = [v_size(1),v_size(2),v_size(3)*v_size(4)];
t10_reshape = reshape(t10,shape1);
t20_reshape = reshape(t20,shape1);
t30_reshape = reshape(t30,shape1);
t40_reshape = reshape(t40,shape1);
iteration_round = 0;
while 1

L0_pre = L0;

[~,L0] = fn_qr(L0,t10_reshape);
disp('corss 1');
disp(L0);

[~,L0] = fn_qr(L0,t20_reshape);
disp('corss 2');
disp(L0);

[~,L0] = fn_qr(L0,t30_reshape);
disp('corss 3');
disp(L0);

[~,L0] = fn_qr(L0,t40_reshape);
disp('corss 4');
disp(L0);

normfactor = max(abs(reshape(L0,[1,numel(L0)])));
L0 = L0./normfactor;
store_normfactor = store_normfactor * normfactor;

disp(['store_normfactor=',num2str(store_normfactor)]);
iteration_round = iteration_round + 1;
disp(['finished round ',num2str(iteration_round)]);
disp(num2str(abs(L0-L0_pre)));

if max(abs(reshape(L0 - L0_pre,[1,numel(L0)])))<1e-10
    break;
end

end

end


function [L,Q] = fn_lq(T1,T2)

t = mcon({T1,T2},{[-1,1,-3],[1,-2]});
vs = size(t);
t = reshape(permute(t,[1,3,2]),vs(1),vs(3)*vs(2));
[Q,R]=qr(t');
L = R';
Q = Q';
Q = Q(1:vs(2),:);
L = L(:,1:vs(1));

v_diag = sign(diag(L));
tran_m = diag(v_diag);
L = L*tran_m;
Q = tran_m*Q;

Q = reshape(Q,[vs(1),vs(3),vs(2)]);
Q = permute(Q,[1,3,2]);
L = reshape(L,vs(2),vs(2));

end


function [Q,R] = fn_qr(T1,T2)

t = mcon({T1,T2},{[-1,1],[1,-2,-3]});
vs = size(t);
t = reshape(permute(t,[1,3,2]),vs(1)*vs(3),vs(2));
[Q,R]=qr(t);
Q = Q(:,1:vs(2));
R = R(1:vs(2),:);

v_diag = sign(diag(R));
tran_m = diag(v_diag);
Q = Q*tran_m;
R = tran_m*R;

Q = reshape(Q,[vs(1),vs(3),vs(2)]);
Q = permute(Q,[1,3,2]);
R = reshape(R,vs(2),vs(2));

end


function res = fn_comb_abcd(A,B,C,D)

res = mcon({D,C,A,B},...
    {[1,4,-1],[2,1,-2],[3,2,-3],[4,3,-4]});

end



function [A,D] = fn_svd_ad(T,truncate)

v_size = size(T);
if length(v_size)<4
   v_size = [v_size,ones(1,4-length(v_size))]; 
end

t = reshape(T,v_size(1)*v_size(2),v_size(3)*v_size(4));
[u,s,v] = svd(t,'econ');
v_diag = diag(s);
tuncate_small = sum(v_diag>1e-12);
disp('the AD decomposition value: ');
disp(num2str(v_diag(1:min(10,length(v_diag)))'));
truncate = min(truncate,tuncate_small);
truncate = min(truncate,tuncate_small);
trun = min([truncate,v_size(1)*v_size(2),v_size(3)*v_size(4)]);
u = u(:,1:trun);
s = s(1:trun,1:trun);
v = v(:,1:trun);
u1 = u*sqrt(s);
v1 = sqrt(s)*v';

A = reshape(u1,[v_size(1),v_size(2),trun]);
D = reshape(v1,[trun,v_size(3),v_size(4)]);
D = permute(D,[2,3,1]);

end

function [B,C] = fn_svd_bc(T,truncate)

v_size = size(T);
if length(v_size)<4
   v_size = [v_size,ones(1,4-length(v_size))]; 
end

t = permute(T,[2,3,4,1]);
t = reshape(t,v_size(2)*v_size(3),v_size(4)*v_size(1));
[u,s,v] = svd(t,'econ');
v_diag = diag(s);
tuncate_small = sum(v_diag>1e-12);
disp('the BC decomposition value: ');
disp(num2str(v_diag(1:min(10,length(v_diag)))'));
truncate = min(truncate,tuncate_small);
trun = min([truncate,v_size(2)*v_size(3),v_size(4)*v_size(1)]);
u = u(:,1:trun);
s = s(1:trun,1:trun);
v = v(:,1:trun);
u1 = u*sqrt(s);
v1 = sqrt(s)*v';

B = reshape(u1,[v_size(2),v_size(3),trun]);
C = reshape(v1,[trun,v_size(4),v_size(1)]);
C = permute(C,[2,3,1]);

end


function res = mcon(cTensor,cIndex)

numTen = length(cIndex);
numIndex = 0;
numIndexOut = 0;
for i=1:numTen

    numIndex = numIndex + length(cIndex{i});
    numIndexOut = numIndexOut + sum(cIndex{i}<0); 
%     disp(numIndexOut);
%     disp(numIndex);
    
end

numIndexIn = (numIndex - numIndexOut)/2;

v_order(netcon(cIndex,0,2,1,1)) = 1:numIndexIn;
c_order = cell(1,numTen);
for i=1:numTen
    
   tv = cIndex{i};
   for j = 1:length(tv)
       if tv(j)>0
           tv(j) = v_order(tv(j));
       end
   end
   
   c_order{i} = tv;
   
end

res = ncon(cTensor,c_order);

end