%% prepare the original T matrix
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

lambda = 0.001;

p1 = 1/(1-lambda);
p2 = lambda/(lambda-1);

T1 = p1*m4p + p2*m4id;
T2 = p1*m4p2 + p2*m4id;

tm = mcon({T1,T1,T2,T2},...
    {[-1,-2,1,4],[1,-3,-4,2],[3,2,-5,-6],[-8,4,3,-7]});

[A1,D1] = fn_svd_ad(T1,4);
[B1,C1] = fn_svd_bc(T1,4);
[A2,D2] = fn_svd_ad(T2,4);
[B2,C2] = fn_svd_bc(T2,4);

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

if length(size(T))<4
    
   T = reshape(T,[size(T),ones(1,4-length(size(T)))]); 
    
end

T = fn_full_dim(T,4,4);


%%
n = 4;
nD = 8;

T01 = ttttt1;
T02 = ttttt2;

% T01 = T;
% T02 = T;

% Z = mcon({T01,T02,T01,T02,conj(T01),conj(T02),conj(T01),conj(T02)},...
%     {[1,2,10,9],[10,3,4,11],[12,11,5,6],[8,9,12,7],...
%     [1,2,14,13],[14,3,4,15],[16,15,5,6],[8,13,16,7]});

Z = mcon({T01,T02,T01,T02},...
    {[1,2,6,5],[6,3,1,7],[4,5,8,2],[8,7,4,3]});

disp(['inital Z = ',num2str(Z)]);

Num_iteration = 10;
v_norm_factor = ones(1,Num_iteration);
iteration = 1;
store_norm_factor = 1;
while 1

T1 = permute(T01,[4,3,1,2]);
T2 = permute(T02,[1,4,2,3]);
T3 = permute(T01,[2,1,3,4]);
T4 = permute(T02,[3,2,4,1]);
op = {T1,T2,T3,T4};
Z = mcon({op{1},op{2},op{3},op{4}},...
    {[5,6,1,3],[6,7,4,1],[7,8,2,4],[8,5,3,2]});
Z = store_norm_factor^4 * Z; 
disp(['normalized T gives Z = ',num2str(Z)]);

op = fn_filter_op(op);

Z = mcon({op{1},op{2},op{3},op{4}},...
    {[5,6,1,3],[6,7,4,1],[7,8,2,4],[8,5,3,2]});
Z = store_norm_factor^4 * Z; 
disp(['normalized T gives Z = ',num2str(Z)]);

norm_op = fn_norm_op(op);
t = norm_op^(-1/8);
op{1} = t.*op{1}; op{2} = t.*op{2};
op{3} = t.*op{3}; op{4} = t.*op{4};
v_norm_factor(iteration) = 1/t;
store_norm_factor = store_norm_factor^2 * 1/t;

Z = mcon({op{1},op{2},op{3},op{4}},...
    {[5,6,1,3],[6,7,4,1],[7,8,2,4],[8,5,3,2]});
Z = store_norm_factor^4 * Z; 
disp(['normalized T gives Z = ',num2str(Z)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
truncate = 8;
[B1,C1] = fn_svd_bc(permute(op{1},[3,4,2,1]),truncate);
[A2,D2] = fn_svd_ad(permute(op{2},[1,3,4,2]),truncate);
[B3,C3] = fn_svd_bc(permute(op{3},[2,1,3,4]),truncate);
[A4,D4] = fn_svd_ad(permute(op{4},[4,2,1,3]),truncate);

tm4 = mcon({C1,B1,A2,D2,B3,C3,D4,A4},...
    {[8,-1,1],[-2,2,1],[2,-3,3],[-4,4,3],...
    [4,-5,5],[-6,6,5],[6,-7,7],[-8,8,7]});

tm5 = mcon({op{1},op{2},op{3},op{4}},{[4,1,-1,-2],[1,2,-3,-4],...
    [2,3,-5,-6],[3,4,-7,-8]});

disp(max(abs(reshape(tm4-tm5,[1,numel(tm4)]))));

% tran_n_nD = zeros(8,8);
% tran_n_nD(logical(eye(min(n,n)))) = 1;

mps2{1} = permute(C1,[1,3,2]);
mps2{2} = permute(B1,[3,2,1]);
mps2{3} = permute(A2,[1,3,2]);
mps2{4} = permute(D2,[3,2,1]);
mps2{5} = permute(B3,[1,3,2]);
mps2{6} = permute(C3,[3,2,1]);
mps2{7} = permute(D4,[1,3,2]);
mps2{8} = permute(A4,[3,2,1]);

% 
% for i=1:8
%    mps2{i} = fn_full_dim(mps2{i},4,8);
% end

tm6 = mcon({mps2{1},mps2{2},mps2{3},mps2{4},...
    mps2{5},mps2{6},mps2{7},mps2{8}},...
    {[1,2,-1],[2,3,-2],[3,4,-3],[4,5,-4],...
    [5,6,-5],[6,7,-6],[7,8,-7],[8,1,-8]});

disp(max(abs(reshape(tm6-tm5,[1,numel(tm6)]))));
% disp(max(abs(reshape(tm6-tm4,[1,numel(tm6)]))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if iteration > 1
mps = mps2;
else
mps = fn_createrandommps_periodic(8,nD,size(op{1},1));
end

norm_mps = fn_norm_mps(mps);
t = norm_mps^(-1/16);
for i=1:8
   mps{i} = t.*mps{i}; 
end

temp = fn_norm_op(op);
temp = temp + fn_norm_mps(mps);
temp = temp - fn_norm_op_mps(op,mps);
temp = temp - conj(fn_norm_op_mps(op,mps));
disp(temp);

point = 3;
tensor_N = fn_cal_N(mps,point);
tensor_W = fn_cal_W(mps,op,point);

term1 = fn_norm_op(op);
term2 = mcon({tensor_N,mps{point},conj(mps{point})},...
    {[1,2,3,4],[1,2,5],[3,4,5]});
term3 = mcon({tensor_W,conj(mps{point})},...
    {[1,2,3],[2,3,1]});
term4 = mcon({conj(tensor_W),mps{point}},...
    {[1,2,3],[2,3,1]});

delta_op_mps = term1 + term2 - term3 - term4;
disp(['delta_op_mps = ',num2str(delta_op_mps)]);

point = 1;
direction = -1;
while 1
    
    disp(['point = ', num2str(point)]);
    
    tensor_N = fn_cal_N(mps,point);
    tensor_W = fn_cal_W(mps,op,point);
    vs_N = size(tensor_N);
    vs_W = size(tensor_W);   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mN = reshape(tensor_N,vs_N(1)*vs_N(2),vs_N(3)*vs_N(4));
    mN = (mN + mN')/2;    
    mW = reshape(tensor_W,vs_W(1),vs_W(2)*vs_W(3));  
    [u,s] = eig(mN); 
    X = mW * u * diag(1./diag(s)) * u';
    norm_X = trace(X * mN * X');    
    X = X./sqrt(norm_X);   
    X = permute(reshape(X,vs_W),[2,3,1]);
    mps{point} = X;
    
    term2 = mcon({tensor_N,X,conj(X)},...
        {[1,2,3,4],[1,2,5],[3,4,5]});
    term3 = mcon({tensor_W,conj(X)},...
        {[1,2,3],[2,3,1]});
    term4 = mcon({conj(tensor_W),X},...
        {[1,2,3],[2,3,1]});
    
    delta_op_mps = term1 + term2 - term3 - term4;
    disp(['delta_op_mps = ',num2str(delta_op_mps)]);
    
    if delta_op_mps < 2e-4
        break;
    end
    
    mps0 = mps;
    mps = fn_canon_mps(mps,point);
    
    if point<8
        point = point+1;
    else
        point = 1;
    end
    
%     if point == 8 || point == 1
%         direction = direction * (-1);   
%     end    
%     point = point + direction;
    
%     point = floor(rand()*8)+1;
    
end

T01 = mcon({mps{4},mps{1},mps{8},mps{5}},...
    {[-1,1,2],[3,-2,2],[-3,3,4],[1,-4,4]});

T02 = mcon({mps{7},mps{6},mps{3},mps{2}},...
    {[2,-1,1],[-2,2,3],[4,-3,3],[-4,4,1]});

Z2 = mcon(mps,...
    {[8,1,9],[1,2,11],[2,3,12],[3,4,9],...
    [4,5,10],[5,6,12],[6,7,11],[7,8,10]});

Z2 = Z2 * store_norm_factor^4; disp(Z2);

tm1 = mcon(mps,...
    {[1,2,-1],[2,3,-2],[3,4,-3],[4,5,-4],...
    [5,6,-5],[6,7,-6],[7,8,-7],[8,1,-8]});
tm2 = mcon(op,...
    {[1,2,-1,-2],[2,3,-3,-4],[3,4,-5,-6],[4,1,-7,-8]});

disp(max(abs(reshape(tm1-tm2,[1,numel(tm1)]))));

Z = mcon({T01,T02},{[1,2,3,4],[3,4,1,2]});
disp(Z);

%%%%%%%%%%%%%
iteration = iteration + 1;

if iteration > Num_iteration
    break;
end

end

%%%%%%%%%%%%%%%  END SCRIPTS   %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function mps = fn_canon_mps(mps,point)

tensor = mps{point};
v_size = size(tensor);
tensor = permute(tensor,[1,3,2]);
m_ten = reshape(tensor,v_size(1)*v_size(3),v_size(2));
[u,s,v] = svd(m_ten,'econ');
tensor = reshape(u,[v_size(1),v_size(3),v_size(2)]);
tensor = permute(tensor,[1,3,2]);
mps{point} = tensor;

if point<8
    point_next = point+1;
else
    point_next = 1;
end

ten_tran = s*v';

tensor = mps{point_next};
tensor = ncon({ten_tran,tensor},{[-1,1],[1,-2,-3]});
mps{point_next} = tensor;

end



function T = fn_full_dim(T,Num,d)

if Num < length(size(T))
    error('Num should be larger than dimension of T');
end

if length(size(T))<Num
    v1 = [size(T),ones(1,Num - length(size(T)))];
else
    v1 = size(T);
end
v2 = d*ones(1,Num);
v_ind = max([v1;v2]);

for i=1:Num
   
    d_in = v1(i);
    d_out = v_ind(i);   
    tm = zeros(d_in,d_out);
    tm(logical(eye(min(d_in,d_out)))) = 1;
   
    T = fn_contract(T,Num,i,tm,2,1);
    T = permute(T,[1:(i-1),Num,i:(Num-1)]);
    
end

end


function op = fn_filter_op(op)

T1 = op{1};
T2 = op{2};
T3 = op{3};
T4 = op{4};

[R41,norm_R41] = fn_filter_R(T1,T2,T3,T4);
[L41,norm_L41] = fn_filter_L(T1,T2,T3,T4);
[R12,norm_R12] = fn_filter_R(T2,T3,T4,T1);
[L12,norm_L12] = fn_filter_L(T2,T3,T4,T1);
[R23,norm_R23] = fn_filter_R(T3,T4,T1,T2);
[L23,norm_L23] = fn_filter_L(T3,T4,T1,T2);
[R34,norm_R34] = fn_filter_R(T4,T1,T2,T3);
[L34,norm_L34] = fn_filter_L(T4,T1,T2,T3);

[u,s,v] = svd(L41*R41);
PR41 = R41 * v * diag(power(diag(s),-1/2));
PL41 = diag(power(diag(s),-1/2)) * u' * L41;

[u,s,v] = svd(L12*R12);
PR12 = R12 * v * diag(power(diag(s),-1/2));
PL12 = diag(power(diag(s),-1/2)) * u' * L12;

[u,s,v] = svd(L23*R23);
PR23 = R23 * v * diag(power(diag(s),-1/2));
PL23 = diag(power(diag(s),-1/2)) * u' * L23;

[u,s,v] = svd(L34*R34);
PR34 = R34 * v * diag(power(diag(s),-1/2));
PL34 = diag(power(diag(s),-1/2)) * u' * L34;

PL41 = permute(PL41,[2,1]);
PL12 = permute(PL12,[2,1]);
PL23 = permute(PL23,[2,1]);
PL34 = permute(PL34,[2,1]);

T1 = mcon({T1,PL41,PR12,PR34,PL23},...
    {[1,2,3,4],[1,-1],[2,-2],[3,-3],[4,-4]});
T2 = mcon({T2,PL12,PR23,PR41,PL34},...
    {[1,2,3,4],[1,-1],[2,-2],[3,-3],[4,-4]});
T3 = mcon({T3,PL23,PR34,PR12,PL41},...
    {[1,2,3,4],[1,-1],[2,-2],[3,-3],[4,-4]});
T4 = mcon({T4,PL34,PR41,PR23,PL12},...
    {[1,2,3,4],[1,-1],[2,-2],[3,-3],[4,-4]});

op = {T1,T2,T3,T4};

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


function res = fn_8leg_norm(psi)

res = mcon({psi,conj(psi)},...
    {[1,2,3,4,5,6,7,8],[1,2,3,4,5,6,7,8]});

end

function res = fn_norm_op(op)

tt = cell(1,8);
tt{1} = op{1}; tt{2} = op{2}; tt{3} = op{3}; tt{4} = op{4};
tt{5} = conj(op{1}); tt{6} = conj(op{2});
tt{7} = conj(op{3}); tt{8} = conj(op{4});

tc = cell(1,8);
tc{1} = [4,1,9,10];
tc{2} = [1,2,11,12];
tc{3} = [2,3,13,14];
tc{4} = [3,4,15,16];

tc{5} = [8,5,9,10];
tc{6} = [5,6,11,12];
tc{7} = [6,7,13,14];
tc{8} = [7,8,15,16];

res = mcon(tt,tc);
end

function res = fn_norm_mps(mps)

tt = cell(1,16);
for i=1:8
   tt{i} = mps{i};
   tt{i+8} = conj(mps{i});
end

tc = cell(1,16);
tc{1} = [8,1,17];
tc{2} = [1,2,18];
tc{3} = [2,3,19];
tc{4} = [3,4,20];
tc{5} = [4,5,21];
tc{6} = [5,6,22];
tc{7} = [6,7,23];
tc{8} = [7,8,24];

tc{9} = [16,9,17];
tc{10} = [9,10,18];
tc{11} = [10,11,19];
tc{12} = [11,12,20];
tc{13} = [12,13,21];
tc{14} = [13,14,22];
tc{15} = [14,15,23];
tc{16} = [15,16,24];

res = mcon(tt,tc);

end

function res = fn_norm_op_mps(op,mps)

tt = cell(1,12);
for i=1:8
   tt{i} = mps{i};
end
for i=1:4
   tt{i+8} = op{i}; 
end

tc = cell(1,12);
tc{1} = [8,1,13];
tc{2} = [1,2,14];
tc{3} = [2,3,15];
tc{4} = [3,4,16];
tc{5} = [4,5,17];
tc{6} = [5,6,18];
tc{7} = [6,7,19];
tc{8} = [7,8,20];

tc{9} = [12,9,13,14];
tc{10} = [9,10,15,16];
tc{11} = [10,11,17,18];
tc{12} = [11,12,19,20];

res = mcon(tt,tc);

end

function tensor = fn_cal_W(input_t,input_T,point)

c_ten_t = cell(1,7);
c_ten_T = cell(1,4);

if mod(point,2)==1
    
    for i=1:7
        c_ten_t{i} = conj(input_t{mod(i+point-1,8)+1});
%         disp(i);
%         disp(mod(i+point-1,8)+1);
    end
    
    point_T = (point + 1)/2;
    for i=1:4
        c_ten_T{i} = input_T{mod(i+point_T-2,4)+1};
%         disp(i);
%         disp(mod(i+point_T-2,4)+1);
    end
    
    tc = cell(1,11);
    tc{1} = [-3,1,11];
    tc{2} = [1,2,12];
    tc{3} = [2,3,13];
    tc{4} = [3,4,14];
    tc{5} = [4,5,15];
    tc{6} = [5,6,16];
    tc{7} = [6,-2,17];
    
    tc{8} = [10,7,-1,11];
    tc{9} = [7,8,12,13];
    tc{10} = [8,9,14,15];
    tc{11} = [9,10,16,17];
    
    tensor = mcon([c_ten_t,c_ten_T],tc);
    
else
    
    for i=1:7
        c_ten_t{i} = conj(input_t{mod(i+point-1,8)+1});
%         disp(i);
%         disp(mod(i+point-1,8)+1);
    end
    
    point_T = point/2;
    for i=1:4
        c_ten_T{i} = input_T{mod(i+point_T-2,4)+1};
%         disp(i);
%         disp(mod(i+point_T-2,4)+1);
    end
    
    tc = cell(1,11);
    tc{1} = [-3,1,12];
    tc{2} = [1,2,13];
    tc{3} = [2,3,14];
    tc{4} = [3,4,15];
    tc{5} = [4,5,16];
    tc{6} = [5,6,17];
    tc{7} = [6,-2,11];
    
    tc{8} = [10,7,11,-1];
    tc{9} = [7,8,12,13];
    tc{10} = [8,9,14,15];
    tc{11} = [9,10,16,17];
    
    tensor = mcon([c_ten_t,c_ten_T],tc);
    
end

end

function tensor = fn_cal_N(input_tensor,point)

c_ten = cell(1,14);
for i=1:7
    c_ten{i} = input_tensor{mod(i+point-1,8)+1};
    c_ten{i+7} = conj(c_ten{i});
    %     disp(i);
    %     disp(mod(i+point-1,8)+1);
end

tc = cell(1,14);
tc{1} = [-2,1,13];
tc{2} = [1,2,14];
tc{3} = [2,3,15];
tc{4} = [3,4,16];
tc{5} = [4,5,17];
tc{6} = [5,6,18];
tc{7} = [6,-1,19];

tc{8} = [-4,7,13];
tc{9} = [7,8,14];
tc{10} = [8,9,15];
tc{11} = [9,10,16];
tc{12} = [10,11,17];
tc{13} = [11,12,18];
tc{14} = [12,-3,19];

tensor = mcon(c_ten,tc);

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
