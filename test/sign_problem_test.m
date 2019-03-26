sx = [0,1;1,0];
sy = [0,-1j;1j,0];
sz = [1,0;0,-1];
sI = eye(2);

sx1 = kron(sx,sI);
sx2 = kron(sI,sx);
sy1 = kron(sy,sI);
sy2 = kron(sI,sy);
sz1 = kron(sz,sI);
sz2 = kron(sI,sz);
sII = eye(4);
mp = zeros(4,4);
mp(1,1) = 1;
mp(4,4) = 1;
mp(2,3) = 1;
mp(3,2) = 1;
mp2 = zeros(4,4);
mp2(1,1) = 1;
mp2(4,4) = 1;
mp2(1,4) = 1;
mp2(4,1) = 1;

a = 0.1;

op = sII - a*mp;
disp(op);

disp(sx1'*sy2'*op*sy2*sx1);
disp(sx1'*op*sx1);
disp(sx1'*sx1'*sy2'*op*sy2*sx1*sx1);
disp(sy2'*op*sy2);

p1 = permute(reshape(mp,[2,2,2,2]),[1,2,3,4]);
p2 = permute(reshape(mp2,[2,2,2,2]),[1,2,3,4]);
id = permute(reshape(sII,[2,2,2,2]),[1,2,3,4]);

% Y1*P1*Y2 = P1
t = ncon({p1,sy,sy},{[-1,1,2,-4],[-2,1],[2,-3]});
fprintf('check Y1*P1*Y2 = P1 : ');
disp(max(reshape(t-p1,[1,numel(t)])));

t = ncon({sy,sy},{[-1,-3],[-2,-4]});
fprintf(' check Y1*Y2 = P1 - P2 : ');
disp(max(reshape(t-p1+p2,[1,numel(t)])));

t = ncon({sy,sy,p2},{[-1,1],[-2,2],[1,2,-3,-4]});
fprintf(' check Y1*Y2*P2 =  - P2 : ');
disp(max(reshape(t+p2,[1,numel(t)])));

t = ncon({sy,sy,p2},{[1,-3],[2,-4],[-1,-2,1,2]});
fprintf(' check P2*Y1*Y2 =  - P2 : ');
disp(max(reshape(t+p2,[1,numel(t)])));

lambda = rand();
theta = lambda/(1-lambda);

R1 = 1/(1-lambda)*p1 - lambda/(1-lambda)*id;
R2 = 1/(1-lambda)*p2 - lambda/(1-lambda)*id;

B = rand();
mB = expm(B*[1,0;0,-1]);

Z = mcon({R1,R1,R1,R1,R2,R2,R2,R2,mB,mB,mB,mB},...
    {[1,10,2,11],[2,13,3,14],[3,16,4,17],[4,19,1,20],...
    [5,11,6,9],[6,14,7,12],[7,17,8,15],[8,20,5,18],...
    [9,10],[12,13],[15,16],[18,19]});
disp(Z);

R1 = p1 + theta*p2;
R2 = theta*p1 + p2;
mB1 = mB;
mB2 = sy'*mB*sy;
Z = mcon({R1,R1,R1,R1,R2,R2,R2,R2,mB2,mB1,mB2,mB1},...
    {[1,10,2,11],[2,13,3,14],[3,16,4,17],[4,19,1,20],...
    [5,11,6,9],[6,14,7,12],[7,17,8,15],[8,20,5,18],...
    [9,10],[12,13],[15,16],[18,19]});
disp(Z);

