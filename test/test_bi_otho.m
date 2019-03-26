mpsA = fn_createrandommps(5,4,2);
mpsB = fn_createrandommps(5,4,2);

[mps1,mps2,cs] = fn_prepare_mps_ns(mpsA,mpsB);


tm = ncon({mps1{3},mps1{4},mps1{5},...
    conj(mps2{3}),conj(mps2{4}),conj(mps2{5})},...
    {[-1,4,1],[4,5,2],[5,8,3],[-2,6,1],[6,7,2],[7,8,3]});
disp(max(abs(reshape(tm - cs{2},[1,numel(tm)]))));

tm = ncon({mps1{2},mps1{3},mps1{4},mps1{5},...
    conj(mps2{2}),conj(mps2{3}),conj(mps2{4}),conj(mps2{5})},...
    {[-1,6,1],[6,7,2],[7,8,3],[8,5,4],...
    [-2,9,1],[9,10,2],[10,11,3],[11,5,4]});
disp(max(abs(reshape(tm - cs{1},[1,numel(tm)]))));

A = rand(2,2,2)+rand(2,2,2)*1j;
B = rand(2,2,2) + rand(2,2,2)*1j;
s = diag(rand(1,2));
[UA,UB,A_out,B_out,AB_out] = fn_prepare_onesite_lr_ns(A,B,s);

tm1 = ncon({A,conj(B)},{[-1,-2,1],[-3,-4,1]});
tm2 = ncon({UA,conj(UB),A_out,conj(B_out)},...
    {[-1,2,1],[-3,3,1],[2,-2],[3,-4]});
disp(max(abs(reshape(tm1-tm2,[1,numel(tm1)]))));

ncon({UA,conj(UA)},{[1,-1,2],[1,-2,2]})
ncon({UB,conj(UB)},{[1,-1,2],[1,-2,2]})

tm1 = ncon({UA,A_out},{[-1,1,-3],[1,-2]});
disp(max(abs(reshape(tm1-A,[1,numel(tm1)]))));
tm1 = ncon({UB,B_out},{[-1,1,-3],[1,-2]});
disp(max(abs(reshape(tm1-B,[1,numel(tm1)]))));

ncon({UA,conj(UB)},{[1,-1,2],[1,-2,2]})
ncon({UA,conj(UB),s},{[1,-1,3],[2,-2,3],[1,2]}) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mpsA = fn_createrandommps(5,4,2);
mpsB = fn_createrandommps(5,4,2);

[mps1,mps2,cs] = fn_prepare_mps_lr_ns(mpsA,mpsB);

ncon({mps1{1},conj(mps2{1})},{[1,-1,2],[1,-2,2]})

ncon({mps1{1},mps1{2},mps1{3},mps1{4},...
    conj(mps2{1}),conj(mps2{2}),conj(mps2{3}),conj(mps2{4})},...
    {[7,1,8],[1,2,9],[2,3,10],[3,-1,11],...
    [7,4,8],[4,5,9],[5,6,10],[6,-2,11]}) - cs{5}






