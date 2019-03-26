function res = fn_mps_to_bigmatrix(mps)

Num = length(mps);

index = cell(1,Num);

size_up = 1;
size_down = 1;
for k=2:Num-1
   size_up = size_up*size(mps{k},2);
   size_down = size_down*size(mps{k},4);
   index{k} = [k-1,-k,k,-(k+Num)]; 
end

index{1} = [-2*Num-1,-1,1,-(Num+1)];
size_up = size_up*size(mps{1},2);
size_down = size_down*size(mps{1},4);
index{Num} = [Num-1,-Num,-(2*Num+2),-2*Num];
size_up = size_up*size(mps{Num},2);
size_down = size_down*size(mps{Num},4);

size_left = size(mps{1},1);
size_right = size(mps{Num},3);

res = mcon(mps,index);
res = reshape(res,[size_up,size_down,size_left,size_right]);

end