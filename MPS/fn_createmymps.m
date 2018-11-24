function [mps]=fn_createmymps(N,d)

left =eye(d);
left = reshape(left,[1,d,d]);

right = eye(d);
right = reshape(right,[d,1,d]);

type1 = right;
type2 = left;

mps=cell(1,N);
mps{1}=left;
mps{N}=right;

for i=2:(N-1)
    if mod(i,2) == 1
    mps{i}=type1;
    else
        mps{i} = type2;
    end
end

end