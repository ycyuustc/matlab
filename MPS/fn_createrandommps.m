function [mps]=fn_createrandommps(N,D,d)

mps=cell(1,N);
mps{1}=randn(1,D,d)/sqrt(D);
mps{N}=randn(D,1,d)/sqrt(D);

for i=2:(N-1)
    mps{i}=randn(D,D,d)/sqrt(D);
end

end