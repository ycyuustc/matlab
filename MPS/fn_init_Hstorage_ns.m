function [Hstorage]=fn_init_Hstorage_ns(mpsA,mpsB,hset)

d=size(hset{1,1},1);
[M,N]=size(hset);
Hstorage=cell(M,N+1);

for m=1:M
    Hstorage{m,1}=1;
    Hstorage{m,N+1}=1;
end

for j=N:-1:2
    for m=1:M
%         h=reshape(hset{m,j},[1,1,d,d]);
        h = hset{m,j};
        Hstorage{m,j}=fn_updateCright(Hstorage{m,j+1},mpsB{j},h,mpsA{j});
    end
end

end