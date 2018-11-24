Num=4;
Num_term=5;
d=2;
D=4;

mps=fn_createrandommps(Num,D,d);
mps=fn_prepare_mps(mps);

hset=cell(Num_term,Num);

for i=1:Num_term
    for j=1:Num
        
        hset{i,j}=rand(d,d);
    
    end
end

Hstorage=fn_init_Hstorage(mps,hset);

Hleft=Hstorage(:,1);
Hright=Hstorage(:,2);
hsetj=hset(:,j);