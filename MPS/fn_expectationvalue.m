function [e,n]=fn_expectationvalue(mps,hset)

[M,N]=size(hset);
d=size(mps{1},3);

e=0;
for m=1:M
    em=1;
    for j=N:-1:1
        h=hset{m,j};
        h=reshape(h,[1,1,d,d]);
        em=fn_updateCright(em,mps{j},h,mps{j});
    end
    e=e+em;
end

n=1;
X=eye(d);
X=reshape(X,[1,1,d,d]);
for j=N:-1:1
    n=fn_updateCright(n,mps{j},X,mps{j});
end

e=e/n;

end