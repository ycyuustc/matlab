mid3 = eye(3);

c0 = ncon({mid3,mid3},{[-1,-3],[-2,-4]});
c1 = ncon({mid3,mid3},{[-1,-4],[-2,-3]});
c2 = ncon({mid3,mid3},{[-1,-2],[-3,-4]});

fn9to3 = @(m) reshape(permute(m,[2,1,4,3]),9,9);

m0 = fn9to3(c0);
m1 = fn9to3(c1);
m2 = fn9to3(c2);
