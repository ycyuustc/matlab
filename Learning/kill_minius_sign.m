mI = eye(4);
m1 = [1,0,0,0;0,0,1,0;0,1,0,0;0,0,0,1];
m2 = [1,0,0,1;0,0,0,0;0,0,0,0;1,0,0,1];
mz = [1,0,0,0;0,-1,0,0;0,0,-1,0;0,0,0,1];

cI = reshape(mI,[2,2,2,2]);
c1 = reshape(m1,[2,2,2,2]);
c2 = reshape(m2,[2,2,2,2]);

y = [0,-1j;1j,0];
x = [0,1;1,0];
z = [1,0;0,-1];

fn4to2 = @(m) reshape(permute(m,[2,1,4,3]),4,4);

fn_add_1 = @(m,y) ncon({m,y,y},{[-1,1,2,-4],[-2,1],[2,-3]});
fn_add = @(m,y1,y2) ncon({m,y1,y2},{[-1,1,2,-4],[-2,1],[2,-3]});

fn_add_2 = @(m,y) ncon({m,y,y},{[1,2,-3,-4],[-1,1],[-2,2]});
fn_add_3 = @(m,y) ncon({m,y,y},{[2,-2,-3,1],[-1,2],[1,-4]});
fn_add_4 = @(m,y) ncon({m,y,y},{[-1,-2,1,2],[1,-3],[2,-4]});
fn_add_heng = @(m,y) ncon({m,y,y},{[1,-2,3,-4],[-1,1],[3,-3]});
fn_add_heng2 = @(m,y1,y2) ncon({m,y1,y2},{[1,-2,3,-4],[-1,1],[3,-3]});
fn_add_shu = @(m,y) ncon({m,y,y},{[-1,2,-3,4],[-2,2],[4,-4]});
fn_add_shu2 = @(m,y1,y2) ncon({m,y1,y2},{[-1,2,-3,4],[-2,2],[4,-4]});

c1add1 = fn_add_1(c1,y);  
c1add2 = fn_add_2(c1,y);
c1add3 = fn_add_3(c1,y);
c1add4 = fn_add_4(c1,y);
m1add1 = fn4to2(c1add1);
m1add2 = fn4to2(c1add2);
m1add3 = fn4to2(c1add3);
m1add4 = fn4to2(c1add4);

c2add1 = fn_add_1(c2,y);
c2add2 = fn_add_2(c2,y);
c2add3 = fn_add_3(c2,y);
c2add4 = fn_add_4(c2,y);
m2add1 = fn4to2(c2add1);
m2add2 = fn4to2(c2add2);
m2add3 = fn4to2(c2add3);
m2add4 = fn4to2(c2add4);

cIadd1 = fn_add_1(cI,y);
cIadd2 = fn_add_2(cI,y);
cIadd3 = fn_add_3(cI,y);
cIadd4 = fn_add_4(cI,y);
mIadd1 = fn4to2(cIadd1);
mIadd2 = fn4to2(cIadd2);
mIadd3 = fn4to2(cIadd3);
mIadd4 = fn4to2(cIadd4);



