function [mR,mL,vd] = fn_eig(M, Keeped)
%FN_EIG 此处显示有关此函数的摘要
%   此处显示详细说明

keep = min(size(M,1),Keeped);
[u,d,v] = eig(M);
v = conj(v);
[d,index] = sort(diag(d),'descend');
u = u(:,index);
v = v(:,index);

u = u(:,1:keep);
v = v(:,1:keep);
d = d(1:keep);

point = 1;
while point<=keep-1
    
    ep1 = d(point);
    ep2 = d(point+1);
   
    if abs(imag(ep1))>1e-10 && abs(imag(ep1)+imag(ep2))<1e-10
        ur = real(u(:,point));
        ui = imag(u(:,point));
        vr = real(v(:,point));
        vi = imag(v(:,point));
        if abs(transpose(vr)*ur)<abs(transpose(vr)*ui)
            t = ur;
            ur = ui;
            ui = t;
        end
        [ur,vr] = fn_normalize(ur,vr);
        ui = ui - transpose(vr)*ui*ur;
        vi = vi - transpose(vi)*ur*vr;
        [ui,vi] = fn_normalize(ui,vi);
        u(:,point) = ur;
        u(:,point+1) = ui;
        v(:,point) = vr;
        v(:,point+1) = vi;
        point = point + 2;
    else
        [t1,t2] = fn_normalize(u(:,point),v(:,point));
        u(:,point) = t1;
        v(:,point) = t2;
        point = point + 1;
    end
         
end

if  point == keep
   ur = real(u(:,keep));
   vr = real(v(:,keep));
   [ur,vr] = fn_normalize(ur,vr);
   u(:,keep) = ur;
   v(:,keep) = vr;
end

mR = u;
mL = v;
vd = d;

end

function [v1,v2] = fn_normalize(v1,v2)

r = transpose(v2)*v1;
if r < 0
    v2 = -v2;
end

v1 = v1/sqrt(abs(r));
v2 = v2/sqrt(abs(r));

end