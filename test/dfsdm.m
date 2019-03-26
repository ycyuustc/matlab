vk1 = linspace(kk1-0.00000001,kk1+0.00000001,1001);
vk2 = linspace(kk2-0.00000001,kk2+0.00000001,1001);

[mk1,mk2] = meshgrid(vk1,vk2);

% vk1 = linspace(kk1-0.00001,kk1+0.00001,1000);
% vk2 = linspace(kk2-0.00001,kk2+0.00001,1000);
% 
% [mk1,mk2] = meshgrid(vk1,vk2);
% 
tm = abs(f5(mk1,mk2,g,L,d))+abs(f5(mk2,mk1,g,L,d));
[ind1,ind2] = find(tm==min(min(tm)));
disp(abs(f5(kk1,kk2,g,L,d))+abs(f5(kk2,kk1,g,L,d)));
disp(min(min(tm)));