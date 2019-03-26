mH = rand(3,3);
[U,D] = eig(mH);

m_left_vector = inv(U);
m_right_vector =U;

mL = m_left_vector;
mR = m_right_vector;

for i=1:3
    for j=1:3
        disp([i,j]);
        disp(mL(i,:)*mR(:,j));
        disp(mL(i,:)*mH*mR(:,j));
    end
end

tm = 0;
for i=1:3
   tm = tm + mR(:,i)*mL(i,:);
end
disp(tm);

T = U'*U;

bx1_R = [1;0;0];
bx1_L = [1,0,0];

aR = mL*bx1_R;
aL = bx1_L*mR;

disp(aL);
disp(aR'*T);