[result,m_pers,m3]=canonical_presentation_2(3);

mg1=m3(:,:,1);
mg2=m3(:,:,2);
mg3=m3(:,:,3);
mg4=m3(:,:,4);
mg5=m3(:,:,5);
mg6=m3(:,:,6);

% �����ʾ�Ѿ���mg1��mg6�����ˡ�
% mg3Ϊ�Ի���mg4Ϊ�û���������Ԫ����һ������Ԫ��

[v,d]=eig(mg4);

v1=v(:,1);
v2=mg3*v1;
%  v1��v2�ųɵĶ�ά�ռ���Ա�ʾ���Ⱥ

tm=eye(6);
tm(:,1)=v1;
tm(:,2)=v2;
[mQ,mR]=qr(tm);
% mQ�ĵ�һ�к͵ڶ���Ϊv1��v2.
% mQ��������ͬ���У�����һ����������

mg1=mQ'*mg1*mQ;
mg2=mQ'*mg2*mQ;
mg3=mQ'*mg3*mQ;
mg4=mQ'*mg4*mQ;
mg5=mQ'*mg5*mQ;
mg6=mQ'*mg6*mQ;

mgp1=zeros(2,2,6);
mgp1(:,:,1)=mg1(1:2,1:2);
mgp1(:,:,2)=mg2(1:2,1:2);
mgp1(:,:,3)=mg3(1:2,1:2);
mgp1(:,:,4)=mg4(1:2,1:2);
mgp1(:,:,5)=mg5(1:2,1:2);
mgp1(:,:,6)=mg6(1:2,1:2);
% ��һ���ʾ�Ѿ��ó�

mga1=mg1(3:6,3:6);
mga2=mg2(3:6,3:6);
mga3=mg3(3:6,3:6);
mga4=mg4(3:6,3:6);
mga5=mg5(3:6,3:6);
mga6=mg6(3:6,3:6);
%  �������������һ�ν�ά���ղŷ�����Ķ�ά��ʾ����
%  ��׼��ʾ��

[v,d]=eig(mga4);
v1=v(:,1);
% v1��Ӧ������ֵΪ1

v2=mga3*v1;

vp=v1+v2;
vm=v1-v2;

v1=vp/norm(vp);
v2=vm/norm(vm);
tm=eye(4);
tm(:,1)=v1;
tm(:,2)=v2;

[mQ,mR]=qr(tm);

mga1=mQ'*mga1*mQ;
mga2=mQ'*mga2*mQ;
mga3=mQ'*mga3*mQ;
mga4=mQ'*mga4*mQ;
mga5=mQ'*mga5*mQ;
mga6=mQ'*mga6*mQ;

mgp2=zeros(2,2,6);
mgp2(:,:,1)=mga1(1:2,1:2);
mgp2(:,:,2)=mga2(1:2,1:2);
mgp2(:,:,3)=mga3(1:2,1:2);
mgp2(:,:,4)=mga4(1:2,1:2);
mgp2(:,:,5)=mga5(1:2,1:2);
mgp2(:,:,6)=mga6(1:2,1:2);

mgb1=mga1(3:4,3:4);
mgb2=mga2(3:4,3:4);
mgb3=mga3(3:4,3:4);
mgb4=mga4(3:4,3:4);
mgb5=mga5(3:4,3:4);
mgb6=mga6(3:4,3:4);

mgp3=zeros(2,2,6);
mgp3(:,:,1)=mgb1;
mgp3(:,:,2)=mgb2;
mgp3(:,:,3)=mgb3;
mgp3(:,:,4)=mgb4;
mgp3(:,:,5)=mgb5;
mgp3(:,:,6)=mgb6;

% �������Ƿֽ��������ʾ
% ���ҷ������������������ά�ı�ʾ�������ɶ�һά��ʾ��ֱ�͡�
% �����������Ǳ�׼��ʾ��





