% ����ʱ��΢�ֶԲ����
clear all
clc

%% ��������
% ������
earthRadius         =   6378.137 * 1000;                % ����뾶��m 
gravParam           =   3.986004418E14;                 % ��������������m3/s2
% �ο�Բ�������
refHeight           =   400 * 1000;                     % �ο�����߶ȣ�m
refRadius           =   earthRadius + refHeight;        % �ο�����볤�ᣬm
refOmega            =   sqrt(gravParam / (refRadius*refRadius*refRadius)); % �ο�������ٶȣ�rad/s
% ׷�ٺ����������ݺ�������ʼ�������
Xp0 = [-6*1000; -16*1000; 4*1000; -9; 13.6; 0];      % �ռ�Բ�Ʒ�״̬
% Xp0 = [-10*1000; -50*1000; 4*1000; 0; 0; 0];
Xe0 = [0; 0; 0; 0; 0; 0];
% �Բ߳�ʼ״̬����
X0 = Xe0 - Xp0;
% ΢�ֶԲ�ʱ��ͻ��ֲ���
gameT               =   60 * 90;                         % s
gameT               =   2*2*pi/refOmega;                  % s,һ���Ʒ�����
stepT               =   5;                              % s

%% ����λ
% ��������λ����Ҫ����
scales.length       =   refRadius;                      % ��������λ���ο�����뾶
scales.speed        =   sqrt(gravParam/scales.length);  % �ٶ�����λ����������λ��Ӧ��Բ����ٶȴ�С
% ��������λ���ɻ�������λ����
scales.time         =   scales.length/scales.speed;     % ʱ������λ������/�ٶ�
scales.acceleration =   scales.speed/scales.time;       % ���ٶ�����λ���ٶ�/ʱ��
scales.omega        =   1/scales.time;                  % ���ٶ�����λ��1/ʱ��
scales.gravparam    =   scales.acceleration*scales.length^2;    % ������������λ�����ٶ�*����^2
% ����֮��ĳ�����ȫ�ֱ������������ļ���ʹ��
global CONSTANTS
CONSTANTS.refOmega  =   refOmega / scales.omega;        % ����֮��Ĳο�������ٶ�
CONSTANTS.gameT     =   gameT / scales.time;            % ����֮��ĶԲ�ʱ��
CONSTANTS.stepT     =   stepT / scales.time;            % ����֮��Ļ��ֲ���
CONSTANTS.X0 = X0;
CONSTANTS.X0(1:3) = X0(1:3) / scales.length;
CONSTANTS.X0(4:6) = X0(4:6) / scales.speed;
CONSTANTS.Xp0 = Xp0;
CONSTANTS.Xp0(1:3) = Xp0(1:3) / scales.length;
CONSTANTS.Xp0(4:6) = Xp0(4:6) / scales.speed;
CONSTANTS.Xe0 = Xe0;
CONSTANTS.Xe0(1:3) = Xe0(1:3) / scales.length;
CONSTANTS.Xe0(4:6) = Xe0(4:6) / scales.speed;

%% CW�����е�����ϵ������
A = zeros(6,6);
A(1,4) = 1;
A(2,5) = 1;
A(3,6) = 1;
A(4,1) = 3*CONSTANTS.refOmega*CONSTANTS.refOmega;
A(4,5) = 2*CONSTANTS.refOmega;
A(5,4) = -2*CONSTANTS.refOmega;
A(6,3) = -CONSTANTS.refOmega*CONSTANTS.refOmega;
B = zeros(6,3);
B(4,1) = 1;
B(5,2) = 1;
B(6,3) = 1;

%% ΢�ֶԲ�֧�������еļ�����������
% ״̬����Ӧ����Q�����Է�Ϊ���������һ��ֻ���������������룬һ��ͬʱ�������ٶ�Ҫ��
Q = zeros(6,6);
Q(1,1) = 1;
Q(2,2) = 1;
Q(3,3) = 1;
Q = eye(6);        % �����������Ƿ���Ҫ�ٶȴﵽһ��
% ��������Ӧ����
R1 = eye(3);
R2 = sqrt(2) * R1;

%% �����迨�᷽�����
R = inv(inv(R1) - inv(R2));
S = zeros(6,3);
E = eye(6);
P = care(A, B, Q, R, S, E);

%% ���ֻ�öԲ߹켣����Ӧ�Ŀ�����
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
tspan = [0: CONSTANTS.stepT: CONSTANTS.gameT];
% ���ֵõ��Բ߹켣
[TSPAN, X] = ode45(@(t,X)accInfCW(t, X, P, B, R1, R2), tspan, CONSTANTS.X0, options);
% ��һ����ÿ�����
len = length(TSPAN);
Up = zeros(len,3);
Ue = zeros(len,3);
for i = 1 : len
    Up(i,:) = inv(R1) * B' * P * X(i,:)';
    Ue(i,:) = inv(R2) * B' * P * X(i,:)';
end

%% ��Ծ��롢����ٶȴ�С
DR = TSPAN;
DV = TSPAN;
for i = 1 : len
    DR(i) = norm(X(i,1:3));
    DV(i) = norm(X(i,4:6));
end

%% ��׷��������������״̬���ֱ�������
% ��Ҫ���������õĶԲ߿��������л��ֻ�ȡ�����ȶ����Ҫ��ȫ�ֱ�����ʹ����΢�ַ������ܹ����ʵ�������
tPVec = TSPAN;
tEVec = TSPAN;
aPxVec = TSPAN; aPyVec = TSPAN; aPzVec = TSPAN; aExVec = TSPAN; aEyVec = TSPAN; aEzVec = TSPAN;
for i = 1 : len
    aPxVec(i) = Up(i,1); aPyVec(i) = Up(i,2); aPzVec(i) = Up(i,3);
    aExVec(i) = Ue(i,1); aEyVec(i) = Ue(i,2); aEzVec(i) = Ue(i,3);
end
% ���ֻ������������״̬�켣
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
[tPSpan, Xp] = ode45(@(t,X)accCW(t, X, tPVec, aPxVec, aPyVec, aPzVec), tPVec, CONSTANTS.Xp0, options);
[tESpan, Xe] = ode45(@(t,X)accCW(t, X, tEVec, aExVec, aEyVec, aEzVec), tEVec, CONSTANTS.Xe0, options);

%% �ֱ�õ�׷����������������һ�����������Ĺ켣
% ���׷�����������ƣ������Ĺ켣
[tPSpan_N, Xp_N] = ode45(@(t,X)accCW(t, X, tPVec, 0*aPxVec, 0*aPyVec, 0*aPzVec), tPVec, CONSTANTS.Xp0, options);
% ����������������ƣ�׷����������
[tESpan_N, Xe_N] = ode45(@(t,X)accCW(t, X, tEVec, 0*aExVec, 0*aEyVec, 0*aEzVec), tEVec, CONSTANTS.Xe0, options);

%% ���һ����������֮ǰ,�õ��Բ�֧����������ʱ��ı仯����
J = TSPAN;
J_C = TSPAN;
for i = 1 : len
    % ÿ��ʱ��˲ʱ��Ӧ�Ĵ�����֧������
    J_C(i) = X(i,:)*Q*X(i,:)' + Up(i,:)*R1*Up(i,:)' - Ue(i,:)*R2*Ue(i,:)';
end
J(1) = 0;
for i = 2 : len
    % ������ֵ���ֵ���ʽ
    t_temp = TSPAN(1:i);
    J_temp = J_C(1:i);
    J(i) = trapz(t_temp, J_temp);
end


%% ��״̬���Ϳ������������һ�����õ����һ��֮���״̬���Ϳ�����
T = TSPAN;
x = zeros(len,1);
y = zeros(len,1);
z = zeros(len,1);
vx = zeros(len,1);
vy = zeros(len,1);
vz = zeros(len,1);
upx = zeros(len,1);
upy = zeros(len,1);
upz = zeros(len,1);
uex = zeros(len,1);
uey = zeros(len,1);
uez = zeros(len,1);
up = zeros(len,1);
ue = zeros(len,1);
for i = 1 : len
    T(i) = T(i) * scales.time;
    X(i,1:3) = X(i,1:3) * scales.length;
    X(i,4:6) = X(i,4:6) * scales.speed;
    Up(i) = Up(i) * scales.acceleration;
    Ue(i) = Ue(i) * scales.acceleration;
    x(i) = X(i,1); y(i) = X(i,2); z(i) = X(i,3);
    vx(i) = X(i,4); vy(i) = X(i,5); vz(i) = X(i,6);
    upx(i) = Up(i,1); upy(i) = Up(i,2); upz(i) = Up(i,3);
    uex(i) = Ue(i,1); uey(i) = Ue(i,2); uez(i) = Ue(i,3);
    up(i) = norm(Up(i));
    ue(i) = norm(Ue(i));
end

%% ��׷�ٺ����������ݺ�������״̬�����й�һ��
for i = 1 : len
    tPSpan(i) = tPSpan(i) * scales.time;
    Xp(i,1:3) = Xp(i,1:3) * scales.length;
    Xp(i,4:6) = Xp(i,4:6) * scales.speed;
    tESpan(i) = tESpan(i) * scales.time;
    Xe(i,1:3) = Xe(i,1:3) * scales.length;
    Xe(i,4:6) = Xe(i,4:6) * scales.speed;
    % �������Ƶ�׷�����켣
    Xp_N(i,1:3) = Xp_N(i,1:3) * scales.length;
    Xp_N(i,4:6) = Xp_N(i,4:6) * scales.speed;
    % �������Ƶ��������켣
    Xe_N(i,1:3) = Xe_N(i,1:3) * scales.length;
    Xe_N(i,4:6) = Xe_N(i,4:6) * scales.speed;
end

%% ����ͼ��
%% �Բ�״̬�仯
% ��ά�ռ�켣
figure;
plot3(x, y, z);
% �Բ߸���״̬�����仯
figure;
plot(T, x/1000);
title('x');
figure;
plot(T,y/1000);
title('y');
figure;
plot(T,z/1000);
title('z');
figure;
plot(T,vx);
title('vx');
figure;
plot(T,vy);
title('vy');
figure;
plot(T,vz);
title('vz');
% DR
figure;
plot(T, DR*scales.length/1000);
title('DR');
%DV
figure;
plot(T, DV*scales.speed);
title('DV');

%% ׷������������״̬���仯 
% ��ά�켣
figure;
hold all;
plot3(Xp(1,1)/1000, Xp(1,2)/1000, Xp(1,3)/1000,'s');
plot3(Xe(1,1)/1000, Xe(1,2)/1000, Xe(1,3)/1000,'o');
legend('׷�ٺ�����', '���ݺ�����');
plot3(Xp(:,1)/1000, Xp(:,2)/1000, Xp(:,3)/1000);
plot3(Xe(:,1)/1000, Xe(:,2)/1000, Xe(:,3)/1000);
title('׷�ٺ����������ݺ�������ά�ռ�켣');
% ׷����������״̬�����仯
figure;
hold all;
plot(T, Xp(:,1)/1000);
plot(T, Xe(:,1)/1000);
legend('׷�ٺ�����', '���ݺ�����');
title('x');
figure;
hold all;
plot(T, Xp(:,2)/1000);
plot(T, Xe(:,2)/1000);
legend('׷�ٺ�����', '���ݺ�����');
title('y');
figure;
hold all;
plot(T, Xp(:,3)/1000);
plot(T, Xe(:,3)/1000);
legend('׷�ٺ�����', '���ݺ�����');
title('z');
figure;
hold all;
plot(T, Xp(:,4));
plot(T, Xe(:,4));
legend('׷�ٺ�����', '���ݺ�����');
title('vx');
figure;
hold all;
plot(T, Xp(:,5));
plot(T, Xe(:,5));
legend('׷�ٺ�����', '���ݺ�����');
title('vy');
figure;
hold all;
plot(T, Xp(:,6));
plot(T, Xe(:,6));
legend('׷�ٺ�����', '���ݺ�����');
title('vz');

%% ׷������������������
% ���Ե����������ϵĿ������仯
figure;
hold all;
plot(T, upx);
plot(T, upy);
plot(T, upz);
legend('a_{p,x}','a_{p,y}','a_{p,z}');
title('u_p');
figure;
hold all;
plot(T, uex);
plot(T, uey);
plot(T, uez);
legend('a_{e,x}','a_{e,y}','a_{e,z}');
title('u_e');
% �Ա�����׷������������������
figure;
hold all;
plot(T, upx);
plot(T, uex);
legend('׷�ٺ�����', '���ݺ�����');
title('a_x');
figure;
hold all;
plot(T, upy);
plot(T, uey);
legend('׷�ٺ�����', '���ݺ�����');
title('a_y');
figure;
hold all;
plot(T, upz);
plot(T, uez);
legend('׷�ٺ�����', '���ݺ�����');
title('a_z');
figure;
hold all;
plot(T, up);
plot(T, ue);
legend('׷�ٺ�����', '���ݺ�����');
title('a');

%% ׷���������Ͳ����������״̬���仯�Ƚ�
% ��ά�켣
figure;
hold all;
plot3(Xp(:,1)/1000, Xp(:,2)/1000, Xp(:,3)/1000);
plot3(Xp_N(:,1)/1000, Xp_N(:,2)/1000, Xp_N(:,3)/1000);
legend('����', '������');
plot3(Xp(1,1)/1000, Xp(1,2)/1000, Xp(1,3)/1000, 's');
title('׷������ά�ռ�켣');
% % ״̬�����仯
% figure;
% hold all;
% plot(T, Xp(:,1));
% plot(T, Xp_N(:,1));
% legend('����', '������');
% title('x');
% figure;
% hold all;
% plot(T, Xp(:,2));
% plot(T, Xp_N(:,2));
% legend('����', '������');
% title('y');
% figure;
% hold all;
% plot(T, Xp(:,3));
% plot(T, Xp_N(:,3));
% legend('����', '������');
% title('z');
% figure;
% hold all;
% plot(T, Xp(:,4));
% plot(T, Xp_N(:,4));
% legend('����', '������');
% title('vx');
% figure;
% hold all;
% plot(T, Xp(:,5));
% plot(T, Xp_N(:,5));
% legend('����', '������');
% title('vy');
% figure;
% hold all;
% plot(T, Xp(:,6));
% plot(T, Xp_N(:,6));
% legend('����', '������');
% title('vz');

%% �����������Ͳ����������״̬���仯�Ƚ�
% ��ά�켣
figure;
hold all;
plot3(Xe(:,1)/1000, Xe(:,2)/1000, Xe(:,3)/1000);
plot3(Xe_N(:,1)/1000, Xe_N(:,2)/1000, Xe_N(:,3)/1000);
legend('����', '������');
plot3(Xe(1,1)/1000, Xe(1,2)/1000, Xe(1,3)/1000,'o');
title('��������ά�ռ�켣');
% % ״̬�����仯
% figure;
% hold all;
% plot(T, Xe(:,1));
% plot(T, Xe_N(:,1));
% legend('����', '������');
% title('x');
% figure;
% hold all;
% plot(T, Xe(:,2));
% plot(T, Xe_N(:,2));
% legend('����', '������');
% title('y');
% figure;
% hold all;
% plot(T, Xe(:,3));
% plot(T, Xe_N(:,3));
% legend('����', '������');
% title('z');
% figure;
% hold all;
% plot(T, Xe(:,4));
% plot(T, Xe_N(:,4));
% legend('����', '������');
% title('vx');
% figure;
% hold all;
% plot(T, Xe(:,5));
% plot(T, Xe_N(:,5));
% legend('����', '������');
% title('vy');
% figure;
% hold all;
% plot(T, Xe(:,6));
% plot(T, Xe_N(:,6));
% legend('����', '������');
% title('vz');

%% ׷��������������������������µı仯����
X_PM_EN = Xe_N - Xp;
DR_PM_EN = TSPAN;
DV_PM_EN = TSPAN;
for i = 1 : len
    DR_PM_EN(i) = norm(X_PM_EN(i,1:3));
    DV_PM_EN(i) = norm(X_PM_EN(i,4:6));
end
figure;
plot(T, DR_PM_EN/1000);
title('DR_PM_EN');
figure;
plot(T, DV_PM_EN);
title('DV_PM_EN');

%% ֧��������ʱ��仯����
figure;
plot(TSPAN*scales.time, J);
title('֧������');

