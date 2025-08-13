% 无限时域微分对策求解
clear all
clc

%% 参数配置
% 物理常数
earthRadius         =   6378.137 * 1000;                % 地球半径，m 
gravParam           =   3.986004418E14;                 % 地球引力常数，m3/s2
% 参考圆轨道参数
refHeight           =   400 * 1000;                     % 参考轨道高度，m
refRadius           =   earthRadius + refHeight;        % 参考轨道半长轴，m
refOmega            =   sqrt(gravParam / (refRadius*refRadius*refRadius)); % 参考轨道角速度，rad/s
% 追踪航天器与逃逸航天器初始轨道参数
Xp0 = [-6*1000; -16*1000; 4*1000; -9; 13.6; 0];      % 空间圆绕飞状态
% Xp0 = [-10*1000; -50*1000; 4*1000; 0; 0; 0];
Xe0 = [0; 0; 0; 0; 0; 0];
% 对策初始状态参数
X0 = Xe0 - Xp0;
% 微分对策时间和积分步长
gameT               =   60 * 90;                         % s
gameT               =   2*2*pi/refOmega;                  % s,一个绕飞周期
stepT               =   5;                              % s

%% 正则单位
% 基本正则单位，需要定义
scales.length       =   refRadius;                      % 长度正则单位，参考轨道半径
scales.speed        =   sqrt(gravParam/scales.length);  % 速度正则单位，长度正则单位对应的圆轨道速度大小
% 其它正则单位，由基本正则单位导出
scales.time         =   scales.length/scales.speed;     % 时间正则单位，长度/速度
scales.acceleration =   scales.speed/scales.time;       % 加速度正则单位，速度/时间
scales.omega        =   1/scales.time;                  % 角速度正则单位，1/时间
scales.gravparam    =   scales.acceleration*scales.length^2;    % 引力常数正则单位，加速度*长度^2
% 正则化之后的常数，全局变量，在所有文件中使用
global CONSTANTS
CONSTANTS.refOmega  =   refOmega / scales.omega;        % 正则化之后的参考轨道角速度
CONSTANTS.gameT     =   gameT / scales.time;            % 正则化之后的对策时间
CONSTANTS.stepT     =   stepT / scales.time;            % 正则化之后的积分步长
CONSTANTS.X0 = X0;
CONSTANTS.X0(1:3) = X0(1:3) / scales.length;
CONSTANTS.X0(4:6) = X0(4:6) / scales.speed;
CONSTANTS.Xp0 = Xp0;
CONSTANTS.Xp0(1:3) = Xp0(1:3) / scales.length;
CONSTANTS.Xp0(4:6) = Xp0(4:6) / scales.speed;
CONSTANTS.Xe0 = Xe0;
CONSTANTS.Xe0(1:3) = Xe0(1:3) / scales.length;
CONSTANTS.Xe0(4:6) = Xe0(4:6) / scales.speed;

%% CW方程中的两个系数矩阵
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

%% 微分对策支付函数中的几个矩阵配置
% 状态量相应矩阵Q，可以分为两种情况，一种只考虑两航天器距离，一种同时还考虑速度要求
Q = zeros(6,6);
Q(1,1) = 1;
Q(2,2) = 1;
Q(3,3) = 1;
Q = eye(6);        % 这里来控制是否还需要速度达到一致
% 控制量相应矩阵
R1 = eye(3);
R2 = sqrt(2) * R1;

%% 代数黎卡提方程求解
R = inv(inv(R1) - inv(R2));
S = zeros(6,3);
E = eye(6);
P = care(A, B, Q, R, S, E);

%% 积分获得对策轨迹及对应的控制律
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
tspan = [0: CONSTANTS.stepT: CONSTANTS.gameT];
% 积分得到对策轨迹
[TSPAN, X] = ode45(@(t,X)accInfCW(t, X, P, B, R1, R2), tspan, CONSTANTS.X0, options);
% 进一步获得控制律
len = length(TSPAN);
Up = zeros(len,3);
Ue = zeros(len,3);
for i = 1 : len
    Up(i,:) = inv(R1) * B' * P * X(i,:)';
    Ue(i,:) = inv(R2) * B' * P * X(i,:)';
end

%% 相对距离、相对速度大小
DR = TSPAN;
DV = TSPAN;
for i = 1 : len
    DR(i) = norm(X(i,1:3));
    DV(i) = norm(X(i,4:6));
end

%% 对追踪器和逃逸器的状态量分别进行求解
% 需要根据上面获得的对策控制量进行积分获取，首先定义必要的全局变量，使得在微分方程中能够访问到控制量
tPVec = TSPAN;
tEVec = TSPAN;
aPxVec = TSPAN; aPyVec = TSPAN; aPzVec = TSPAN; aExVec = TSPAN; aEyVec = TSPAN; aEzVec = TSPAN;
for i = 1 : len
    aPxVec(i) = Up(i,1); aPyVec(i) = Up(i,2); aPzVec(i) = Up(i,3);
    aExVec(i) = Ue(i,1); aEyVec(i) = Ue(i,2); aEzVec(i) = Ue(i,3);
end
% 积分获得两航天器的状态轨迹
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
[tPSpan, Xp] = ode45(@(t,X)accCW(t, X, tPVec, aPxVec, aPyVec, aPzVec), tPVec, CONSTANTS.Xp0, options);
[tESpan, Xe] = ode45(@(t,X)accCW(t, X, tEVec, aExVec, aEyVec, aEzVec), tEVec, CONSTANTS.Xe0, options);

%% 分别得到追踪器和逃逸器单独一方不做机动的轨迹
% 如果追踪器不做控制，求他的轨迹
[tPSpan_N, Xp_N] = ode45(@(t,X)accCW(t, X, tPVec, 0*aPxVec, 0*aPyVec, 0*aPzVec), tPVec, CONSTANTS.Xp0, options);
% 如果逃逸器不做控制，追踪器做控制
[tESpan_N, Xe_N] = ode45(@(t,X)accCW(t, X, tEVec, 0*aExVec, 0*aEyVec, 0*aEzVec), tEVec, CONSTANTS.Xe0, options);

%% 逆归一化各物理量之前,得到对策支付函数关于时间的变化规律
J = TSPAN;
J_C = TSPAN;
for i = 1 : len
    % 每个时刻瞬时对应的待积分支付函数
    J_C(i) = X(i,:)*Q*X(i,:)' + Up(i,:)*R1*Up(i,:)' - Ue(i,:)*R2*Ue(i,:)';
end
J(1) = 0;
for i = 2 : len
    % 采用数值积分的形式
    t_temp = TSPAN(1:i);
    J_temp = J_C(1:i);
    J(i) = trapz(t_temp, J_temp);
end


%% 对状态量和控制量进行逆归一化，得到逆归一化之后的状态量和控制量
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

%% 对追踪航天器和逃逸航天器的状态量进行归一化
for i = 1 : len
    tPSpan(i) = tPSpan(i) * scales.time;
    Xp(i,1:3) = Xp(i,1:3) * scales.length;
    Xp(i,4:6) = Xp(i,4:6) * scales.speed;
    tESpan(i) = tESpan(i) * scales.time;
    Xe(i,1:3) = Xe(i,1:3) * scales.length;
    Xe(i,4:6) = Xe(i,4:6) * scales.speed;
    % 不做控制的追踪器轨迹
    Xp_N(i,1:3) = Xp_N(i,1:3) * scales.length;
    Xp_N(i,4:6) = Xp_N(i,4:6) * scales.speed;
    % 不做控制的逃逸器轨迹
    Xe_N(i,1:3) = Xe_N(i,1:3) * scales.length;
    Xe_N(i,4:6) = Xe_N(i,4:6) * scales.speed;
end

%% 绘制图形
%% 对策状态变化
% 三维空间轨迹
figure;
plot3(x, y, z);
% 对策各个状态分量变化
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

%% 追踪器和逃逸器状态量变化 
% 三维轨迹
figure;
hold all;
plot3(Xp(1,1)/1000, Xp(1,2)/1000, Xp(1,3)/1000,'s');
plot3(Xe(1,1)/1000, Xe(1,2)/1000, Xe(1,3)/1000,'o');
legend('追踪航天器', '逃逸航天器');
plot3(Xp(:,1)/1000, Xp(:,2)/1000, Xp(:,3)/1000);
plot3(Xe(:,1)/1000, Xe(:,2)/1000, Xe(:,3)/1000);
title('追踪航天器和逃逸航天器三维空间轨迹');
% 追逃两航天器状态分量变化
figure;
hold all;
plot(T, Xp(:,1)/1000);
plot(T, Xe(:,1)/1000);
legend('追踪航天器', '逃逸航天器');
title('x');
figure;
hold all;
plot(T, Xp(:,2)/1000);
plot(T, Xe(:,2)/1000);
legend('追踪航天器', '逃逸航天器');
title('y');
figure;
hold all;
plot(T, Xp(:,3)/1000);
plot(T, Xe(:,3)/1000);
legend('追踪航天器', '逃逸航天器');
title('z');
figure;
hold all;
plot(T, Xp(:,4));
plot(T, Xe(:,4));
legend('追踪航天器', '逃逸航天器');
title('vx');
figure;
hold all;
plot(T, Xp(:,5));
plot(T, Xe(:,5));
legend('追踪航天器', '逃逸航天器');
title('vy');
figure;
hold all;
plot(T, Xp(:,6));
plot(T, Xe(:,6));
legend('追踪航天器', '逃逸航天器');
title('vz');

%% 追踪器和逃逸器控制量
% 各自的三个方向上的控制量变化
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
% 对比来看追踪器和逃逸器控制量
figure;
hold all;
plot(T, upx);
plot(T, uex);
legend('追踪航天器', '逃逸航天器');
title('a_x');
figure;
hold all;
plot(T, upy);
plot(T, uey);
legend('追踪航天器', '逃逸航天器');
title('a_y');
figure;
hold all;
plot(T, upz);
plot(T, uez);
legend('追踪航天器', '逃逸航天器');
title('a_z');
figure;
hold all;
plot(T, up);
plot(T, ue);
legend('追踪航天器', '逃逸航天器');
title('a');

%% 追踪器机动和不机动情况下状态量变化比较
% 三维轨迹
figure;
hold all;
plot3(Xp(:,1)/1000, Xp(:,2)/1000, Xp(:,3)/1000);
plot3(Xp_N(:,1)/1000, Xp_N(:,2)/1000, Xp_N(:,3)/1000);
legend('机动', '不机动');
plot3(Xp(1,1)/1000, Xp(1,2)/1000, Xp(1,3)/1000, 's');
title('追踪器三维空间轨迹');
% % 状态分量变化
% figure;
% hold all;
% plot(T, Xp(:,1));
% plot(T, Xp_N(:,1));
% legend('机动', '不机动');
% title('x');
% figure;
% hold all;
% plot(T, Xp(:,2));
% plot(T, Xp_N(:,2));
% legend('机动', '不机动');
% title('y');
% figure;
% hold all;
% plot(T, Xp(:,3));
% plot(T, Xp_N(:,3));
% legend('机动', '不机动');
% title('z');
% figure;
% hold all;
% plot(T, Xp(:,4));
% plot(T, Xp_N(:,4));
% legend('机动', '不机动');
% title('vx');
% figure;
% hold all;
% plot(T, Xp(:,5));
% plot(T, Xp_N(:,5));
% legend('机动', '不机动');
% title('vy');
% figure;
% hold all;
% plot(T, Xp(:,6));
% plot(T, Xp_N(:,6));
% legend('机动', '不机动');
% title('vz');

%% 逃逸器机动和不机动情况下状态量变化比较
% 三维轨迹
figure;
hold all;
plot3(Xe(:,1)/1000, Xe(:,2)/1000, Xe(:,3)/1000);
plot3(Xe_N(:,1)/1000, Xe_N(:,2)/1000, Xe_N(:,3)/1000);
legend('机动', '不机动');
plot3(Xe(1,1)/1000, Xe(1,2)/1000, Xe(1,3)/1000,'o');
title('逃逸器三维空间轨迹');
% % 状态分量变化
% figure;
% hold all;
% plot(T, Xe(:,1));
% plot(T, Xe_N(:,1));
% legend('机动', '不机动');
% title('x');
% figure;
% hold all;
% plot(T, Xe(:,2));
% plot(T, Xe_N(:,2));
% legend('机动', '不机动');
% title('y');
% figure;
% hold all;
% plot(T, Xe(:,3));
% plot(T, Xe_N(:,3));
% legend('机动', '不机动');
% title('z');
% figure;
% hold all;
% plot(T, Xe(:,4));
% plot(T, Xe_N(:,4));
% legend('机动', '不机动');
% title('vx');
% figure;
% hold all;
% plot(T, Xe(:,5));
% plot(T, Xe_N(:,5));
% legend('机动', '不机动');
% title('vy');
% figure;
% hold all;
% plot(T, Xe(:,6));
% plot(T, Xe_N(:,6));
% legend('机动', '不机动');
% title('vz');

%% 追踪器机动，逃逸器不机动情况下的变化规律
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

%% 支付函数随时间变化规律
figure;
plot(TSPAN*scales.time, J);
title('支付函数');

