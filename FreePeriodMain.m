% 生存型微分对策求解，关键是构造两点边值问题并进行求解
clear all
clc

%% 参数配置
% 物理常数
earthRadius         =   6378.137 * 1000;                % 地球半径，m 
gravParam           =   3.986004418E14;                 % 地球引力常数，m3/s2
% 参考圆轨道参数
refHeight           =   36000 * 1000;                     % 参考轨道高度，m
refRadius           =   earthRadius + refHeight;        % 参考轨道半长轴，m
refOmega            =   sqrt(gravParam / (refRadius*refRadius*refRadius)); % 参考轨道角速度，rad/s
% 追踪航天器与逃逸航天器初始轨道参数
Xp = [-6*1000; -16*1000; 4*1000; -9; 13.6; 0];
Xe = [0; 0; 0; 0; 0; 0];
% 追踪器与逃逸器的最大推力加速度
g0 = 9.78;                      % 海平面重力加速度
TpMax = 0.06;
TeMax = 0.04;
% 对策初始状态参数
X0 = Xe - Xp;
% 微分对策积分步长
stepT               =   10;                              % s
% 两航天器初始相对速度和初始相对距离
r0 = norm([X0(1), X0(2), X0(3)]);
v0 = norm([X0(4), X0(5), X0(6)]);
% 简单估计得到追逃时间
a0 = TpMax - TeMax;
estT = (sqrt(v0*v0+2*a0*r0) - v0) / a0;

%% 正则单位
% 基本正则单位，需要定义
scales.length       =   1;%refRadius;                      % 长度正则单位，参考轨道半径
scales.speed        =   1;%sqrt(gravParam/scales.length);  % 速度正则单位，长度正则单位对应的圆轨道速度大小
% 其它正则单位，由基本正则单位导出
scales.time         =   scales.length/scales.speed;     % 时间正则单位，长度/速度
scales.acceleration =   scales.speed/scales.time;       % 加速度正则单位，速度/时间
scales.omega        =   1/scales.time;                  % 角速度正则单位，1/时间
scales.gravparam    =   scales.acceleration*scales.length^2;    % 引力常数正则单位，加速度*长度^2
% 正则化之后的常数，全局变量，在所有文件中使用
global CONSTANTS
CONSTANTS.refOmega  =   refOmega / scales.omega;        % 正则化之后的参考轨道角速度
CONSTANTS.TpMax = TpMax / scales.acceleration;
CONSTANTS.TeMax = TeMax / scales.acceleration;          % 正则化之后的两航天器最大推力加速度
CONSTANTS.stepT     =   stepT / scales.time;            % 正则化之后的积分步长
CONSTANTS.X0 = X0;
CONSTANTS.X0(1:3) = X0(1:3) / scales.length;
CONSTANTS.X0(4:6) = X0(4:6) / scales.speed;
CONSTANTS.Xp = Xp;
CONSTANTS.Xp(1:3) = Xp(1:3) / scales.length;
CONSTANTS.Xp(4:6) = Xp(4:6) / scales.speed;
CONSTANTS.Xe = Xe;
CONSTANTS.Xe(1:3) = Xe(1:3) / scales.length;
CONSTANTS.Xe(4:6) = Xe(4:6) / scales.speed;
CONSTANTS.estT = estT / scales.time;

%% 求解两边边值问题，本质上是求解一个非线性方程组
%% 采用进化算法获得初值
% 调用GA


% 调用DE


% 调用PSO
npars = 4;          % 设计变量有4个
lb = [600/scales.time -10 -10 -10];
ub = [1200/scales.time 10 10 10];
optPSO.npop = npars*15;         optPSO.niter = 300;
optPSO.Display = 1;          optPSO.Plot = 0;
optPSO.Vectorized = 0;     optPSO.Parallel = 0;
optPSO.StallIterLimit = 100;
[yopt, fval, exitflag, output] = ppso(@psoFunc,npars,lb,ub, optPSO);

% % 自变量为追逃时间tf，以及协态量终值中v1,v2,v3，需要给出初值
Y0 = yopt;
options = optimoptions('fsolve','Display','iter', 'TolFun', 1e-10, 'TolX', 1e-10, 'MaxFunEvals', 1000); % Option to display output
[Y,fval] = fsolve(@endGame, Y0, options) % Call solver




