function dX = accInfCW(t, X, P, B, R1, R2)
% 无限时域微分对策中构造CW方程

% 载入常数
global CONSTANTS
omg = CONSTANTS.refOmega;

% CW方程，控制量直接为推力加速度，采用反馈控制律获得
up = inv(R1) * B' * P * X;
ue = inv(R2) * B' * P * X;

% 状态微分方程
dX = zeros(6,1);
dX(1) = X(4);
dX(2) = X(5);
dX(3) = X(6);
dX(4) = 2*omg*X(5) + 3*omg*omg*X(1) + ue(1) - up(1);
dX(5) = -2*omg*X(4) + ue(2) - up(2);
dX(6) = -omg*omg*X(3) + ue(3) - up(3);
