function dX = angCW(t, X, Tp, Te, lambda0)
% 生存型微分对策中构造CW方程

% 载入常数
global CONSTANTS
omg = CONSTANTS.refOmega;

% 由时刻t得到当前时刻的协态量
PhiMtx = lambdaTranMtx(omg, t);
lbd = PhiMtx * lambda0;

% 由该时刻协态量计算追逃双方控制量
r = norm([lbd(4), lbd(5), lbd(6)]);
if abs(r) < 1e-6
    alpha = 0;
    beta = 0;
else
    beta = asin(lbd(6)/r);
    calpha = lbd(4) / (norm([lbd(4), lbd(5)]));
    salpha = lbd(5) / (norm([lbd(4), lbd(5)]));
    alpha = atan2(salpha, calpha);
end

% CW方程，控制量为推力方向角

% 状态微分方程
dX = zeros(6,1);
dX(1) = X(4);
dX(2) = X(5);
dX(3) = X(6);
dX(4) = 2*omg*X(5) + 3*omg*omg*X(1) + (Te-Tp) * cos(alpha)*cos(beta);
dX(5) = -2*omg*X(4) + (Te - Tp) * sin(alpha)*cos(beta);
dX(6) = -omg*omg*X(3) + (Te - Tp) * sin(beta);
