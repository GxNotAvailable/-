function dX = accCW(t, X, tVec, axVec, ayVec, azVec)
% 固定逗留期微分对策中构造CW方程

% 载入常数
global CONSTANTS
omg = CONSTANTS.refOmega;

% CW方程，控制量直接为推力加速度，采用插值获得
ux =  interp1(tVec,axVec,t,'spline'); 
uy =  interp1(tVec,ayVec,t,'spline'); 
uz =  interp1(tVec,azVec,t,'spline'); 

% 状态微分方程
dX = zeros(6,1);
dX(1) = X(4);
dX(2) = X(5);
dX(3) = X(6);
dX(4) = 2*omg*X(5) + 3*omg*omg*X(1) + ux;
dX(5) = -2*omg*X(4) + uy;
dX(6) = -omg*omg*X(3) + uz;
