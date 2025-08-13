function F = endGame(y)
% 生存型微分对策中构造两边边值问题进行求解
% 自变量为追逃时间tf，协态量终值中的v1,v2,v3
% 等式约束为终端位置约束、终端哈密顿量横截条件约束

global CONSTANTS
omg = CONSTANTS.refOmega;

% 自变量中的追逃结束时间和协态量终值
tf = y(1);
lbd_f = zeros(6,1);
lbd_f(1) = y(2);
lbd_f(2) = y(3);
lbd_f(3) = y(4);

% 由协态量终值获得协态量初始值
PhiMtx_t02tf = lambdaTranMtx(omg, tf);      % 初始到终端协态量状态转移矩阵
if abs(det(PhiMtx_t02tf)) < 1e-10
    F = zeros(4,1);
    F(1) = 1e4;
    F(2) = 1e4;
    F(3) = 1e4;
    F(4) = 1e4;
    return;
else
    lbd_0 = inv(PhiMtx_t02tf) * lbd_f;               % t0时刻的协态量
    % 基于上述得到的协态量初始值，积分获得对策状态终值
    tspan = [0: CONSTANTS.stepT: tf];
    options_ode = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
    [TSPAN, XSPAN] = ode45(@(t,X)angCW(t, X, CONSTANTS.TpMax, CONSTANTS.TeMax, lbd_0), tspan, CONSTANTS.X0, options_ode);
    Xf = XSPAN(end,:);

    % 计算终端哈密顿量
    Hf = lbd_f(1)*Xf(4) + lbd_f(2)*Xf(5) + lbd_f(3)*Xf(6);

    % 构造终端等式约束
    F = zeros(4,1);
    F(1) = Xf(1);
    F(2) = Xf(2);
    F(3) = Xf(3);
    F(4) = Hf + 1;
end