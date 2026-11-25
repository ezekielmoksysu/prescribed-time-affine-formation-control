function [Rout,Pout] = formation_lissajous_sixstate
% =========================================================================
% 三车全六状态 UGV：BTG 预设时间跟踪 + FAS 全驱逆 + 指定时间 ESO
% -------------------------------------------------------------------------
% 单工况版本（不再包含无补偿/anti-windup/Yang 等对比工况）：
%   1) 三辆 UGV 用完整六状态模型（X,Y,phi,vx,vy,r + delta 一阶滤波）；
%   2) 分布式指定时间观测器估计虚拟领航者 (q0,v0)；
%   3) 每辆车用 BTG 预设时间律 + FAS 全驱逆得到 u_c=[Fx;Fyf]；
%   4) 物理限幅作用于 Fx 与前轮转角 delta（tanh 光滑饱和）；
%   5) 降阶 ESO 估计集总失配 Λ=BΔu+δ_i，并在控制律中前馈抵消。
%
% 状态布局（P.nS = 4N + 11N = 45）：
%   x(1:2N)      : qh，三车对领航者位置的估计
%   x(2N+1:4N)   : vh，三车对领航者速度的估计
%   之后每车 11 维：
%     b+1,b+2,b+3    : X,Y,phi
%     b+4,b+5,b+6    : vx,vy,r
%     b+7            : delta（实际前轮转角，一阶滤波）
%     b+8:b+9        : z1（ESO 内部状态，估计控制点速度 qdot）
%     b+10:b+11      : z2（ESO 估计的集总失配 Λ）
%
% 想改什么看哪里：
%   - 图/权重/领航者/仿射机动 : params()
%   - 初始偏差               : initial_state()
%   - BTG 增益 Tc/eps/kc/lamc: params() 中 "BTG 跟踪"
%   - 观测器 Te/gq/gv/iota  : params() 中 "领航者观测器"
%   - ESO Ta/b/k1/k2        : params() 中 "饱和失配观测器"
%   - 物理限幅 Fxmax/dmax    : params() 中 "UGV"
%   - 未知扰动幅值/频率       : rhs() 中 dvx_d/dvy_d 两行
%   - 仿真步长/时长           : params() 中 P.dt / P.T
%
% 输出：
%   formation_sixstate_lissajous_sat.png    三车轨迹
%   formation_lissajous_estimation.png      分布式估计误差（瞬态+包络）
%   formation_lissajous_tracking.png        动力学跟踪误差（瞬态+包络）
%   formation_lissajous_eso.png             ESO 重构误差（瞬态+包络）
%   formation_lissajous_inputs.png          指令/实际输入与饱和
%   formation_lissajous_metrics.mat         全部采集信号（P.save_stride 抽样）
%
% 用法：
%   formation_lissajous_sixstate
% =========================================================================
clc;                                    % 清空命令行
close all;                              % 关闭所有图窗
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','fulyactBC')); % 车形绘制

P = params();                           % 读取全部参数（结构体 P）
[t,X] = rk4(@(t,x) rhs(t,x,P),0:P.dt:P.T,initial_state(P),P.save_stride); % RK4 积分闭环（按 stride 记录）
R = unpack(X,t,P);                      % 从轨迹提取指标
print_metrics(R);                       % 打印终端指标
plot_trajectory(R,P);                   % 三车轨迹图
plot_estimation(R,P);                   % 分布式估计误差（瞬态+稳态）
plot_tracking(R,P);                     % 动力学跟踪误差（瞬态+稳态）
plot_eso(R,P);                          % ESO 重构误差（瞬态+稳态）
plot_inputs(R,P);                       % 指令/实际输入与饱和
plot_control_compare(R,P);               % 限幅前指令与实际输入的物理量对比
save_metrics(R,P);                      % 保存采集数据到 .mat

% ---- 把所有信号缓存到 MATLAB 基础工作区 ----
SIM=R;                                  % 结构体：字段即各信号
assignin('base','SIM',SIM);             % 结构体整体
assignin('base','P',P);                 % 参数
assignin('base','X',X);                 % 原始状态轨迹
fn=fieldnames(R);                       % R 的所有信号字段
for ii=1:numel(fn)
    assignin('base',fn{ii},R.(fn{ii}));  % 每个字段单独成变量
end
assignin('base','t',R.t);               % 时间轴
fprintf('signals cached in workspace: SIM, P, X, t, and %d R-fields\n',numel(fn));

if nargout>0, Rout=R; Pout=P; end       % 供诊断脚本取结果
end

% =========================================================================
function P = params
% 所有可调参数集中在这里；改参数只改这个函数
P.N=3;                                  % 车辆数
P.A=[0 4 3; 2 0 .5; 3 1 0];             % 邻接矩阵：A(i,j)=车 i 从车 j 收信息的权重
P.sigma=[1;0;0];                        % pinned 权重：只有 1 号车直接连领航者
P.L=diag(sum(P.A,2))-P.A;               % 拉普拉斯矩阵 L = D - A
P.H=P.L+diag(P.sigma);                  % 组合矩阵 H = L + diag(sigma)

P.x0=-38; P.vx_path=1.2; P.A_sin=5; P.w_sin=0.10; % 正弦轨迹
P.s0=0; P.s1=1; P.Tpath=96;
P.poff=[0 -4 -4; 0 -4 4];                % 编队偏置（随参考航向一起旋转）
P.Tramp=10;                               % 参考平滑过渡时间（reference governor）
P.e0=[2.5-4 1.0; -1.2-4 -1.8-4; -1.5-4 0.6+5]'; % 初始编队偏差（2xN）

% 领航者观测器（本文指定时间观测器）
P.Te=8; P.qe=1; P.epse=1.0; P.aa_max=0; P.aa_t1=8.2; P.aa_t2=9.2; P.aa_w=0.15;            % 截止时间、幂次、冻结参数
P.g0=5; P.iota0=0.5; P.eg=0.5; P.cg=5; P.ei=0.5; P.ci=5; P.eps_s=1e-3; P.rot_gain=1; % 自适应增益初值/增长率、tanh 边界层

% BTG 跟踪
P.Tc=10; P.q=2; P.epsi=1.5;             % 跟踪截止时间与增益参数（增大 εi 降低冻结增益）
P.lamc=1; P.kc=2; P.alpha=8;            % 常数增益、比例增益、滤波极点

% 饱和失配观测器
P.Ta=0.5; P.b=0.05; P.k1=4; P.k2=4;     % 冻结时间、正则化常数、增益

% UGV 物理参数
P.m=1500; P.Jz=1500; P.lf=1.2; P.lr=1.4;                   % 质量、转动惯量、前后轴距
P.Cf=5e4; P.Cr=5e4; P.Kg=0.5; P.xi=0.5; P.taud=0.02;       % 侧偏刚度、驱动分配、位移增益、转向滤波
P.Fxmax=3000; P.dmax=0.65; P.vxeps=0.05;                    % 纵向力限幅、转角限幅、低速保护

P.dt=5e-5; P.T=96;                                          % 固定增益；步长兼顾数值稳定与运行时间
P.save_stride=40;                                          % 记录/后处理抽样步长（40=每 2 ms）
P.nV=11;                                                   % 每车 11 个状态
P.nS=7*P.N + P.nV*P.N;                                     % 总状态数：估计器 4N + 自适应增益 3N + 车辆 11N
P.q0=@(t) liss_q(t,P);                                     % 领航者位置
P.v0=@(t) liss_v(t,P);                                     % 领航者速度
P.a0=@(t) liss_a(t,P);                                     % 领航者加速度
end

% =========================================================================
function q = liss_q(t,P)
q = [P.x0 + P.vx_path*t; P.A_sin*sin(P.w_sin*t)];
end

function v = liss_v(t,P)
v = [P.vx_path; P.A_sin*P.w_sin*cos(P.w_sin*t)];
end

function a = liss_a(t,P)
a = [0; -P.A_sin*P.w_sin^2*sin(P.w_sin*t)];
end

function [phi,phid,phidd] = ref_heading(t,P)
% 参考航向及其一、二阶导（由领航者速度方向定义；二阶导用中心差分）
v=P.v0(t); a=P.a0(t);
D=v(1)^2+v(2)^2; N=v(1)*a(2)-v(2)*a(1);
phi=atan2(v(2),v(1)); phid=N/D;
h=1e-3;
v1=P.v0(t+h); a1=P.a0(t+h); v2=P.v0(t-h); a2=P.a0(t-h);
N1=v1(1)*a1(2)-v1(2)*a1(1); D1=v1(1)^2+v1(2)^2;
N2=v2(1)*a2(2)-v2(2)*a2(1); D2=v2(1)^2+v2(2)^2;
phidd=((N1/D1)-(N2/D2))/(2*h);
end

function [s,sd,sdd] = scale_profile(t,P)
% 平滑缩放脉冲（对应 fullyactBC1210.m 的 aa：8.2-9.2 s 先放大再收回）
w=P.aa_w;
x1=(t-P.aa_t1)/w; x2=(t-P.aa_t2)/w;
S1=1/(1+exp(-x1)); S2=1/(1+exp(-x2));
s=P.aa_max*(S1-S2);
sd=P.aa_max*(S1*(1-S1)-S2*(1-S2))/w;
sdd=P.aa_max*(S1*(1-S1)*(1-2*S1)-S2*(1-S2)*(1-2*S2))/w^2;
end

function [w,wd,wdd] = smooth_ramp(t,P)
% 参考引导权重：五次 smoothstep，1→0，端点一二阶导为零
T=P.Tramp;
if t>=T
    w=0; wd=0; wdd=0; return;
end
u=t/T;
w =1-(6*u^5-15*u^4+10*u^3);
wd =-(30*u^4-60*u^3+30*u^2)/T;
wdd=-(120*u^3-180*u^2+60*u)/T^2;
end

function [A,Ad,Add] = af_all(t,P)
% 仿射矩阵 A_f=(1+aa(t))R(phi_d)，与 fullyactBC1210.m 的 Aff('A',...) 一致
[phi,phid,phidd] = ref_heading(t,P);
phi=P.rot_gain*phi; phid=P.rot_gain*phid; phidd=P.rot_gain*phidd; % 阻尼编队旋转
[s,sd,sdd] = scale_profile(t,P);                        % 缩放脉冲及其导数
R = [cos(phi),-sin(phi); sin(phi),cos(phi)];
Rp = [-sin(phi),-cos(phi); cos(phi),-sin(phi)];         % dR/dphi
Rpp = [-cos(phi),sin(phi); -sin(phi),-cos(phi)];        % d2R/dphi2
A = (1+s)*R;
Ad = sd*R + (1+s)*phid*Rp;
Add = sdd*R + 2*sd*phid*Rp + (1+s)*(phidd*Rp + phid^2*Rpp);
end

% =========================================================================
function x0 = initial_state(P)
% 组装初始状态；状态顺序与文件头说明一致
N=P.N; t0=0;                                               % 车辆数、初始时刻
q0=P.q0(t0); v0=P.v0(t0); [A0,Ad0,~]=af_all(t0,P);         % 领航者初始位置/速度/仿射矩阵
qh=[q0+[0.5;-0.3], q0+[-0.4;0.4], q0+[0.35;0.5]];           % 三车对领航者位置的初始估计（拉开但可达）
vh=zeros(2,N);                                              % 三车对领航者速度的初始估计（全为 0）
e0=P.e0;                                                % 初始编队偏差（取自参数）
d0=zeros(2,N);                                              % 三车初始速度偏差（全为 0）

x0=zeros(P.nS,1);                                          % 预分配总状态向量
x0(1:2*N)=reshape(qh,[],1);                                % 填入位置估计 qh
x0(2*N+1:4*N)=reshape(vh,[],1);                            % 填入速度估计 vh
x0(4*N+1:5*N)=P.g0*ones(N,1);                              % 自适应位置增益初值
x0(5*N+1:6*N)=P.g0*ones(N,1);                              % 自适应速度增益初值
x0(6*N+1:7*N)=P.iota0*ones(N,1);                           % 自适应 sign 增益初值
for i=1:N                                                  % 逐车填 11 维状态
    b=7*N+(i-1)*P.nV;                                      % 第 i 车状态的起始下标偏移
    qr=qh(:,i)+A0*P.poff(:,i);                             % 第 i 车的参考位置 q_{i,r}
    q=qr+e0(:,i);                                          % 第 i 车的初始控制点位置 q
    qd=zeros(2,1);                                         % 初始速度全为 0
    [phi,~,~]=ref_heading(t0,P);                           % 初始航向对齐领航者切线
    X=q(1)-P.xi*cos(phi); Y=q(2)-P.xi*sin(phi);            % 由 q 反推物理位置 X,Y
    vb=[0;0];                                              % 车体系初始速度也为 0
    x0(b+1)=X; x0(b+2)=Y; x0(b+3)=phi;                     % 填 X,Y,phi
    x0(b+4)=vb(1); x0(b+5)=vb(2); x0(b+6)=0;               % 填 vx,vy,r（初始 r=0）
    x0(b+7)=0;                                             % 前轮转角初值 0
    x0(b+8:b+9)=qd;                                        % ESO 的 z1 初值取 qdot(0)
    x0(b+10:b+11)=[0;0];                                   % ESO 的 z2 初值取 0
end
end

% =========================================================================
function dx = rhs(t,x,P)
% 闭环右端：领航者观测器 + BTG 跟踪 + FAS 逆 + 物理限幅 + 车辆 + ESO
N=P.N;                                                     % 车辆数
qh=reshape(x(1:2*N),2,N); vh=reshape(x(2*N+1:4*N),2,N);    % 取出估计器状态
gq=x(4*N+1:5*N); gv=x(5*N+1:6*N); iot=x(6*N+1:7*N);        % 取出逐车自适应增益
q0=P.q0(t); v0=P.v0(t); a0=P.a0(t);                        % 当前领航者位置/速度/加速度
[A,Ad,Add]=af_all(t,P);                                    % 当前仿射矩阵及其导数（一次算完）
[w,wd,wdd]=smooth_ramp(t,P);                                % 参考平滑引导权重

% ---- 领航者观测器（逐车自适应增益）----
muE = btg_mu(t,P.Te,P.epse,P.qe);                          % 观测器 BTG 增益 mu_e(t)
dqh=zeros(2,N); dvh=zeros(2,N);                            % 预分配估计器导数
dgq=zeros(N,1); dgv=zeros(N,1); diot=zeros(N,1);           % 预分配增益自适应律
for i=1:N                                                  % 逐车计算观测器
    eq=zeros(2,1); ev=zeros(2,1);                          % 第 i 车的位置/速度一致性误差
    for j=1:N                                              % 遍历邻居
        if P.A(i,j)~=0                                     % 若 j 是 i 的邻居
            eq=eq+P.A(i,j)*(qh(:,i)-qh(:,j));              % 累加位置邻居差
            ev=ev+P.A(i,j)*(vh(:,i)-vh(:,j));              % 累加速度邻居差
        end
    end
    eq=eq+P.sigma(i)*(qh(:,i)-q0);                         % 加 pinned 位置误差
    ev=ev+P.sigma(i)*(vh(:,i)-v0);                         % 加 pinned 速度误差
    dqh(:,i)=vh(:,i)-gq(i)*muE*eq;                         % 位置观测器：qhat_dot
    dvh(:,i)=-gv(i)*muE*ev-iot(i)*tanh(ev/P.eps_s);        % 速度观测器：自适应增益 + tanh 边界层
    if t < P.Te-P.epse                                     % 冻结前更新增益
        dgq(i)=muE*(P.eg+P.cg*norm(eq)^2);                 % 位置增益自适应律
        dgv(i)=muE*(P.eg+P.cg*norm(ev)^2);                 % 速度增益自适应律
        diot(i)=muE*(P.ei+P.ci*norm(ev));                  % sign 增益自适应律
    end                                                    % 冻结后导数为零
end

% ---- BTG 跟踪参数 ----
mu=btg_mu(t,P.Tc,P.epsi,P.q);                              % 跟踪增益 mu(t)
gamma=P.lamc+btg_mudot(t,P.Tc,P.epsi,P.q)/mu;              % gamma = lamc + mudot/mu
g=P.kc*mu+(P.q/P.Tc)*mu^(1/P.q)+gamma;                     % 合成增益 g(t)

dx=zeros(P.nS,1);                                          % 预分配总导数
dx(1:2*N)=reshape(dqh,[],1);                               % 填估计器位置导数
dx(2*N+1:4*N)=reshape(dvh,[],1);                           % 填估计器速度导数
dx(4*N+1:5*N)=dgq;                                         % 填位置增益导数
dx(5*N+1:6*N)=dgv;                                         % 填速度增益导数
dx(6*N+1:7*N)=diot;                                        % 填 sign 增益导数

for i=1:N                                                  % 逐车计算控制与动力学
    b=7*N+(i-1)*P.nV;                                      % 第 i 车状态起始偏移
    X=x(b+1); Y=x(b+2); phi=x(b+3);                        % 取出 X,Y,phi
    vx=x(b+4); vy=x(b+5); r=x(b+6); delta=x(b+7);          % 取出 vx,vy,r,delta
    z1=x(b+8:b+9); z2=x(b+10:b+11);                        % 取出 ESO 状态
    z1=min(max(z1,-1e3),1e3); z2=min(max(z2,-50),50);       % ESO 状态投影（工程保护）
    cs=cos(phi); sn=sin(phi);                              % 预计算 cos/sin(phi)

    q=[X+P.xi*cs; Y+P.xi*sn];                              % FAS 控制点 q=[X+xi cos, Y+xi sin]
    qd=[vx*cs-vy*sn-P.xi*r*sn; vx*sn+vy*cs+P.xi*r*cs];     % 控制点速度 qdot

    qr=qh(:,i)+A*P.poff(:,i);                              % 第 i 车参考位置 q_{i,r}
    qdr=vh(:,i)+Ad*P.poff(:,i);                            % 参考速度 qdot_{i,r}
    qddr=a0+Add*P.poff(:,i);                               % 参考加速度
    qr=qr+w*P.e0(:,i); qdr=qdr+wd*P.e0(:,i); qddr=qddr+wdd*P.e0(:,i); % 平滑吸收初偏

    eta=q-qr; etad=qd-qdr; E=etad+P.alpha*eta;             % 协同误差 eta、其导数、滤波误差 E

    [B,H]=vehicle_BH(P,X,Y,phi,vx,vy,r,delta);             % 计算 FAS 的 B,H（用实际转角）
    vxs=max(vx,P.vxeps); beta_r=(vy+P.lf*r)/vxs;           % 低速保护后的 vx 与后轴侧偏角

    eddot=-P.alpha*etad-g*E;                               % eta_ddot* = -alpha etad - g E
    uc=safe_solve(B,eddot+qddr-H-z2);                      % FAS 逆：含 z2 前馈
    uc(1)=2*P.Fxmax*tanh(uc(1)/(2*P.Fxmax));            % 纵向指令上限（防止 ESO 积分饱和）
    uc(2)=2*P.Cf*P.dmax*tanh(uc(2)/(2*P.Cf*P.dmax));    % 侧向指令上限

    % ---- 转向前置保护：把指令转角投影到 0.9*dmax ----
    dsafe=0.9*P.dmax;
    dcg=uc(2)/P.Cf+beta_r;
    dcg=dsafe*tanh(dcg/dsafe);
    uc(2)=P.Cf*(dcg-beta_r);

    % ---- 物理执行器饱和（tanh 光滑限幅）----
    dc=uc(2)/P.Cf+beta_r;                                  % 指令对应的转角
    Fxa=P.Fxmax*tanh(uc(1)/P.Fxmax);                       % 实际纵向力（光滑饱和）
    da=P.dmax*tanh(dc/P.dmax);                             % 实际转角（光滑饱和）
    Fyfa=P.Cf*(da-beta_r);                                 % 由实际转角反算实际 Fyf
    ua=[Fxa;Fyfa]; Du=ua-uc;                               % 实际输入与失配
    dvx_d=0.35*sin(0.8*t+0.4*(i-1));                       % 未知纵向扰动（控制器不知道）
    dvy_d=0.25*cos(0.6*t+0.3*(i-1));                       % 未知侧向扰动（控制器不知道）
    dq_d=[dvx_d*cs-dvy_d*sn; dvx_d*sn+dvy_d*cs];           % 旋转到 FAS 坐标的未知扰动
    Lam=B*Du+dq_d;                                         % 真实集总失配 = 饱和失配 + 未知扰动

    % ---- 车辆六状态 ----
    Fyf=P.Cf*(delta-beta_r);                               % 用实际转角算侧向力
    Fyr=-P.Cr*(vy-P.lr*r)/vxs;                             % 后轴侧向力
    cd=cos(delta); sd=sin(delta);                          % 预计算 cos/sin(delta)
    dx(b+1)=vx*cs-vy*sn;                                   % Xdot
    dx(b+2)=vx*sn+vy*cs;                                   % Ydot
    dx(b+3)=r;                                             % phidot
    dx(b+4)=(Fxa*((1-P.Kg)+P.Kg*cd)-Fyf*sd)/P.m+vy*r+dvx_d;% vxdot（含未知纵向扰动）
    dx(b+5)=(P.Kg*Fxa*sd+Fyf*cd+Fyr)/P.m-vx*r+dvy_d;       % vydot（含未知侧向扰动）
    dx(b+6)=(P.lf*(P.Kg*Fxa*sd+Fyf*cd)-P.lr*Fyr)/P.Jz;     % rdot
    dx(b+7)=(da-delta)/P.taud;                             % 转角一阶滤波

    % ---- ESO（本文：rho 在 T_a 后冻结）----
    if t<P.Ta                                              % 冻结前
        rho=1/(P.Ta-t+P.b);                                % rho(t)=1/(Ta-t+b)
    else                                                   % 冻结后
        rho=1/P.b;                                         % rho=1/b
    end
    dx(b+8:b+9)=z2+P.k1*rho*(qd-z1)+B*uc+H;                % z1dot
    dx(b+10:b+11)=P.k2*rho^2*(qd-z1);                      % z2dot
end
end

% =========================================================================
function [B,H]=vehicle_BH(P,X,Y,phi,vx,vy,r,delta)
% 由六状态算 FAS 的输入矩阵 B 与非线性项 H（含论文修正后的 +Gamma sin 项）
vxs=max(vx,P.vxeps);                                       % 低速保护
Fyr=-P.Cr*(vy-P.lr*r)/vxs;                                 % 后轴侧向力
Gamma=P.xi*P.lr*Fyr/P.Jz;                                  % Gamma 项
cs=cos(phi); sn=sin(phi); cd=cos(delta); sd=sin(delta);    % 预计算三角函数
B=(1/P.m)*[(1-P.Kg)*cs+P.Kg*cos(delta+phi), -sin(delta+phi); % B 第一行第一项/第二项
           (1-P.Kg)*sn+P.Kg*sin(delta+phi),  cos(delta+phi)];% B 第二行第一项/第二项
B=B+(P.xi*P.lf/P.Jz)*[-sn*P.Kg*sd, -sn*cd;                 % 叠加 xi*lf/Jz 修正项（第一行）
                        cs*P.Kg*sd,  cs*cd];               % 叠加 xi*lf/Jz 修正项（第二行）
H=[-Fyr*sn/P.m+Gamma*sn-P.xi*r^2*cs;                       % H 第一行
    Fyr*cs/P.m-Gamma*cs-P.xi*r^2*sn];                      % H 第二行
end

% =========================================================================
function mu=btg_mu(t,T,eps,q)
% BTG 增益：t<T-eps 时发散，之后冻结为 (T/eps)^q
if t<T-eps                                                % 冻结前
    mu=(T/(T-t))^q;                                        % mu(t)
else                                                       % 冻结后
    mu=(T/eps)^q;                                          % mu_max
end
end

function gmm=btg_mudot(t,T,eps,q)
% BTG 增益的对数导数 mudot/mu（冻结前 0.3 s 内平滑降到 0）
tf=T-eps; d=0.3;                                           % 冻结时刻与过渡宽度
if t<tf-d                                                  % 远未冻结
    gmm=q/(T-t);                                           % q/(T-t)
elseif t<tf                                                % 平滑过渡
    gmm=q/(T-t)*max(0,(tf-t)/d);                           % 线性降到 0
else                                                       % 冻结后
    gmm=0;                                                 % 导数为 0
end
end

% =========================================================================
function R=unpack(X,t,P)
% 从轨迹重新算指标；注意这里复算的控制量必须与 rhs 完全一致
N=P.N; n=numel(t);                                        % 车辆数、时间点数
R.t=t;                                                    % 存时间
% ---- 轨迹与参考 ----
R.q=zeros(2,N,n); R.qr=zeros(2,N,n); R.eta=zeros(2,N,n);   % q、参考、位置误差
R.qd=zeros(2,N,n); R.qdr=zeros(2,N,n); R.etad=zeros(2,N,n);% 速度、参考速度、速度误差
R.E=zeros(2,N,n);                                          % 滤波误差
% ---- 领航者估计与估计误差 ----
R.qh=zeros(2,N,n); R.vh=zeros(2,N,n);                      % 三车对领航者的估计
R.q0=zeros(2,n); R.v0=zeros(2,n);                          % 真实领航者状态
R.est_q=zeros(N,n); R.est_v=zeros(N,n);                    % 逐车位置/速度估计误差范数
R.est_q_max=zeros(1,n); R.est_v_max=zeros(1,n);            % 三车最大估计误差
R.gq=zeros(N,n); R.gv=zeros(N,n); R.iota=zeros(N,n);       % 自适应增益轨迹
% ---- ESO 相关 ----
R.Lam=zeros(2,N,n); R.z1=zeros(2,N,n); R.z2=zeros(2,N,n);  % 真失配、ESO 状态
R.e2=zeros(2,N,n);                                         % z2-Lam
R.eso_err=zeros(N,n); R.eso_err_max=zeros(1,n);            % 逐车/最大 ESO 残差
R.z1_err=zeros(N,n); R.z1_err_max=zeros(1,n);              % qd-z1 逐车/最大
% ---- 输入与饱和 ----
R.uc=zeros(2,N,n); R.ua=zeros(2,N,n); R.Du=zeros(2,N,n);   % 指令/实际/失配
R.uc_raw=zeros(2,N,n);                                     % 限幅前的原始指令（FAS 逆输出）
R.fx=zeros(1,N,n); R.delta=zeros(1,N,n); R.dcmd=zeros(N,n);% 实际 Fx、实际/指令转角
R.dcmd_raw=zeros(N,n);                                     % 限幅前的原始指令转角
R.satF=zeros(N,n); R.satD=zeros(N,n);                      % 逐车饱和标志
R.sat=zeros(1,n);                                          % 任一通道饱和
% ---- 跟踪误差 ----
R.track_ref=zeros(N,n); R.track_true=zeros(N,n);           % 相对估计/真实参考的位置误差
R.track_v=zeros(N,n); R.E_err=zeros(N,n);                  % 速度误差、滤波误差范数
R.form=zeros(1,n); R.lamnorm=zeros(1,n);                   % 最大真实编队误差、最大失配范数
R.e2norm=zeros(1,n); R.z2norm=zeros(1,n);                  % 最大残差/估计范数
R.err_obs=zeros(1,n);                                      % 最大领航者估计误差
for k=1:n                                                  % 遍历每个时刻
    tk=t(k); x=X(:,k);                                     % 当前时刻与状态
    q0=P.q0(tk); v0=P.v0(tk); [A,Ad,Add]=af_all(tk,P);     % 当前领航者/仿射矩阵（一次算完）
    [w,wd,wdd]=smooth_ramp(tk,P);                           % 参考平滑引导权重
    qh=reshape(x(1:2*N),2,N);                              % 估计器位置
    vh=reshape(x(2*N+1:4*N),2,N);                          % 估计器速度
    gq=x(4*N+1:5*N); gv=x(5*N+1:6*N); iot=x(6*N+1:7*N);    % 自适应增益
    R.q0(:,k)=q0; R.v0(:,k)=v0;                            % 存领航者位置与速度
    R.gq(:,k)=gq; R.gv(:,k)=gv; R.iota(:,k)=iot;           % 存逐车增益轨迹
    for i=1:N                                              % 逐车
        b=7*N+(i-1)*P.nV;                                  % 状态偏移
        Xi=x(b+1); Yi=x(b+2); phi=x(b+3);                  % 取 X,Y,phi
        vx=x(b+4); vy=x(b+5); r=x(b+6); delta=x(b+7);      % 取 vx,vy,r,delta
        z1=x(b+8:b+9); z2=x(b+10:b+11);                    % 取 z1,z2
        z1=min(max(z1,-1e3),1e3); z2=min(max(z2,-50),50);   % 同步状态投影
        cs=cos(phi); sn=sin(phi);                          % 三角函数
        q=[Xi+P.xi*cs; Yi+P.xi*sn];                        % 控制点位置
        qd=[vx*cs-vy*sn-P.xi*r*sn; vx*sn+vy*cs+P.xi*r*cs]; % 控制点速度
        qr=qh(:,i)+A*P.poff(:,i);                          % 参考位置
        qr=qr+w*P.e0(:,i);                                 % 平滑吸收初偏
        % ---- 轨迹、参考、跟踪误差 ----
        R.q(:,i,k)=q; R.qr(:,i,k)=qr; R.eta(:,i,k)=q-qr;   % 存 q/参考/位置误差
        R.track_ref(i,k)=norm(q-qr);                       % 相对估计参考的跟踪误差
        R.track_true(i,k)=norm(q-(q0+A*P.poff(:,i)));      % 相对真实领航者的跟踪误差
        % ---- 分布式估计误差 ----
        R.qh(:,i,k)=qh(:,i); R.vh(:,i,k)=vh(:,i);          % 存三车估计
        R.est_q(i,k)=norm(qh(:,i)-q0);                     % 位置估计误差
        R.est_v(i,k)=norm(vh(:,i)-v0);                     % 速度估计误差
        vxs=max(vx,P.vxeps); beta_r=(vy+P.lf*r)/vxs;       % 低速保护与侧偏角
        [B,H]=vehicle_BH(P,Xi,Yi,phi,vx,vy,r,delta);       % 复算 B,H
        % 复算控制量（用于统计）
        mu=btg_mu(tk,P.Tc,P.epsi,P.q);                     % 复算 mu
        gamma=P.lamc+btg_mudot(tk,P.Tc,P.epsi,P.q)/mu;     % 复算 gamma
        g=P.kc*mu+(P.q/P.Tc)*mu^(1/P.q)+gamma;             % 复算 g
        qdr=vh(:,i)+Ad*P.poff(:,i);                        % 复算参考速度
        qddr=P.a0(tk)+Add*P.poff(:,i);                     % 参考加速度
        qdr=qdr+wd*P.e0(:,i); qddr=qddr+wdd*P.e0(:,i);     % 平滑吸收初偏
        etad=qd-qdr; E=etad+P.alpha*(q-qr);                % 复算 etad 与 E
        R.qd(:,i,k)=qd; R.qdr(:,i,k)=qdr; R.etad(:,i,k)=etad; R.E(:,i,k)=E; % 存速度/参考/误差
        R.track_v(i,k)=norm(etad); R.E_err(i,k)=norm(E);   % 速度跟踪误差、滤波误差范数
        eddot=-P.alpha*etad-g*E;                           % 复算 eta_ddot*
        uc=safe_solve(B,eddot+qddr-H-z2);                  % 复算 FAS 逆
        R.uc_raw(:,i,k)=uc;                                % 存限幅前的原始指令
        uc(1)=2*P.Fxmax*tanh(uc(1)/(2*P.Fxmax));        % 指令上限（与 rhs 一致）
        uc(2)=2*P.Cf*P.dmax*tanh(uc(2)/(2*P.Cf*P.dmax));% 指令上限（与 rhs 一致）
        dsafe=0.9*P.dmax;                                  % 转向前置保护
        dcg=uc(2)/P.Cf+beta_r;
        dcg=dsafe*tanh(dcg/dsafe);
        uc(2)=P.Cf*(dcg-beta_r);
        dc=uc(2)/P.Cf+beta_r;                              % 转角
        Fxa=P.Fxmax*tanh(uc(1)/P.Fxmax);                   % 实际 Fx
        da=P.dmax*tanh(dc/P.dmax);                         % 实际转角
        Fyfa=P.Cf*(da-beta_r);                             % 实际 Fyf
        ua=[Fxa;Fyfa]; Du=ua-uc;                           % 实际输入与失配
        dvx_d=0.35*sin(0.8*tk+0.4*(i-1));                  % 复算未知纵向扰动
        dvy_d=0.25*cos(0.6*tk+0.3*(i-1));                  % 复算未知侧向扰动
        dq_d=[dvx_d*cs-dvy_d*sn; dvx_d*sn+dvy_d*cs];       % 旋转到 FAS 坐标
        Lam=B*Du+dq_d;                                     % 真实集总失配
        R.uc(:,i,k)=uc; R.ua(:,i,k)=ua; R.Du(:,i,k)=Du;    % 存指令/实际/失配
        R.Lam(:,i,k)=Lam; R.z1(:,i,k)=z1; R.z2(:,i,k)=z2;  % 存真失配与 ESO 状态
        R.e2(:,i,k)=z2-Lam;                                % ESO 残差
        R.eso_err(i,k)=norm(z2-Lam);                       % 逐车 ESO 残差范数
        R.z1_err(i,k)=norm(qd-z1);                         % ESO 速度重构误差
        R.fx(:,i,k)=Fxa; R.delta(:,i,k)=delta;             % 存实际 Fx 与转角
        R.dcmd(i,k)=dc;                                    % 存指令转角
        R.dcmd_raw(i,k)=R.uc_raw(2,i,k)/P.Cf+beta_r;       % 存限幅前的原始指令转角
        R.phi(i,k)=phi; R.beta(i,k)=beta_r;                % 存航向与 beta_r
        R.satF(i,k)=abs(uc(1))>0.8*P.Fxmax;                % 纵向力饱和标志（80% 阈值）
        R.satD(i,k)=abs(dc)>0.8*P.dmax;                    % 转角饱和标志
    end
    % ---- 汇总本时刻的最大值 ----
    R.form(k)=max(R.track_true(:,k));                      % 最大真实跟踪误差
    R.lamnorm(k)=max(vecnorm(R.Lam(:,:,k),2,1));           % 最大失配范数
    R.est_q_max(k)=max(R.est_q(:,k));                      % 最大位置估计误差
    R.est_v_max(k)=max(R.est_v(:,k));                      % 最大速度估计误差
    R.eso_err_max(k)=max(R.eso_err(:,k));                  % 最大 ESO 残差
    R.z1_err_max(k)=max(R.z1_err(:,k));                    % 最大速度重构误差
    R.z2norm(k)=max(vecnorm(R.z2(:,:,k),2,1));             % 最大 z2 范数
    R.e2norm(k)=R.eso_err_max(k);                          % 兼容旧字段
    R.err_obs(k)=R.est_q_max(k);                           % 兼容旧字段
    R.sat(k)=any(R.satF(:,k) | R.satD(:,k));               % 任一通道饱和
end
end

% =========================================================================
function print_metrics(R)
% 打印终端指标
last=find(R.sat>0,1,'last');                               % 最后一次出现饱和的时刻索引
if isempty(last)                                           % 如果从未饱和
    tlast=0;                                               % 记 0
else                                                       % 否则
    tlast=R.t(last);                                       % 记最后饱和时刻
end
idx=find(R.t>tlast & R.form<0.01,1,'first');               % 最后饱和后首次误差<1cm
if isempty(idx), rec=NaN; else rec=R.t(idx); end           % 恢复时刻（无则 NaN）
fprintf('form: max %.4f m, final %.4f m; sat fraction %.1f%%, last sat %.3f s, recovery %.3f s\n', ... % 编队误差/饱和摘要
    max(R.form),R.form(end),100*mean(R.sat),tlast,rec-tlast); % 对应数值
fprintf('max ||Lam|| %.3f m/s^2; max |Fx| %.1f N; max |delta| %.3f rad\n', ... % 失配/输入摘要
    max(R.lamnorm),max(abs(R.fx),[],'all'),max(abs(R.delta),[],'all')); % 对应数值
fprintf('distributed estimation: |q| max/final %.3e/%.3e m; |v| max/final %.3e/%.3e m/s\n', ... % 分布式估计误差
    max(R.est_q_max),R.est_q_max(end),max(R.est_v_max),R.est_v_max(end)); % 位置/速度
fprintf('ESO residual: at Ta %.3e, max %.3e, final %.3e m/s^2; |qd-z1| max %.3e\n', ... % ESO 残差
    interp1(R.t,R.eso_err_max,0.5),max(R.eso_err_max),R.eso_err_max(end),max(R.z1_err_max)); % 残差与速度重构
fprintf('adaptive observer gains (final): gq = [%.2f %.2f %.2f], gv = [%.2f %.2f %.2f], iota = [%.2f %.2f %.2f]\n', ...
    R.gq(:,end),R.gv(:,end),R.iota(:,end));                % 冻结后的自适应增益
fprintf('dynamics tracking: true max/final %.3e/%.3e m; ref max/final %.3e/%.3e m; |E| max %.3e\n', ... % 动力学跟踪
    max(R.track_true,[],'all'),max(R.track_true(:,end)), ... % 真实参考误差
    max(R.track_ref,[],'all'),max(R.track_ref(:,end)), ...   % 估计参考误差
    max(R.E_err,[],'all'));                                  % 滤波误差
fprintf('max |uc| = %.3e\n',max(abs(R.uc),[],'all'));      % 最大指令
fprintf('control: pre-limit max |Fx| %.3e N, |Fyf| %.3e N, |delta| %.3f rad\n', ... % 限幅前原始指令
    max(abs(R.uc_raw(1,:,:)),[],'all'),max(abs(R.uc_raw(2,:,:)),[],'all'),max(abs(R.dcmd_raw),[],'all')); % 三通道
fprintf('control: applied  max |Fx| %.3f N, |Fyf| %.3e N, |delta| %.3f rad\n', ... % 实际施加输入
    max(abs(R.fx),[],'all'),max(abs(R.ua(2,:,:)),[],'all'),max(abs(R.delta),[],'all')); % 三通道
end

% =========================================================================
function plot_trajectory(R,P)
% 三车全六状态轨迹图（车头用真实航向）
tt=R.t;
figure('Color','w','Position',[100 100 1000 700]); hold on; box on; grid off;
c1=[0.99 0.49 0.00]; c2=[0.17 0.62 0.47]; c3=[0.88 0.17 0.56];
h0=plot(R.q0(1,:),R.q0(2,:),'--','Color','k','LineWidth',3);
h1=plot(squeeze(R.q(1,1,:)),squeeze(R.q(2,1,:)),'-','Color',c1,'LineWidth',2);
h2=plot(squeeze(R.q(1,2,:)),squeeze(R.q(2,2,:)),'-','Color',c2,'LineWidth',2);
h3=plot(squeeze(R.q(1,3,:)),squeeze(R.q(2,3,:)),'-','Color',c3,'LineWidth',2);
tss=[0 .25 .5 .75 1]*P.T;                            % 快照 0,25,50,75,96 s
hqdel=[];                                                  % 前轮转角箭头句柄
for tq=tss
    [~,k]=min(abs(tt-tq)); Q=R.q(:,:,k); lq=R.q0(:,k);
    plot([Q(1,1) Q(1,2) Q(1,3) Q(1,1)],[Q(2,1) Q(2,2) Q(2,3) Q(2,1)], ...
        ':','Color','k','LineWidth',1.5);
    for i=1:3
        psi=R.phi(i,k);                                % 车体真实航向
        delta=R.delta(1,i,k);                          % 实际前轮转角
        Carplot_Ezekiel(Q(1,i),Q(2,i),psi,delta,'r','k',1.3);
        fx=Q(1,i)+1.3*P.lf*cos(psi); fy=Q(2,i)+1.3*P.lf*sin(psi); % 前轴位置
        hqi=quiver(fx,fy,5.2*cos(psi+delta),5.2*sin(psi+delta),0, ...
            'Color','b','LineWidth',2.0,'MaxHeadSize',1.0,'AutoScale','off');
        if isempty(hqdel), hqdel=hqi; end               % 记首个句柄供图例
    end
    text(lq(1),lq(2)+0.5,sprintf('#t=%.1f s',tq),'FontSize',12, ...
        'HorizontalAlignment','center','FontName','Times');
end
axis equal;
allx=[R.q0(1,:), reshape(R.q(1,:,:),1,[])];
ally=[R.q0(2,:), reshape(R.q(2,:,:),1,[])];
xlim([min(allx)-20 max(allx)+20]); ylim([min(ally)-10 max(ally)+10]);
set(gca,'FontSize',15); xlabel('x (m)'); ylabel('y (m)');
legend([h0 h1 h2 h3 hqdel],{'Virtual leader','Vehicle 1','Vehicle 2','Vehicle 3', ...
    'Front-wheel angle \delta_f'},'Location','northeast');
exportgraphics(gcf,'formation_sixstate_lissajous_sat.png','Resolution',200);
fprintf('saved: formation_sixstate_lissajous_sat.png\n');
end

% =========================================================================
function plot_estimation(R,P)
% 分布式观测器估计误差：左列 0-3 s 瞬态，右列全时段 RMS 包络
cc=[0.99 0.49 0.00;0.17 0.62 0.47;0.88 0.17 0.56];         % 三车配色
w=max(1,round(0.02/P.dt));                                 % 20 ms 滑动 RMS 窗
tw=transient_window(R,P); i1=R.t<=tw;                                                 % 瞬态区间索引
tt=R.t(i1);
figure('Color','w','Position',[50 50 1150 760]);
sgtitle('Distributed leader-observer estimation errors','FontSize',12,'FontWeight','bold');

subplot(4,1,1);                                            % 瞬态：最大误差
plot(tt,R.est_q_max(i1)+eps,'b-','LineWidth',1.4); hold on;
plot(tt,R.est_v_max(i1)+eps,'r--','LineWidth',1.4);
xline(P.Te,'k:','T_e'); grid on; xlabel('t [s]'); ylabel('max_i norm');
legend({'$\|\hat{q}_i-q_0\|$ [m]','$\|\hat{v}_i-v_0\|$ [m/s]'},'Location','best','Interpreter','latex');
title(sprintf('(a) Transient (0-%g s)',tw)); xlim([0 tw]);

subplot(4,1,2);                                            % 全时段：RMS 包络
plot(R.t,rms_env(R.est_q_max,w)+eps,'b-','LineWidth',1.3); hold on;
plot(R.t,rms_env(R.est_v_max,w)+eps,'r--','LineWidth',1.3);
xline(P.Te,'k:','T_e'); grid on; xlabel('t [s]'); ylabel('RMS envelope');
legend({'position [m]','velocity [m/s]'},'Location','best');
title('(b) Envelope (0-96 s)');

subplot(4,1,3);                                            % 逐车位置误差
for i=1:P.N, plot(tt,R.est_q(i,i1)+eps,'Color',cc(i,:),'LineWidth',1.3); hold on; end
xline(P.Te,'k:','T_e'); grid on; xlabel('t [s]'); ylabel('$\|\hat{q}_i-q_0\|$ [m]','Interpreter','latex');
legend({'veh 1','veh 2','veh 3'},'Location','best'); title('(c) Position, per vehicle');
xlim([0 tw]);

subplot(4,1,4);                                            % 逐车速度误差
for i=1:P.N, plot(tt,R.est_v(i,i1)+eps,'Color',cc(i,:),'LineWidth',1.3); hold on; end
xline(P.Te,'k:','T_e'); grid on; xlabel('t [s]'); ylabel('$\|\hat{v}_i-v_0\|$ [m/s]','Interpreter','latex');
legend({'veh 1','veh 2','veh 3'},'Location','best'); title('(d) Velocity, per vehicle');
xlim([0 tw]);

exportgraphics(gcf,'formation_lissajous_estimation.png','Resolution',200);
fprintf('saved: formation_lissajous_estimation.png\n');
end

% =========================================================================
function plot_tracking(R,P)
% 六状态动力学跟踪误差：左列 0-3 s 瞬态，右列全时段 RMS 包络
cc=[0.99 0.49 0.00;0.17 0.62 0.47;0.88 0.17 0.56];         % 三车配色
w=max(1,round(0.02/P.dt));                                 % 20 ms 滑动 RMS 窗
tw=transient_window(R,P); i1=R.t<=tw; tt=R.t(i1);                                     % 瞬态区间
figure('Color','w','Position',[50 50 1150 760]);
sgtitle('Six-state dynamics tracking errors','FontSize',12,'FontWeight','bold');

subplot(4,1,1);                                            % 瞬态：最大真实跟踪误差
plot(tt,max(R.track_true(:,i1),[],1)+eps,'b-','LineWidth',1.4); hold on;
xline(P.Tc,'k:','T_c'); grid on; xlabel('t [s]'); ylabel('max_i||q_i-q_0-Ap_i|| [m]');
title(sprintf('(a) Transient (0-%g s)',tw)); xlim([0 tw]);

subplot(4,1,2);                                            % 全时段：RMS 包络
plot(R.t,rms_env(max(R.track_true,[],1),w)+eps,'b-','LineWidth',1.3); hold on;
xline(P.Tc,'k:','T_c'); grid on; xlabel('t [s]'); ylabel('RMS envelope [m]');
title('(b) Envelope (0-96 s)');

subplot(4,1,3);                                            % 逐车位置跟踪误差
for i=1:P.N, plot(tt,R.track_true(i,i1)+eps,'Color',cc(i,:),'LineWidth',1.3); hold on; end
xline(P.Tc,'k:','T_c'); grid on; xlabel('t [s]'); ylabel('||q_i-q_0-Ap_i|| [m]');
legend({'veh 1','veh 2','veh 3'},'Location','best'); title('(c) Position, per vehicle');
xlim([0 tw]);

subplot(4,1,4);                                            % 逐车速度跟踪误差
for i=1:P.N, plot(tt,R.track_v(i,i1)+eps,'Color',cc(i,:),'LineWidth',1.3); hold on; end
xline(P.Tc,'k:','T_c'); grid on; xlabel('t [s]'); ylabel('$\|\dot{q}_i-\dot{q}_{i,r}\|$ [m/s]','Interpreter','latex');
legend({'veh 1','veh 2','veh 3'},'Location','best'); title('(d) Velocity, per vehicle');
xlim([0 tw]);

exportgraphics(gcf,'formation_lissajous_tracking.png','Resolution',200);
fprintf('saved: formation_lissajous_tracking.png\n');
end

% =========================================================================
function plot_eso(R,P)
% ESO 失配重构：左列 0-3 s 重构过程，右列残差瞬态与全时段包络
cc=[0.99 0.49 0.00;0.17 0.62 0.47;0.88 0.17 0.56];         % 三车配色
w=max(1,round(0.02/P.dt));                                 % 20 ms 滑动 RMS 窗
tw=transient_window(R,P); i1=R.t<=tw; tt=R.t(i1);                                     % 瞬态区间
figure('Color','w','Position',[50 50 1150 760]);
sgtitle('Prescribed-time ESO mismatch reconstruction','FontSize',12,'FontWeight','bold');

subplot(4,1,1);                                            % 1 号车 X 分量
plot(tt,squeeze(R.Lam(1,1,i1)),'k-','LineWidth',1.2); hold on;
plot(tt,squeeze(R.z2(1,1,i1)),'b--','LineWidth',1.2);
xline(P.Ta,'k:','T_a'); grid on; xlabel('t [s]'); ylabel('\Lambda_x, z_{2,x} [m/s^2]');
legend({'\Lambda_{1,x}','z_{2,1,x}'},'Location','best'); title('(a) Vehicle 1, X component');
xlim([0 tw]);

subplot(4,1,2);                                            % 1 号车 Y 分量
plot(tt,squeeze(R.Lam(2,1,i1)),'k-','LineWidth',1.2); hold on;
plot(tt,squeeze(R.z2(2,1,i1)),'b--','LineWidth',1.2);
xline(P.Ta,'k:','T_a'); grid on; xlabel('t [s]'); ylabel('\Lambda_y, z_{2,y} [m/s^2]');
legend({'\Lambda_{1,y}','z_{2,1,y}'},'Location','best'); title('(b) Vehicle 1, Y component');
xlim([0 tw]);

subplot(4,1,3);                                            % 残差瞬态
for i=1:P.N, plot(tt,R.eso_err(i,i1)+eps,'Color',cc(i,:),'LineWidth',1.3); hold on; end
xline(P.Ta,'k:','T_a'); grid on; xlabel('t [s]'); ylabel('||z_{2,i}-\Lambda_i|| [m/s^2]');
legend({'veh 1','veh 2','veh 3'},'Location','best'); title(sprintf('(c) Residual (0-%g s)',tw));
xlim([0 tw]);

subplot(4,1,4);                                            % 残差全时段包络
plot(R.t,rms_env(R.eso_err_max,w)+eps,'b-','LineWidth',1.3); hold on;
plot(R.t,rms_env(R.z1_err_max,w)+eps,'r--','LineWidth',1.3);
xline(P.Ta,'k:','T_a'); grid on; xlabel('t [s]'); ylabel('RMS envelope');
legend({'$\|z_2-\Lambda\|$ [m/s$^2$]','$\|\dot{q}-z_1\|$ [m/s]'},'Location','best','Interpreter','latex');
title('(d) Envelope (0-96 s)');

exportgraphics(gcf,'formation_lissajous_eso.png','Resolution',200);
fprintf('saved: formation_lissajous_eso.png\n');
end

% =========================================================================
function plot_inputs(R,P)
% 指令与实际输入（0-3 s）及饱和占比
tw=transient_window(R,P); i1=R.t<=tw; tt=R.t(i1);                                     % 瞬态区间
fxc=squeeze(max(abs(R.uc(1,:,i1)),[],2))';                 % 三车最大指令
fxa=squeeze(max(abs(R.ua(1,:,i1)),[],2))';                 % 三车最大实际
fyc=squeeze(max(abs(R.uc(2,:,i1)),[],2))';                 % 侧向最大指令
fya=squeeze(max(abs(R.ua(2,:,i1)),[],2))';                 % 侧向最大实际
dcmd=max(abs(R.dcmd(:,i1)),[],1);                          % 最大指令转角（1xk）
dact=squeeze(max(abs(R.delta(:,:,i1)),[],2))';             % 最大实际转角（1xk）
figure('Color','w','Position',[50 50 1150 760]);
sgtitle(sprintf('Commanded and applied inputs (0-%g s)',tw),'FontSize',12,'FontWeight','bold');

subplot(4,1,1);                                            % Fx
plot(tt,fxc/P.Fxmax,'b-','LineWidth',1.3); hold on;
plot(tt,fxa/P.Fxmax,'r--','LineWidth',1.2);
grid on; xlabel('t [s]'); ylabel('max_i |\cdot|/F_{x,max}');
legend({'u_c','u_a'},'Location','best'); title('(a) Driving force'); xlim([0 tw]);

subplot(4,1,2);                                            % Fyf
plot(tt,fyc/(P.Cf*P.dmax),'b-','LineWidth',1.3); hold on;
plot(tt,fya/(P.Cf*P.dmax),'r--','LineWidth',1.2);
grid on; xlabel('t [s]'); ylabel('max_i |\cdot|/(C_f\delta_{max})');
legend({'u_c','u_a'},'Location','best'); title('(b) Lateral force'); xlim([0 tw]);

subplot(4,1,3);                                            % 转角
plot(tt,dcmd/P.dmax,'b-','LineWidth',1.3); hold on;
plot(tt,dact/P.dmax,'r--','LineWidth',1.2);
grid on; xlabel('t [s]'); ylabel('max_i |\delta|/\delta_{max}');
legend({'\delta_c','\delta_a'},'Location','best'); title('(c) Steering angle'); xlim([0 tw]);

subplot(4,1,4);                                            % 饱和标志
plot(tt,100*max(R.satF(:,i1),[],1),'b-','LineWidth',1.2); hold on;
plot(tt,100*max(R.satD(:,i1),[],1),'r--','LineWidth',1.2);
grid on; xlabel('t [s]'); ylabel('percent');
legend({'F_x near limit','\delta near limit'},'Location','best');
title(sprintf('(d) Saturation flags (total %.1f%%/%.1f%%)', ...
    100*mean(R.satF(:)),100*mean(R.satD(:)))); xlim([0 tw]);

exportgraphics(gcf,'formation_lissajous_inputs.png','Resolution',200);
fprintf('saved: formation_lissajous_inputs.png\n');
end

% =========================================================================
function plot_control_compare(R,P)
% 限幅前原始指令与实际施加输入的对比（物理量纲，0-3 s，逐车）
N=P.N; tw=min(3,P.T); i1=R.t<=tw; tt=R.t(i1);              % 只看限幅发生的瞬态段
figure('Color','w','Position',[40 40 1500 900]);
tiledlayout(3,3,'Padding','compact','TileSpacing','compact');

for ch=1:3                                                % 三个通道：Fx、Fyf、delta
    for i=1:N                                             % 逐车
        nexttile;
        switch ch                                         % 取出该通道的原始/实际曲线
            case 1
                yraw=squeeze(R.uc_raw(1,i,i1));           % 限幅前纵向力
                yact=squeeze(R.ua(1,i,i1));               % 实际纵向力
                ylab='F_x [N]';
            case 2
                yraw=squeeze(R.uc_raw(2,i,i1));           % 限幅前侧向力
                yact=squeeze(R.ua(2,i,i1));               % 实际侧向力
                ylab='F_{yf} [N]';
            case 3
                yraw=squeeze(R.dcmd_raw(i,i1));           % 限幅前指令转角
                yact=squeeze(R.delta(1,i,i1));            % 实际前轮转角
                ylab='\delta_f [rad]';
        end
        plot(tt,yraw,'b-','LineWidth',1.2); hold on;      % 限幅前（细实线）
        plot(tt,yact,'r--','LineWidth',1.5);              % 实际施加（粗虚线）
        grid on; xlabel('t [s]'); ylabel(ylab);
        if ch==1 && i==1
            legend({'pre-limit u_c','applied u_a'},'Location','best');
        end
        title(sprintf('Vehicle %d',i));
        xlim([0 tw]);
    end
end
sgtitle('Pre-limit command vs applied input (physical units)');
exportgraphics(gcf,'formation_control_compare.png','Resolution',200);
fprintf('saved: formation_control_compare.png\n');
end

% =========================================================================
function tw=transient_window(R,P)
% 自适应瞬态窗口：真实编队误差首次进入 1 cm 后留 50% 余量
idx=find(R.form<1e-2,1,'first');
if isempty(idx), ts=P.T; else ts=R.t(idx); end
tw=min(P.T,max(3,1.5*ts));
end

% =========================================================================
function y=rms_env(x,w)
% 滑动 RMS 包络：消除范数穿零造成的对数十倍尖坑
y=sqrt(movmean(x(:).'.^2,w));
end

% =========================================================================
function save_metrics(R,P)
% 保存所有采集信号到 .mat（数据已按 P.save_stride 记录）
n=numel(R.t); idx=1:n;                                     % 不再二次抽样
S=struct();                                               % 保存结构体
fn=fieldnames(R);                                         % 遍历 R 的字段
for ii=1:numel(fn)
    v=R.(fn{ii});                                         % 当前字段
    if ~isnumeric(v), continue; end                       % 非数值跳过
    if isvector(v) && numel(v)==n                         % 时间序列（行/列）
        S.(fn{ii})=v(idx);
    elseif ndims(v)==3 && size(v,3)==n                    % 2xN xn 数组
        S.(fn{ii})=v(:,:,idx);
    elseif ismatrix(v) && size(v,2)==n                    % Nxn 数组
        S.(fn{ii})=v(:,idx);
    else
        S.(fn{ii})=v;                                     % 标量或其他
    end
end
Q=rmfield(P,{'q0','v0','a0'});                            % 去掉函数句柄
S.params=Q;                                               % 存参数
save('formation_lissajous_metrics.mat','S','-v7.3');      % 保存
fprintf('saved: formation_lissajous_metrics.mat (%d samples)\n',numel(idx));
end

% =========================================================================
function [t,X]=rk4(f,tspan,x0,stride)
% 定步长 RK4：tspan 为 0:dt:T；只按 stride 记录状态以省内存/后处理
if nargin<4 || isempty(stride), stride=1; end
h=tspan(2)-tspan(1); n=numel(tspan);                       % 步长与总步数
t=tspan(1:stride:n); m=numel(t);                           % 记录时间点
X=zeros(numel(x0),m); X(:,1)=x0;                           % 预分配输出
j=1; xk=x0;                                                % 记录指针与当前状态
for k=1:n-1                                                % 逐步推进（全步长积分）
    tk=tspan(k);                                           % 当前时刻
    k1=f(tk,xk);                                           % 第 1 个斜率
    k2=f(tk+h/2,xk+h/2*k1);                                % 第 2 个斜率
    k3=f(tk+h/2,xk+h/2*k2);                                % 第 3 个斜率
    k4=f(tk+h,xk+h*k3);                                    % 第 4 个斜率
    xk=xk+h/6*(k1+2*k2+2*k3+k4);                           % RK4 更新
    if mod(k,stride)==0                                    % 到记录点
        j=j+1; X(:,j)=xk;                                  % 存入输出
    end
end
end

function [h1,th1,h2,h3]=Carplot_Ezekiel(x,y,psi,phy,box_color,wheel_color,sc)
% 复用 fullyactBC 的车形绘制；sc 为车形缩放（只影响画图）
if nargin<7, sc=1; end
h1 =  PlotCarbox(x, y, psi, 'Color', box_color, LineWidth=2);
h2 = PlotCarWheels(x, y, psi,phy,'k', 'EdgeColor',wheel_color,'Linewidth',3);
h1.XData = x + sc*(h1.XData - x);
h1.YData = y + sc*(h1.YData - y);
for kk=1:numel(h2)
    h2{kk}.XData = x + sc*(h2{kk}.XData - x);
    h2{kk}.YData = y + sc*(h2{kk}.YData - y);
end
h3 = plot(x , y, Color='b',LineStyle='-',LineWidth=4);
end

function uc = safe_solve(B,rhs)
% 2x2 显式解；仅当 B 近奇异时退回截断 SVD
d = B(1,1)*B(2,2) - B(1,2)*B(2,1);                         % 行列式
if abs(d) > 1e-8*max(1,norm(B,'fro')^2)                   % 良态：用伴随矩阵
    uc = [ B(2,2)*rhs(1) - B(1,2)*rhs(2);
          -B(2,1)*rhs(1) + B(1,1)*rhs(2)]/d;
else                                                       % 近奇异：保留原伪逆
    uc = pinv(B,1e-6)*rhs;
end
end
