%%%%%%%%%%%%%%%%%% 金融输电权市场的收入充裕度问题仿真 %%%%%%%%%%%%%%%%%%%%%%%
%{
Version 1.1
作者: 曾 鹏骁
单位: 华南理工大学 电力学院
日期: 2021/01/05
---------------------------------------------------------------------------
本程序基于IEEE-30节点标准测试系统进行设计，基本功能包含：
1. 搭建同时可行性测试模型，执行金融输电权拍卖
2. 搭建日前市场模型，执行日前市场出清与结算
3. 计算收入充裕度n ...Based on the IEEE 30 Bus System 
补充说明：
1. 需要确保已安装MATPOWER
2. 在同时可行性测试模型与日前市场模型中，均考虑预想事故约束
3. 采用基态拓扑和预想事故“线路36停运”时，线路35的潮流下限约束达界，其对应的影
   子价格为 25.8302￥/MW
4. 设置了4个仿真场景讨论收入充裕度问题
   - case 1: 基态
   - case 2: 线路非计划停运，收入充裕度不足（取消第277行的注释）
   - case 3: 线路非计划停运+全局减容因子，恢复收入充裕度至100%（取消第227、155行的注释）
   - case 4: 线路非计划停运+最优线路开断，恢复收入充裕度至100%（取消第227、230行的注释）
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PART 0. 初始化

clear;
clc;
define_constants;
mpc         = loadcase('case30');                                          % 加载IEEE-30节点标准测试系统
T           = 24;                                                          % 将日前市场划分为24个结算时段

% 统计发电机组、母线、以及支路数目
GEN_N       = size(mpc.gen,1);                                             
BUS_N       = size(mpc.bus,1);
BRANCH_N    = size(mpc.branch,1);

% 定义平衡节点
Bus_slack   = 6;

% 对各定价节点进行编号
LoadZone    = 0;
G_1         = 1;
G_2         = 2; 
G_22        = 22; 
G_27        = 27;
G_23        = 23;
G_13        = 13;

% 预设节点边际电价矩阵
LMP.sys     = zeros(1,T);
LMP.bus     = zeros(BUS_N,T);

% 预设收入充裕度相关参数
FTR.CongestSurplus ...                                                     % 日前市场阻塞盈余
            = zeros(1,T);
CongestPayoff ...                                                          % 金融输电权阻塞收益
            = zeros(1,T);
RevAdequacy = zeros(1,T);                                                  % 各结算时段的收入充裕度
TotRevAdequacy ...                                                         % 配置平衡账户后的总收入充裕度
            = NaN;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PART 1. 同时可行性测试模型
%{  
说明：
1. 基于同时可行性测试模型发行义务型基荷金融输电权；
2. 考虑显著影响电网阻塞状况的10个“N-1”预想事故；
%}
%%   1.1：模型参数配置
fprintf('FTR-PRIMARY MARKET --------------------------------------------');
fprintf('\n');

% 定义同时可行性测试模型拓扑
Topology_sft = ones(BRANCH_N,1);

% 为同时可行性测试模型设置负荷分布因子
LDF.gs       = mpc.bus(:,PD) ./ sum(mpc.bus(:,PD));
LDF.sft      = LDF.gs;

%%   1.2：投标指令的预处理 ------------------------------------------

%           量     价差   源节点    汇节点
% BidOrder = [80   -1.7748   G_1    LoadZone; 
%             80   -1.7748   G_2    LoadZone;
%             50   -1.7748   G_22   LoadZone;
%             55   24.0554   G_27   LoadZone;
%             30   -1.7748   G_23   LoadZone;
%             40   -1.7748   G_13   LoadZone];
BidOrder   = [80   1.7748   G_1    LoadZone;
              80  -1.7748   G_2    LoadZone;
              50  -1.7748   G_22   LoadZone;
              55  24.0554   G_27   LoadZone;
              30  -1.7748   G_23   LoadZone;
              40  -1.7748   G_13   LoadZone];
BIDORDER_N = size(BidOrder,1);
Cost_ftr   = BidOrder(:,2);

% 引入线性映射矩阵'Mapping'，存储金融输电权的源、汇节点位置信息 
Map_ftr    = zeros(BIDORDER_N,BUS_N);                                         
for i = 1:BIDORDER_N
    
    % 存储源节点的位置信息
    if BidOrder(i,3) == LoadZone
        Map_ftr(i,:) = Map_sft(i,:) + LDF.sft';
    else
        Map_ftr(i,BidOrder(i,3)) = Map_ftr(i,BidOrder(i,3)) + 1;
    end
    
    % 存储汇节点的位置信息
    if BidOrder(i,4) == LoadZone
        Map_ftr(i,:) = Map_ftr(i,:) - LDF.sft';
    else
        Map_ftr(i,BidOrder(i,4)) = Map_ftr(i,BidOrder(i,3)) - 1;
    end
    
end

%%   1.3: 电网络约束的表达

% 计算潮流转移分布因子矩阵 -------------------------------------------------

% 设置预想事故场景
%“N-1”预想事故             开断的线路编号                                  开断的线路路径（源节点 -> 汇节点）
Contingency_sft          = [12 ...                                         % Bus 6  -> Bus 10                 
                            36 ...                                         % Bus 28 -> Bus 27
                            14 ...                                         % Bus 9  -> Bus 10
                            15 ...                                         % Bus 4  -> Bus 12
                            26 ...                                         % Bus 10 -> Bus 17
                            32 ...                                         % Bus 23 -> Bus 24
                            17 ...                                         % Bus 12 -> Bus 14
                            1  ...                                         % Bus 1  -> Bus 2
                            27 ...                                         % Bus 10 -> Bus 21
                            7  ...                                         % Bus 4  -> Bus 6
                            ];
CONTINGENCY_N            = size(Contingency_sft,2);                        % 统计预想事故的个数
PTDF_sft                 = zeros((CONTINGENCY_N + 1) * BRANCH_N,BUS_N);    % 预分配内存给潮流转移分布因子矩阵

% N-0 基态运方安全约束
mpc.branch(:,BR_STATUS)  = Topology_sft;
PTDF_sft(1:BRANCH_N,:)   = makePTDF(mpc,Bus_slack);

% N-1 预想事故安全约束
for i = 1:CONTINGENCY_N
    mpc.branch(Contingency_sft(i),BR_STATUS) = 0;                          % 开断“N-1”预想事故对应的线路
    PTDF_sft((BRANCH_N * i + 1) : BRANCH_N * (i + 1),:) ...
                         = makePTDF(mpc,Bus_slack);
    mpc.branch(:,BR_STATUS) ...
                         = Topology_sft;                                   % 将拓扑恢复至基态
end
      
% 线路潮流上、下限计算 -----------------------------------------------------
%{  
说明：为保障收入充裕度，可设置小于1的线路全局减容因子
%}

% 设置线路全局减容因子
% DerateFactor_sft         = 1;
DerateFactor_sft         = 0.88;                                           

% 设置环流占用输电容量的预测值
LoopFlow_sft             = zeros(BRANCH_N,1);

% 折算线路潮流上、下限
Line_max_sft             =   DerateFactor_sft * mpc.branch(:,RATE_A) - LoopFlow_sft;
Line_min_sft             = - DerateFactor_sft * mpc.branch(:,RATE_A) - LoopFlow_sft;

% FTR头寸约束 -------------------------------------------------------------
Flow_ftr_max             = BidOrder(:,1);                                  % 金融输电权的出清头寸不超过其投标头寸

%%   1.4： 金融输电权一级市场出清

% 目标函数 ----------------------------------------------------------------
% 备注：以按报价支付的拍卖收益最大化为目标函数，求解线性规划问题
f       = - Cost_ftr;

% 约束条件 ----------------------------------------------------------------

% 系统功率平衡约束
Aeq     = ones(BUS_N,1)' * Map_ftr' ;
beq     = 0;

% 线路潮流上、下限约束 + 金融输电权投标头寸约束
A       = [ PTDF_sft * Map_ftr';...
          - PTDF_sft * Map_ftr';...
            eye(BIDORDER_N,BIDORDER_N) ];
Line_max_sft_extend ...
        = zeros(BRANCH_N * (CONTINGENCY_N + 1),1);
Line_min_sft_extend ...
        = zeros(BRANCH_N * (CONTINGENCY_N + 1),1);
for i   = 1 : (CONTINGENCY_N + 1)
    Line_max_sft_extend( BRANCH_N * (i - 1) + 1 : BRANCH_N * i) ...
        = Line_max_sft;
    Line_min_sft_extend( BRANCH_N * (i - 1) + 1 : BRANCH_N * i) ...
        = Line_min_sft;
end
b       = [ Line_max_sft_extend;...
          - Line_min_sft_extend;...
            Flow_ftr_max ];
lb      = zeros(1,BIDORDER_N);
ub      = BidOrder(:,1)';                                                  

% 调用求解器求解 -----------------------------------------------------------
options = optimoptions('linprog','Algorithm','dual-simplex','Display','off','ConstraintTolerance',1e-9);
tic
[x,fval,exitflag,output,lambda] = linprog(f,A,b,Aeq,beq,lb,ub,options);
toc
if exitflag == 1 || exitflag == 0
    fprintf('FTR primary market clearing completed');
else
    fprintf('FTR primary market clearing uncompleted');
end
fprintf('\n\n');

% 出清结果存储 ------------------------------------------------------------
FTR.quantity         = x;
FTR.spreads.sft      = BidOrder(:,2);
FTR.path.source      = BidOrder(:,3);
FTR.path.sink        = BidOrder(:,4);
FTR.AuctRev          = - fval * T;
SP_sft.upper         = ...
    lambda.ineqlin(1 : size(PTDF_sft,1));                                  % 线路潮流上限约束的影子价格
SP_sft.lower         = ...
    lambda.ineqlin(size(PTDF_sft,1) + 1 : 2 * size(PTDF_sft,1));           % 线路潮流下限约束的影子价格

%%   1.5： 出清结果展示

figure(1);
bar(1:BIDORDER_N,[BidOrder(:,1) FTR.quantity]);
axis([0  (BIDORDER_N + 1)  0  100]);
grid on;
legend('FTR Bidding Quantity','FTR Clearing Quantity');
xlabel('FTR Bidding Order');
ylabel('FTR Quantity (MW)');
title('The Results of FTR Primary Market','FontWeight','bold');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PART 2. 日前市场模型
%{

%} 
%%   2.1：模型参数配置
fprintf('DAY-AHEAD MARKET ----------------------------------------------');
fprintf('\n');

% 网络参数配置 -------------------------------------------------------------

% 定义日前市场模型拓扑
Topology_da          = ones(BRANCH_N,1);       


% 线路非计划停运
Topology_da(18)      = 0;                                                  % 开断18号线路：12 -> 15
fprintf('Outage: Line 12->15');
fprintf('\n');

% % 线路非计划停运 + 最优线路开断
% Topology_da(32)      = 0;                                                  % 开断32号线路：23 -> 24
% mpc.branch(:,BR_STATUS) ...
%                      = Topology_da;
% fprintf('Optimal Transmission Switching: Line 23->24');
% fprintf('\n');

% 设置预想事故场景
Contingency_da          = Contingency_sft;                                 % 日前市场模型中，设置与同时可行性测试模型一致的预想事故场景
                 
% 计算潮流转移分布因子矩阵
% PTDF_da              = PTDF_sft;
PTDF_da                 = zeros(((CONTINGENCY_N + 1) * BRANCH_N),BUS_N);   % 预分配内存给潮流转移分布因子矩阵

% N-0 基态运方安全约束
mpc.branch(:,BR_STATUS) = Topology_da;
ptdf_da_gs              = makePTDF(mpc,Bus_slack);
PTDF_da(1:BRANCH_N,:)   = ptdf_da_gs;

% N-1 预想事故安全约束
for i = 1:CONTINGENCY_N
    mpc.branch(Contingency_da(i),BR_STATUS) = 0;                           % 开断“N-1”预想事故对应的线路
    PTDF_da((BRANCH_N * i + 1) : BRANCH_N * (i + 1),:) ...
                     = makePTDF(mpc,Bus_slack);
    mpc.branch(:,BR_STATUS) ...
                     = Topology_da;                                        % 将拓扑恢复至基态
end

% 线路潮流上、下限计算 -----------------------------------------------------

% 设置线路全局减容因子
DerateFactor_da      = 1;

% 设置环流占用输电容量的预测值
LoopFlow             = zeros(BUS_N,T);
% LoopFlow(5,:)        = randn(1,T);
LoopFlow(5,:)        = zeros(1,T);
LoopFlow(30,:)       = - LoopFlow(5,:);

% 折算线路潮流上、下限
Line_max_da          =   DerateFactor_da .* mpc.branch(:,RATE_A) * ones(1,T) ...
                         - ptdf_da_gs * LoopFlow;
Line_min_da          = - DerateFactor_da .* mpc.branch(:,RATE_A) * ones(1,T) ...
                         - ptdf_da_gs * LoopFlow;

% 负荷参数配置 -------------------------------------------------------------

% 设置系统负荷曲线（参考广东统调负荷曲线特征进行设置）
% 备注：考虑“N-1”的情况下，负荷率最大不应超过1.23，否则无解
LoadRatio            = [0.97 0.92 0.88 0.84 0.81 0.85 ...
                        1.01 1.13 1.16 1.19 1.23 1.14 ...
                        1.12 1.18 1.21 1.22 1.20 1.13 ...
                        1.16 1.18 1.17 1.14 1.10 1.02];
Load.sys             = sum(mpc.bus(:,PD)) * LoadRatio;

% 绘出系统负荷曲线
figure(2);
bar(1:24, Load.sys);
axis normal;
axis([0 25 100 250]);
set(gca,'XTick',1:1:24) ;
set(gca,'PlotBoxAspectRatio',[2 1 1]);
grid on;
xlabel('Hours'); ylabel('System Load (MW)');
title('The Total Load Curve','FontWeight','bold');

% 设置节点负荷
LDF.da               = LDF.sft;                                            % 为日前市场模型设置负荷分布因子
Load.bus             = zeros(BUS_N,T);
for t = 1:T
    Load.bus(:,t)    = Load.sys(t) * LDF.da;
end

% 发电机组参数配置 ---------------------------------------------------------

% 构造矩阵'Map_gen'，表达发电机组与母线间的关联关系
Map_gen              = zeros(GEN_N,BUS_N);
for i = 1:GEN_N
    Map_gen(i,mpc.gen(i,GEN_BUS)) = 1;
end
% Map_gen(:,Bus_slack) = - ones(Gen_N,1);

% 设置发电机组成本参数
%{
备注: 
1. Total Cost = a * x^2 + b * x + c
2. 为将价格单位表达为“元/MWh”，系数a/b/c均乘以100
%}
Gen.Cost.a   = 100 * mpc.gencost(:,5);
% Gen.Cost.a   = 100 * [0.0200; 0.0175; 0.0725; 0.0083; 0.0250; 0.0250];
Gen.Cost.b   = 100 * mpc.gencost(:,6);
Gen.Cost.c   = 100 * mpc.gencost(:,7);

%%   2.2：日前市场出清

% 目标函数 ----------------------------------------------------------------
% 备注：以总发电成本最小化为目标，求解二次规划问题 min ( 1/2 * x' * H * x + f' * x )

% 目标函数二次项
H             = zeros(BUS_N,BUS_N);
for i = 1:GEN_N
    H(mpc.gen(i,GEN_BUS),mpc.gen(i,GEN_BUS)) ...
              = 2 * Gen.Cost.a(i);
end

% 目标函数一次项
f             = Map_gen' * Gen.Cost.b;

% 约束条件 -----------------------------------------------------------------

% 系统功率平衡约束
Aeq           = ones(1,BUS_N);
beq           = NaN;

% 支路潮流上、下限约束
A             = [ PTDF_da; ...
                - PTDF_da ];
b             = NaN * ones(size(A,1),1);
Line_max_da_extend ...
              = zeros(BRANCH_N * (CONTINGENCY_N + 1),T);
Line_min_da_extend ...
              = zeros(BRANCH_N * (CONTINGENCY_N + 1),T);
for i   = 1 : (CONTINGENCY_N + 1)
    Line_max_da_extend( BRANCH_N * (i - 1) + 1 : BRANCH_N * i,:) ...
              = Line_max_da;
    Line_min_da_extend( BRANCH_N * (i - 1) + 1 : BRANCH_N * i,:) ...
              = Line_min_da;
end
      
% 发电机组出力上、下限约束
%{
备注：假设机组出力可降至零，因此不存在上抬（uplift）费用;
%}
lb            = zeros(BUS_N,1);
ub            = Map_gen' * mpc.gen(:,PMAX);

% 调用求解器求解 -----------------------------------------------------------
% x0            = Load.sys(1) * mpc.gen(:,PMAX) ./ sum(mpc.gen(:,PMAX)); 
Gen.Output    = zeros(GEN_N,T);
Gen.TotCost   = zeros(1,T);
SP_da.upper   = zeros(size(PTDF_da,1),T);                                  % 线路潮流上限约束的影子价格
SP_da.lower   = zeros(size(PTDF_da,1),T);                                  % 线路潮流下限约束的影子价格
exitflag      = zeros(1,T);
tic
bar = waitbar(0,'Loading Data for Day-ahead Market ...');                  % waitbar显示进度条
for t = 1:T
    
    % 修改系统功率平衡约束RHS
    beq       =  Load.sys(t);
       
    % 修改支路潮流上、下限约束RHS
    b         = [ PTDF_da * Load.bus(:,t) + Line_max_da_extend(:,t);
                - PTDF_da * Load.bus(:,t) - Line_min_da_extend(:,t) ];
    options   = optimoptions('quadprog','Algorithm','interior-point-convex','Display','off', ...
        'MaxIterations',200,'OptimalityTolerance',1e-12,'StepTolerance',1e-18,'ConstraintTolerance',1e-19);
    [x,fval,exitflag(t),output,lambda] ...
              = quadprog(H,f,A,b,Aeq,beq,lb,ub,[],options);
    
    if exitflag(t) < 0
        break;
    end
    
    % 出清结果存储
    Gen.Output(:,t)  = Map_gen * x;
    Gen.TotCost(:,t) = fval;
    LMP.sys(t)       = - lambda.eqlin;                                     % 平衡节点处的节点边际电价，即所谓“系统电能量价格”
    SP_da.upper(:,t) = lambda.ineqlin(1 : size(PTDF_da,1));
    SP_da.lower(:,t) = lambda.ineqlin(size(PTDF_da,1) + 1 : end);
    
    % 显示进度条
    str              = ...
        ['Day-ahead Market Clearing...',num2str(100*t/T),'%'];             % 以百分比形式显示处理进程
    waitbar(t/T,bar,str)                                                   % 更新进度条bar，配合bar使用
    
end

close(bar)
clearvars bar;
toc

if min(exitflag) < 0
    fprintf('Day-ahead market clearing uncompleted');
else
    fprintf('Day-ahead market clearing completed');
end
fprintf('\n\n');

%%   2.3：日前市场结算

% 预分配内存 ---------------------------------------------------------------
Gen.SpotRev           = zeros(GEN_N,T);
Gen.TotSpotRev        = zeros(1,T);
Load.TotPayout        = zeros(1,T);
FTR.spreads.da        = zeros(BIDORDER_N,T);
FTR.CongestPayoff     = zeros(BIDORDER_N,T);
FTR.TotCongestPayoff  = zeros(1,T);

% 求节点边际电价 -----------------------------------------------------------
for t = 1:T
    LMP.bus(:,t)      = LMP.sys(t) - ...                                   % 求t时段的节点边际电价
        PTDF_da' * (SP_da.upper(:,t) - SP_da.lower(:,t));
    LMP.LoadZone(t)   = LMP.bus(:,t)' * Load.bus(:,t) ./ Load.sys(t);      % 求负荷区t时段的加权平均节点边际电价
end

% 对各市场参与者作结算 -----------------------------------------------------

% 发电侧 + 负荷侧
for t = 1:T
    Gen.SpotRev(:,t)  = Gen.Output(:,t) .* Map_gen * LMP.bus(:,t);         % 求各机组t时段的现货收益
    Gen.TotSpotRev(t) = sum(Gen.SpotRev(:,t));                             % 求发电侧t时段的现货收益
    Load.TotPayout(t) = LMP.LoadZone(t) * Load.sys(t);                     % 求负荷侧t时段的现货支出
end

% 金融输电权
for t = 1:T
    FTR.spreads.da(:,t)     ...
                      = - Map_ftr * LMP.bus(:,t);                          % 求金融输电权在现货市场中的结算价差
    FTR.CongestPayoff(:,t)  ...
                      = FTR.spreads.da(:,t) .* FTR.quantity;               % 各路径金融输电权t时段的阻塞收益
    FTR.TotCongestPayoff(t) ...
                      = sum(FTR.CongestPayoff(:,t));                       % 全部金融输电权t时段的总阻塞收益
end

% 求收入充裕度 ------------------------------------------------------------
FTR.CongestSurplus    = Load.TotPayout - Gen.TotSpotRev;                   % 求市场运营机构t时段收入的阻塞盈余
for t = 1:T
    if (-1e-4 < FTR.TotCongestPayoff(t)) && ...
            (FTR.TotCongestPayoff(t) < 1e-4)                               % 避免总阻塞收益太小时，产生数值计算问题
        RevAdequacy(t) = 1;
    else
        RevAdequacy(t) = FTR.CongestSurplus(t) ./ FTR.TotCongestPayoff(t); 
    end
end
TotRevAdequacy         = sum(FTR.CongestSurplus) ./ sum(FTR.TotCongestPayoff);

%% PART 3. Post-Processing

% 市场基本情况简报 ---------------------------------------------------------
fprintf('RESULTS REPORT ------------------------------------------------');
fprintf('\n');
fprintf('The total cost paid to GENCOs is ￥ %.2f %\n',sum(Gen.TotSpotRev));
fprintf('\n');
fprintf('The total revenue collected from LSEs is ￥ %.2f %\n',sum(Load.TotPayout));
fprintf('\n');
fprintf('The total revenue adequacy of FTR market is %.2f %% \n', 100 * TotRevAdequacy);

% 绘出各场景下的节点边际电价走势 -------------------------------------------
figure(3);
plot(1:24,LMP.LoadZone,'--o','color','g','linewidth',1.5);
legend('LMP ( Ground-State )');
xlabel('Hours');
ylabel('LMP ( ￥/MWh )');
title('The LMP Demonstration','FontWeight','bold');
grid on;
set(gca,'XTick',1:1:24);
set(gcf,'unit','centimeters','position',[20 10 20 10]);

% 绘出阻塞最严重时段下的最优潮流情况 ----------------------------------------

OPF_CongestHour        = PTDF_da(1:BRANCH_N,:) ...
                         * (Map_gen' * Gen.Output(:,11) - Load.bus(:,11));
figure(4);
bar(1:BRANCH_N,OPF_CongestHour);
xlabel('Line');
ylabel('OPF ( MW )');
title('The OPF Distribution at the Most Serious Congested Hour','FontWeight','bold');
grid on;
set(gca,'XTick',0:5:(BRANCH_N + 1));
set(gcf,'unit','centimeters','position',[20 10 20 10]);

% 绘图展示各时段的不平衡资金 -----------------------------------------------
figure(5);
bar(1:24, FTR.CongestSurplus - FTR.TotCongestPayoff);
set(gca,'XTick',1:1:24) ;
% set(gca,'PlotBoxAspectRatio',[2 1 1]);
xlabel('Hours');
ylabel('Amount of Money ( ￥ )');
title('The Unbalanced Funds','FontWeight','bold');
grid on;

%% PART 4. 结束
load train;
sound(y,Fs);

% % Save Results
% filename = 'testdata.xlsx';