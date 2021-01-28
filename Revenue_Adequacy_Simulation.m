%%%%%%%%%%%%%%%%%% �������Ȩ�г��������ԣ��������� %%%%%%%%%%%%%%%%%%%%%%%
%{
Version 1.1
����: �� ����
��λ: ��������ѧ ����ѧԺ
����: 2021/01/05
---------------------------------------------------------------------------
���������IEEE-30�ڵ��׼����ϵͳ������ƣ��������ܰ�����
1. �ͬʱ�����Բ���ģ�ͣ�ִ�н������Ȩ����
2. ���ǰ�г�ģ�ͣ�ִ����ǰ�г����������
3. ���������ԣ��n ...Based on the IEEE 30 Bus System 
����˵����
1. ��Ҫȷ���Ѱ�װMATPOWER
2. ��ͬʱ�����Բ���ģ������ǰ�г�ģ���У�������Ԥ���¹�Լ��
3. ���û�̬���˺�Ԥ���¹ʡ���·36ͣ�ˡ�ʱ����·35�ĳ�������Լ����磬���Ӧ��Ӱ
   �Ӽ۸�Ϊ 25.8302��/MW
4. ������4�����泡�����������ԣ������
   - case 1: ��̬
   - case 2: ��·�Ǽƻ�ͣ�ˣ������ԣ�Ȳ��㣨ȡ����277�е�ע�ͣ�
   - case 3: ��·�Ǽƻ�ͣ��+ȫ�ּ������ӣ��ָ������ԣ����100%��ȡ����227��155�е�ע�ͣ�
   - case 4: ��·�Ǽƻ�ͣ��+������·���ϣ��ָ������ԣ����100%��ȡ����227��230�е�ע�ͣ�
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PART 0. ��ʼ��

clear;
clc;
define_constants;
mpc         = loadcase('case30');                                          % ����IEEE-30�ڵ��׼����ϵͳ
T           = 24;                                                          % ����ǰ�г�����Ϊ24������ʱ��

% ͳ�Ʒ�����顢ĸ�ߡ��Լ�֧·��Ŀ
GEN_N       = size(mpc.gen,1);                                             
BUS_N       = size(mpc.bus,1);
BRANCH_N    = size(mpc.branch,1);

% ����ƽ��ڵ�
Bus_slack   = 6;

% �Ը����۽ڵ���б��
LoadZone    = 0;
G_1         = 1;
G_2         = 2; 
G_22        = 22; 
G_27        = 27;
G_23        = 23;
G_13        = 13;

% Ԥ��ڵ�߼ʵ�۾���
LMP.sys     = zeros(1,T);
LMP.bus     = zeros(BUS_N,T);

% Ԥ�������ԣ����ز���
FTR.CongestSurplus ...                                                     % ��ǰ�г�����ӯ��
            = zeros(1,T);
CongestPayoff ...                                                          % �������Ȩ��������
            = zeros(1,T);
RevAdequacy = zeros(1,T);                                                  % ������ʱ�ε������ԣ��
TotRevAdequacy ...                                                         % ����ƽ���˻�����������ԣ��
            = NaN;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PART 1. ͬʱ�����Բ���ģ��
%{  
˵����
1. ����ͬʱ�����Բ���ģ�ͷ��������ͻ��ɽ������Ȩ��
2. ��������Ӱ���������״����10����N-1��Ԥ���¹ʣ�
%}
%%   1.1��ģ�Ͳ�������
fprintf('FTR-PRIMARY MARKET --------------------------------------------');
fprintf('\n');

% ����ͬʱ�����Բ���ģ������
Topology_sft = ones(BRANCH_N,1);

% Ϊͬʱ�����Բ���ģ�����ø��ɷֲ�����
LDF.gs       = mpc.bus(:,PD) ./ sum(mpc.bus(:,PD));
LDF.sft      = LDF.gs;

%%   1.2��Ͷ��ָ���Ԥ���� ------------------------------------------

%           ��     �۲�   Դ�ڵ�    ��ڵ�
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

% ��������ӳ�����'Mapping'���洢�������Ȩ��Դ����ڵ�λ����Ϣ 
Map_ftr    = zeros(BIDORDER_N,BUS_N);                                         
for i = 1:BIDORDER_N
    
    % �洢Դ�ڵ��λ����Ϣ
    if BidOrder(i,3) == LoadZone
        Map_ftr(i,:) = Map_sft(i,:) + LDF.sft';
    else
        Map_ftr(i,BidOrder(i,3)) = Map_ftr(i,BidOrder(i,3)) + 1;
    end
    
    % �洢��ڵ��λ����Ϣ
    if BidOrder(i,4) == LoadZone
        Map_ftr(i,:) = Map_ftr(i,:) - LDF.sft';
    else
        Map_ftr(i,BidOrder(i,4)) = Map_ftr(i,BidOrder(i,3)) - 1;
    end
    
end

%%   1.3: ������Լ���ı��

% ���㳱��ת�Ʒֲ����Ӿ��� -------------------------------------------------

% ����Ԥ���¹ʳ���
%��N-1��Ԥ���¹�             ���ϵ���·���                                  ���ϵ���··����Դ�ڵ� -> ��ڵ㣩
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
CONTINGENCY_N            = size(Contingency_sft,2);                        % ͳ��Ԥ���¹ʵĸ���
PTDF_sft                 = zeros((CONTINGENCY_N + 1) * BRANCH_N,BUS_N);    % Ԥ�����ڴ������ת�Ʒֲ����Ӿ���

% N-0 ��̬�˷���ȫԼ��
mpc.branch(:,BR_STATUS)  = Topology_sft;
PTDF_sft(1:BRANCH_N,:)   = makePTDF(mpc,Bus_slack);

% N-1 Ԥ���¹ʰ�ȫԼ��
for i = 1:CONTINGENCY_N
    mpc.branch(Contingency_sft(i),BR_STATUS) = 0;                          % ���ϡ�N-1��Ԥ���¹ʶ�Ӧ����·
    PTDF_sft((BRANCH_N * i + 1) : BRANCH_N * (i + 1),:) ...
                         = makePTDF(mpc,Bus_slack);
    mpc.branch(:,BR_STATUS) ...
                         = Topology_sft;                                   % �����˻ָ�����̬
end
      
% ��·�����ϡ����޼��� -----------------------------------------------------
%{  
˵����Ϊ���������ԣ�ȣ�������С��1����·ȫ�ּ�������
%}

% ������·ȫ�ּ�������
% DerateFactor_sft         = 1;
DerateFactor_sft         = 0.88;                                           

% ���û���ռ�����������Ԥ��ֵ
LoopFlow_sft             = zeros(BRANCH_N,1);

% ������·�����ϡ�����
Line_max_sft             =   DerateFactor_sft * mpc.branch(:,RATE_A) - LoopFlow_sft;
Line_min_sft             = - DerateFactor_sft * mpc.branch(:,RATE_A) - LoopFlow_sft;

% FTRͷ��Լ�� -------------------------------------------------------------
Flow_ftr_max             = BidOrder(:,1);                                  % �������Ȩ�ĳ���ͷ�粻������Ͷ��ͷ��

%%   1.4�� �������Ȩһ���г�����

% Ŀ�꺯�� ----------------------------------------------------------------
% ��ע���԰�����֧���������������ΪĿ�꺯����������Թ滮����
f       = - Cost_ftr;

% Լ������ ----------------------------------------------------------------

% ϵͳ����ƽ��Լ��
Aeq     = ones(BUS_N,1)' * Map_ftr' ;
beq     = 0;

% ��·�����ϡ�����Լ�� + �������ȨͶ��ͷ��Լ��
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

% ������������ -----------------------------------------------------------
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

% �������洢 ------------------------------------------------------------
FTR.quantity         = x;
FTR.spreads.sft      = BidOrder(:,2);
FTR.path.source      = BidOrder(:,3);
FTR.path.sink        = BidOrder(:,4);
FTR.AuctRev          = - fval * T;
SP_sft.upper         = ...
    lambda.ineqlin(1 : size(PTDF_sft,1));                                  % ��·��������Լ����Ӱ�Ӽ۸�
SP_sft.lower         = ...
    lambda.ineqlin(size(PTDF_sft,1) + 1 : 2 * size(PTDF_sft,1));           % ��·��������Լ����Ӱ�Ӽ۸�

%%   1.5�� ������չʾ

figure(1);
bar(1:BIDORDER_N,[BidOrder(:,1) FTR.quantity]);
axis([0  (BIDORDER_N + 1)  0  100]);
grid on;
legend('FTR Bidding Quantity','FTR Clearing Quantity');
xlabel('FTR Bidding Order');
ylabel('FTR Quantity (MW)');
title('The Results of FTR Primary Market','FontWeight','bold');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PART 2. ��ǰ�г�ģ��
%{

%} 
%%   2.1��ģ�Ͳ�������
fprintf('DAY-AHEAD MARKET ----------------------------------------------');
fprintf('\n');

% ����������� -------------------------------------------------------------

% ������ǰ�г�ģ������
Topology_da          = ones(BRANCH_N,1);       


% ��·�Ǽƻ�ͣ��
Topology_da(18)      = 0;                                                  % ����18����·��12 -> 15
fprintf('Outage: Line 12->15');
fprintf('\n');

% % ��·�Ǽƻ�ͣ�� + ������·����
% Topology_da(32)      = 0;                                                  % ����32����·��23 -> 24
% mpc.branch(:,BR_STATUS) ...
%                      = Topology_da;
% fprintf('Optimal Transmission Switching: Line 23->24');
% fprintf('\n');

% ����Ԥ���¹ʳ���
Contingency_da          = Contingency_sft;                                 % ��ǰ�г�ģ���У�������ͬʱ�����Բ���ģ��һ�µ�Ԥ���¹ʳ���
                 
% ���㳱��ת�Ʒֲ����Ӿ���
% PTDF_da              = PTDF_sft;
PTDF_da                 = zeros(((CONTINGENCY_N + 1) * BRANCH_N),BUS_N);   % Ԥ�����ڴ������ת�Ʒֲ����Ӿ���

% N-0 ��̬�˷���ȫԼ��
mpc.branch(:,BR_STATUS) = Topology_da;
ptdf_da_gs              = makePTDF(mpc,Bus_slack);
PTDF_da(1:BRANCH_N,:)   = ptdf_da_gs;

% N-1 Ԥ���¹ʰ�ȫԼ��
for i = 1:CONTINGENCY_N
    mpc.branch(Contingency_da(i),BR_STATUS) = 0;                           % ���ϡ�N-1��Ԥ���¹ʶ�Ӧ����·
    PTDF_da((BRANCH_N * i + 1) : BRANCH_N * (i + 1),:) ...
                     = makePTDF(mpc,Bus_slack);
    mpc.branch(:,BR_STATUS) ...
                     = Topology_da;                                        % �����˻ָ�����̬
end

% ��·�����ϡ����޼��� -----------------------------------------------------

% ������·ȫ�ּ�������
DerateFactor_da      = 1;

% ���û���ռ�����������Ԥ��ֵ
LoopFlow             = zeros(BUS_N,T);
% LoopFlow(5,:)        = randn(1,T);
LoopFlow(5,:)        = zeros(1,T);
LoopFlow(30,:)       = - LoopFlow(5,:);

% ������·�����ϡ�����
Line_max_da          =   DerateFactor_da .* mpc.branch(:,RATE_A) * ones(1,T) ...
                         - ptdf_da_gs * LoopFlow;
Line_min_da          = - DerateFactor_da .* mpc.branch(:,RATE_A) * ones(1,T) ...
                         - ptdf_da_gs * LoopFlow;

% ���ɲ������� -------------------------------------------------------------

% ����ϵͳ�������ߣ��ο��㶫ͳ���������������������ã�
% ��ע�����ǡ�N-1��������£����������Ӧ����1.23�������޽�
LoadRatio            = [0.97 0.92 0.88 0.84 0.81 0.85 ...
                        1.01 1.13 1.16 1.19 1.23 1.14 ...
                        1.12 1.18 1.21 1.22 1.20 1.13 ...
                        1.16 1.18 1.17 1.14 1.10 1.02];
Load.sys             = sum(mpc.bus(:,PD)) * LoadRatio;

% ���ϵͳ��������
figure(2);
bar(1:24, Load.sys);
axis normal;
axis([0 25 100 250]);
set(gca,'XTick',1:1:24) ;
set(gca,'PlotBoxAspectRatio',[2 1 1]);
grid on;
xlabel('Hours'); ylabel('System Load (MW)');
title('The Total Load Curve','FontWeight','bold');

% ���ýڵ㸺��
LDF.da               = LDF.sft;                                            % Ϊ��ǰ�г�ģ�����ø��ɷֲ�����
Load.bus             = zeros(BUS_N,T);
for t = 1:T
    Load.bus(:,t)    = Load.sys(t) * LDF.da;
end

% �������������� ---------------------------------------------------------

% �������'Map_gen'����﷢�������ĸ�߼�Ĺ�����ϵ
Map_gen              = zeros(GEN_N,BUS_N);
for i = 1:GEN_N
    Map_gen(i,mpc.gen(i,GEN_BUS)) = 1;
end
% Map_gen(:,Bus_slack) = - ones(Gen_N,1);

% ���÷������ɱ�����
%{
��ע: 
1. Total Cost = a * x^2 + b * x + c
2. Ϊ���۸�λ���Ϊ��Ԫ/MWh����ϵ��a/b/c������100
%}
Gen.Cost.a   = 100 * mpc.gencost(:,5);
% Gen.Cost.a   = 100 * [0.0200; 0.0175; 0.0725; 0.0083; 0.0250; 0.0250];
Gen.Cost.b   = 100 * mpc.gencost(:,6);
Gen.Cost.c   = 100 * mpc.gencost(:,7);

%%   2.2����ǰ�г�����

% Ŀ�꺯�� ----------------------------------------------------------------
% ��ע�����ܷ���ɱ���С��ΪĿ�꣬�����ι滮���� min ( 1/2 * x' * H * x + f' * x )

% Ŀ�꺯��������
H             = zeros(BUS_N,BUS_N);
for i = 1:GEN_N
    H(mpc.gen(i,GEN_BUS),mpc.gen(i,GEN_BUS)) ...
              = 2 * Gen.Cost.a(i);
end

% Ŀ�꺯��һ����
f             = Map_gen' * Gen.Cost.b;

% Լ������ -----------------------------------------------------------------

% ϵͳ����ƽ��Լ��
Aeq           = ones(1,BUS_N);
beq           = NaN;

% ֧·�����ϡ�����Լ��
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
      
% �����������ϡ�����Լ��
%{
��ע�������������ɽ����㣬��˲�������̧��uplift������;
%}
lb            = zeros(BUS_N,1);
ub            = Map_gen' * mpc.gen(:,PMAX);

% ������������ -----------------------------------------------------------
% x0            = Load.sys(1) * mpc.gen(:,PMAX) ./ sum(mpc.gen(:,PMAX)); 
Gen.Output    = zeros(GEN_N,T);
Gen.TotCost   = zeros(1,T);
SP_da.upper   = zeros(size(PTDF_da,1),T);                                  % ��·��������Լ����Ӱ�Ӽ۸�
SP_da.lower   = zeros(size(PTDF_da,1),T);                                  % ��·��������Լ����Ӱ�Ӽ۸�
exitflag      = zeros(1,T);
tic
bar = waitbar(0,'Loading Data for Day-ahead Market ...');                  % waitbar��ʾ������
for t = 1:T
    
    % �޸�ϵͳ����ƽ��Լ��RHS
    beq       =  Load.sys(t);
       
    % �޸�֧·�����ϡ�����Լ��RHS
    b         = [ PTDF_da * Load.bus(:,t) + Line_max_da_extend(:,t);
                - PTDF_da * Load.bus(:,t) - Line_min_da_extend(:,t) ];
    options   = optimoptions('quadprog','Algorithm','interior-point-convex','Display','off', ...
        'MaxIterations',200,'OptimalityTolerance',1e-12,'StepTolerance',1e-18,'ConstraintTolerance',1e-19);
    [x,fval,exitflag(t),output,lambda] ...
              = quadprog(H,f,A,b,Aeq,beq,lb,ub,[],options);
    
    if exitflag(t) < 0
        break;
    end
    
    % �������洢
    Gen.Output(:,t)  = Map_gen * x;
    Gen.TotCost(:,t) = fval;
    LMP.sys(t)       = - lambda.eqlin;                                     % ƽ��ڵ㴦�Ľڵ�߼ʵ�ۣ�����ν��ϵͳ�������۸�
    SP_da.upper(:,t) = lambda.ineqlin(1 : size(PTDF_da,1));
    SP_da.lower(:,t) = lambda.ineqlin(size(PTDF_da,1) + 1 : end);
    
    % ��ʾ������
    str              = ...
        ['Day-ahead Market Clearing...',num2str(100*t/T),'%'];             % �԰ٷֱ���ʽ��ʾ�������
    waitbar(t/T,bar,str)                                                   % ���½�����bar�����barʹ��
    
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

%%   2.3����ǰ�г�����

% Ԥ�����ڴ� ---------------------------------------------------------------
Gen.SpotRev           = zeros(GEN_N,T);
Gen.TotSpotRev        = zeros(1,T);
Load.TotPayout        = zeros(1,T);
FTR.spreads.da        = zeros(BIDORDER_N,T);
FTR.CongestPayoff     = zeros(BIDORDER_N,T);
FTR.TotCongestPayoff  = zeros(1,T);

% ��ڵ�߼ʵ�� -----------------------------------------------------------
for t = 1:T
    LMP.bus(:,t)      = LMP.sys(t) - ...                                   % ��tʱ�εĽڵ�߼ʵ��
        PTDF_da' * (SP_da.upper(:,t) - SP_da.lower(:,t));
    LMP.LoadZone(t)   = LMP.bus(:,t)' * Load.bus(:,t) ./ Load.sys(t);      % �󸺺���tʱ�εļ�Ȩƽ���ڵ�߼ʵ��
end

% �Ը��г������������� -----------------------------------------------------

% ����� + ���ɲ�
for t = 1:T
    Gen.SpotRev(:,t)  = Gen.Output(:,t) .* Map_gen * LMP.bus(:,t);         % �������tʱ�ε��ֻ�����
    Gen.TotSpotRev(t) = sum(Gen.SpotRev(:,t));                             % �󷢵��tʱ�ε��ֻ�����
    Load.TotPayout(t) = LMP.LoadZone(t) * Load.sys(t);                     % �󸺺ɲ�tʱ�ε��ֻ�֧��
end

% �������Ȩ
for t = 1:T
    FTR.spreads.da(:,t)     ...
                      = - Map_ftr * LMP.bus(:,t);                          % ��������Ȩ���ֻ��г��еĽ���۲�
    FTR.CongestPayoff(:,t)  ...
                      = FTR.spreads.da(:,t) .* FTR.quantity;               % ��·���������Ȩtʱ�ε���������
    FTR.TotCongestPayoff(t) ...
                      = sum(FTR.CongestPayoff(:,t));                       % ȫ���������Ȩtʱ�ε�����������
end

% �������ԣ�� ------------------------------------------------------------
FTR.CongestSurplus    = Load.TotPayout - Gen.TotSpotRev;                   % ���г���Ӫ����tʱ�����������ӯ��
for t = 1:T
    if (-1e-4 < FTR.TotCongestPayoff(t)) && ...
            (FTR.TotCongestPayoff(t) < 1e-4)                               % ��������������̫Сʱ��������ֵ��������
        RevAdequacy(t) = 1;
    else
        RevAdequacy(t) = FTR.CongestSurplus(t) ./ FTR.TotCongestPayoff(t); 
    end
end
TotRevAdequacy         = sum(FTR.CongestSurplus) ./ sum(FTR.TotCongestPayoff);

%% PART 3. Post-Processing

% �г���������� ---------------------------------------------------------
fprintf('RESULTS REPORT ------------------------------------------------');
fprintf('\n');
fprintf('The total cost paid to GENCOs is �� %.2f %\n',sum(Gen.TotSpotRev));
fprintf('\n');
fprintf('The total revenue collected from LSEs is �� %.2f %\n',sum(Load.TotPayout));
fprintf('\n');
fprintf('The total revenue adequacy of FTR market is %.2f %% \n', 100 * TotRevAdequacy);

% ����������µĽڵ�߼ʵ������ -------------------------------------------
figure(3);
plot(1:24,LMP.LoadZone,'--o','color','g','linewidth',1.5);
legend('LMP ( Ground-State )');
xlabel('Hours');
ylabel('LMP ( ��/MWh )');
title('The LMP Demonstration','FontWeight','bold');
grid on;
set(gca,'XTick',1:1:24);
set(gcf,'unit','centimeters','position',[20 10 20 10]);

% �������������ʱ���µ����ų������ ----------------------------------------

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

% ��ͼչʾ��ʱ�εĲ�ƽ���ʽ� -----------------------------------------------
figure(5);
bar(1:24, FTR.CongestSurplus - FTR.TotCongestPayoff);
set(gca,'XTick',1:1:24) ;
% set(gca,'PlotBoxAspectRatio',[2 1 1]);
xlabel('Hours');
ylabel('Amount of Money ( �� )');
title('The Unbalanced Funds','FontWeight','bold');
grid on;

%% PART 4. ����
load train;
sound(y,Fs);

% % Save Results
% filename = 'testdata.xlsx';