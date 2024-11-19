load('light.mat');
load('light_0.mat');
load('wuhan_unicom_mobility.mat');
load('wuhan2019.mat');
[factor_delta,factor] = distributeL(wuhan_outof,wuhan_into);
[factor_delta2,factor2] = distributeL2(wuhan2019);
light_daily = [light(:,1).*factor_delta{1},light(:,2).*factor_delta{2},light(:,3).*factor_delta{3}];
light0_daily = [light_0(:,1).*factor_delta2{1},light_0(:,2).*factor_delta2{2},light_0(:,3).*factor_delta2{3}];

global L_0;
L_0=light_daily;

doly=length(light0_daily);
loop_interval = 10;
ntl_struc_col = [253,192,134;
    190,174,212;
    127,201,127;
    244,202,228;
    56,108,176]/255;
figure('color','w');
tiledlayout('flow','Padding','compact','TileSpacing','compact');
nexttile
for i=1:size(light_daily,1)
    plot(light_daily(i,:),'-','LineWidth',1.5,'Color',ntl_struc_col(i,:)); hold on
end
legend({'RES','COM','IND','TSP','PUB'},'Location','southeast','NumColumns',5);
xlim([0 90]);
ylabel('dNTL 2020');
set(gca,'FontSize',12);
nexttile
for i=1:size(light_daily,1)
    plot(light0_daily(i,:),'--','LineWidth',1.5,'Color',ntl_struc_col(i,:)); hold on
end
legend({'RES','COM','IND','TSP','PUB'},'Location','southeast','NumColumns',5);
xlim([0 90]);
xlabel('DAY');
ylabel('dNTL 2019');
set(gca,'FontSize',12);

%% calibration or certain parameters
%% params
%SEIR
params_setting = zeros(1,20);
params_setting(1)=0.177;%beta0;_____
params_setting(2)=1/14;%alpha1;
params_setting(3)=1;%symptomatic proportion
params_setting(4)=1/10;%gamma1;
params_setting(5)=1/20;%gamma20;
params_setting(6)=0.05;%p_no_infectious;_____
params_setting(7)=23;%L_beta_cfr;_____
params_setting(8)=38;%L_p;_____
params_setting(9)=55;%L_g;
params_setting(10)=3;%k_beta;_____
params_setting(11)=1;%k_p;_____
params_setting(12)=10;%k_g;
params_setting(13)=0.04;%cfr;
params_setting(14)=0.9999;%proportion;
params_setting(15)=1.002;%pSE;
thetaup = 1;
params_setting(16) = 302.37994*thetaup;%theta_r;
params_setting(17) = 69.85718554*thetaup;%theta_c;
params_setting(18) = 18.68914715*thetaup;%theta_i; 
params_setting(19) = 3.216516185*thetaup;%theta_t;
params_setting(20) = 116.1426018*thetaup;%theta_p; 

flag = 1;
if flag==1
    global params_certainty;
    [I_ci,R_ci,Iconf, Rconf,params_all,I_mean,R_mean] = NSGA();
    params = mean(params_all,1);
    params_certainty = params_all;
elseif flag==0
    params = params_setting;
end

%NTL
res = light_daily(1,:);
com = light_daily(2,:);
ind = light_daily(3,:);
tsp = light_daily(4,:);
pub = light_daily(5,:);

res0 = light0_daily(1,:);
com0 = light0_daily(2,:);
ind0 = light0_daily(3,:);
tsp0 = light0_daily(4,:);
pub0 = light0_daily(5,:);

%% SEIR model
T=1:90;
lockdown_day = 23;
% Real world model: 2020 model
L=[res;com;ind;tsp;pub];
type_all = [1 1 1 1];
type_open = [0 0 0 0];
type_res = [1 0 0 0];
type_com = [0 1 0 0];
type_ind = [0 0 1 0];
type_pub = [0 0 0 1];
type_cip = [0 1 1 1];
type_rip = [1 0 1 1];
type_rcp = [1 1 0 1];
type_rci = [1 1 1 0];
type_ci = [0 1 1 0];
type_cp = [0 1 0 1];
type_ip = [0 0 1 1];

if flag==0
    [I_hat2020_medium,R_hat2020_medium,D2020_medium]=seir_real(T,params,L,lockdown_day,type_all);
else
    [I_hat2020_medium,R_hat2020_medium,D2020_medium]=seir_real(T,mean(params_all,1),L,lockdown_day,type_all);
end

%% Scenario1: 2019 no_beta + 2019 NTL
res0_1 = [res(1:lockdown_day-1),res0(lockdown_day:end)]; 
com0_1 = [com(1:lockdown_day-1),com0(lockdown_day:end)]; 
ind0_1 = [ind(1:lockdown_day-1),ind0(lockdown_day:end)];
tsp0_1 = [tsp(1:lockdown_day-1),tsp0(lockdown_day:end)]; 
pub0_1 = [pub(1:lockdown_day-1),pub0(lockdown_day:end)]; 
L=[res0_1;com0_1;ind0_1;tsp0_1;pub0_1];
[I_hat2019_nobeta_medium,R_hat2019_nobeta_medium,D2019_nobeta_medium]=seir_real(T,params,L,lockdown_day,type_open);

%% Scenario2: 2019 no_beta + 2019 NTL + LU=0
%res=0
res1_1 = res; res1_1(res0~=0) = 0; res1_1(1:lockdown_day-1)=res(1:lockdown_day-1);
L=[res1_1;com0_1;ind0_1;tsp0_1;pub0_1];
[I_hat2019_nobeta_res_medium,R_hat2019_nobeta_res_medium,D2019_nobeta_res_medium]=seir(T,params,L,lockdown_day,type_res);
%com=0
com1_1 = com; com1_1(com0~=0) = 0; com1_1(1:lockdown_day-1)=com(1:lockdown_day-1);
L=[res0_1;com1_1;ind0_1;tsp0_1;pub0_1];
[I_hat2019_nobeta_com_medium,R_hat2019_nobeta_com_medium,D2019_nobeta_com_medium]=seir(T,params,L,lockdown_day,type_com);
%ind=0
ind1_1 = ind; ind1_1(ind0~=0) = 0; ind1_1(1:lockdown_day-1)=ind(1:lockdown_day-1);
L=[res0_1;com0_1;ind1_1;tsp0_1;pub0_1];
[I_hat2019_nobeta_ind_medium,R_hat2019_nobeta_ind_medium,D2019_nobeta_ind_medium]=seir(T,params,L,lockdown_day,type_ind);
%pub=0
pub1_1 = pub; pub1_1(pub0~=0) = 0; pub1_1(1:lockdown_day-1)=pub(1:lockdown_day-1);
L=[res0_1;com0_1;ind0_1;tsp0_1;pub1_1];
[I_hat2019_nobeta_pub_medium,R_hat2019_nobeta_pub_medium,D2019_nobeta_pub_medium]=seir(T,params,L,lockdown_day,type_pub);
%rcip=0
L=[res1_1;com1_1;ind1_1;tsp0_1;pub1_1];
[I_hat2019_nobeta_rcip_medium,R_hat2019_nobeta_rcip_medium,D2019_nobeta_rcip_medium]=seir(T,params,L,lockdown_day,type_all);
%ci=0
L=[res0_1;com1_1;ind1_1;tsp0_1;pub0_1];
type_ci = [0 1 1 0];
[I_hat2019_nobeta_ci_medium,R_hat2019_nobeta_ci_medium,D2019_nobeta_ci_medium]= seir(T,params,L,lockdown_day,type_ci);
%cp=0
L=[res0_1;com1_1;ind0_1;tsp0_1;pub1_1];
type_cp = [0 1 0 1];
[I_hat2019_nobeta_cp_medium,R_hat2019_nobeta_cp_medium,D2019_nobeta_cp_medium]= seir(T,params,L,lockdown_day,type_cp);
%ip=0
L=[res0_1;com0_1;ind1_1;tsp0_1;pub1_1];
type_ip = [0 0 1 1];
[I_hat2019_nobeta_ip_medium,R_hat2019_nobeta_ip_medium,D2019_nobeta_ip_medium]= seir(T,params,L,lockdown_day,type_ip);
%cip=0
L=[res0_1;com1_1;ind1_1;tsp0_1;pub1_1];
type_cip = [0 1 1 1];
[I_hat2019_nobeta_cip_medium,R_hat2019_nobeta_cip_medium,D2019_nobeta_cip_medium]= seir(T,params,L,lockdown_day,type_cip);

%% Scenario3: 2019 no_beta + 2019 NTL + LU=max(2020,2019)
%res=0
res2_1 = res1_1; res2_1(lockdown_day:end)=max(res(lockdown_day:end),res0(lockdown_day:end));
L=[res2_1;com1_1;ind1_1;tsp0_1;pub1_1];
[I_hat2019_2020_res_medium,R_hat2019_2020_res_medium,D2019_2020_res_medium]=seir(T,params,L,lockdown_day,type_cip);
%com=0
com2_1 = com1_1; com2_1(lockdown_day:end)=max(com(lockdown_day:end),com0(lockdown_day:end));
L=[res2_1;com2_1;ind1_1;tsp0_1;pub1_1];
[I_hat2019_2020_com_medium,R_hat2019_2020_com_medium,D2019_2020_com_medium]=seir(T,params,L,lockdown_day,type_ip);
%ind=0
ind2_1 = ind1_1; ind2_1(lockdown_day:end)=max(ind(lockdown_day:end),ind0(lockdown_day:end));
L=[res2_1;com1_1;ind2_1;tsp0_1;pub1_1];
[I_hat2019_2020_ind_medium,R_hat2019_2020_ind_medium,D2019_2020_ind_medium]=seir(T,params,L,lockdown_day,type_cp);
%pub=0
pub2_1 = pub1_1; pub2_1(lockdown_day:end)=max(pub(lockdown_day:end),pub0(lockdown_day:end));
L=[res2_1;com1_1;ind1_1;tsp0_1;pub2_1];
[I_hat2019_2020_pub_medium,R_hat2019_2020_pub_medium,D2019_2020_pub_medium]=seir(T,params,L,lockdown_day,type_ci);
%rcip=0
L=[res2_1;com2_1;ind2_1;tsp0_1;pub2_1];
[I_hat2019_2020_rcip_medium,R_hat2019_2020_rcip_medium,D2019_2020_rcip_medium]=seir(T,params,L,lockdown_day,type_open);
%rci=0
L=[res2_1;com2_1;ind2_1;tsp0_1;pub1_1];
[I_hat2019_2020_rci_medium,R_hat2019_2020_rci_medium,D2019_2020_rci_medium]=seir(T,params,L,lockdown_day,type_pub);
%rcp=0
L=[res2_1;com2_1;ind1_1;tsp0_1;pub2_1];
[I_hat2019_2020_rcp_medium,R_hat2019_2020_rcp_medium,D2019_2020_rcp_medium]= seir(T,params,L,lockdown_day,type_ind);
%rip=0
L=[res2_1;com1_1;ind2_1;tsp0_1;pub2_1];
[I_hat2019_2020_rip_medium,R_hat2019_2020_rip_medium,D2019_2020_rip_medium]=seir(T,params,L,lockdown_day,type_com);

%border control
L=[res1_1;com1_1;ind1_1;tsp0_1;pub1_1];
[I_nobeta,R_nobeta]=seir_real(T,params,L,lockdown_day,type_open);

%% Visualization
%real world
Tbl         = readtable('cases.xlsx');
EstmInci    = Tbl.infection;
R_real      = Tbl.recovery;
size_font   = 15;
f=figure('color','w');
tiledlayout('flow','Padding','tight','TileSpacing','compact');
set(f,'position',[100 100 2000 320]);
% pre\post\opening: Onset
nexttile
data_period = datetime(2020,1,1):datetime(2020,3,30);
T=data_period;
xconf = [T T(end:-1:1)];
plot(T,EstmInci(1:90),'k.','MarkerSize',10,'DatetimeTickFormat','MMM dd');hold on;
plot(T,I_hat2020_medium,'b-','LineWidth',1.5,'DatetimeTickFormat','MMM dd');grid on;
plot(T,I_nobeta,'Color',[238,155,0]/255,'LineWidth',1.5,'DatetimeTickFormat','MMM dd');
plot(T,I_hat2019_nobeta_medium,'Color',[174,32,18]/255,'LineWidth',1.5,'DatetimeTickFormat','MMM dd');
fill([lockdown_day-1 95 95 lockdown_day-1],[0 0 3000 3000],'k','FaceAlpha',.1,'EdgeAlpha',0)
ylabel('# of Onset','FontSize',size_font);
xlim([T(1),T(end)]);
ylim([0 3000]);
xlabel('DAY','FontSize',size_font);
set(gca,'FontSize',size_font);
set(gca,'LineWidth',1.5);
% pre\post\opening: Removal
nexttile
plot(T,R_real(1:90),'k.','MarkerSize',10,'DatetimeTickFormat','MMM dd');hold on;
plot(T,R_hat2020_medium,'b-','LineWidth',1.5,'DatetimeTickFormat','MMM dd');grid on;
plot(T,R_nobeta,'Color',[238,155,0]/255,'LineWidth',1.5,'DatetimeTickFormat','MMM dd');
plot(T,R_hat2019_nobeta_medium,'Color',[174,32,18]/255,'LineWidth',1.5,'DatetimeTickFormat','MMM dd');
fill([lockdown_day-1 95 95 lockdown_day-1],[0 0 3000 3000],'k','FaceAlpha',.1,'EdgeAlpha',0)
ylabel('# of Removal','FontSize',size_font);
ylim([0 3000]);
xlim([T(1),T(end)]);
xlabel('DAY','FontSize',size_font);
set(gca,'FontSize',size_font);
set(gca,'LineWidth',1.5);
% Pre-NTL: res\com\pub
base = I_hat2019_2020_rcip_medium;
nexttile
I_diff = zeros(90,9);
for i=1:90
    I_diff(i,:) = [round(sum(base(1:i)))-round(sum(I_hat2019_nobeta_res_medium(1:i))),...
        round(sum(base(1:i)))-round(sum(I_hat2019_nobeta_rcip_medium(1:i))),...
        round(sum(base(1:i)))-round(sum(I_hat2019_nobeta_com_medium(1:i))),...
        round(sum(base(1:i)))-round(sum(I_hat2019_nobeta_ind_medium(1:i))),...
        round(sum(base(1:i)))-round(sum(I_hat2019_nobeta_pub_medium(1:i))),...
        round(sum(base(1:i)))-round(sum(I_hat2019_nobeta_ci_medium(1:i))),...
        round(sum(base(1:i)))-round(sum(I_hat2019_nobeta_cp_medium(1:i))),...
        round(sum(base(1:i)))-round(sum(I_hat2019_nobeta_ip_medium(1:i))),...
        round(sum(base(1:i)))-round(sum(I_hat2019_nobeta_cip_medium(1:i)))
        ];
end
bar(2:5,I_diff(90,2:5),'FaceColor',[27,158,119]/255,'EdgeColor',[27,158,119]/255,'LineWidth',1.5);
ylim([0,1.2*max(I_diff(90,2:5))]);hold on;
%xticklabels({'X COM','X IND','X PUB','X RCIP'});
%xlabel('Lockdown Certain Land Use','FontSize',size_font);
xticklabels({'','','',''});
ylabel('# of Infection Averted','FontSize',size_font);
set(gca,'FontSize',size_font);
set(gca,'LineWidth',1.5);
nexttile
% death_base = D2020_medium;
base2 = I_hat2019_nobeta_rcip_medium;
I_diff2 = zeros(90,8);
for i=1:90
    I_diff2(i,:) = [round(sum(base2(1:i)))-round(sum(I_hat2019_2020_res_medium(1:i))),...
        round(sum(base2(1:i)))-round(sum(I_hat2019_2020_com_medium(1:i))),...
        round(sum(base2(1:i)))-round(sum(I_hat2019_2020_ind_medium(1:i))),...
        round(sum(base2(1:i)))-round(sum(I_hat2019_2020_pub_medium(1:i))),...
        round(sum(base2(1:i)))-round(sum(I_hat2019_2020_rcip_medium(1:i))),...
        round(sum(base2(1:i)))-round(sum(I_hat2019_2020_rci_medium(1:i))),...
        round(sum(base2(1:i)))-round(sum(I_hat2019_2020_rcp_medium(1:i))),...
        round(sum(base2(1:i)))-round(sum(I_hat2019_2020_rip_medium(1:i)))
        ];
end
bar(1:4,-I_diff2(90,1:4),'FaceColor',[217,95,2]/255,'EdgeColor',[217,95,2]/255,'LineWidth',1.5);
ylim([0,1.2*max(-I_diff2(90,1:4))]);hold on;
xticklabels({'V RES','V COM+RES','V IND+RES','V PUB+RES','V RCIP'});
xlabel('Activate Certain Land Use','FontSize',size_font);
xticklabels({'','','',''});
ylabel('# of Excess Infection','FontSize',size_font);
set(gca,'FontSize',size_font);
set(gca,'LineWidth',1.5);