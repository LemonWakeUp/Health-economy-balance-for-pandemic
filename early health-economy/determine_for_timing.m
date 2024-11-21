%% health
% lockdown period mobility data
city_light = readtable("NTL2.xlsx");
city_light0 = readtable("NTL_Before2.xlsx");
city_outof = readtable("Mobility Distribution.xlsx");
city_into = readtable("Mobility Distribution.xlsx");
city_light = table2array(city_light(1:5,1:3));
city_light0 = table2array(city_light0(1:5,1:3));
city_outof = table2array(city_outof(2,2:93));
city_into = table2array(city_into(6,2:93));
city2 = readtable("Mobility.xlsx");
city2 = table2array(city2(1,2:93));

[factor_delta,factor] = distributeL(city_outof,city_into);
[factor_delta2,factor2] = distributeL2(city2);
light_daily = [city_light(:,1).*factor_delta{1},city_light(:,2).*factor_delta{2},city_light(:,3).*factor_delta{3}];
light0_daily = [city_light0(:,1).*factor_delta2{1},city_light0(:,2).*factor_delta2{2},city_light0(:,3).*factor_delta2{3}];

global L_0;
L_0=light_daily;

city_light_normal = readtable("NTL_normal.xlsx");
city_mobility_normal = readtable("Mobility_normal.xlsx");
city_light_normal = table2array([city_light_normal(4,2:13);city_light_normal(8,2:13);city_light_normal(12,2:13);city_light_normal(16,2:13);city_light_normal(20,2:13)]);
city_mobility_normal = table2array(city_mobility_normal(4,3:end));
[factor_delta_nor,factor_nor] = distributeL_normal(city_mobility_normal);
light_daily_nor = [city_light_normal(:,1).*factor_delta_nor{1},city_light_normal(:,2).*factor_delta_nor{2},city_light_normal(:,3).*factor_delta_nor{3},city_light_normal(:,4).*factor_delta_nor{4},city_light_normal(:,5).*factor_delta_nor{5},city_light_normal(:,6).*factor_delta_nor{6},city_light_normal(:,7).*factor_delta_nor{7},city_light_normal(:,8).*factor_delta_nor{8},city_light_normal(:,9).*factor_delta_nor{9},city_light_normal(:,10).*factor_delta_nor{10},city_light_normal(:,11).*factor_delta_nor{11},city_light_normal(:,12).*factor_delta_nor{12}];

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

res_n = light_daily_nor(1,:);
com_n = light_daily_nor(2,:);
ind_n = light_daily_nor(3,:);
tsp_n = light_daily_nor(4,:);
pub_n = light_daily_nor(5,:);

%% economy
load('damage_department.mat');
load('params.mat');
load('Shanghai_final2.mat');
load('Shanghai_mid2.mat');
load('Shanghai_total2.mat');
mid_produce = Shanghai_mid2;
final_produce = Shanghai_final2(:,1);
total_produce = Shanghai_total2';
import = sum(Shanghai_final2(:,2:3),2)';
export = zeros(42,1);
error = zeros(42,1);

%% general
%intermediate comsuption
mid_produce_IandI = sum(sum(mid_produce(2:27,2:27),1),2);
mid_produce_IandC = sum(sum(mid_produce(2:27,28:32),1),2);
mid_produce_IandP = sum(sum(mid_produce(2:27,33:42),1),2);

mid_produce_CandI = sum(sum(mid_produce(28:32,2:27),1),2);
mid_produce_CandC = sum(sum(mid_produce(28:32,28:32),1),2);
mid_produce_CandP = sum(sum(mid_produce(28:32,33:42),1),2);

mid_produce_PandI = sum(sum(mid_produce(33:42,2:27),1),2);
mid_produce_PandC = sum(sum(mid_produce(33:42,28:32),1),2);
mid_produce_PandP = sum(sum(mid_produce(33:42,33:42),1),2);

final_produce = sum(final_produce,2)+export+error-import';
final_produce_category = [sum(final_produce(2:27)),sum(final_produce(28:32)),sum(final_produce(33:42))];

total_produce_category = [sum(total_produce(2:27)),sum(total_produce(28:32)),sum(total_produce(33:42))];
mid_produce_category = [mid_produce_IandI,mid_produce_IandC,mid_produce_IandP;mid_produce_CandI,mid_produce_CandC,mid_produce_CandP;mid_produce_PandI,mid_produce_PandC,mid_produce_PandP];

%real damage in category
damage_category=[sum(total_produce(2:26)),sum(total_produce(27)),sum(total_produce(28)),sum(total_produce(30)),sum(total_produce(31)),sum(total_produce(32)),sum(total_produce(33)),sum(total_produce(34)),sum(total_produce(35:36)),sum(total_produce(37)),sum(total_produce(38)),sum(total_produce(39)),sum(total_produce(40)),sum(total_produce(41))].*damage_department;
damage = -[sum(damage_category(1:2))/sum(total_produce_category(1)),sum(damage_category(3:6))/sum(total_produce_category(2)),sum(damage_category(7:14))/sum(total_produce_category(3))];
damage = 1-power((1-damage),1/64);

a = mid_produce_category*inv(diag(total_produce_category));
a_p = (mid_produce_category'*inv(diag(total_produce_category)))';

%lockdown对总生产值的冲击
Q0=total_produce_category;
Y0=Q0;

%% SEIR model
T=1:(365+28);
lockdown_day = [0,14,28,42,56,70,84,98,112,126,140,154,168,196,224,280,336]+29;
lockdown_day2 = [0,14,28,42,56,70,84,98,112,126,140,154,168,196,224,280,336]+1;

%% Scenario1: 0_1:open 0_2:lockdown
res0_1 = [max(res0,res),res_n(93:end),res_n(1:28)]; 
com0_1 = [max(com0,com),com_n(93:end),com_n(1:28)]; 
ind0_1 = [max(ind0,ind),ind_n(93:end),ind_n(1:28)]; 
tsp0_1 = [max(tsp0,tsp),tsp_n(93:end),tsp_n(1:28)]; 
pub0_1 = [max(pub0,pub),pub_n(93:end),pub_n(1:28)]; 
L=[res0_1;com0_1;ind0_1;tsp0_1;pub0_1];
type_open=[0 0 0 0];
[I_hat2019_nobeta_medium,R_hat2019_nobeta_medium,D2019_nobeta_medium]=seir_timing(T,params,L,1,type_open);
res0_2 = zeros(1,365+28);
com0_2 = zeros(1,365+28);
ind0_2 = zeros(1,365+28);
tsp0_2 = zeros(1,365+28);
pub0_2 = zeros(1,365+28);

L_jud = [0,1,0,0,0;
    0,0,1,0,0;
    0,0,0,0,1;
    0,1,1,0,0;
    0,1,0,0,1;
    0,0,1,0,1;
    0,1,1,0,1;
    1,1,1,0,1];

iter=8;
n=1;
I_hat2019_nobeta_lock = zeros(1,365+28);
R_hat2019_nobeta_lock = zeros(1,365+28);
D2019_nobeta_lock = zeros(1,365+28);
loss_type = zeros(1,1000);
direct_loss_type= zeros(1,1000);
damage_cur = zeros(3,365+28);
time_min = 85;
res1_1 = res0_1;
com1_1 = com0_1;
ind1_1 = ind0_1;
tsp1_1 = tsp0_1;
pub1_1 = pub0_1;
type = zeros(1,4);
time_min = [];
time_max = [];
time_end = [];
max_produce = [1.007,1.007,1.007];
t_rec=30;

for lock_type = 1:8
    type = zeros(1,4);
    for lock_timing=1:length(lockdown_day)-1
        res1_1 = res0_1;
        com1_1 = com0_1;
        ind1_1 = ind0_1;
        tsp1_1 = tsp0_1;
        pub1_1 = pub0_1;
        time_min(n) = lockdown_day2(lock_timing);
        damage_cur = zeros(3,365+28);
        if L_jud(lock_type,1)==1
            type(1,1) = 1;
            res1_1(lockdown_day(1):lockdown_day(lock_timing))  = 0;
        end
        if L_jud(lock_type,2)==1
            type(1,2) = 1;
            com1_1(lockdown_day(1):lockdown_day(lock_timing))  = 0;
            damage_cur(2,1:lockdown_day2(lock_timing)) = damage(2);
        end
        if L_jud(lock_type,3)==1
            type(1,3) = 1;
            ind1_1(lockdown_day(1):lockdown_day(lock_timing))  = 0;
            damage_cur(1,1:lockdown_day2(lock_timing)) = damage(1);
        end
        if L_jud(lock_type,5)==1
            type(1,4) = 1;
            pub1_1(lockdown_day(1):lockdown_day(lock_timing))  = 0;
            damage_cur(3,1:lockdown_day2(lock_timing)) = damage(3);
        end
        L=[res1_1;com1_1;ind1_1;tsp0_1;pub1_1];
        time_max(n) = time_min(n)+365;
        [I_hat2019_nobeta_lock(n,:),R_hat2019_nobeta_lock(n,:),D2019_nobeta_lock(n,:),]=seir_timing(T,params,L,lock_timing,type);
        [loss_type(n,:),time_end(n)]=ARIO(damage_cur,Y0,Q0,a,final_produce_category,time_min(n),time_max(n),max_produce,t_rec);
        direct_loss_type(n,:)=direct_loss(damage_cur,Y0,Q0,a,time_min(n),time_end(n),max_produce,t_rec);
        n=n+1;
    end
end

%%full relaxation
L=[res0_1;com0_1;ind0_1;tsp0_1;pub0_1];
type=[0 0 0 0];
lock_timing=1;
[I_open,R_open,D_open]=seir_timing(T,params,L,lock_timing,type);