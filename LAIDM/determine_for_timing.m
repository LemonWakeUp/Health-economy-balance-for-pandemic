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

params_setting(1)=0.455737;
params_setting(2)=0.229054;
params_setting(3)=0.092325;
params_setting(4)=0.215386;
params_setting(5)=0.116980;
params_setting(6)=0.049888;
params_setting(7)=52.258318;
params_setting(8)=29.114688;
params_setting(9)=55.705386;
params_setting(10)=3.104142;
params_setting(11)=3.179818;
params_setting(12)=5.125632;
params_setting(13)=0.001286;
params_setting(14)=0.3;
params_setting(15)=0.999;
theta=1;
params_setting(16)=38.623846*theta;
params_setting(17)=79.826101*theta;
params_setting(18)=77.159742*theta;
params_setting(19)=38.831442*theta;
params_setting(20)=43.086164*theta;

params = params_setting;

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
res0_1 = [max(res0,res),res_n(93:end),res_n(1:28)]; %res0_1(1:23) = res(1+8:23+8);
com0_1 = [max(com0,com),com_n(93:end),com_n(1:28)]; %com0_1(1:23) = com(1+8:23+8);
ind0_1 = [max(ind0,ind),ind_n(93:end),ind_n(1:28)]; %ind0_1(1:23) = ind(1+8:23+8);
tsp0_1 = [max(tsp0,tsp),tsp_n(93:end),tsp_n(1:28)]; %tsp0_1(1:23) = tsp(1+8:23+8);
pub0_1 = [max(pub0,pub),pub_n(93:end),pub_n(1:28)]; %pub0_1(1:23) = pub(1+8:23+8);
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

%% function
function [onset,R_hat,D,S,E,I,R]=seir_timing(T,params,L,lock_timing,type)
population  = 24758900;
S           = zeros(1,length(T));
E           = zeros(1,length(T));
I           = zeros(1,length(T));
R           = zeros(1,length(T));
D           = zeros(1,length(T));
I(1,1)      = 11;
E(1,:)      = 10*I(1,:);
S(1,:)      = population - E(1,:) - I(1,:);
%params
k_1             = 0.1;
k_2             = 2;
k_3             = 1;
beta0           = params(1);
alpha1          = params(2);
symp_pro         = params(3);
gamma1          = params(4);
gamma20         = params(5);
p_no_infectious = params(6);
L_beta_cfr      = params(7);
L_p             = params(8);
L_g             = params(9);
k_beta          = params(10);
k_p             = params(11);
k_g             = params(12);
cfr             = params(13);
proportion      = params(14);
pSE             = params(15);
proportionS     = pSE*proportion;
proportionE     = (1-pSE)*proportion;
proportionR     = 1-proportion;
theta_r         = params(16);
theta_c         = params(17);
theta_i         = params(18);
theta_t         = params(19);
theta_p         = params(20);
res             = L(1,:);
com             = L(2,:);
ind             = L(3,:);
tsp             = L(4,:);
pub             = L(5,:);

lockday = [0,14,28,42,56,70,84,98,112,126,140,154,168,196,224,280,336]+28;

load('xishu.mat','xishu');
load('NTL_weight_SH.mat','NTL_weight_SH');
NTL_pro = sum(NTL_weight_SH.*type'.*xishu')./sum(NTL_weight_SH.*xishu');

% a = 1./(L_beta_cfr*1);
% L_beta_cfr = 1./(a.*NTL_pro);

gamma2_model    = zeros(size(T));
E_hat           = zeros(size(T));
R_hat           = zeros(size(T));
R_cum           = zeros(size(T));

start_date = 2; 
end_date = length(T);
for doy = start_date : end_date
    gamma2_model(doy)=gamma20*(1-exp(-(T(doy)/L_g).^k_g));
    if doy~=0 && lockday(1)<doy && doy <= lockday(lock_timing) && NTL_pro~=0
        beta =  min(beta0*exp(-(T(doy)/L_beta_cfr).^k_beta),beta0);
        beta = beta0-(beta0-beta)*NTL_pro;
    else
        beta = beta0*(k_3-k_1/2+k_1/2*cos(k_2*pi/365*T(doy)));
        %  beta = beta0*(k_3-k_1*sin(k_2*pi/length(T)*T(doy)));
    end
    % proportionS     = S(doy-1)/population;
    % proportionE     = E(doy-1)/population;
    % proportionR     = R(doy-1)/population;
    if doy>0
        p_infectious = 1-p_no_infectious*exp(-(T(doy)/L_p).^k_p);
        S(doy) = S(doy-1) - beta*p_infectious*(I(doy-1)+E(doy-1))/population*S(doy-1) +...
            proportionS*(theta_r*res(doy)+theta_c*com(doy)+theta_i*ind(doy)+theta_t*tsp(doy)+theta_p*pub(doy));
        E(doy) = E(doy-1) + beta*p_infectious*(I(doy-1)+E(doy-1))/population*S(doy-1) - (symp_pro*alpha1 + (1-symp_pro)*gamma1)*E(doy-1)+...
            proportionE*(theta_r*res(doy)+theta_c*com(doy)+theta_i*ind(doy)+theta_t*tsp(doy)+theta_p*pub(doy));
        E_hat(doy)=beta*p_infectious*(I(doy-1)+E(doy-1))/population*S(doy-1);
        I(doy) = I(doy-1) + symp_pro*alpha1*E(doy-1) - gamma20*I(doy-1);
        D(doy) = gamma20*I(doy-1)*cfr;
        R(doy) = R(doy-1) + (1-symp_pro)*gamma1*E(doy-1) + (gamma20*I(doy-1) - D(doy)) +...
            proportionR*(theta_r*res(doy)+theta_c*com(doy)+theta_i*ind(doy)+theta_t*tsp(doy)+theta_p*pub(doy));
        population = population - D(doy)+...
            proportionS*(theta_r*res(doy)+theta_c*com(doy)+theta_i*ind(doy)+theta_t*tsp(doy)+theta_p*pub(doy))+...
            proportionE*(theta_r*res(doy)+theta_c*com(doy)+theta_i*ind(doy)+theta_t*tsp(doy)+theta_p*pub(doy))+...
            proportionR*(theta_r*res(doy)+theta_c*com(doy)+theta_i*ind(doy)+theta_t*tsp(doy)+theta_p*pub(doy));
        R_hat(doy) = gamma20*I(doy-1) - D(doy);
        %+(1-symp_pro)*gamma1*E(doy-1);
        if R_hat(doy)<0
            R_hat(doy)=0;
        elseif R_hat(doy)<=0
            R_hat(doy) = 0;
        end
    end
    R_cum(doy)=R_cum(doy-1)+R_hat(doy);
    checkState = [S(doy),E(doy),I(doy),D(doy),R(doy),population];
    checkState(checkState(:)<0)=0;
    S(doy)=checkState(1);E(doy)=checkState(2);I(doy)=checkState(3);D(doy)=checkState(4);R(doy)=checkState(5);
end
%onset   = alpha1*E+(1-symp_pro)*E_hat;
onset = symp_pro*alpha1*E+(1-symp_pro)*E_hat;
Ia      = gamma1*E;
end


function [factor_delta,factor] = distributeL(city_outof,city_into)
% in&out
data_period = datetime(2022,3,1):datetime(2022,5,31);
% population+in-out
city_delta_pop = zeros(size(data_period,2));
for day_ = 1:size(city_delta_pop,2)
     city_delta_pop(day_)=city_into(day_)-city_outof(day_);
end
% delta distribution factor
fir_delta = city_delta_pop(1:31);
sec_delta = city_delta_pop(32:61);
thir_delta = city_delta_pop(62:92);
fir_factor_delta = fir_delta/sum(fir_delta);
sec_factor_delta = sec_delta/sum(sec_delta);
thir_factor_delta = thir_delta/sum(thir_delta);
factor = [fir_factor_delta,sec_factor_delta,thir_factor_delta];
factor_delta = {fir_factor_delta,sec_factor_delta,thir_factor_delta};
end

function [factor_delta,factor] = distributeL2(city_delta_pop)
% in&out
data_period = datetime(2022,3,1):datetime(2022,5,31);
% delta distribution factor
fir_delta = city_delta_pop(1:31);
sec_delta = city_delta_pop(32:61);
thir_delta = city_delta_pop(62:92);
fir_factor_delta = fir_delta/sum(fir_delta);
sec_factor_delta = sec_delta/sum(sec_delta);
thir_factor_delta = thir_delta/sum(thir_delta);
factor = [fir_factor_delta,sec_factor_delta,thir_factor_delta];
factor_delta = {fir_factor_delta,sec_factor_delta,thir_factor_delta};
end


function [factor_delta,factor] = distributeL_normal(city_delta_pop)
% in&out
data_period = datetime(2021,3,1):datetime(2022,2,28);
% delta distribution factor
thir_delta = city_delta_pop(data_period>=datetime(2021,3,1)&data_period<=datetime(2021,3,31));
four_delta = city_delta_pop(data_period>=datetime(2021,4,1)&data_period<=datetime(2021,4,30));
fif_delta = city_delta_pop(data_period>=datetime(2021,5,1)&data_period<=datetime(2021,5,31));
six_delta = city_delta_pop(data_period>=datetime(2021,6,1)&data_period<=datetime(2021,6,30));
seven_delta = city_delta_pop(data_period>=datetime(2021,7,1)&data_period<=datetime(2021,7,31));
eight_delta = city_delta_pop(data_period>=datetime(2021,8,1)&data_period<=datetime(2021,8,31));
nine_delta = city_delta_pop(data_period>=datetime(2021,9,1)&data_period<=datetime(2021,9,30));
ten_delta = city_delta_pop(data_period>=datetime(2021,10,1)&data_period<=datetime(2021,10,31));
ele_delta = city_delta_pop(data_period>=datetime(2021,11,1)&data_period<=datetime(2021,11,30));
twe_delta = city_delta_pop(data_period>=datetime(2021,12,1)&data_period<=datetime(2021,12,31));
fir_delta = city_delta_pop(data_period>=datetime(2022,1,1)&data_period<=datetime(2022,1,31));
sec_delta = city_delta_pop(data_period>=datetime(2022,2,1)&data_period<=datetime(2022,2,28));
thir_factor_delta = thir_delta/sum(thir_delta);
four_factor_delta = four_delta/sum(four_delta);
fif_factor_delta = fif_delta/sum(fif_delta);
six_factor_delta = six_delta/sum(six_delta);
seven_factor_delta = seven_delta/sum(seven_delta);
eight_factor_delta = eight_delta/sum(eight_delta);
nine_factor_delta = nine_delta/sum(nine_delta);
ten_factor_delta = ten_delta/sum(ten_delta);
ele_factor_delta = ele_delta/sum(ele_delta);
twe_factor_delta = twe_delta/sum(twe_delta);
fir_factor_delta = fir_delta/sum(fir_delta);
sec_factor_delta = sec_delta/sum(sec_delta);
factor = [thir_factor_delta,four_factor_delta,fif_factor_delta,six_factor_delta,seven_factor_delta,eight_factor_delta,nine_factor_delta,ten_factor_delta,ele_factor_delta,twe_factor_delta,fir_factor_delta,sec_factor_delta];
factor_delta = {thir_factor_delta,four_factor_delta,fif_factor_delta,six_factor_delta,seven_factor_delta,eight_factor_delta,nine_factor_delta,ten_factor_delta,ele_factor_delta,twe_factor_delta,fir_factor_delta,sec_factor_delta};
end

function curve=filtered_curvature(x,y,N,h,judge)
x_mean = mean(x);
y_mean = mean(y);
x_diff = x-x_mean;
y_diff = y-y_mean;
Suuu = sum(power(x_diff,3));
Svvv = sum(power(y_diff,3));
Suu = sum(power(x_diff,2));
Svv = sum(power(y_diff,2));
Suuv = sum(power(x_diff,2).*y_diff);
Suvv = sum(power(y_diff,2).*x_diff);
Suv = sum(x_diff.*y_diff);
center_x = (Suuv*Suv-Suuu*Svv-Suvv*Svv+Suv*Svvv)/2*(Suv*Suv-Suu*Svv)+x_mean;
center_y = (-Suu*Suuv+Suuu*Suv+Suv*Suvv-Suu*Svvv)/2*(Suv*Suv-Suu*Svv)+y_mean;
center = [center_x,center_y];
curve = 1/sqrt(sum(power(x-center_x,2)+power(y-center_y,2))/N);
x_middle = x((N+1)/2);
y_middle = y((N+1)/2);
theta = asin(sqrt(power(x_middle-center_x,2)+power(y_middle-center_y,2))/(x_middle-center_x))/pi*180;
if judge==1
    curve=0;
else
    if ~(y_middle<=h && theta<=0 && theta>=-90)
        curve=0;
    end
end
end
