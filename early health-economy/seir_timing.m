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

