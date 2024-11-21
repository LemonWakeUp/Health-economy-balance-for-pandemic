function [onset,R_hat,D,S,E,I,R]=seir_real(T,params,L,lockday,type)
population  = 9319700;
S           = zeros(1,length(T));
E           = zeros(1,length(T));
I           = zeros(1,length(T));
R           = zeros(1,length(T));
D           = zeros(1,length(T));
I(1,1)      = 66;
E(1,:)      = 10*I(1,:);
S(1,:)      = population - E(1,:) - I(1,:);
%params
k_1             = 0.1;
k_2             = 2;
k_3             = 1;
beta0           = params(1);
alpha1          = params(2);
symp            = params(3);
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

load('xishu.mat','xishu');
load('NTL_weight_WH.mat','NTL_weight_WH');
NTL_pro = sum(NTL_weight_WH.*type')./sum(NTL_weight_WH);
% a = 1./(L_beta_cfr*1);
% L_beta_cfr = 1./(a.*NTL_pro);

gamma2_model    = zeros(size(T));
R_hat           = zeros(size(T));
R_cum           = zeros(size(T));
start_date = 2; 
end_date = 90;
for doy = start_date : end_date
    gamma2_model(doy)=gamma20*(1-exp(-(T(doy)/L_g).^k_g));
    beta = beta0*(k_3-k_1/2+k_1/2*cos(k_2*pi/365*T(doy)));
    if doy~=0 && lockday<=doy
        beta_min =  min(beta*exp(-(T(doy)/L_beta_cfr).^k_beta),beta);
        beta = beta-(beta-beta_min)*NTL_pro;
    end
    if doy>0
        p_infectious = 1-p_no_infectious*exp(-(T(doy)/L_p).^k_p);
        S(doy) = S(doy-1) - beta*p_infectious*(I(doy-1)+E(doy-1))/population*S(doy-1) +...
            proportionS*(theta_r*res(doy)+theta_c*com(doy)+theta_i*ind(doy)+theta_t*tsp(doy)+theta_p*pub(doy));
        E(doy) = E(doy-1) + beta*p_infectious*(I(doy-1)+E(doy-1))/population*S(doy-1) - (symp*alpha1 + (1-symp)*gamma1)*E(doy-1)+...
            proportionE*(theta_r*res(doy)+theta_c*com(doy)+theta_i*ind(doy)+theta_t*tsp(doy)+theta_p*pub(doy));
        I(doy) = I(doy-1) + symp*alpha1*E(doy-1) - gamma2_model(doy)*I(doy-1);
        D(doy) = gamma2_model(doy)*I(doy-1)*cfr;
        R(doy) = R(doy-1) + (1-symp)*gamma1*E(doy-1) + (gamma2_model(doy)*I(doy-1) - D(doy));
        population = population - D(doy)+...
            proportionS*(theta_r*res(doy)+theta_c*com(doy)+theta_i*ind(doy)+theta_t*tsp(doy)+theta_p*pub(doy))+...
            proportionE*(theta_r*res(doy)+theta_c*com(doy)+theta_i*ind(doy)+theta_t*tsp(doy)+theta_p*pub(doy))+...
            proportionR*(theta_r*res(doy)+theta_c*com(doy)+theta_i*ind(doy)+theta_t*tsp(doy)+theta_p*pub(doy));
        R_hat(doy) = gamma2_model(doy)*I(doy-1)-D(doy);
        %+(1-symp_pro)*gamma1*E(doy-1);
        if R_hat(doy)<0
            R_hat(doy)=0;
        elseif R_hat(doy)<=0
            R_hat(doy) = 0;
        end
    end
    R_cum(doy) = R_cum(doy-1)+R_hat(doy);
    checkState = [S(doy),E(doy),I(doy),D(doy),R(doy),population];
    checkState(checkState(:)<0)=0;
    S(doy)=checkState(1);E(doy)=checkState(2);I(doy)=checkState(3);D(doy)=checkState(4);R(doy)=checkState(5);
end
onset   = symp*alpha1*E;
end
