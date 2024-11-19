function [factor_delta2,factor2] = distributeL2(wuhan2019)
data_period = 1:90;
delta_ = wuhan2019;
jan_delta = delta_(data_period<32);
feb_delta = delta_(31<data_period & data_period<60);
mar_delta = delta_(59<data_period);
jan_factor_delta = jan_delta';
feb_factor_delta = feb_delta';
mar_factor_delta = mar_delta';
factor2 =delta_;
factor_delta2 = {jan_factor_delta,feb_factor_delta,mar_factor_delta};
end