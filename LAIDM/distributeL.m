function [factor_delta,factor] = distributeL(wuhan_outof,wuhan_into)
% in&out
data_period = datetime(2020,1,1):datetime(2020,3,23);
% population+in-out
wuhan_delta_pop = zeros(1,size(data_period,2));
for day_ = 1:size(wuhan_delta_pop,2)
     wuhan_delta_pop(1,day_)=-sum(wuhan_outof(day_,:))+sum(wuhan_into(:,day_));
end
% delta distribution factor
data_period = datetime(2020,1,1):datetime(2020,3,31);
delta_ = [wuhan_delta_pop,wuhan_delta_pop(1,77:83)];
jan_delta = delta_(data_period<datetime(2020,2,1));
feb_delta = delta_(datetime(2020,1,31)<data_period & data_period<datetime(2020,2,29));
mar_delta = delta_(datetime(2020,2,28)<data_period& data_period<datetime(2020,3,31));
jan_factor_delta = jan_delta/sum(jan_delta);
feb_factor_delta = feb_delta/sum(feb_delta);
mar_factor_delta = mar_delta/sum(mar_delta);
factor = delta_;
factor_delta = {jan_factor_delta,feb_factor_delta,mar_factor_delta};
end