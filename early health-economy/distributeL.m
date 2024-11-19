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

