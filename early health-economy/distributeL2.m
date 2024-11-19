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
