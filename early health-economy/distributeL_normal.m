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
