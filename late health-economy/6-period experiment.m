n=1; 
time_min = 168; 
Type_timing = zeros(1,6); 
time_end = [];
type = zeros(6,4); 

res1_1 = res0_1;
com1_1 = com0_1;
ind1_1 = ind0_1;
tsp1_1 = tsp0_1;
pub1_1 = pub0_1;

%setting state of each six period
for fir=1:iter
    res1_1(lockdown_day(1):lockdown_day(2))  = res0_1(lockdown_day(1):lockdown_day(2));
    com1_1(lockdown_day(1):lockdown_day(2))  = com0_1(lockdown_day(1):lockdown_day(2));
    ind1_1(lockdown_day(1):lockdown_day(2))  = ind0_1(lockdown_day(1):lockdown_day(2));
    pub1_1(lockdown_day(1):lockdown_day(2))  = pub0_1(lockdown_day(1):lockdown_day(2));
    damage_cur(:,lockdown_day2(1):lockdown_day2(2)) = 0;
    type(1,:) =0;
    Type_timing(n,1)=fir;
    if fir == 1 
        time_min = lockdown_day2(1);
    end
    if L_jud(fir,1)==1
        res1_1(lockdown_day(1):lockdown_day(2))  = 0;
        type(1,1) = 1;
    end
    if L_jud(fir,2)==1
        com1_1(lockdown_day(1):lockdown_day(2))  = 0;
        type(1,2) = 1;
        damage_cur(2,lockdown_day2(1):lockdown_day2(2)) = damage(2);
    end
    if L_jud(fir,3)==1
        ind1_1(lockdown_day(1):lockdown_day(2))  = 0;
        type(1,3) = 1;
        damage_cur(1,lockdown_day2(1):lockdown_day2(2)) = damage(1);
    end
    if L_jud(fir,5)==1
        pub1_1(lockdown_day(1):lockdown_day(2))  = 0;
        type(1,4) = 1;
        damage_cur(3,lockdown_day2(1):lockdown_day2(2)) = damage(3);
    end
    for sec=1:iter
        res1_1(lockdown_day(2):lockdown_day(3))  = res0_1(lockdown_day(2):lockdown_day(3));
        com1_1(lockdown_day(2):lockdown_day(3))  = com0_1(lockdown_day(2):lockdown_day(3));
        ind1_1(lockdown_day(2):lockdown_day(3))  = ind0_1(lockdown_day(2):lockdown_day(3));
        pub1_1(lockdown_day(2):lockdown_day(3))  = pub0_1(lockdown_day(2):lockdown_day(3));
        damage_cur(:,lockdown_day2(2):lockdown_day2(3)) = 0;
        type(2,:) =0;
        Type_timing(n,2)=sec;
        Type_timing(n,1)=fir;
        if sec == 1 && fir~=1
            time_min = lockdown_day2(2);
        end
        if L_jud(sec,1)==1
            type(2,1) = 1;
            res1_1(lockdown_day(2):lockdown_day(3))  = 0;
        end
        if L_jud(sec,2)==1
            com1_1(lockdown_day(2):lockdown_day(3))  = 0;
            type(2,2) = 1;
            damage_cur(2,lockdown_day2(2):lockdown_day2(3)) = damage(2);
        end
        if L_jud(sec,3)==1
            ind1_1(lockdown_day(2):lockdown_day(3))  = 0;
            type(2,3) = 1;
            damage_cur(1,lockdown_day2(2):lockdown_day2(3)) = damage(1);
        end
        if L_jud(sec,5)==1
            pub1_1(lockdown_day(2):lockdown_day(3))  = 0;
            type(2,4) = 1;
            damage_cur(3,lockdown_day2(2):lockdown_day2(3)) = damage(3);
        end
        for thir=1:iter
            res1_1(lockdown_day(3):lockdown_day(4))  = res0_1(lockdown_day(3):lockdown_day(4));
            com1_1(lockdown_day(3):lockdown_day(4))  = com0_1(lockdown_day(3):lockdown_day(4));
            ind1_1(lockdown_day(3):lockdown_day(4))  = ind0_1(lockdown_day(3):lockdown_day(4));
            pub1_1(lockdown_day(3):lockdown_day(4))  = pub0_1(lockdown_day(3):lockdown_day(4));
            damage_cur(:,lockdown_day2(3):lockdown_day2(4)) = 0;
            type(3,:) =0;
            Type_timing(n,2)=sec;
            Type_timing(n,1)=fir;
            Type_timing(n,3)=thir;
            if thir == 1 && sec~=1
                time_min = lockdown_day2(3);
            end
            if L_jud(thir,1)==1
                res1_1(lockdown_day(3):lockdown_day(4))  = 0;
                type(3,1) = 1;
            end
            if L_jud(thir,2)==1
                com1_1(lockdown_day(3):lockdown_day(4))  = 0;
                type(3,2) = 1;
                damage_cur(2,lockdown_day2(3):lockdown_day2(4)) = damage(2);
            end
            if L_jud(thir,3)==1
                ind1_1(lockdown_day(3):lockdown_day(4))  = 0;
                type(3,3) = 1;
                damage_cur(1,lockdown_day2(3):lockdown_day2(4)) = damage(1);
            end
            if L_jud(thir,5)==1
                pub1_1(lockdown_day(3):lockdown_day(4))  = 0;
                type(3,4) = 1;
                damage_cur(3,lockdown_day2(3):lockdown_day2(4)) = damage(3);
            end
            for four=1:iter
                res1_1(lockdown_day(4):lockdown_day(5))  = res0_1(lockdown_day(4):lockdown_day(5));
                com1_1(lockdown_day(4):lockdown_day(5))  = com0_1(lockdown_day(4):lockdown_day(5));
                ind1_1(lockdown_day(4):lockdown_day(5))  = ind0_1(lockdown_day(4):lockdown_day(5));
                pub1_1(lockdown_day(4):lockdown_day(5))  = pub0_1(lockdown_day(4):lockdown_day(5));
                damage_cur(:,lockdown_day2(4):lockdown_day2(5)) = 0;
                type(4,:) =0;
                Type_timing(n,4)=four;
                Type_timing(n,2)=sec;
                Type_timing(n,1)=fir;
                Type_timing(n,3)=thir;
                if four == 1 && thir~=1
                    time_min = lockdown_day2(4);
                end
                if L_jud(four,1)==1
                    res1_1(lockdown_day(4):lockdown_day(5))  = 0;
                    type(4,1) = 1;
                end
                if L_jud(four,2)==1
                    com1_1(lockdown_day(4):lockdown_day(5))  = 0;
                    type(4,2) = 1;
                    damage_cur(2,lockdown_day2(4):lockdown_day2(5)) = damage(2);
                end
                if L_jud(four,3)==1
                    ind1_1(lockdown_day(4):lockdown_day(5))  = 0;
                    type(4,3) = 1;
                    damage_cur(1,lockdown_day2(4):lockdown_day2(5)) = damage(1);
                end
                if L_jud(four,5)==1
                    pub1_1(lockdown_day(4):lockdown_day(5))  = 0;
                    type(4,4) = 1;
                    damage_cur(3,lockdown_day2(4):lockdown_day2(5)) = damage(3);
                end
                for fif=1:iter
                    res1_1(lockdown_day(5):lockdown_day(6))  = res0_1(lockdown_day(5):lockdown_day(6));
                    com1_1(lockdown_day(5):lockdown_day(6))  = com0_1(lockdown_day(5):lockdown_day(6));
                    ind1_1(lockdown_day(5):lockdown_day(6))  = ind0_1(lockdown_day(5):lockdown_day(6));
                    pub1_1(lockdown_day(5):lockdown_day(6))  = pub0_1(lockdown_day(5):lockdown_day(6));
                    damage_cur(:,lockdown_day2(5):lockdown_day2(6)) = 0;
                    type(5,:) =0;
                    Type_timing(n,5)=fif;
                    Type_timing(n,4)=four;
                    Type_timing(n,2)=sec;
                    Type_timing(n,1)=fir;
                    Type_timing(n,3)=thir;
                    if fif == 1 && four~=1
                        time_min = lockdown_day2(5);
                    end
                    if L_jud(fif,1)==1
                        res1_1(lockdown_day(5):lockdown_day(6))  = 0;
                        type(5,1) = 1;
                    end
                    if L_jud(fif,2)==1
                        com1_1(lockdown_day(5):lockdown_day(6))  = 0;
                        type(5,2) = 1;
                        damage_cur(2,lockdown_day2(5):lockdown_day2(6)) = damage(2);
                    end
                    if L_jud(fif,3)==1
                        ind1_1(lockdown_day(5):lockdown_day(6))  = 0;
                        type(5,3) = 1;
                        damage_cur(1,lockdown_day2(5):lockdown_day2(6)) = damage(1);
                    end
                    if L_jud(fif,5)==1
                        pub1_1(lockdown_day(5):lockdown_day(6))  = 0;
                        type(5,4) = 1;
                        damage_cur(3,lockdown_day2(5):lockdown_day2(6)) = damage(3);
                    end
                    for six=1:iter
                        res1_1(lockdown_day(6):lockdown_day(7))  = res0_1(lockdown_day(6):lockdown_day(7));
                        com1_1(lockdown_day(6):lockdown_day(7))  = com0_1(lockdown_day(6):lockdown_day(7));
                        ind1_1(lockdown_day(6):lockdown_day(7))  = ind0_1(lockdown_day(6):lockdown_day(7));
                        pub1_1(lockdown_day(6):lockdown_day(7))  = pub0_1(lockdown_day(6):lockdown_day(7));
                        damage_cur(:,lockdown_day2(6):lockdown_day2(7)) = 0;
                        type(6,:) =0;
                        Type_timing(n,6)=six;
                        Type_timing(n,5)=fif;
                        Type_timing(n,4)=four;
                        Type_timing(n,2)=sec;
                        Type_timing(n,1)=fir;
                        Type_timing(n,3)=thir;
                        if six == 1 && fif~=1
                            time_min = lockdown_day2(6);
                        end
                        if L_jud(six,1)==1
                            res1_1(lockdown_day(6):lockdown_day(7))  = 0;
                            type(6,1) = 1;
                        end
                        if L_jud(six,2)==1
                            com1_1(lockdown_day(6):lockdown_day(7))  = 0;
                            type(6,2) = 1;
                            damage_cur(2,lockdown_day2(6):lockdown_day2(7)) = damage(2);
                        end
                        if L_jud(six,3)==1
                            ind1_1(lockdown_day(6):lockdown_day(7))  = 0;
                            type(6,3) = 1;
                            damage_cur(1,lockdown_day2(6):lockdown_day2(7)) = damage(1);
                        end
                        if L_jud(six,5)==1
                            pub1_1(lockdown_day(6):lockdown_day(7))  = 0;
                            type(6,4) = 1;
                            damage_cur(3,lockdown_day2(6):lockdown_day2(7)) = damage(3);
                        end
                        L=[res1_1;com1_1;ind1_1;tsp0_1;pub1_1];
                        T=1:393;
                        [I_hat2019_nobeta_lock(n,:),R_hat2019_nobeta_lock(n,:),D2019_nobeta_lock(n,:),S_exist(n,:),E_exist(n,:),I_exist(n,:),R_exist(n,:)]=seir_timing(T,params,L,type);
                        [loss_type(n,1:1000)]=ARIO(damage_cur,Y0,Q0,a,final_produce_category,time_min,time_min+365,max_produce,t_rec);
                        n=n+1;
                    end
                end
            end
        end
    end
end