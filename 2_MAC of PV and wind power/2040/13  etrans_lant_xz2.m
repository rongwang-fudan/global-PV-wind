tic
clear;
load('H:\global-PV-wind\ANS\choo_type_8_2040xz.mat'); % choo
load('H:\global-PV-wind\ANS\etrans_cou1_num_1023_pro2_testt_8_2070.mat')% 实际有效发电量, TWh/year,etrans_cou1_num1
load('H:\global-PV-wind\ANS\etrans_cou_num_1023_pro2_testt_8_2070.mat')% 全部发电量, TWh/year,etrans_cou_num1
etrans_cou1_num1_xz = etrans_cou1_num1*0;
etrans_cou_num1_xz = etrans_cou_num1*0;
[m,n]=find(choo==1);
etrans_cou1_num1_xz(m,:,:) = etrans_cou1_num1(m,:,:);
etrans_cou_num1_xz(m,:,:) = etrans_cou_num1(m,:,:);
clear etrans_cou1_num1
clear etrans_cou_num1

load('H:\global-PV-wind\ANS\etrans_cou1_num_1023_pro2_testt_interUHV_8_2070.mat')% 实际有效发电量, TWh/year,etrans_cou1_num1
load('H:\global-PV-wind\ANS\etrans_cou_num_1023_pro2_testt_interUHV_8_2070.mat')% 全部发电量, TWh/year,etrans_cou_num1
[m,n]=find(choo==2);
etrans_cou1_num1_xz(m,:,:) = etrans_cou1_num(m,:,:);
etrans_cou_num1_xz(m,:,:) = etrans_cou_num(m,:,:);
clear etrans_cou1_num
clear etrans_cou_num

load('H:\global-PV-wind\ANS\etrans_cou1_num_1023_pro2_testt_sto_8_2070.mat')% 实际有效发电量, TWh/year,etrans_cou1_num1
load('H:\global-PV-wind\ANS\etrans_cou_num_1023_pro2_testt_sto_8_2070.mat')% 全部发电量, TWh/year,etrans_cou_num1
[m,n]=find(choo==3);
etrans_cou1_num1_xz(m,:,:) = etrans_cou1_num(m,:,:);
etrans_cou_num1_xz(m,:,:) = etrans_cou_num(m,:,:);
clear etrans_cou1_num
clear etrans_cou_num

load('H:\global-PV-wind\ANS\etrans_cou1_num_1023_pro2_testt_UHV_8_2070.mat')% 实际有效发电量, TWh/year,etrans_cou1_num1
load('H:\global-PV-wind\ANS\etrans_cou_num_1023_pro2_testt_UHV_8_2070.mat')% 全部发电量, TWh/year,etrans_cou_num1
[m,n]=find(choo==4);
etrans_cou1_num1_xz(m,:,:) = etrans_cou1_num1(m,:,:);
etrans_cou_num1_xz(m,:,:) = etrans_cou_num1(m,:,:);
clear etrans_cou1_num1
clear etrans_cou_num1

load('H:\global-PV-wind\ANS\etrans_cou1_num_1023_pro2_testt_UHV_interUHV_8_2070.mat')% 实际有效发电量, TWh/year,etrans_cou1_num1
load('H:\global-PV-wind\ANS\etrans_cou_num_1023_pro2_testt_UHV_interUHV_8_2070.mat')% 全部发电量, TWh/year,etrans_cou_num1
[m,n]=find(choo==5);
etrans_cou1_num1_xz(m,:,:) = etrans_cou1_num(m,:,:);
etrans_cou_num1_xz(m,:,:) = etrans_cou_num(m,:,:);
clear etrans_cou1_num
clear etrans_cou_num

load('H:\global-PV-wind\ANS\etrans_cou1_num_1023_pro2_testt_sto_interUHV_8_2070.mat')% 实际有效发电量, TWh/year,etrans_cou1_num1
load('H:\global-PV-wind\ANS\etrans_cou_num_1023_pro2_testt_sto_interUHV_8_2070.mat')% 全部发电量, TWh/year,etrans_cou_num1
[m,n]=find(choo==6);
etrans_cou1_num1_xz(m,:,:) = etrans_cou1_num(m,:,:);
etrans_cou_num1_xz(m,:,:) = etrans_cou_num(m,:,:);
clear etrans_cou1_num
clear etrans_cou_num

load('H:\global-PV-wind\ANS\etrans_cou1_num_1023_pro2_testt_UHV_sto_8_2070.mat')% 实际有效发电量, TWh/year,etrans_cou1_num1
load('H:\global-PV-wind\ANS\etrans_cou_num_1023_pro2_testt_UHV_sto_8_2070.mat')% 全部发电量, TWh/year,etrans_cou_num1
[m,n]=find(choo==7);
etrans_cou1_num1_xz(m,:,:) = etrans_cou1_num1(m,:,:);
etrans_cou_num1_xz(m,:,:) = etrans_cou_num1(m,:,:);
clear etrans_cou1_num1
clear etrans_cou_num1

load('H:\global-PV-wind\ANS\etrans_cou1_num_1023_pro2_testt_UHV_sto_interUHV_8_2070.mat')% 实际有效发电量, TWh/year,etrans_cou1_num
load('H:\global-PV-wind\ANS\etrans_cou_num_1023_pro2_testt_UHV_sto_interUHV_8_2070.mat')% 全部发电量, TWh/year,etrans_cou_num
[m,n]=find(choo==8);
etrans_cou1_num1_xz(m,:,:) = etrans_cou1_num(m,:,:);
etrans_cou_num1_xz(m,:,:) = etrans_cou_num(m,:,:);
clear etrans_cou1_num
clear etrans_cou_num

etrans_cou1_num = etrans_cou1_num1_xz;
etrans_cou_num = etrans_cou_num1_xz;
save('H:\global-PV-wind\ANS\etrans_cou1_num_1023_pro2_8_2040.mat','etrans_cou1_num','-v7.3')% 实际有效发电量, TWh/year
save('H:\global-PV-wind\ANS\etrans_cou_num_1023_pro2_8_2040.mat','etrans_cou_num','-v7.3')% 全部发电量, TWh/year

