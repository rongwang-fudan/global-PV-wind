tic
clear;
load('H:\world\code\choo_type_8_2040xz.mat'); % choo
load('H:\world\code\etrans_cou1_num_1023_pro2_testt_8_2070.mat')% 实际有效发电量, TWh/year,etrans_cou1_num1
load('H:\world\code\etrans_cou_num_1023_pro2_testt_8_2070.mat')% 全部发电量, TWh/year,etrans_cou_num1
etrans_cou1_num1_xz = etrans_cou1_num1*0;
etrans_cou_num1_xz = etrans_cou_num1*0;
[m,n]=find(choo==1);
etrans_cou1_num1_xz(m,:,:) = etrans_cou1_num1(m,:,:);
etrans_cou_num1_xz(m,:,:) = etrans_cou_num1(m,:,:);

load('H:\Global PV and wind\ANS\etrans_cou1_num_1023_pro2_testt_UHV_8_2070.mat')% 实际有效发电量, TWh/year,etrans_cou1_num1
load('H:\Global PV and wind\ANS\etrans_cou_num_1023_pro2_testt_UHV_8_2070.mat')% 全部发电量, TWh/year,etrans_cou_num1
[m,n]=find(choo==2);
etrans_cou1_num1_xz(m,:,:) = etrans_cou1_num1(m,:,:);
etrans_cou_num1_xz(m,:,:) = etrans_cou_num1(m,:,:);

load('H:\Global PV and wind\ANS\etrans_cou1_num_1023_pro2_testt_UHV_sto_8_2070.mat')% 实际有效发电量, TWh/year,etrans_cou1_num1
load('H:\Global PV and wind\ANS\etrans_cou_num_1023_pro2_testt_UHV_sto_8_2070.mat')% 全部发电量, TWh/year,etrans_cou_num1
[m,n]=find(choo==3);
etrans_cou1_num1_xz(m,:,:) = etrans_cou1_num1(m,:,:);
etrans_cou_num1_xz(m,:,:) = etrans_cou_num1(m,:,:);

load('H:\Global PV and wind\ANS\etrans_cou1_num_1023_pro2_testt_UHV_sto_interUHV_8_2070.mat')% 实际有效发电量, TWh/year,etrans_cou1_num
load('H:\Global PV and wind\ANS\etrans_cou_num_1023_pro2_testt_UHV_sto_interUHV_8_2070.mat')% 全部发电量, TWh/year,etrans_cou_num
[m,n]=find(choo==4);
etrans_cou1_num1_xz(m,:,:) = etrans_cou1_num(m,:,:);
etrans_cou_num1_xz(m,:,:) = etrans_cou_num(m,:,:);

etrans_cou1_num = etrans_cou1_num1_xz;
etrans_cou_num = etrans_cou_num1_xz;
save('H:\Global PV and wind\ANS\etrans_cou1_num_1023_pro2_8_2040.mat','etrans_cou1_num','-v7.3')% 实际有效发电量, TWh/year
save('H:\Global PV and wind\ANS\etrans_cou_num_1023_pro2_8_2040.mat','etrans_cou_num','-v7.3')% 全部发电量, TWh/year

