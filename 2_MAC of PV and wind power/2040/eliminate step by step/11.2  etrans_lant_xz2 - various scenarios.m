tic
clear;
load('H:\Global PV and wind\ANS\choo_type_8_2040xz_CaseB.mat'); % choo
load('H:\Global PV and wind\ANS\etrans_cou1_num_1023_pro2_testt_8_2070.mat')% 实际有效发电量, TWh/year,etrans_cou1_num1
load('H:\Global PV and wind\ANS\etrans_cou_num_1023_pro2_testt_8_2070.mat')% 全部发电量, TWh/year,etrans_cou_num1
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

etrans_cou1_num = etrans_cou1_num1_xz;
etrans_cou_num = etrans_cou_num1_xz;
save('H:\Global PV and wind\ANS\etrans_cou1_num_1023_pro2_8_2040_nosto.mat','etrans_cou1_num','-v7.3')% 实际有效发电量, TWh/year
save('H:\Global PV and wind\ANS\etrans_cou_num_1023_pro2_8_2040_nosto.mat','etrans_cou_num','-v7.3')% 全部发电量, TWh/year

%%
tic
clear;
load('H:\Global PV and wind\ANS\choo_type_8_2040xz.mat'); % choo
load('H:\Global PV and wind\ANS\etrans_cou1_num_1023_pro2_testt_8_2070.mat')% 实际有效发电量, TWh/year,etrans_cou1_num1
load('H:\Global PV and wind\ANS\etrans_cou_num_1023_pro2_testt_8_2070.mat')% 全部发电量, TWh/year,etrans_cou_num1
etrans_cou1_num1_xz = etrans_cou1_num1*0;
etrans_cou_num1_xz = etrans_cou_num1*0;
[m,n]=find(choo<=4);
etrans_cou1_num1_xz(m,:,:) = etrans_cou1_num1(m,:,:);
etrans_cou_num1_xz(m,:,:) = etrans_cou_num1(m,:,:);

etrans_cou1_num = etrans_cou1_num1_xz;
etrans_cou_num = etrans_cou_num1_xz;
save('H:\Global PV and wind\ANS\etrans_cou1_num_1023_pro2_8_2040_nosto_noUHV.mat','etrans_cou1_num','-v7.3')% 实际有效发电量, TWh/year
save('H:\Global PV and wind\ANS\etrans_cou_num_1023_pro2_8_2040_nosto_noUHV.mat','etrans_cou_num','-v7.3')% 全部发电量, TWh/year

