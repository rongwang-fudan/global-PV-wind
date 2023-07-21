% 实际各国learning rate，当为nan和inf时和全球统一
tic
clear;
load('H:\Global PV and wind\ANS\index_mineral_pv_county0811_2_all_CNfirst_5%_inilow_pro2_8_2070.mat')
load('H:\Global PV and wind\ANS\index_mineral_ons_county0811_2_all_CNfirst_5%_inilow_pro2_8_2070.mat')
load('H:\Global PV and wind\ANS\index_mineral_off_county0811_2_all_CNfirst_5%_inilow_pro2_8_2070.mat')

index_mineral_pv_time2 = index_mineral_pv;
index_mineral_ons_time2 = index_mineral_ons;
index_mineral_off_time2 = index_mineral_off;

save('H:\Global PV and wind\ANS\index_mineral_pv_time2_county0811_all_cou_5%_LRglobal_inilow_pro2_8_2070.mat','index_mineral_pv_time2') % 按照成本排序&考虑建厂时间后保留的PV电厂原始序号
save('H:\Global PV and wind\ANS\index_mineral_ons_time2_county0811_all_cou_5%_LRglobal_inilow_pro2_8_2070.mat','index_mineral_ons_time2') % 按照成本排序后&考虑建厂时间后保留的onshorewind电厂原始序号
save('H:\Global PV and wind\ANS\index_mineral_off_time2_county0811_all_cou_5%_LRglobal_inilow_pro2_8_2070.mat','index_mineral_off_time2') % 按照成本排序后&考虑建厂时间后保留的offshorewind电厂原始序号
