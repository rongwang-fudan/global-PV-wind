tic
clear;
load('H:\global-PV-wind\ANS\Ph_PVwind2070_8_2s_6_2a.mat')  % million $,Ph_PVwind2070
load('H:\global-PV-wind\ANS\Inv_CCS2070_c_8_2s_6_2a.mat')  % million $,CCS投资
load('H:\global-PV-wind\ANS\Inv_others2070_8_2s_6_2a.mat') % million $，+附带fossil fuel成本，总投资
load('H:\global-PV-wind\ANS\Inv_others_e2070_8_2s_6_2a.mat')  % million $,抵消的fossil fuel的成本
Inv_others_c = Inv_others-Inv_others_e;
load('H:\global-PV-wind\ANS\Inv_PVwind2070_noe_8_2s_6_2a.mat')  % million $,净投资
Ph_PVwinda = cumsum(Ph_PVwind5);
Inv_cleana = Inv_PVwind2070_noe+Inv_CCS2070_c+Inv_others_c; % 净投资
clear Inv_others
clear Inv_others_e
clear Inv_PVwind2070_noe
clear Inv_CCS2070_c
clear Inv_others_c
clear Inv_PVwinda
clear Inv_CCS2070a
clear Inv_othersa
[B,IX] = sort([diff(sum(Inv_cleana,2)/10)./diff(Ph_PVwinda(:,10))]);
A = diff(Inv_cleana(:,10));
B = diff(Ph_PVwinda(:,10));
% A(A>0)=0;
[m,n]=find(A>0);
save('H:\global-PV-wind\ANS\m_2070_8_2s_6_2a.mat','m','-v7.3')  % TWh/yr


load('H:\global-PV-wind\ANS\Ph_PVwind2040_8_2s_6_2a.mat')  % million $,Ph_PVwind2070
load('H:\global-PV-wind\ANS\Ph_CCS2040_8_2s_6_2a.mat') % Ph_CCS5，TWh/yr
CCSa = diff(Ph_CCS5(:,10));
load('H:\global-PV-wind\ANS\Inv_CCS2040_c_8_2s_6_2a.mat')  
% trillion $,Inv_clean，CCS投资+CCS附带的Fossile fuel成本
load('H:\global-PV-wind\ANS\Inv_others2040_8_2s_6_2a.mat')  % trillion $,Inv_clean
load('H:\global-PV-wind\ANS\Inv_others_e2040_8_2s_6_2a.mat')  % trillion $,Inv_clean
Inv_others_c = Inv_others-Inv_others_e;
load('H:\global-PV-wind\ANS\Inv_PVwind2040_noe_8_2s_6_2a.mat')  % million $,Inv_PVwind2070
Ph_PVwinda = cumsum(Ph_PVwind5);
Inv_cleana = Inv_PVwind2070_noe+Inv_CCS2070_c+Inv_others_c; % 净投资
clear Inv_others
clear Inv_others_e
clear Inv_PVwind2070_noe
clear Inv_CCS2070_c
clear Inv_others_c
clear Inv_PVwinda
clear Inv_CCS2070a
clear Inv_othersa
[B,IX] = sort([diff(sum(Inv_cleana,2)/10)./diff(Ph_PVwinda(:,10))]);
save('H:\global-PV-wind\ANS\IX_netcost_2040_8_2s_6_2a.mat','IX','-v7.3')  % 
A = diff(Inv_cleana(:,10));
B = diff(Ph_PVwinda(:,10));
[m,n]=find(A>0);
save('H:\global-PV-wind\ANS\m_2040_8_2s_6_2a.mat','m','-v7.3')  % TWh/yr

