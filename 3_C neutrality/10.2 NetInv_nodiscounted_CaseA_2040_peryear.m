%%
tic
clear;
load('H:\Global PV and wind\ANS\Inv_CCS2040_c_cou0_8_2s_6_2a.mat')  % Inv_CCS2070_c_cou0
load('H:\Global PV and wind\ANS\Inv_others2040_c_cou0_8_2s_6_2a.mat')  % Inv_others2070_c_cou0
Inv_others2070_c_cou0_dis = Inv_others2070_c_cou0;
Inv_others2070_c_cou0_dis2 = [zeros(size(Inv_others2070_c_cou0_dis,1),1) Inv_others2070_c_cou0_dis];
Inv_CCS2070_c_cou0_dis = Inv_CCS2070_c_cou0;
Inv_CCS2070_c_cou0_dis2 = [zeros(size(Inv_CCS2070_c_cou0_dis,1),1) Inv_CCS2070_c_cou0_dis];
    Inv_CCS2070_c_cou0_dis2a = Inv_CCS2070_c_cou0_dis2;
    Inv_others2070_c_cou0_dis2a = Inv_others2070_c_cou0_dis2;
save('H:\Global PV and wind\ANS\Inv_CCS2040_c_cou0_nodis2a_peryear.mat','Inv_CCS2070_c_cou0_dis2a','-v7.3')  % TWh/year
save('H:\Global PV and wind\ANS\Inv_others2040_c_cou0_nodis2a_peryear.mat','Inv_others2070_c_cou0_dis2a','-v7.3')  % TWh/year

%%
tic
clear;
nn = 1;
for i = 1:192
    for j = 1:192
        if i ~=j
            index_ij_s(nn,1) = i;
            index_ij_s(nn,2) = j;
            nn = nn+1;
        end
    end
end

% 没有PV和wind时各国各年的CCS和other renewable的discounted to 2019之后的投资
load('H:\Global PV and wind\ANS\Inv_CCS2040_c_cou0_nodis2a_peryear.mat')  % Inv_CCS2070_c_cou0
load('H:\Global PV and wind\ANS\Inv_others2040_c_cou0_nodis2a_peryear.mat')  % Inv_others2070_c_cou0_dis2a
Inv0_cou2040 = Inv_CCS2070_c_cou0_dis2a+Inv_others2070_c_cou0_dis2a;
clear Inv_CCS2070_c_cou0_dis2a
clear Inv_others2070_c_cou0_dis2a

% 没有国际合作时的CCS,other renewable和PV&wind的discounted to 2019之后的投资的变化
load('H:\Global PV and wind\ANS\Inv_CCS_dif_noconnect_cou_dis2a_8_2040.mat') % Inv_CCS_dif_noconnect_cou_dis2a
load('H:\Global PV and wind\ANS\Inv_others_dif_noconnect_cou_dis2a_8_2040.mat') % Inv_others_dif_noconnect_cou_dis2a
aa = Inv_CCS_dif_noconnect_cou_dis2a+Inv_others_dif_noconnect_cou_dis2a;
a1 = zeros(192,1);
for i = 1:10
    a1(:,i+1) = aa(:,i)*2-a1(:,i);
end
discount = 0.03;
rr2(1,1) = 0;
rr2(2,1) = 1/(1+discount)^(25-20);
rr2(3,1) = 1/(1+discount)^(30-20);
rr2(4,1) = 1/(1+discount)^(35-20);
rr2(5,1) = 1/(1+discount)^(40-20);
rr2(6,1) = 1/(1+discount)^(45-20);
rr2(7,1) = 1/(1+discount)^(50-20);
rr2(8,1) = 1/(1+discount)^(55-20);
rr2(9,1) = 1/(1+discount)^(60-20);
rr2(10,1) = 1/(1+discount)^(65-20);
rr2(11,1) = 1/(1+discount)^(70-20);
a2 = a1./repmat(rr2',[192,1]);
a2(find(isnan(a2)==1))=0;
% for i = 1:5
%     Inv1_cou2040(:,i) = sum(a2(:,i:i+1),2)/2;
% end
Inv1_cou2040 = a2;
clear Inv_CCS_dif_noconnect_cou_dis2a
clear Inv_others_dif_noconnect_cou_dis2a

%
load('H:\Global PV and wind\ANS\Inv_CCSall_dis2a_8_2040.mat')  %Inv_CCSall_dis2a
load('H:\Global PV and wind\ANS\Inv_othersall_dis2a_8_2040.mat')  %Inv_othersall_dis2a
% sum(Inv_CCS_dif_noconnect_cou_dis2a)-Inv_CCSall_dis2a(1,:)
% sum(Inv_others_dif_noconnect_cou_dis2a)-Inv_othersall_dis2a(1,:)
Inv_CCSall_dis2a(1,:) = [];
Inv_othersall_dis2a(1,:) = [];
Inv_CCSall_dis2a_reduct_cou = zeros(192,10);
Inv_othersall_dis2a_reduct_cou = zeros(192,10);
for i = 1:192
    [m,n]=find(index_ij_s(:,2)==i);
    if size(m,1)~=1
    Inv_CCSall_dis2a_reduct_cou(i,:) = sum(Inv_CCSall_dis2a(m,:));
    Inv_othersall_dis2a_reduct_cou(i,:) = sum(Inv_othersall_dis2a(m,:));
    else
    Inv_CCSall_dis2a_reduct_cou(i,:) = Inv_CCSall_dis2a(m,:);
    Inv_othersall_dis2a_reduct_cou(i,:) = Inv_othersall_dis2a(m,:);
    end
end
clear Inv_CCSall_dis2a
clear Inv_othersall_dis2a
aa = Inv_CCSall_dis2a_reduct_cou+Inv_othersall_dis2a_reduct_cou;
a1 = zeros(192,1);
for i = 1:10
    a1(:,i+1) = aa(:,i)*2-a1(:,i);
end
a2 = a1./repmat(rr2',[192,1]);
a2(find(isnan(a2)==1))=0;
% for i = 1:5
%     Inv2_cou2040(:,i) = sum(a2(:,i:i+1),2)/2;
% end
Inv2_cou2040 = a2;
clear Inv_CCSall_dis2a_reduct_cou
clear Inv_othersall_dis2a_reduct_cou

Inv_CCS_others_dis2a_cou = Inv0_cou2040+Inv1_cou2040+Inv2_cou2040;
clear Inv0_cou2040
clear Inv1_cou2040
clear Inv2_cou2040


%% 
load('H:\Global PV and wind\ANS\Inv_PWall_dis2a_2_2040CaseA_cou.mat')  % Inv_PWall_dis2a
load('H:\Global PV and wind\ANS\Inv_PWtoCCS_dis2a_2_2040CaseA_cou.mat')  % Inv_PWtoCCS_dis2a
aa = Inv_PWall_dis2a+Inv_PWtoCCS_dis2a;
a1 = zeros(192,1);
for i = 1:10
    a1(:,i+1) = aa(:,i)*2-a1(:,i);
end
discount = 0.03;
a2 = a1./repmat(rr2',[192,1]);
a2(find(isnan(a2)==1))=0;
% for i = 1:5
%     Inv_PW(:,i) = sum(a2(:,i:i+1),2)/2;
% end
Inv_PW = a2;
clear Inv_PWall_dis2a
clear Inv_PWtoCCS_dis2a

%%
rmb2us=6.8967; % RMB to USD2019
rmb2us2020=6.8996; % RMB to USD2019
rrr1 = rmb2us2020./rmb2us;
rrr = rrr1;
Inv_netcost_2040Cneutrality_cou = (Inv_CCS_others_dis2a_cou+Inv_PW)*rrr;
% sum(Inv_netcost_2040Cneutrality_cou)/10^6*10;
NetInv = Inv_netcost_2040Cneutrality_cou/10^6*5;
save('H:\Global PV and wind\ANS\NetInv_nodiscounted_3_2040peryear.mat','NetInv'); 
% trillion discounted net cost, trillion 2020$
% 2020:10:1070
