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

rmb2us=6.8967; % RMB to USD2019
rmb2us2020=6.8996; % RMB to USD2019
rrr1 = rmb2us2020./rmb2us;
rrr = rrr1;

load('H:\Global PV and wind\ANS\Inv_CCSall_dis2a_7_2040.mat')  %Inv_CCSall_dis2a
load('H:\Global PV and wind\ANS\Inv_othersall_dis2a_7_2040.mat')  %Inv_othersall_dis2a
load('H:\Global PV and wind\ANS\Inv_PWall_dis2a_6_2040xz.mat')  %Inv_PWall_dis2a
load('H:\Global PV and wind\ANS\Inv_PWtoCCS_dis2a_6_2040xz.mat')  %Inv_PWall_dis2a
load('H:\Global PV and wind\ANS\Inv_CCS0_dis2a_7_2040.mat')  %Inv_CCS0_dis2a
load('H:\Global PV and wind\ANS\Inv_others_c0_dis2a_7_2040.mat')  %Inv_others_c0_dis2a
aa = Inv_CCSall_dis2a+Inv_othersall_dis2a+Inv_PWall_dis2a+Inv_PWtoCCS_dis2a;
a1 = zeros(size(aa,1),1);
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
a2 = a1./repmat(rr2',[size(aa,1),1]);
a2(find(isnan(a2)==1))=0;
% for i = 1:5
%     Inv_dif2040(:,i) = sum(a2(:,i:i+1),2)/2;
% end
Inv_dif2040 = a2;

% 没有PV和wind时各国各年的CCS和other renewable的discounted to 2019之后的投资
load('H:\Global PV and wind\ANS\Inv_CCS2040_c_cou0_dis2a.mat')  % Inv_CCS2070_c_cou0
load('H:\Global PV and wind\ANS\Inv_others2040_c_cou0_dis2a.mat')  % Inv_others2070_c_cou0_dis2a
aa = Inv_CCS2070_c_cou0_dis2a+Inv_others2070_c_cou0_dis2a;
a1 = zeros(192,1);
for i = 1:10
    a1(:,i+1) = aa(:,i)*2-a1(:,i);
end
a2 = a1./repmat(rr2',[192,1]);
a2(find(isnan(a2)==1))=0;
% for i = 1:5
%     Inv0_cou2040(:,i) = sum(a2(:,i:i+1),2)/2;
% end
Inv0_cou2040 = a2;


% 没有国际合作时的CCS,other renewable和PV&wind的discounted to 2019之后的投资的变化
load('H:\Global PV and wind\ANS\Inv_CCS_dif_noconnect_cou_dis2a_8_2040.mat') % Inv_CCS_dif_noconnect_cou_dis2a
load('H:\Global PV and wind\ANS\Inv_others_dif_noconnect_cou_dis2a_8_2040.mat') % Inv_others_dif_noconnect_cou_dis2a
load('H:\Global PV and wind\ANS\inv_cou2040_dis2a_noconnect.mat')  % inv_cou2040_dis2a
aa = Inv_CCS_dif_noconnect_cou_dis2a+Inv_others_dif_noconnect_cou_dis2a+inv_cou2040_dis2a;
a1 = zeros(192,1);
for i = 1:10
    a1(:,i+1) = aa(:,i)*2-a1(:,i);
end
a2 = a1./repmat(rr2',[192,1]);
a2(find(isnan(a2)==1))=0;
% for i = 1:5
%     Inv_dif2040_noconnect_cou(:,i) = sum(a2(:,i:i+1),2)/2;
% end
Inv_dif2040_noconnect_cou = a2;

Inv0_coudif = Inv0_cou2040*rrr/10^6*5;
Inv_dif_noconnect_cou = Inv_dif2040_noconnect_cou*rrr/10^6*5;
A = Inv_dif2040(2:end,:)*rrr;
A(A>0)=0;

Cost_reduct_cou = zeros(192,11);
for i = 1:192
    [m,n]=find(index_ij_s(:,2)==i);
    Cost_reduct_cou(i,:) = sum(A(m,:))/10^6*5;
end
NetInv = Inv0_coudif+Inv_dif_noconnect_cou+Cost_reduct_cou;
sum(NetInv)
save('H:\Global PV and wind\ANS\NetInv_nodiscounted_6_2040peryear.mat','NetInv'); % trillion discounted net cost, trillion 2020$
