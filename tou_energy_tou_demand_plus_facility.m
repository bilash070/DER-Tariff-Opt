clear all;
close all;
clc;

load=xlsread('rate structure ieua data.xlsx','Sheet1','C2:C2977');
solar=xlsread('rate structure ieua data.xlsx','Sheet1','D2:D2977');

opt_load=load; %declaring optimal load
n=96; %declaring number of timestpes for each optimization
del_t=1/4; %time delta
d=length(load)/n; %number of days

%% tou energy charge array
OFF1=0.05727*ones(31,1);
MID1=0.07566*ones(16,1);
ON=0.10258*ones(25,1);
MID2=0.07566*ones(20,1);
OFF2=0.05727*ones(4,1);
alpha=[OFF1;MID1;ON;MID2;OFF2];

%% tou demand charge matrix
beta_fac=19.02;%facilities related demand charge
beta_MID_val=4.17;
beta_ON_val=21.73;

% beta_OFF=zeros(n);
% for i=1:31
%     beta_OFF(i,i)=beta_OFF_val;
% end
% for i=73:n
%     beta_OFF(i,i)=beta_OFF_val;
% end

beta_MID=zeros(n);
for i=32:47
    beta_MID(i,i)=beta_MID_val;
end
for i=73:92
    beta_MID(i,i)=beta_MID_val;
end

beta_ON=zeros(n);
for i=48:72
    beta_ON(i,i)=beta_ON_val;
end

eta_plus=0.96; %charging efficiency
eta_minus=0.96; %discharging efficiency
Emax=288; %SOC upper limit
Emin=64;%SOC lower limit
E_init=160;%initial state of charge
P_B_plus_max=120; %charging power limit
P_B_minus_max=120; %discharging power limit


%% optimization for one month
for j=1:d
    l1=(j-1)*n+1;
    l2=j*n;
    P_L=load(l1:l2);
    P_S=solar(l1:l2);
    
    %% optimization for a day
    cvx_begin
        variables P_G(n) P_SL(n) P_SB(n) P_B_plus(n) P_B_minus(n)
        expression E_B(n)
        E_B(1)=E_init;
        for t=2:n
                E_B(t)=E_B(t-1)+del_t*(P_B_plus(t-1)-P_B_minus(t-1));
        end
        
        minimize(alpha'*P_G*del_t+max(beta_MID*P_G)+max(beta_ON*P_G)+beta_fac*max(P_G))
        subject to
        for t=1:n
                E_B(t)>=Emin;
                E_B(t)<=Emax;
                P_B_plus(t)>=0;
                P_B_plus(t)<=P_B_plus_max;
                P_B_minus(t)>=0;
                P_B_minus(t)<=P_B_minus_max;
                P_SL(t)+P_B_plus(t)/eta_plus==P_S(t);
                P_SL(t)+P_G(t)+P_B_minus(t)*eta_minus==P_L(t);
                P_SL(t)>=0;
                P_SB(t)>=0;
        end
    cvx_end
    
    opt_load(l1:l2)=P_G;
    E_init=E_B(n);
end

alpha_month=repmat(alpha,d,1);
beta_ON_month=diag(repmat(diag(beta_ON),d,1));
beta_MID_month=diag(repmat(diag(beta_MID),d,1));
% beta_OFF_month=diag(repmat(diag(beta_OFF),d,1));

unopt_cost_1=alpha_month'*load*del_t+max(beta_MID_month*load)+max(beta_ON_month*load)+beta_fac*max(load)
unopt_cost_2=alpha_month'*(load-solar)*del_t+max(beta_MID_month*(load-solar))+max(beta_ON_month*(load-solar))+beta_fac*max(load-solar)
opt_cost=alpha_month'*opt_load*del_t+max(beta_MID_month*opt_load)+max(beta_ON_month*opt_load)+beta_fac*max(opt_load)
savings_1=unopt_cost_1-unopt_cost_2
savings_2=unopt_cost_2-opt_cost

load_with_solar=load-solar;
t_ax=(1:length(load))/n;
plot(t_ax,load,t_ax,opt_load)
xlabel('Time(days)')
ylabel('Power(kW)')
legend('Load w/o Optimization','Load with Optimization','Location','southeast')
title('Building 3 Load with and without Optimization:Rate Structure Type C')
   
