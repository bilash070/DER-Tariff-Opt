clear all;
close all;
clc;

load=xlsread('rate structure cajalco data.xlsx','Sheet1','C2:C2881');
solar=xlsread('rate structure cajalco data.xlsx','Sheet1','D2:D2881');

opt_load=load; %declaring optimal load
n=96; %declaring number of timestpes for each optimization
del_t=1/4; %time delta
d=length(load)/n; %number of days

%% tou energy charge array
OFF1=0.07637*ones(31,1);
MID1=0.13837*ones(16,1);
ON=0.3397*ones(25,1);
MID2=0.13837*ones(20,1);
OFF2=0.07637*ones(4,1);
alpha=[OFF1;MID1;ON;MID2;OFF2];

beta=11.3;%demand charge

eta_plus=0.96; %charging efficiency
eta_minus=0.96; %discharging efficiency
Emax=450; %SOC upper limit
Emin=100;%SOC lower limit
E_init=250;%initial state of charge
P_B_plus_max=200; %charging power limit
P_B_minus_max=200; %discharging power limit


%% optimization for one month
for j=1:d
    l1=(j-1)*n+1;
    l2=j*n;
    P_L=load(l1:l2);
    P_S=solar(l1:l2);
    
    %% optimization for a day
    cvx_begin
        variables P_G(n) P_SL(n) P_B_plus(n) P_B_minus(n)
        expression E_B(n)
        E_B(1)=E_init;
        for t=2:n
                E_B(t)=E_B(t-1)+del_t*(P_B_plus(t-1)-P_B_minus(t-1));
        end
        
        minimize(alpha'*P_G*del_t+beta*max(P_G))
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
%                 P_SB(t)>=0;
        end
    cvx_end
    
    opt_load(l1:l2)=P_G;
    E_init=E_B(n);
end

alpha_month=repmat(alpha,d,1);
% beta_ON_month=diag(repmat(diag(beta_ON),d,1));
% beta_MID_month=diag(repmat(diag(beta_MID),d,1));
% beta_OFF_month=diag(repmat(diag(beta_OFF),d,1));

load_with_solar=load-solar;
unopt_cost_1=alpha_month'*load*del_t+beta*max(load)
unopt_cost_2=alpha_month'*(load-solar)*del_t+beta*max(load-solar)
opt_cost=alpha_month'*opt_load*del_t+beta*max(opt_load)
savings_1=unopt_cost_1-unopt_cost_2
savings_2=unopt_cost_2-opt_cost

t_ax=(1:length(load))/n;
plot(t_ax,load,t_ax,opt_load)
xlabel('Time(days)')
ylabel('Power(kW)')
legend('Load w/o Optimization','Load with Optimization','Location','southeast')
title('Building 1 Load with and without Optimization:Rate Structure Type A')
   
