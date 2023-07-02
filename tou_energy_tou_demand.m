clear all;
close all;
clc;

load=xlsread('rate structure 1200 data.xlsx','Sheet1','C2:C2977');
solar=xlsread('rate structure 1200 data.xlsx','Sheet1','E2:E2977');

opt_load=load; %declaring optimal load
n=96; %declaring number of timestpes for each optimization
del_t=1/4; %time delta
d=length(load)/n; %number of days

%% tou energy charge array
OFF1=0.0755*ones(31,1);
MID1=0.0874*ones(16,1);
ON=0.1079*ones(25,1);
MID2=0.0874*ones(20,1);
OFF2=0.0755*ones(4,1);
alpha=[OFF1;MID1;ON;MID2;OFF2];

%% tou demand charge matrix
beta_OFF_val=1.53;
beta_MID_val=3.13;
beta_ON_val=7.06;

beta_OFF=zeros(n);
for i=1:31
    beta_OFF(i,i)=beta_OFF_val;
end
for i=93:n
    beta_OFF(i,i)=beta_OFF_val;
end

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
Emax=450; %SOC upper limit
Emin=100;%SOC lower limit
E_init=250;%initial state of charge
P_B_plus_max=100; %charging power limit
P_B_minus_max=100; %discharging power limit


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
        
        minimize(alpha'*P_G*del_t+max(beta_OFF*P_G)+max(beta_MID*P_G)+max(beta_ON*P_G))
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
beta_ON_month=diag(repmat(diag(beta_ON),d,1));
beta_MID_month=diag(repmat(diag(beta_MID),d,1));
beta_OFF_month=diag(repmat(diag(beta_OFF),d,1));

unopt_cost_1=alpha_month'*load*del_t+max(beta_OFF_month*load)+max(beta_MID_month*load)+max(beta_ON_month*load)
unopt_cost_2=alpha_month'*(load-solar)*del_t+max(beta_OFF_month*(load-solar))+max(beta_MID_month*(load-solar))+max(beta_ON_month*(load-solar))
opt_cost=alpha_month'*opt_load*del_t+max(beta_OFF_month*opt_load)+max(beta_MID_month*opt_load)+max(beta_ON_month*opt_load)
savings_1=unopt_cost_1-unopt_cost_2
savings_2=unopt_cost_2-opt_cost

load_with_solar=load-solar;
t_ax=(1:length(load))/n;
plot(t_ax,load,t_ax,opt_load)
xlabel('Time(days)')
ylabel('Power(kW)')
legend('Load w/o Optimization','Load with Optimization','Location','southeast')
title('Building 2 Load with and without Optimization:Rate Structure Type B')

load_mat=reshape(load,96,[]);
opt_load_mat=reshape(opt_load,96,[]);
load_avg=mean(load_mat,2);
opt_load_avg=mean(opt_load_mat,2);

hr=(0:n-1)/4;
figure
plot(hr,load_avg,hr,opt_load_avg)
legend('avg_load','avg_opt_load')
   

   
