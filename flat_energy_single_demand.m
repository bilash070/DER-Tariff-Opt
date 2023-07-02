clear all;
close all;
clc;

load=xlsread('rate structure admin data.xlsx','Sheet1','D2:D2881');
solar=xlsread('rate structure admin data.xlsx','Sheet1','E2:E2881');

opt_load=load; %declaring optimal load
n=96; %declaring number of timestpes for each optimization
del_t=1/4; %time delta
alpha=0.139; %flat energy charge
beta=10.58; %single demand charge
eta_plus=0.96; %charging efficiency
eta_minus=0.96; %discharging efficiency
Emax=450; %SOC upper limit
Emin=100;%SOC lower limit
E_init=250;%initial state of charge
P_B_plus_max=100; %charging power limit
P_B_minus_max=100; %discharging power limit
d=length(load)/n; %number of days



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
        
        minimize(alpha*del_t*sum(P_G)+beta*max(P_G))
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

unopt_cost_1=alpha*del_t*sum(load)+beta*max(load)
unopt_cost_2=alpha*del_t*sum(load-solar)+beta*max(load-solar)
opt_cost=alpha*del_t*sum(opt_load)+beta*max(opt_load)
savings_1=unopt_cost_1-unopt_cost_2
savings_1=unopt_cost_2-opt_cost

load_with_solar=load-solar;
t_ax=(1:length(load))/n;
plot(t_ax,load,t_ax,opt_load)
xlabel('Time(days)')
ylabel('Power(kW)')
legend('Load w/o Optimization','Load with Optimization','Location','southeast')
title('Building 4 Load with and without Optimization: Rate Structure Type E')

load_mat=reshape(load,96,[]);
opt_load_mat=reshape(opt_load,96,[]);
load_avg=mean(load_mat,2);
opt_load_avg=mean(opt_load_mat,2);

hr=(0:n-1)/4;
figure
plot(hr,load_avg,hr,opt_load_avg)
legend('avg_load','avg_opt_load')
   
