%% Step 1 抄TI
% 电机参数和控制器参数
clear all;clc
K = 162; % motor parameter 
L = 0.6e-3; % inductor in current loop
B=2126; % Bandwidth of Current Loop

% K = 4395; % motor parameter 
% L = 7.986e-3; % inductor in current loop
% B=2126; % Bandwidth of Current Loop

% TI parameter
delta = 4;
Kiss = (1.28)/(L*delta^2)
Kps = (1.28)/(L*delta*K)
% PDFF parameter 4-5-4.42
Kpf = 0.15; % 
Kpr = 0.15; %
Kis = Kiss.*Kps;

% Transfer function
s=tf('s');
% PI Transfer function
sys_open_PI = (B*K*Kps)*(s+Kiss)/(s*s)/(s+B);
sys_close_PI = feedback(sys_open_PI,1);
sys_load_PI = tf([-K -K*B 0],[1 B K*Kps*B K*Kps*Kiss*B]); % 负载扰动传函
% PDFF Transfer function
sys_open_PDFF = tf([B*K*Kpr B*K*Kis],[1 B B*K*(Kpf-Kpr) 0]);
sys_close_PDFF = tf([B*K*Kpr B*K*Kis],[1 B B*K*Kpf B*K*Kis]);  
sys_load_PDFF = tf([-K -K*B 0],[1 B K*Kpf*B K*Kis*B]); % 负载扰动传函

%% Step 2 RootLoci找Kpf（Kds）
figure(1)
delta_vector = [1.5 3 4 6 10 20];
num_delta_vector = length(delta_vector);

for i = 1:1:num_delta_vector
    Kis = (1.28)/(L.*delta_vector(i)^2).*(1.28)/(L.*delta_vector(i).*K);
    num_Kds = [K*B 0];
    den_Kds = [1 B 0 Kis*K*B];
    sys_rootloci_Kds = tf(num_Kds,den_Kds);
    subplot(2,3,i)
    rlocus(sys_rootloci_Kds);
    plot_num = i;
    title("Damping Factor" + " " +"\delta = "+ num2str(delta_vector(plot_num)));
    % grid on;
end
Kis = Kiss.*Kps; % 重新设置Kis MLGBD 这个循环害我总优化成PI
%% Step 3 图解法找Kpr（Kfs）
figure(2)
% syms Kpr_need;
% solve(B*K*Kpr_need^3-B*K*Kpf*Kpr_need^2+B*Kis*Kpr_need-Kis^2 == 0,Kpr_need)
Kpr_need_plot = 0:0.01:Kpf;
y_Value = B*K*Kpr_need_plot.^3-B*K*Kpf*Kpr_need_plot.^2+B*Kis*Kpr_need_plot-Kis.^2;

k = find(abs(y_Value-0)<=0.01e+06);% 求解Kpr3次曲线和y=0的交点,有bug，还没有解决，因为y_Value是个离散的量，我们要找的交点可能刚好从负跳到了正
y_Equal = B*K*Kpr_need_plot(k).^3-B*K*Kpf*Kpr_need_plot(k).^2+B*Kis*Kpr_need_plot(k)-Kis.^2;

plot(Kpr_need_plot,B*K*Kpr_need_plot.^3-B*K*Kpf*Kpr_need_plot.^2+B*Kis*Kpr_need_plot-Kis.^2,'LineStyle','-','LineWidth',1.5);hold on;
plot(Kpr_need_plot(k),y_Equal,'*')
plot(Kpr_need_plot,0*Kpr_need_plot,'LineStyle','-','LineWidth',1.5); hold off;

title('寻找Kpr之旅,Kpr-K');xlabel('K_{pr}'); ylabel('Value');
xlim([0 Kpf]);
legend('Value','相交点','0 curve');
%% Step 4 Draw the plot and SEE YOUR NICE WORK
% % 设置画图板
figure(3)

% 画关于Kpf的Root Loci
subplot(2,3,1);
num_Kds = [K*B 0];
den_Kds = [1 B 0 Kis*K*B];
sys_rootloci_Kds = tf(num_Kds,den_Kds);
rlocus(sys_rootloci_Kds)
grid on;

% 阶跃响应
subplot(2,3,2);
[y_PDFF,tOut_PDFF] = step(sys_close_PDFF,0.04);
plot(tOut_PDFF,y_PDFF,'LineStyle','-','LineWidth',1.5);hold on;
% 设定Kpr前的响应
[y_PDFF_Kpf,tOut_PDFF_Kpf] = step(tf([B*K*Kps B*K*Kis],[1 B B*K*Kpf B*K*Kis]),0.04);
plot(tOut_PDFF_Kpf,y_PDFF_Kpf,'LineStyle','-','LineWidth',1.5);hold on;
% 监测0.8K到1.2K的系统robust
[y_PDFF_KL,tOut_PDFF_KL] = step(tf([B*0.8*K*Kpr B*0.8*K*Kis],[1 B B*0.8*K*Kpf B*0.8*K*Kis]),0.04);
plot(tOut_PDFF_KL,y_PDFF_KL,'LineStyle','-','LineWidth',1.5);hold on;
[y_PDFF_KH,tOut_PDFF_KH] = step(tf([B*1.2*K*Kpr B*1.2*K*Kis],[1 B B*1.2*K*Kpf B*1.2*K*Kis]),0.04);
plot(tOut_PDFF_KH,y_PDFF_KH,'LineStyle','-','LineWidth',1.5);hold on;
% TI的响应
[y_PI,tOut_PI] = step(sys_close_PI,0.04);
plot(tOut_PI,y_PI,'LineStyle','-','LineWidth',1.5);hold on;
grid on;
title('Step Response');xlabel('Time'); ylabel('Amplitude');
legend('PDFF','PDFF_{NoKprOptimal}','PDFF_{0.8K}','PDFF_{1.2K}','PI');
hold off;

% 阶跃响应性能指标
S_PDFF = stepinfo(y_PDFF,tOut_PDFF,1)
S_PI = stepinfo(y_PI,tOut_PI,1)

% 开环传函的Bode
subplot(2,3,3);
bode(sys_open_PDFF,sys_open_PI,{1e-1,1e+6});
grid on; hold on;
h = findobj(gcf,'type','line');
set(h,'linewidth',1.5);
legend('PDFF_{Open}','PI_{Open}');
hold off;

% 计算Bode图性能指标
[Gm_PDFF,Pm_PDFF,Wcg_PDFF,Wcp_PDFF]=margin(sys_open_PDFF)
[Gm_PI,Pm_PI,Wcg_PI,Wcp_PI]=margin(sys_open_PI) 

% 闭环传函的Bode
subplot(2,3,4);
bode(sys_close_PDFF,sys_close_PI,{1e-1,1e+6});
grid on; hold on;
h = findobj(gcf,'type','line');
set(h,'linewidth',1.5);
legend('PDFF_{Close}','PI_{Close}');
hold off;

% 抗扰响应（扰动传函的阶跃响应
subplot(2,3,5);

t_Load = linspace(0,0.1);
Load = 0.5*square(t_Load,50)+0.5;% randn(size(t_Load))/10; % 生成峰值为1，占空比为50%的方波信号

L_PDFF = lsim(sys_load_PDFF,Load,t_Load); 
L_PI = lsim(sys_load_PI,Load,t_Load); 
plot(t_Load,L_PDFF,'LineStyle','-','LineWidth',1.5), hold on;
plot(t_Load,L_PI,'LineStyle','-','LineWidth',1.5), hold on;
grid on, axis([0 0.1 -1.5 0.25]);
title('Disturbanec Response')
xlabel('Time'); ylabel('Amplitude');
legend('L_{PDFF}','L_{PI}');
hold off;

% 扰动Bode
subplot(2,3,6);
bode(sys_load_PDFF,sys_load_PI)
grid on; hold on;
h = findobj(gcf,'type','line');
set(h,'linewidth',1.5);
legend('PDFF_{Disturbance}','PI_{Disturbance}');
hold off;

% pzmap
figure(4)
subplot(2,2,1)
pzmap(sys_open_PDFF);legend('pzmap_{PDFF_{Open}}');
subplot(2,2,2)
pzmap(sys_close_PDFF); legend('pzmap_{PDFF_{Close}}');
subplot(2,2,3)
pzmap(tf([B*0.8*K*Kpr B*0.8*K*Kis],[1 B B*0.8*K*(Kpf-Kpr) 0])); legend('pzmap_{PDFF_{Oepn}0.8K}');
subplot(2,2,4)
pzmap(tf([B*1.2*K*Kpr B*1.2*K*Kis],[1 B B*1.2*K*(Kpf-Kpr) 0])); legend('pzmap_{PDFF_{Oepn}1.2K}');


% t = 0:0.001:6;
% sys_close_PDFF = sys_close_PDFF/s;
% y_PDFF = step(sys_close_PDFF,t);
% plot(t,y_PDFF,'LineStyle','-','LineWidth',1.5);hold on;
% plot(t,t,'LineStyle','-','LineWidth',1.5);hold off;
% legend('Proposal','Input');

figure(5)
t_ramp = linspace(0,1,1000);
t_constant = linspace(0.2,1,800);
t_constant_compensate = linspace(0,0,200);
t_constant = [t_constant_compensate t_constant]
t_NB = linspace(1,1,800);
t_NB = [t_constant_compensate,t_NB];

ramp = 2500 .* t_ramp - 2500 .* t_constant + 1000 + 500 .* t_NB % 设置的输入信号
% plot(t_ramp,ramp)

ramp_PDFF = lsim(sys_close_PDFF,ramp,t_ramp);
ramp_PI = lsim(sys_close_PI,ramp,t_ramp);
plot(t_ramp,ramp_PDFF,'LineStyle','-','LineWidth',1.5), hold on;
plot(t_ramp,ramp_PI,'LineStyle','-','LineWidth',1.5), hold on;
plot(t_ramp,ramp,'LineStyle','-','LineWidth',1.5), hold off;
grid on
% title('Ramp Response');


set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 12, 5]); % 设置图窗大小
ax = gca;
exportgraphics(ax,"C:\Users\Wu\Desktop\Fullpaper_ECCE-Asia\Simulation Result/Ramp Response of Proposal Under delta = " + num2str(delta) + ".jpg","Resolution",600);
hold off;







