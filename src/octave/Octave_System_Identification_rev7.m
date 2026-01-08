%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%     SAFEHX SYSTEM IDENTIFICATION HEADER                         %%%%%%%%
%%%%%%%     BY: ALEXANDRE HENRIQUE COSTA ROSSI                          %%%%%%%%
%%%%%%%     REV7                                                        %%%%%%%%
%%%%%%%     Data: 21/05/2024                                            %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



clc;
close all;
clear all;

pkg load control  % Carrega o pacote de controle
pkg load io       % Carrega o pacote para ler arquivos csv
pkg load signal
pkg load image

f = 10; %frequencia de aquisicao
sampleTime = 1/f;


%------------------------     Leitura de dados      ----------------------------


M1 = dlmread("Resposta_Valvula_2.csv"); %dados do CSV com o experimento de identificação da resposta da válvula. Os dados são no formato resposta pressao; resposta vazao; ref_inv
M2 = dlmread("Resposta_Inversor_2.csv"); %dados do CSV com o experimento de identificação da resposta do inversor. Os dados são no formato resposta pressao; resposta vazao; ref_

n1 = length(M1)-2;
x1 = 0:sampleTime:sampleTime*(n1-2); %vetor de tempos da valvula
n2 = length(M2)-2;
x2 = 0:sampleTime:sampleTime*(n2-2); %vetor de tempos do inversor


%---------------------     Processamento de dados      -------------------------


corte_inicio1 = 100; %corte do inicio do teste que demora para começar a coleta de dados definitiva

corte_inicio2 = 1100;
corte_fim2 = 1000;

resp_p_valv = M1(2:(length(M1)-2),1)';
resp_p_valv =resp_p_valv(corte_inicio1:length(resp_p_valv));
y1_barra = mean(resp_p_valv);
resp_v_valv = M1(2:(length(M1)-2),2)';
resp_v_valv =resp_v_valv(corte_inicio1:length(resp_v_valv));
y2_barra = mean(resp_v_valv);

%Variáveis incrementais pressão valvula
%x_barra = 3.208;
x1_barra = mean(resp_p_valv(1:15));
resp_p_valv = resp_p_valv-x1_barra;

%Variáveis incrementais vazão valvula
%x2_barra = 0.763;
x2_barra = mean(resp_v_valv(1:15));

resp_v_valv = resp_v_valv - x2_barra;
ref_valv = M1(2:(length(M1)-2),3)';

ref_valv =ref_valv(corte_inicio1:length(ref_valv));

%Variáveis incrementais referência válvula
%u1_barra = 25;
u1_barra = mean(ref_valv(1:15));

ref_valv = ref_valv - u1_barra;

x1 = 0:sampleTime:sampleTime*(n1-corte_inicio1-1);

resp_p_inv = M2(2:(length(M2)-2),1)';
resp_p_inv =resp_p_inv(corte_inicio2:(length(resp_p_inv)-corte_fim2));

x3_barra = mean(resp_p_inv(1:50));
##x3_barra = 0;


%Variáveis incrementais Pressão x Inversor
resp_p_inv = resp_p_inv - x3_barra;

resp_v_inv = M2(2:(length(M2)-2),2)';
resp_v_inv =resp_v_inv(corte_inicio2:(length(resp_v_inv)-corte_fim2));
ref_inv = M2(2:(length(M2)-2),3)';

x4_barra = mean(resp_v_inv(1:50));
##x4_barra = 0;

%Variáveis incrementais Pressão x Inversor
resp_v_inv = resp_v_inv - x4_barra;

%Variáveis incrementais referência inversor
u2_barra = mean(ref_inv(1:50));
##u2_barra = 0;

ref_inv = ref_inv - u2_barra;

%ref_inv = ref_inv - u1_barra;
ref_inv =ref_inv(corte_inicio2:(length(ref_inv)-corte_fim2));


x2 = 0:sampleTime:sampleTime*(n2-corte_inicio2-corte_fim2-1);

%------------------------      Identificação       -----------------------------

%% P - Valv --------------------------------------------------------------------

##dsys_p_valv_ss = moen4(iddata([resp_p_valv',resp_v_valv'], ref_valv', sampleTime),25);
vetor_estaveis1 = zeros(1,1);
best_i1 = 0;
r_square_p_valv_best = 0;

##for i = 1:80
##  dsys_p_valv_ss = moen4(iddata([resp_p_valv'], ref_valv', sampleTime),i);
##  csys_p_valv_ss = d2c(dsys_p_valv_ss);
##  [sys_p_valv_num, sys_p_valv_den] = ss2tf(csys_p_valv_ss);
##  csys_p_valv_tf = tf(sys_p_valv_num, sys_p_valv_den);
##  [p,z] = pzmap(csys_p_valv_tf);
##  j = 1;
##
##  if all(real(p)<0) vetor_estaveis1(end+1) = i;
##      dsys_p_valv_ss = moen4(iddata([resp_p_valv'], ref_valv', sampleTime),i);
##      csys_p_valv_ss = d2c(dsys_p_valv_ss);
##      [sys_p_valv_num, sys_p_valv_den] = ss2tf(csys_p_valv_ss);
##      csys_p_valv_tf = tf(sys_p_valv_num, sys_p_valv_den);
##      [Ysim_dsys_p_valv_ss,Tsim_dsys_p_valv_ss] = lsim(dsys_p_valv_ss, [repmat(ref_valv(1),1000,1);ref_valv'], sampleTime);
##      Ysim_dsys_p_valv_ss_aux = zeros(1,length(x1));
##      Ysim_dsys_p_valv_ss_aux(1:length(x1)) = Ysim_dsys_p_valv_ss((length(Ysim_dsys_p_valv_ss)-length(x1)+1):length(Ysim_dsys_p_valv_ss));
##      Ysim_dsys_p_valv_ss_aux2(1:length(Ysim_dsys_p_valv_ss_aux)) = (Ysim_dsys_p_valv_ss_aux(1:length(Ysim_dsys_p_valv_ss_aux)));
##
##      r_square_p_valv = coefficient_of_determination(resp_p_valv ,Ysim_dsys_p_valv_ss_aux2);
##      if (r_square_p_valv > r_square_p_valv_best)
##        r_square_p_valv_best = r_square_p_valv;
##        best_i1 = i;
##      endif
##    endif
##endfor
##

best_i1 = 5;

dsys_p_valv_ss = moen4(iddata([resp_p_valv'], ref_valv', sampleTime),best_i1);
csys_p_valv_ss = d2c(dsys_p_valv_ss);
[sys_p_valv_num, sys_p_valv_den] = ss2tf(csys_p_valv_ss);
csys_p_valv_tf = tf(sys_p_valv_num, sys_p_valv_den);

##dat = iddata(resp_p_valv', ref_valv', sampleTime);
##[dsys_p_valv_ss,x0] = arx(dat, 'na', 1, 'nb', 1, 'tsam', sampleTime);
##
##csys_p_valv_ss = d2c(dsys_p_valv_ss);
##[sys_p_valv_num, sys_p_valv_den] = ss2tf(csys_p_valv_ss);
##csys_p_valv_tf = tf(sys_p_valv_num, sys_p_valv_den);
##csys_p_valv_tf;
##pole(csys_p_valv_tf);


##[dsys_p_valv_ss, K0] = mi.oe([ref_valv'], [resp_p_valv',resp_v_valv'], 4, 1, sampleTime);
##dsys_p_valv_ss = ss(dsys_p_valv_ss.A,dsys_p_valv_ss.B,dsys_p_valv_ss.C, dsys_p_valv_ss.D, sampleTime);
##csys_p_valv_ss = d2c(dsys_p_valv_ss);
##[sys_p_valv_num, sys_p_valv_den] = ss2tf(csys_p_valv_ss);
##csys_p_valv_tf = tf(sys_p_valv_num, sys_p_valv_den);


##----------------------------------------------------------------------------------

%% V - Valv --------------------------------------------------------------------

vetor_estaveis2 = zeros(1,1);
best_i2 = 0;
r_square_v_valv_best = 0;

##for i = 1:21
##  dsys_v_valv_ss = moen4(iddata([resp_v_valv'], ref_valv', sampleTime),i); %melhor resultado com n=20
##  csys_v_valv_ss = d2c(dsys_v_valv_ss);
##  [sys_v_valv_num, sys_v_valv_den] = ss2tf(csys_v_valv_ss);
##  csys_v_valv_tf = tf(sys_v_valv_num, sys_v_valv_den);
##  j = 1;
##  [p,z] = pzmap(csys_p_valv_tf);
##  if all(real(p)<0) vetor_estaveis2(end+1) = i;
##    [Ysim_dsys_v_valv_ss,Tsim_dsys_v_valv_ss] = lsim(dsys_v_valv_ss, [repmat(ref_valv(1),1000,1);ref_valv'], sampleTime);
##    Ysim_dsys_v_valv_ss_aux = zeros(1,length(x1));
##    Ysim_dsys_v_valv_ss_aux(1:length(x1)) = Ysim_dsys_v_valv_ss((length(Ysim_dsys_v_valv_ss)-length(x1)+1):length(Ysim_dsys_v_valv_ss));
##    Ysim_dsys_v_valv_ss_aux2(1:length(Ysim_dsys_v_valv_ss_aux)) = (Ysim_dsys_v_valv_ss_aux(1:length(Ysim_dsys_v_valv_ss_aux)));
##
##    r_square_v_valv = coefficient_of_determination(resp_v_valv ,Ysim_dsys_v_valv_ss_aux2);
##    if (r_square_v_valv > r_square_v_valv_best)
##        r_square_v_valv_best = r_square_v_valv;
##        best_i2 = i;
##      endif
##  endif
##
##
##endfor
##

%best_i2 = 22
best_i2 = 5

dsys_v_valv_ss = moen4(iddata([resp_v_valv'], ref_valv', sampleTime),best_i2); %melhor resultado com melhor i2
csys_v_valv_ss = d2c(dsys_v_valv_ss);
[sys_v_valv_num, sys_v_valv_den] = ss2tf(csys_v_valv_ss);
csys_v_valv_tf = tf(sys_v_valv_num, sys_v_valv_den);

%REDUCE MODEL ---------------------------------------------------------------------------------------------------

csys_v_valv_tf_reduced = tf(sys_v_valv_num(2:end), sys_v_valv_den);
dsys_v_valv_tf_reduced = c2d(csys_v_valv_tf_reduced,sampleTime);
csys_v_valv_ss_reduced = ss(csys_v_valv_tf_reduced);
dsys_v_valv_ss_reduced = c2d(csys_v_valv_ss_reduced, sampleTime);
%dsys_v_valv_ss_reduced = c2d(csys_v_valv_ss_reduced, sampleTime);

##  dat = iddata(resp_v_valv', ref_valv', sampleTime);
##  [dsys_v_valv_ss,x0] = arx(dat, 'na', 1, 'nb', 1, 'tsam', sampleTime);
##
##  csys_v_valv_ss = d2c(dsys_v_valv_ss);
##  [sys_v_valv_num, sys_p_valv_den] = ss2tf(csys_v_valv_ss);
##  csys_v_valv_tf = tf(sys_v_valv_num, sys_v_valv_den);
##  csys_v_valv_tf;
##  pole(csys_v_valv_tf);
##[dsys_v_valv_ss, K0] = mi.oe([ref_valv'], [resp_v_valv',resp_p_valv'], 2, 2, sampleTime);
##dsys_v_valv_ss = ss(dsys_v_valv_ss.A,dsys_v_valv_ss.B,dsys_v_valv_ss.C, dsys_v_valv_ss.D, sampleTime);
##csys_v_valv_ss = d2c(dsys_v_valv_ss);
##[sys_v_valv_num, sys_v_valv_den] = ss2tf(csys_v_valv_ss);
##csys_v_valv_tf = tf(sys_v_valv_num, sys_v_valv_den);
##csys_v_valv_tf;
####pole(csys_v_valv_tf);

%% P - Inv --------------------------------------------------------------------


##dsys_p_inv_ss = moen4(iddata([resp_p_inv'], [ref_inv'], sampleTime),1);
##csys_p_inv_ss = d2c(dsys_p_inv_ss);
##[sys_p_inv_num, sys_p_inv_den] = ss2tf(csys_p_inv_ss);
##csys_p_inv_tf = tf(sys_p_inv_num, sys_p_inv_den);

[dsys_p_inv_ss, K0] = mi.oe([ref_inv'], [resp_p_inv'], 2, 1, sampleTime);
dsys_p_inv_ss = ss(dsys_p_inv_ss.A,dsys_p_inv_ss.B,dsys_p_inv_ss.C, dsys_p_inv_ss.D, sampleTime);
csys_p_inv_ss = d2c(dsys_p_inv_ss);
[sys_p_inv_num, sys_p_inv_den] = ss2tf(csys_p_inv_ss);
csys_p_inv_tf = tf(sys_p_inv_num, sys_p_inv_den);


%% v - Inv --------------------------------------------------------------------


##%dsys_v_inv_ss = moen4(iddata([resp_v_inv',resp_p_inv'], [ref_inv'], sampleTime),1);
##dsys_v_inv_ss = moen4(iddata([resp_v_inv'], [ref_inv'], sampleTime),1);
##csys_v_inv_ss = d2c(dsys_p_inv_ss);
##[sys_v_inv_num, sys_v_inv_den] = ss2tf(csys_v_inv_ss);
##csys_v_inv_tf = tf(sys_v_inv_num, sys_v_inv_den);

[dsys_v_inv_ss, K0] = mi.oe([ref_inv'], [resp_v_inv'], 2, 1, sampleTime);
dsys_v_inv_ss = ss(dsys_v_inv_ss.A,dsys_v_inv_ss.B,dsys_v_inv_ss.C, dsys_v_inv_ss.D, sampleTime);
csys_v_inv_ss = d2c(dsys_v_inv_ss);
[sys_v_inv_num, sys_v_inv_den] = ss2tf(csys_v_inv_ss);
csys_v_inv_tf = tf(sys_v_inv_num, sys_v_inv_den);

save data.mat;

##load data.mat;


##----------- Figuras de comparação entre Modelo e Planta real -----------------
figure;

% Sub-plot Pressão X Válvula
subplot (2, 2, 1);
[Ysim_dsys_p_valv_ss,Tsim_dsys_p_valv_ss] = lsim(dsys_p_valv_ss, [repmat(ref_valv(1),1000,1);ref_valv'], sampleTime);
Ysim_dsys_p_valv_ss_aux = zeros(1,length(x1));
Ysim_dsys_p_valv_ss_aux(1:length(x1)) = Ysim_dsys_p_valv_ss((length(Ysim_dsys_p_valv_ss)-length(x1)+1):length(Ysim_dsys_p_valv_ss));
Ysim_dsys_p_valv_ss_aux2(1:length(Ysim_dsys_p_valv_ss_aux)) = (Ysim_dsys_p_valv_ss_aux(1:length(Ysim_dsys_p_valv_ss_aux)));

r_square_p_valv = coefficient_of_determination(resp_p_valv ,Ysim_dsys_p_valv_ss_aux2);

plot(x1,Ysim_dsys_p_valv_ss_aux2, x1,resp_p_valv);
title (['Plant vs Model: Pressure Response to Valve Perturbation S1, R² = ',num2str(r_square_p_valv,2)]);
legend('Model','Plant',"location","southeast");
ylim([-0.3 0.05]);
xlabel('Time [s]');
ylabel('Normalized pressure [bar(g)]');

% Sub-plot Vazão X Válvula
subplot (2, 2, 2);
[Ysim_dsys_v_valv_ss,Tsim_dsys_v_valv_ss] = lsim(dsys_v_valv_ss, [repmat(ref_valv(1),1000,1);ref_valv'], sampleTime);
Ysim_dsys_v_valv_ss_aux = zeros(1,length(x1));
Ysim_dsys_v_valv_ss_aux(1:length(x1)) = Ysim_dsys_v_valv_ss((length(Ysim_dsys_v_valv_ss)-length(x1)+1):length(Ysim_dsys_v_valv_ss));
Ysim_dsys_v_valv_ss_aux2(1:length(Ysim_dsys_v_valv_ss_aux)) = (Ysim_dsys_v_valv_ss_aux(1:length(Ysim_dsys_v_valv_ss_aux)));

r_square_v_valv = coefficient_of_determination(resp_v_valv ,Ysim_dsys_v_valv_ss_aux2);

%print('model_vs_plant_response.png', '-dpng');  % For compatibility if PDF causes issues

%--------------------------------------------- Sistema reduzido -------------------------------------------------------

%dsys_v_valv_ss_reduced
[Ysim_dsys_v_valv_ss_reduced,Tsim_dsys_v_valv_ss_reduced] = lsim(dsys_v_valv_ss_reduced, [repmat(ref_valv(1),1000,1);ref_valv'], sampleTime);
Ysim_dsys_v_valv_ss_reduced_aux = zeros(1,length(x1));
Ysim_dsys_v_valv_ss_reduced_aux(1:length(x1)) = Ysim_dsys_v_valv_ss_reduced((length(Ysim_dsys_v_valv_ss_reduced)-length(x1)+1):length(Ysim_dsys_v_valv_ss_reduced));
Ysim_dsys_v_valv_ss_reduced_aux2(1:length(Ysim_dsys_v_valv_ss_reduced_aux)) = (Ysim_dsys_v_valv_ss_reduced_aux(1:length(Ysim_dsys_v_valv_ss_reduced_aux)));

%plot(x1,Ysim_dsys_v_valv_ss_aux2,x1, Ysim_dsys_v_valv_ss_reduced_aux2, x1,resp_v_valv);
plot(x1,Ysim_dsys_v_valv_ss_aux2, x1,resp_v_valv);
%plot(x1, Ysim_dsys_v_valv_ss_reduced_aux2, x1,resp_v_valv);

title (['Plant vs Model: Flow Response to Valve Perturbation S1, R² = ',num2str(r_square_v_valv,2)]);
legend('Model','Plant',"location","southeast");
ylim([-0.2 0.4]);
xlabel('Time [s]');
ylabel('Normalized flow [m³/h]');

% Sub-plot Pressão X Inversor
subplot (2, 2, 3);
%mean(ref_inv(1:50)+4)
[Ysim_dsys_p_inv_ss,Tsim_dsys_p_inv_ss] = lsim(dsys_p_inv_ss, [repmat(0,1000,1);ref_inv'], sampleTime);
Ysim_dsys_p_inv_ss_aux = zeros(1,length(x2));
Ysim_dsys_p_inv_ss_aux(1:length(x2)) = Ysim_dsys_p_inv_ss((length(Ysim_dsys_p_inv_ss)-length(x2)+1):length(Ysim_dsys_p_inv_ss));
Ysim_dsys_p_inv_ss_aux2(1:length(Ysim_dsys_p_inv_ss_aux)) = (Ysim_dsys_p_inv_ss_aux(1:length(Ysim_dsys_p_inv_ss_aux)));
plot(x2,Ysim_dsys_p_inv_ss_aux2, x2,resp_p_inv);

r_square_p_inv = coefficient_of_determination(resp_p_inv,Ysim_dsys_p_inv_ss_aux2);

title (['Plant vs Model: Flow Response to VFD Perturbation S2, R² = ', num2str(r_square_p_inv,2)]);
legend('Model','Plant',"location","southeast");
ylim([-0.4 0.7]);
xlabel('Time [s]');
ylabel('Normalized flow [m³/h]');

% Sub-plot Vazão X Inversor
subplot (2, 2, 4);
%mean(ref_inv(1:200))
[Ysim_dsys_v_inv_ss,Tsim_dsys_v_inv_ss] = lsim(dsys_v_inv_ss, [repmat(0,1000,1);ref_inv'], sampleTime);
Ysim_dsys_v_inv_ss_aux = zeros(1,length(x2));
Ysim_dsys_v_inv_ss_aux(1:length(x2)) = Ysim_dsys_v_inv_ss((length(Ysim_dsys_v_inv_ss)-length(x2)+1):length(Ysim_dsys_v_inv_ss));
Ysim_dsys_v_inv_ss_aux2(1:length(Ysim_dsys_v_inv_ss_aux)) = (Ysim_dsys_v_inv_ss_aux(1:length(Ysim_dsys_v_inv_ss_aux)));
plot(x2,Ysim_dsys_v_inv_ss_aux2, x2,resp_v_inv);

r_square_v_inv = coefficient_of_determination(resp_v_inv,Ysim_dsys_v_inv_ss_aux2);

title (['Plant vs Model: Pressure Response to VFD Perturbation S2, R² = ',num2str(r_square_v_inv,2)]);
legend('Model','Plant',"location","southeast");
ylim([-0.06 0.1]);
xlabel('Time [s]');
ylabel('Normalized pressure [bar(g)]');

H = [csys_p_valv_tf(1,1), csys_p_inv_tf(1,1); csys_v_valv_tf(1,1), csys_v_inv_tf(1,1)];

##%G11 -> PID da válvula no TIA Portal: Kp = 0.17, Integral Time = 0.4, Derivate action time = 0.8, derivate delay coeff = 0.5,
##%Derivate action weighting 0.5, sample time = 0.5s
##G11 = pid(0.17,0.4,0.8);
##G11 = G11*-1;
%G22 -> PID do inversor no TIA Portal: Kp = 0.0034, Integral Time = 0.1, Derivate action time = 0.035, derivate delay coeff = 0.5,
%Derivate action weighting-- 0.5, sample time = 0.5s
K12 = pid(0.34,0.1,0.035);
##
##G = [G11; G22];
##
##L = series(G,H);
##
##sys_mf = feedback(L,[1,1],"-");
##
##%sys_mf_ss = tf2ss(sys_mf);
####u = repmat([3.5; 1.5],1,10);
##
##t_final = 4000;
##
##tsim = 0:1:t_final;
##u1 = zeros(length(tsim),1);
##u1(tsim>=2) = 1;
##
##u2 = zeros(length(tsim),1);
##u2(tsim>=2) = 1;
##
##u0 = zeros(length(tsim),1);
##u0(tsim>=2) = 1.5;
##
##uM = [u1, u2];
##tM = [tsim; tsim];
##%[test_y,test_t]=lsim(sys_mf,uM,tsim);
##
##figure;
##
##%lsim(sys_mf,uM,1)
##
##
##
##
##%Teste da resposta do controlador para as diferentes malhas
##figure;
##
##subplot (2, 2, 1);
##test_L11 = series(G11, H(1,1));
##test_sysmf11 = feedback(test_L11,1,"-");
##lsim(test_sysmf11,u1,tsim)
##
##subplot (2, 2, 2);
##test_L12 = series(G11, H(1,2));
##test_sysmf12 = feedback(test_L12,1,"-");
##lsim(test_sysmf12,u1,tsim)
##
##subplot (2, 2, 3);
##test_L21 = series(G22, H(2,1));
##test_sysmf21 = feedback(test_L21,1,"-");
##lsim(test_sysmf21,u1,tsim)
##
##subplot (2, 2, 4);
##test_L22 = series(G22, H(2,2));
##test_sysmf22 = feedback(test_L22,1,"-");
##lsim(test_sysmf22,u2,tsim)
##


##
##D12 = -H(1,2)/H(1,1);
##D21 = -csys_v_valv_tf_reduced/H(1,1);


% ------------------------------------------------- Definição do desacoplador ---------------------------------------------------------------

G11 = H(1,1);
G12 = H(1,2);
G21 = H(2,1);
G22 = H(2,2);

matrix_G = [G11, G12; G21, G22];

%Desacoplador ideal matriz D

D11 = minreal(G11*G22/(G11*G22-G12*G21));
D12 = minreal((-G11*G12)/(G11*G22-G12*G21));
D21 = minreal((-G21*G22)/(G11*G22-G12*G21));
D22 = minreal(G11*G22/(G11*G22-G12*G21));


s = tf('s');
[pG12,zG12] = pzmap(G12);
D12 = (s-pG12(1))*(s-pG12(2))/(s*s);
D12 = (s-pG12(1))*(s-pG12(2))/(s*(s+zG12(1)));

[pG21,zG21] = pzmap(G21);
[numG21,denG21,] = tfdata(G21);
numD21 = tf(denG21,1);
auxdenD21 = minreal(tf(numG21,1)/(s-zG21(1)));
[numG21aux,denG21aux] = tfdata(auxdenD21);
denD21aux1 = (-1*cell2mat(numG21aux));
denD21aux2 = tf(denD21aux1,1);

newZeros = [-1.4982 + 10.2312i, -1.4982 - 10.2312i, -1.7633 + 3.5035i, -1.7633 - 3.5035i];

% Calcula os coeficientes do polinômio cujas raízes são newZeros
num_coef = poly(newZeros);

% Cria a função de transferência usando os coeficientes do numerador
% (supondo que o denominador seja 1, ou adapte conforme sua aplicação)
denD21 = tf(num_coef, 1);

D21 = numD21/(denD21*(s-zG21(1)));
matrix_D = [D11, D12; D21, D22];

##% ----------------------------------------------- Para geração do arquivo -----------------------------------------------------------------
##
[numH11, denH11] = tfdata(H(1,1));
[numH12, denH12] = tfdata(H(1,2));
[numH21, denH21] = tfdata(H(2,1));
[numH22, denH22] = tfdata(H(2,2));

[numD11, denD11] = tfdata(D11, "vector");
[numD12, denD12] = tfdata(D12);
[numD21, denD21] = tfdata(D21);
[numD22, denD22] = tfdata(D22);


% Nome do arquivo de saída e caminho
filename = "C:\\Users\\Public\\saida_scilab.sci"; %output file directory and name
fid = fopen(filename, "w");

for i = 1:2
    for j = 1:2
        % Obter numerador e denominador
        [num, den] = tfdata(H(i,j), "vector");

        % Inverter a ordem dos coeficientes
        num_inv = fliplr(num);
        den_inv = fliplr(den);

        % Escrever no arquivo
        fprintf(fid, "numG%d%d = poly([%s],'s','c');\n", i, j, num2str(num_inv, " %g "));
        fprintf(fid, "denG%d%d = poly([%s],'s','c');\n", i, j, num2str(den_inv, " %g "));
    end
end

% Processar D11
[numD11, denD11] = tfdata(D11, "vector");
fprintf(fid, "numD11 = poly([%s],'s','c');\n", num2str(fliplr(numD11), " %g "));
fprintf(fid, "denD11 = poly([%s],'s','c');\n", num2str(fliplr(denD11), " %g "));

% Processar D12
[numD12, denD12] = tfdata(D12, "vector");
fprintf(fid, "numD12 = poly([%s],'s','c');\n", num2str(fliplr(numD12), " %g "));
fprintf(fid, "denD12 = poly([%s],'s','c');\n", num2str(fliplr(denD12), " %g "));

% Processar D21
[numD21, denD21] = tfdata(D21, "vector");
fprintf(fid, "numD21 = poly([%s],'s','c');\n", num2str(fliplr(numD21), " %g "));
fprintf(fid, "denD21 = poly([%s],'s','c');\n", num2str(fliplr(denD21), " %g "));

% Processar D22
[numD22, denD22] = tfdata(D22, "vector");
fprintf(fid, "numD22 = poly([%s],'s','c');\n", num2str(fliplr(numD22), " %g "));
fprintf(fid, "denD22 = poly([%s],'s','c');\n", num2str(fliplr(denD22), " %g "));


% Fechar o arquivo
fclose(fid);

fprintf("Arquivo '%s' gerado com sucesso!\n", filename);


% Set paper size and figure position (in inches)
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 12 6]);  % [left bottom width height]
set(gcf, 'PaperSize', [11 6]);          % match to PaperPosition

% Save figure
%print('model_vs_plant_response.png', '-dpng');  % For compatibility if PDF causes issues
print('model_vs_plant_response.pdf', '-dpdfwrite');
