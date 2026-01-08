loadXcosLibs(); loadScicos();
clear;
exec("C:\Users\Public\get_last_iteration.sci", -1);
//this is used for stability because the simulation often crashes, so a .bat scrypt is used to run the scilab program in loop until finishing all iterations on the dataset. The "ger_last_iteration.sci is used to keep track of progress"

exec("C:\Users\public\saida_scilab.sci")
//this is the octave model output

//Identified model from octave
//numG11 = poly([-14.0799  4.18853  -1.09254  0.0540003  -0.00881303  -0.000216219],'s','c');
//denG11 = poly([152.834  1506.89  459.092  103.285  8.32633  1],'s','c');
//numG12 = poly([0.408094  -0.0265788],'s','c');
//denG12 = poly([4.01789  18.9276  1],'s','c');
//numG21 = poly([35.1111  -8.62702  2.73192  -0.106446  0.0197366  0.000246824],'s','c');
//denG21 = poly([186.79  1133.01  454.134  100.878  8.29138  1],'s','c');
//numG22 = poly([0.0285995  0.000252339],'s','c');
//denG22 = poly([1.75382  11.9588  1],'s','c');
//numD11 = poly([-4.64533e+07  -4.87195e+08  -1.31085e+09  -2.4554e+08  -7.45894e+07  -1.61312e+07  -7.02268e+06  -724942  -158512  -12110.5  -1160.99  -74.0202  -1.52062  -0.00838656],'s','c');
//denD11 = poly([-6.36811e+08  -1.01499e+10  -4.01053e+10  -3.98143e+09  -2.76599e+09  -4.4893e+08  -1.68107e+08  -1.16003e+07  -2.38718e+06  -80601.9  -2204.65  -219.043  84.0848  1],'s','c');
//numD12 = poly([2.89337e+08  3.62302e+09  1.15118e+10  1.65821e+09  5.32365e+08  1.08245e+08  5.49448e+07  3.41305e+06  961314  39188.2  1803.34  140.738  -40.3301  -0.883352],'s','c');
//denD12 = poly([-6.36811e+08  -1.01499e+10  -4.01053e+10  -3.98143e+09  -2.76599e+09  -4.4893e+08  -1.68107e+08  -1.16003e+07  -2.38718e+06  -80601.9  -2204.65  -219.043  84.0848  1],'s','c');
//numD21 = poly([-9.47823e+07  -1.35856e+09  -4.39066e+09  -6.26237e+08  -3.42908e+08  -6.49818e+07  -2.17929e+07  -2.29807e+06  -407914  -28330.4  -2403.78  -135.606  -2.1115  -0.00957364],'s','c');
//denD21 = poly([-6.36811e+08  -1.01499e+10  -4.01053e+10  -3.98143e+09  -2.76599e+09  -4.4893e+08  -1.68107e+08  -1.16003e+07  -2.38718e+06  -80601.9  -2204.65  -219.043  84.0848  1],'s','c');
//numD22 = poly([-4.64533e+07  -4.87195e+08  -1.31085e+09  -2.4554e+08  -7.45894e+07  -1.61312e+07  -7.02268e+06  -724942  -158512  -12110.5  -1160.99  -74.0202  -1.52062  -0.00838656],'s','c');
//denD22 = poly([-6.36811e+08  -1.01499e+10  -4.01053e+10  -3.98143e+09  -2.76599e+09  -4.4893e+08  -1.68107e+08  -1.16003e+07  -2.38718e+06  -80601.9  -2204.65  -219.043  84.0848  1],'s','c');

G11 = numG11/denG11;
G21 = numG12/denG12;
G12 = numG21/denG21;
G22 = numG22/denG22;

//D11 = numD11/denD11;
//D12 = numD12/denD12;
//D21 = numD21/denD21;
//D22 = numD22/denD22;

Matrix_G = [G11, G12; 
            G21, G22];

Matrix_D = [D11, D12; D21, D22];

//Desacoplador simplificado ou Feedfoward
//D12 = -G12/G11;
//D21 = -G21/G22;

//Valores incrementais
x1_barra = 3.208;
x2_barra = 0.7617;
u1_barra = 25;
u2_barra = 57.971;
y1_barra = 3.1194;
y2_barra = 0.9549;

x2_barra_aux = x2_barra+2;

//Ganho proporcional
Kp1 = 0.17;
Kp2 = 0.0034;
//Proportional action weighting
b1 = 1;
b2 = 1;
//Integration time
Ti1 = 0.4;
Ti2 = 0.1;
//Derivative action time
Td1 = 0.8;
Td2 = 0.035;
//derivative delay coeff
a1 = 0.5;
a2 = 0.5;
//derivative action wheighting
c1 = 0.5;
c2 = 0.5;

//Ruído aplicado ao sinal
noise1 = 0.0;
noise2 = 0.0;

// Fatores do desacoplador independente inicial
//Kp3 = 0;
Kp3 = -6; //fator K proporcional do PID_Compact
c3 = 1; //fator adicional do PID_Compact que reduz o efeito do derivador
Td3 = 0.0001; //Derivative action time do PID_Compact
b3 = 0.01; //Proportional action weighting do PID_Compact
a3 = 1; //derivative delay coefficient
Ti3 = 1000; //Integration time

d3 = 1; //fator para zerar efeito do integrador rodar com 0 ou 1


// Xcos file directory
xcos_file = "C:\Users\Public\Model2x2 - PID_Compact - special decoupler rev2.xcos" 
// Carrega o diagrama no Xcos
scs_m = xcosDiagramToScilab(xcos_file);
// Executa a simulação
xcos_simulate(scs_m, 4);


//Fatores que devem entrar na inspeção para o métod o Response surface methodology (RSM):
//Kp3 ;c3; Td3; b3; a3; Ti3

//_________________________________________ Código de varredura

// Geração dos vetores com variações de ±20%
// Obs.: para valores negativos, a variação mantém o sinal. Assim, por exemplo, para Kp3_0 = -6:
//      -20% => -6*0.8 = -4.8  e +20% => -6*1.2 = -7.2.

// Valores iniciais dos parâmetros
//Kp3_0 = -6;
//c3_0  = 1;
//Td3_0 = 0.0001;
//b3_0  = 0.01;
//a3_0  = 1;
//Ti3_0 = 1000;

// Valores iniciais dos parâmetros
//Kp3_0 = -6;
//c3_0  = 1;
//Td3_0 = 0.0001;
//b3_0  = 0.01;
//a3_0  = 1;
//Ti3_0 = 1000;

//This are the parameters generated in the python scrypt as optimal parameters for following on with the optimization
py_param = [-5.39020958e+00  1.24932007e+02  7.00000000e-05  7.50604000e-01  1.42152000e-01  1.28217400e+00];


Kp3_0 = py_param(1);
c3_0  = py_param(2);
Td3_0 = py_param(3);
b3_0  = py_param(4);
a3_0  = py_param(5);
Ti3_0 = py_param(6);

//// Tentar carregar os valores da última execução
//if fileinfo("saida_final.sci") <> [] then
//    exec("saida_final.sci", -1);
//    disp("Valores carregados de saida.sci");
//    Kp3 = Kp3_0   ; //fator K proporcional do PID_Compact
//    c3 = c3_0; //fator adicional do PID_Compact que reduz o efeito do derivador
//    Td3 = Td3_0; //Derivative action time do PID_Compact
//    b3 = b3_0; //Proportional action weighting do PID_Compact
//    a3 = a3_0; //derivative delay coefficient
//    Ti3 = Ti3_0; //Integration time
//else
//    disp("Usando valores padrão");
//end

fator_multi = 0.05;
fator_adi = 7e-5;

f_low = 1-fator_multi;
f_high = 1+fator_multi;

kp_range  = [Kp3_0,  Kp3_0*f_high,  Kp3_0*f_low];
c3_range  = [ c3_0,   c3_0*f_high,  c3_0*f_low];
Td3_range = [Td3_0,  Td3_0+fator_adi,  Td3_0-fator_adi];
b3_range  = [b3_0,   b3_0*f_high,  b3_0*f_low];
a3_range  = [a3_0,   a3_0*f_high,  a3_0*f_low];
Ti3_range = [Ti3_0,  Ti3_0*f_high,  Ti3_0*f_low];


// Nome do arquivo de resultados ___________________________________________________________________________________________________________________________________
filename = "resultados_RSM48.sci";

// Obtém a última iteração salva
last_seq = get_last_iteration(filename);
mprintf("Última iteração encontrada: %d\n", last_seq);

// Abre o arquivo para continuar escrevendo
fileID = mopen(filename, "at");

// Se for a primeira vez, escreve o cabeçalho
if last_seq == 0
    mfprintf(fileID, "// Resultados da Simulação\n");
    mfprintf(fileID, "// Sequencial, Kp3, c3, Td3, b3, a3, Ti3, peak_value\n");
end

// Inicializa o contador sequencial a partir da última iteração + 1
seq = last_seq + 1;

// Inicializa as variáveis para armazenar o pico mínimo e a melhor combinação
min_peak = %inf;  // valor infinito para garantir que qualquer pico será menor
best_combination = [];
min_pico = %inf;  // Inicializa com um valor alto
max_iter = 9;  // Teste com 10 simulações
iter = 0;
index = 1;

// Loop para testar todas as 3^6 combinações
for j = 1:length(c3_range)
  for k = 1:length(Td3_range)
    for l = 1:length(b3_range)
      for i = 1:length(kp_range)
        for m = 1:length(a3_range)
          for n = 1:length(Ti3_range)
            // Atualiza os parâmetros a serem usados na simulação
                        if index > last_seq
                            mprintf("Rodando a partir de: %d\n", index);

                            // Seleciona os valores da combinação atual
                            Kp3 = kp_range(i);
                            c3 = c3_range(j);
                            Td3 = Td3_range(k);
                            b3 = b3_range(l);
                            a3 = a3_range(m);
                            Ti3 = Ti3_range(n);
                            // Se os parâmetros precisam ser passados de alguma forma para o xcos,
                            // certifique-se de que eles estejam definidos no ambiente ou no diagrama.
                            
                            // Executa a simulação com o diagrama definido em xcos_file
                            scs_m = xcosDiagramToScilab(xcos_file);
                            xcos_simulate(scs_m, 4);
                            peak_value = sum((yout.values(599:1800,3)-yout.values(599:1800,4)).^2);
                            mfprintf(fileID, "[%d, %f, %f, %f, %f, %f, %f, %f]\n", seq, Kp3, c3, Td3, b3, a3, Ti3, peak_value);
                            seq = seq + 1;
                        end
                        index = index + 1;
          end
        end
      end
    end
  end
end

//fic = mopen("saida.sci", "wt");
//mfprintf(fic, "Kp3_0 = %f;\nc3_0 = %f;\nTd3_0 = %f;\nb3_0 = %f;\na3_0 = %f;\nTi3_0 = %f;\n", best_combination(1), best_combination(2), best_combination(3), best_combination(4), best_combination(5), best_combination(6));
//mclose(fic);
//disp("Novos valores salvos em saida.sci");

// Exibe o resultado: a combinação de parâmetros e o pico mínimo encontrado
//disp("Melhor combinação encontrada (Kp3, c3, Td3, b3, a3, Ti3):");
//disp(best_combination);
//disp("Pico máximo mínimo encontrado:");
//disp(min_peak);
//porc_overshoot = (max(yout.values(:,3))-max(yout.values(:,4)))/max(yout.values(:,4));
//disp(porc_overshoot);
//fic = mopen("saida.sci", "wt");
//mfprintf(fic, "Kp3_0 = %f;\nc3_0 = %f;\nTd3_0 = %f;\nb3_0 = %f;\na3_0 = %f;\nTi3_0 = %f;\n", best_combination(1), best_combination(2), best_combination(3), best_combination(4), best_combination(5), best_combination(6));
//mclose(fic);
//disp("Novos valores salvos em saida.sci");

//_________________________________________ Print do gráfico

//Kp3 = Kp3_0; //fator K proporcional do PID_Compact
//c3 = c3_0; //fator adicional do PID_Compact que reduz o efeito do derivador
//Td3 = Td3_0; //Derivative action time do PID_Compact
//b3 = b3_0; //Proportional action weighting do PID_Compact
//a3 = a3_0; //derivative delay coefficient
//Ti3 = Ti3_0; //Integration time

//// Carrega o diagrama no Xcos
//scs_m = xcosDiagramToScilab(xcos_file);
//// Executa a simulação
//xcos_simulate(scs_m, 4);

//// Criando a figura
//clf; // Limpa a figura

//// Primeiro subplot - Colunas 1 e 2
//subplot(3,1,1); // 2 linhas, 1 coluna, primeiro gráfico
//plot(yout.time, yout.values(:,1), 'r') // Coluna 1 em vermelho
//plot(yout.time, yout.values(:,2), 'b') // Coluna 2 em azul
//xlabel("Tempo [s]")
//ylabel("Pressão [Bar(g)]")
//title("Resposta da Pressão ao degrau na pressão e vazão")
//legend("Resposta da Pressão", "Setpoint")
//axis1 = gca();
////axis1.data_bounds = [0,0;length(yout.time),5];

//// Segundo subplot - Colunas 3 e 4
//subplot(3,1,2); // 2 linhas, 1 coluna, segundo gráfico
//plot(yout.time, yout.values(:,3), 'g') // Coluna 3 em verde
//plot(yout.time, yout.values(:,4), 'k') // Coluna 4 em preto


//xlabel("Tempo [s]")
//ylabel("Vazão [m³/h]")
//title("Resposta da Vazão ao degrau na pressão e vazão")
//legend("Resposta da Vazão", "Setpoint", "decoup")
//axis2 = gca();
// Ajusta a exibição
//gcf().color_map = jet(4); // Ajuste opcional de cores
//subplot(3,1,3);
//plot(yout.time, yout.values(:,5), 'r') // Coluna 4 em preto
//plot(yout.time, yout.values(:,6), 'b') // Coluna 4 em preto
//plot(yout.time, yout.values(:,7), 'k') // Coluna 4 em preto
//title("Resposta do decoupler")
//legend("Resposta do decoupler D21", "Resposta direta controle do inv", "resposta somada")

//max( yout.values(:,3))
