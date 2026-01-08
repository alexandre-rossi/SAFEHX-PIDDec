function r2 = coefficient_of_determination(vec1, vec2)
  % Verifica se os vetores têm o mesmo tamanho
  if length(vec1) ~= length(vec2)
    error("Os vetores devem ter o mesmo tamanho.");
  end

  % Calcula a média dos valores observados
  mean_vec1 = mean(vec1);

  % Calcula a Soma dos Quadrados Totais (SST)
  sst = sum((vec1 - mean_vec1) .^ 2);

  % Calcula a Soma dos Quadrados dos Resíduos (SSR)
  ssr = sum((vec1 - vec2) .^ 2);

  % Calcula o Coeficiente de Determinação (R²)
  r2 = 1 - (ssr / sst);
end
