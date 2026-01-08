import numpy as np
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import RidgeCV
from scipy.optimize import minimize
import warnings
import os

# === CONFIGURAÇÕES ===
fator_expansao = 0.2  # Fator de expansão dos limites (20%)
mostrar_resultados = True

#file directory for the RSM results derived from the scilab program "RSM_results.sci"
caminho_arquivo =r"C:\Users\Public\RSM_results.sci"


# === FUNÇÃO DE CARREGAMENTO DE ARQUIVO ===
def carregar_dados(caminho_arquivo):
    matriz_resultados = []

    with open(caminho_arquivo, "r") as arquivo:
        for linha in arquivo:
            linha = linha.strip()
            if linha.startswith("//") or linha == "":
                continue
            linha = linha.replace("[", "").replace("]", "").strip()
            valores = list(map(float, linha.split(",")))[1:]
            matriz_resultados.append(valores)

    data = np.array(matriz_resultados)
    X = data[:, :6]
    y = data[:, 6]
    return X, y

# === FUNÇÃO PRINCIPAL DE AJUSTE E OTIMIZAÇÃO ===
def ajustar_modelo_e_otimizar(X, y):
    # Geração de polinômios de segundo grau
    poly = PolynomialFeatures(degree=2, include_bias=True)
    X_poly = poly.fit_transform(X)

    # Verificar condicionamento
    cond_number = np.linalg.cond(X_poly)
    print(f"\n🔎 Número de condicionamento de X_poly: {cond_number:.2e}")
    if cond_number > 1e12:
        warnings.warn("⚠️ A matriz X_poly está mal condicionada. Isso pode afetar a estabilidade da regressão.", RuntimeWarning)

    # Ajuste do modelo com RidgeCV
    alphas = np.logspace(-6, -2, 50)
    ridge_cv_model = RidgeCV(alphas=alphas, cv=None)  # store_cv_results=True incompatível com cv=None
    ridge_cv_model.fit(X_poly, y)

    # Função para predição
    def predict_peak(x):
        x_poly = poly.transform(x.reshape(1, -1))
        return ridge_cv_model.predict(x_poly)[0]

    # Otimização
    bounds = [(X[:, i].min(), X[:, i].max()) for i in range(6)]
    init_guess = X[y.argmin()]
    resultado = minimize(predict_peak, init_guess, bounds=bounds)

    if resultado.success:
        x_otimo = resultado.x
        y_predito = resultado.fun
        print("\n✅ Próximo ponto ótimo previsto:")
        print(f"Parâmetros:")
        print(f"py_param = {x_otimo};")
        print(f"Peak_value predito: {y_predito}")
        print(f"Melhor alpha: {ridge_cv_model.alpha_:.2e}")

        # Verificar fronteira
        nomes_variaveis = ["Kp3", "c3", "Td3", "b3", "a3", "Ti3"]
        na_fronteira = []
        for i, (val, (lb, ub)) in enumerate(zip(x_otimo, bounds)):
            if np.isclose(val, lb) or np.isclose(val, ub):
                na_fronteira.append((nomes_variaveis[i], val, lb, ub))

        if na_fronteira:
            print("\n⚠️ O ponto ótimo está na fronteira nas variáveis:")
            for nome, val, lb, ub in na_fronteira:
                print(f"  - {nome}: {val:.4f} | Limites: [{lb:.4f}, {ub:.4f}]")
        else:
            print("\n📌 O ponto ótimo está dentro da região experimental.")

        # Cálculo do erro (SSR)
        y_pred = ridge_cv_model.predict(X_poly)
        ssr = np.sum((y - y_pred) ** 2)
        print(f"\n📉 Soma dos quadrados dos resíduos (SSR): {ssr:.6f}")

        # Expandir limites para próxima iteração
        novos_bounds = []
        for i in range(6):
            centro = x_otimo[i]
            delta = (bounds[i][1] - bounds[i][0]) * fator_expansao / 2
            novos_bounds.append((centro - delta, centro + delta))

        return x_otimo, y_predito, ssr, novos_bounds
    else:
        print("❌ Otimização falhou:", resultado.message)
        return None, None, None, None

# === EXECUÇÃO ===
arquivos = [
    "C:/Users/Public/RSM_results.sci", #file directory for the RSM results derived from the scilab program "RSM_results.sci"
    #"C:/Users/alero/OneDrive - CLARK TECNOLOGIA QUÍMICA INDUSTRIA E COMERCIO LTDA/Documentos/Automacao Clark/Controle Avancado SAFEHX/Scilab/resultados_RSM48.sci",
    #"C:/Users/alero/OneDrive - CLARK TECNOLOGIA QUÍMICA INDUSTRIA E COMERCIO LTDA/Documentos/Automacao Clark/Controle Avancado SAFEHX/Scilab/resultados_RSM28.sci",
    #"C:/Users/alero/OneDrive - CLARK TECNOLOGIA QUÍMICA INDUSTRIA E COMERCIO LTDA/Documentos/Automacao Clark/Controle Avancado SAFEHX/Scilab/resultados_RSM29.sci",
    #"C:/Users/alero/OneDrive - CLARK TECNOLOGIA QUÍMICA INDUSTRIA E COMERCIO LTDA/Documentos/Automacao Clark/Controle Avancado SAFEHX/Scilab/resultados_RSM30.sci",
]

print("\n=========================")
print("  ANÁLISE DE ITERAÇÕES RSM")
print("=========================\n")

for i, caminho in enumerate(arquivos, 1):
    print(f"\n📁 Iteração {i} — Arquivo: {os.path.basename(caminho)}")
    X, y = carregar_dados(caminho)
    _, _, _, _ = ajustar_modelo_e_otimizar(X, y)
