
#######Ploter un long call

import numpy as np
import matplotlib.pyplot as plt

# -----------------------------
# Fonction payoff vectorisée
# -----------------------------
def option(option_position, option_type, s, k):
    if option_type == "call" or "Call":
        payoff = np.maximum(s - k, 0) if option_position == "long" else -np.maximum(s - k, 0)
    elif option_type == "put" or "Put":
        payoff = np.maximum(k - s, 0) if option_position == "long" else -np.maximum(k - s, 0)
    return payoff

# -----------------------------
# Paramètres de l'option
# -----------------------------
S_min = int(input("Please enter the minimum value of the underlying to display."))
S_max = int(input("Please enter the maximum value of the underlying to display."))
n_points = 10
K = int(input("Please enter the strike price value."))
option_type = str(input("What is the option type? A call or a put ?"))
option_position = str(input("What is the position on the option? Long or short?"))


# -----------------------------
# Générer les prix du sous-jacent
# -----------------------------
S_range = np.linspace(S_min, S_max, n_points)

# -----------------------------
# Calcul des payoffs
# -----------------------------
payoffs = option(option_position, option_type, S_range, K)

# -----------------------------
# Affichage graphique
# -----------------------------
plt.figure(figsize=(10,6))
plt.plot(S_range, payoffs, label=f"{option_position.capitalize()} {option_type.capitalize()}", color="blue")
plt.axvline(K, linestyle="--", color="red", label="Strike")
plt.axhline(0, color="black", linewidth=0.7)
plt.title(f"Payoff d'un {option_position} {option_type}")
plt.xlabel("Prix du sous-jacent S")
plt.ylabel("Payoff")
plt.legend()
plt.grid(True)
plt.show()
