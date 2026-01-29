#pour lancer le programme, il faut aller dans le terminal et écrire "
# cd "C:\Users\vmmou\OneDrive - Audencia\Documents\Coding\Python\Finance de marché\Pricing option"
#puis
# streamlit run pricer_option.py

import streamlit as st
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt

# --------- Black-Scholes pricing model ---------
def black_scholes_price_greeks(S, K, T, r, sigma, option_type='call', position='long'):
    if T <= 0 or sigma <= 0:
        price = max(S - K, 0) if option_type=='call' else max(K - S, 0)
        delta = gamma = vega = rho = theta = 0
        return price, delta, gamma, vega, rho, theta

    d1 = (np.log(S / K) + (r + 0.5 * sigma**2) * T) / (sigma * np.sqrt(T))
    d2 = d1 - sigma * np.sqrt(T)

    if option_type == 'call':
        price = S * norm.cdf(d1) - K * np.exp(-r * T) * norm.cdf(d2)
        delta = norm.cdf(d1)
        theta = -(S * norm.pdf(d1) * sigma)/(2*np.sqrt(T)) - r*K*np.exp(-r*T)*norm.cdf(d2)
        rho = K*T*np.exp(-r*T)*norm.cdf(d2)
    else:  # put
        price = K * np.exp(-r * T) * norm.cdf(-d2) - S * norm.cdf(-d1)
        delta = -norm.cdf(-d1)
        theta = -(S * norm.pdf(d1) * sigma)/(2*np.sqrt(T)) + r*K*np.exp(-r*T)*norm.cdf(-d2)
        rho = -K*T*np.exp(-r*T)*norm.cdf(-d2)

    gamma = norm.pdf(d1)/(S*sigma*np.sqrt(T))
    vega = S * norm.pdf(d1) * np.sqrt(T)

    # short position
    if position == 'short':
        price *= -1
        delta *= -1
        gamma *= -1
        vega *= -1
        rho *= -1
        theta *= -1

    return price, delta, gamma, vega, rho, theta

# --------- UI ---------
st.set_page_config(page_title="European option pricer", layout="centered")
st.title("EUROPEAN OPTIONS PRICER")

st.sidebar.header("Variables")
S = st.sidebar.number_input("Spot price", value=100.0)
K = st.sidebar.number_input("Strike price", value=100.0)
T = st.sidebar.number_input("Maturity (in years)", value=1.0)
r = st.sidebar.number_input("Risk-free rate", value=0.05)
sigma = st.sidebar.number_input("Volatility", value=0.2)

option_choice = st.sidebar.selectbox("Option type", ['Long Call', 'Short Call', 'Long Put', 'Short Put'])

# Déterminer type et position
option_type = 'call' if 'Call' in option_choice else 'put'
position = 'long' if 'Long' in option_choice else 'short'

# --------- Pricing ---------
price, delta, gamma, vega, rho, theta = black_scholes_price_greeks(S, K, T, r, sigma, option_type, position)

#st.subheader("Greques")
#st.write(f"Delta: {delta:.4f}")
#st.write(f"Gamma: {gamma:.4f}")
#st.write(f"Vega: {vega:.4f}")
#st.write(f"Theta: {theta:.4f}")
#st.write(f"Rho: {rho:.4f}")

# --------- Graphique ---------
st.subheader(f" {option_choice}") #displays the name of the option

st.subheader(f"Price = ${price:.4f}") #displays the price of the option


S_range = np.linspace(0.01, 2*S, 200)
option_values = [black_scholes_price_greeks(s, K, T, r, sigma, option_type, position)[0] for s in S_range]
payoff = np.maximum(S_range - K, 0) if option_type=='call' else np.maximum(K - S_range, 0) #put
payoff = payoff if position=='long' else -payoff #short


# --------- Plotting the option ---------
fig, ax = plt.subplots() 

ax.plot(S_range, option_values, label="Option price")
ax.plot(S_range, payoff, '--', label='Payoff at maturity')
ax.axvline(K, linestyle=':', color='gray', label='Strike price')
ax.set_xlabel('Spot price')
ax.set_ylabel('PnL')
ax.set_title(f'{option_choice}')
ax.legend()
ax.grid(True)
st.pyplot(fig) #displays the plot


# --------- Greeks ---------

def colored_value(name, value):
    color = "green" if value >= 0 else "red"
    return f"<span style='color:{color}'>{name}: {value:.4f}</span>"

st.markdown(colored_value("Delta", delta), unsafe_allow_html=True)
st.markdown(colored_value("Gamma", gamma), unsafe_allow_html=True)
st.markdown(colored_value("Vega", vega), unsafe_allow_html=True)
st.markdown(colored_value("Theta", theta), unsafe_allow_html=True)
st.markdown(colored_value("Rho", rho), unsafe_allow_html=True)


st.caption("European options pricer | Victor Mourgues")

