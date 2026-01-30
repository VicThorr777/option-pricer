# to launch the program locally, go in the terminal and write "
# cd "C:\Users\vmmou\OneDrive - Audencia\Documents\Coding\Python\Finance de marché\Pricing option"
# then
# streamlit run pricer_option.py

import streamlit as st
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt

# --------- Black-Scholes pricing model ---------
def black_scholes_price_greeks(S, K, T, r, sigma, d, option_type='call', position='long'):
    if T <= 0 or sigma <= 0:
        price = max(S - K, 0) if option_type=='call' else max(K - S, 0)
        delta = gamma = vega = rho = theta = 0
        return price, delta, gamma, vega, rho, theta

    d1 = (np.log(S / K) + (r - d + 0.5 * sigma**2) * T) / (sigma * np.sqrt(T))
    d2 = d1 - sigma * np.sqrt(T)

    if option_type == 'call':
        price = S * np.exp(-d * T) * norm.cdf(d1) - K * np.exp(-r * T) * norm.cdf(d2)
        delta = norm.cdf(d1) * np.exp(-d * T)
        theta = (-(S * np.exp(-d*T) * norm.pdf(d1) * sigma)/(2*np.sqrt(T)) - r*K*np.exp(-r*T)*norm.cdf(d2) + d*S*np.exp(-d*T)*norm.cdf(d1) )/255 #multiply this number by 1/255 to get an indication of time decay on a daily basis
        rho = K * T * np.exp(-r * T) * norm.cdf(d2) / 100
    else:  # put
        price = K * np.exp(-r * T) * norm.cdf(-d2) - S * np.exp(-d * T) * norm.cdf(-d1)
        delta = -norm.cdf(-d1) * np.exp(-d * T)
        theta = -(S * np.exp(-d*T) * norm.pdf(d1) * sigma)/(2*np.sqrt(T)) + r*K*np.exp(-r*T)*norm.cdf(-d2) - d*S*np.exp(-d*T)*norm.cdf(-d1)
        rho = -K*T*np.exp(-r*T)*norm.cdf(-d2) / 100

    gamma = np.exp(-d*T) * norm.pdf(d1) / (S*sigma*np.sqrt(T))
    vega = S * np.exp(-d*T) * norm.pdf(d1) * np.sqrt(T) / 100

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

option_choice = st.sidebar.selectbox("Option type", ['Long Call', 'Short Call', 'Long Put', 'Short Put'])
S = st.sidebar.number_input("Spot price", value=100.0)
K = st.sidebar.number_input("Strike price", value=100.0)
T = st.sidebar.number_input("Maturity (in years)", value=1.0)
r = st.sidebar.number_input("Risk-free rate", value=0.05)
sigma = st.sidebar.number_input("Volatility", value=0.2)
d = st.sidebar.number_input("Dividend yield", value =0.04)

# Type and option position
option_type = 'call' if 'Call' in option_choice else 'put'
position = 'long' if 'Long' in option_choice else 'short'

# --------- Pricing ---------
price, delta, gamma, vega, rho, theta = black_scholes_price_greeks(S, K, T, r, sigma, d, option_type, position)


# --------- Graph ---------
st.subheader(f" {option_choice}") #displays the name of the option

st.subheader(f"Price = ${abs(price):.2f}") #displays the price of the option

S_range = np.linspace(0, S*2, 200)
option_values = [black_scholes_price_greeks(S, K, T, r, sigma, d, option_type, position)[0] for S in S_range]
payoff = np.maximum(S_range - K, 0) if option_type=='call' else np.maximum(K - S_range, 0) #put
payoff = payoff if position=='long' else -payoff #short


# --------- Plotting the option ---------

view = st.radio(   #buttons to chose from option and the greeks
    "View",
    ["Option", "Delta", "Gamma", "Vega", "Theta", "Rho"]
)

if view == "Option":
    fig, ax = plt.subplots() 
    
    ax.plot(S_range, option_values, color='orange', label='Option price')
    ax.plot(S_range, payoff, '--', label='Payoff at maturity')
    ax.axvline(S, linestyle=':', color='blue', label='Spot price')
    ax.axvline(K, linestyle=':', color='gray', label='Strike price')
    ax.set_xlabel('Spot price')
    ax.set_ylabel('PnL')
    ax.set_title(f'{option_choice}')
    ax.legend()
    ax.grid(True)
    st.pyplot(fig) #displays the plot

elif view == "Delta":
    
    delta_data = [black_scholes_price_greeks(i, K, T, r, sigma, d, option_type, position)[1] for i in S_range]

    fig, ax = plt.subplots()
    
    ax.plot(S_range, delta_data, color='red', label='Delta')
    ax.axvline(K, linestyle=':', color='gray', label='Strike price')
    ax.axvline(S, linestyle=':', color='blue', label='Spot price')
    ax.axhline(0, color='black', linewidth=0.5)
    
    ax.set_xlabel('Spot price')
    ax.set_ylabel('Delta')
    ax.set_title('Delta is ' + str(round(delta, 4)))
    ax.legend()
    ax.grid(True)
    
    st.pyplot(fig)
    
elif view == "Gamma":
        
    gamma_data = [black_scholes_price_greeks(S, K, T, r, sigma, d, option_type, position)[2] for S in S_range]

    fig, ax = plt.subplots()
        
    ax.plot(S_range, gamma_data, color='red', label='Gamma')
    ax.axvline(K, linestyle=':', color='gray', label='Strike price')
    ax.axvline(S, linestyle=':', color='blue', label='Spot price')
    ax.axhline(0, color='black', linewidth=0.5)
        
    ax.set_xlabel('Spot price')
    ax.set_ylabel('Gamma')
    ax.set_title('Gamma is ' + str(round(gamma, 4)))
    ax.legend()
    ax.grid(True)
        
    st.pyplot(fig)
    
elif view == "Vega":
        
    vega_data = [black_scholes_price_greeks(S, K, T, r, sigma, d, option_type, position)[3] for S in S_range]

    fig, ax = plt.subplots()
        
    ax.plot(S_range, vega_data, color='red', label='Vega')
    ax.axvline(K, linestyle=':', color='gray', label='Strike price')
    ax.axvline(S, linestyle=':', color='blue', label='Spot price')
    ax.axhline(0, color='black', linewidth=0.5)
        
    ax.set_xlabel('Spot price')
    ax.set_ylabel('Vega')
    ax.set_title('Vega is ' + str(round(vega, 4)))
    ax.legend()
    ax.grid(True)
        
    st.pyplot(fig)
    
elif view == "Rho":
        
    rho_data = [black_scholes_price_greeks(S, K, T, r, sigma, d, option_type, position)[4] for S in S_range]

    fig, ax = plt.subplots()
        
    ax.plot(S_range, rho_data, color='red', label='Rho')
    ax.axvline(K, linestyle=':', color='gray', label='Strike price')
    ax.axvline(S, linestyle=':', color='blue', label='Spot price')
    ax.axhline(0, color='black', linewidth=0.5)
        
    ax.set_xlabel('Spot price')
    ax.set_ylabel('Rho')
    ax.set_title('Rho is ' + str(round(rho, 4)))
    ax.legend()
    ax.grid(True)
        
    st.pyplot(fig)

elif view == "Theta":
        
    theta_data = [black_scholes_price_greeks(S, K, T, r, sigma, d, option_type, position)[5] for S in S_range]

    fig, ax = plt.subplots()
        
    ax.plot(S_range, theta_data, color='red', label='Theta')
    ax.axvline(K, linestyle=':', color='gray', label='Strike price')
    ax.axvline(S, linestyle=':', color='blue', label='Spot price')
    ax.axhline(0, color='black', linewidth=0.5)
        
    ax.set_xlabel('Spot price')
    ax.set_ylabel('Theta')
    ax.set_title('Theta is ' + str(round(theta, 4)))
    ax.legend()
    ax.grid(True)
        
    st.pyplot(fig)
        
    
# --------- Greeks ---------

def colored_value(name, value):
    color = "green" if value >= 0 else "red"
    return f"<span style='color:{color}'>{name}: {value:.2f}</span>"

st.markdown(colored_value("Delta", delta), unsafe_allow_html=True)
st.markdown(colored_value("Gamma", gamma), unsafe_allow_html=True)
st.markdown(colored_value("Vega", vega), unsafe_allow_html=True)
st.markdown(colored_value("Theta", theta), unsafe_allow_html=True)
st.markdown(colored_value("Rho", rho), unsafe_allow_html=True)

st.caption("European options pricer | Victor Mourgues - victor.mourgues@audencia.com")

