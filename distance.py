#%%
import decays
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    # initial state
    init = {
        'Iso1': 10000,
        'Iso2': 0, 
        'Iso3': 0, 
    }

    Qvec = {
        'Iso2': 0.3, 
        'Iso3': 1.2, 
    }

    Qvec_prime = {
        'Iso2': 0.3,
        'Iso3': 1.2,
    }

    transitions = {
        'Iso1': [
            ('Iso2', 0.1), 
            ('Iso3', 0.9)
            ],
    }

    thalf = {
        'Iso1': 10, 
        'Iso2': None, 
        'Iso3': None,
    }

    tmax = 1000
    Emax = 1.5      # MeV
    xbins = 10
    final_state, original_emits, w_counts, binedges = decays.sim(init, transitions, thalf, Qvec, Qvec_prime, tmax)

# %%
