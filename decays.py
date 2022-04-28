
# from multiprocessing.sharedctypes import Value
import numpy as np
import random
import matplotlib.pyplot as plt
import sys

def tprobs(atom:str, thalf:dict, transitions:dict, dt:float=1) -> dict:
    """Takes possible transitions dict and outputs list of transition 
    probabilities for that the given isotope

    Args:
        atom (str): name of the atom
        thalf (dict): half
        dt (float): time set, default 1 (units of seconds)

    Returns:
        tuple: (tpossible, pvec) possible transitions (new atomic state) and 
                probability for each transition
    """
    if atom not in transitions.keys():
        return (['nodecay'], [1])     # Probability of nodecay is 1 if atom is stable
    
    tpossible = transitions[atom]    # get list of tuples for allowed 
                                        # daughters and associated 
                                        # bfrac of atom

    bfrac = [item[1] for item in tpossible]
    tposs = [item[0] for item in tpossible]
    thalf_total = thalf[atom]       # total half life of atom (seconds)
    lam_total = np.log(2)/thalf_total
    
    lam_partial = [lam_total * item for item in bfrac]
    pvec = [lam * dt for lam in lam_partial]
    if sum(pvec) <= 1:
        prob_nodecay = 1 - sum(pvec)        # probability of no decay
    else:
        return ValueError("Probabilities in pvec sum to a value greater than 1")
    pvec.append(prob_nodecay)
    tposs.append('nodecay')         # the "no decay" atom. If the code sees 
                                    #   this, then knows that the new atom is 
                                    #   the same as the old atom.

    # thalf_values = [thalf[item] for item in tpossible]
    # lams = [np.log(2)/t12 for t12 in thalf_values]  # gives decay constants for half-lives
    # for i in range(len(lams)):
    
    return (tposs, pvec)


def spec(KE:float, Q:float) -> float:

    global Emax

    if 0 <= KE <= Q:
        out = np.sqrt(KE**2 + 2*KE*0.511) * (Q-KE)**2 * (KE+0.511)
    else:
        out = 0

    return out


def emit(atom: str, Qvec: dict, Qvec_prime:dict):
    """_summary_

    Returns:
        _type_: _description_
    """
    
    global Emax 

    Q = Qvec[atom]
    Qprime = Qvec_prime[atom]
    # Qmax = max(Qvec.values())
    # Qmax = Emax

    KE_bins = np.linspace(0, Emax)

    E_spectrum = [spec(KE, Q) for KE in KE_bins]
    E_spectrum_prime = [spec(KE, Qprime) for KE in KE_bins]

    Esum = sum(E_spectrum)
    Esum_prime = sum(E_spectrum_prime)

    E_spectrum = [item/Esum for item in E_spectrum]

    KE_emission = np.random.choice(KE_bins, p=E_spectrum)

    weight = spec(KE_emission, Qprime) / spec(KE_emission, Q)

    return KE_emission, weight

# generate histogram arrays for the 
def reweight(weights:list, emissions:list) -> tuple:
    """Takes parallel lists of weights and original emission energies 
    and creates the reweighted histogram bins and bin edges edges

    Args:
        weights (list): list of weights from Qprime spectra
        emissions (list): kinetic energies of original emissions
        xbins (int): number of bins in histogram

    Returns:
        tuple: (w_count, binedges) weight sums for each bin, 
        bin edges (kinetic energy) list of length xbins+1
    """
    global Emax
    global xbins

    # print(weights)

    binedges = np.linspace(0, Emax, xbins + 1)
    w_counts = np.zeros(xbins)

    for (w, KE) in zip(weights, emissions):
        for i in range(xbins):
            if binedges[i] <= KE < binedges[i + 1]:
                # print(w_counts)
                # print(w)
                # print(type(w_counts))
                # print(type(w))
                w_counts[i] += w
    
    return w_counts, binedges

# state dict should include every possible atom in the sim

def gen(state:dict, transitions:dict, thalf:dict, Qvec:dict, Qvec_prime:dict, dt:float=1 ) -> dict:
    """Generates time evolved state(t + dt) from input state(t) after time step dt 
        
    Args:
        state (dict): dict of number count of each isotope 
        dt (float): time step, default 1 (units of seconds)

    Returns:
        tuple: (new_ensemble, emissions, all_weights) dict of time evolved state(t + dt) 
                and energies of emitted beta particles (units of MeV)
    """

    atoms = list(state.keys())      # which atoms could be in the initial ensemble
    counts = list(state.values())   # number quantity of each atom in the initial ensemble
    ensemble_out = {}
    batch_states = []               # take all strings of new states and count them later
    particle_emissions = []
    all_weights = []
    for (atom, count) in zip(atoms, counts):
        # now loop over all atoms

        # for atom, get list of transition probabilities and the corresponding new atom
        
        tposs, pvec = tprobs(atom, thalf, transitions, dt)
        
        ensemble_out[atom] = 0
        for i in range(count):
            # pick a new atomic state from possible transitions
            atom_new = np.random.choice(tposs, p=pvec)
            if atom_new == 'nodecay':
                atom_new = atom
            else: 
                emission, weight = emit(atom_new, Qvec, Qvec_prime)  # get the kinetic energy
                particle_emissions.append(emission) #    of emitted beta particle
                all_weights.append(weight)
            batch_states.append(atom_new)
        
    for atom in atoms:
        ensemble_out[atom] = batch_states.count(atom)
    # now we have a dict named ensemble_out containing each daughter state
    #   and the number quantity associated to it

    # print(all_weights)

    return ensemble_out, particle_emissions, all_weights


def sim(init: dict, transitions: dict, thalf: dict, Qvec:dict, Qvec_prime:dict, tmax:float, dt:float=1) -> dict:

    # Expects tmax and dt to be in same time units (assume seconds default)

    # Number of generations where no change occurs needed to trigger premature 
    #   stopping of simulation. 
    nochange_threshold = 20

    simtime = 0
    state = init                    # start by setting the state to initial state
    nochange_counter = 0
    emissions_total = []
    weights_total = []
    while simtime < tmax:
        state_old = state
        state, emissions, weights = gen(state, transitions, thalf, Qvec, Qvec_prime, dt)
        emissions_total += emissions
        weights_total = weights_total + weights
        simtime += dt
        if state == state_old:
            nochange_counter += 1
        else:
            nochange_counter = 0

        if nochange_counter == nochange_threshold:
            early_stop_time = simtime
            print('Simulation stopped early at', early_stop_time, 'generations.')
            print('Maximum generations:', tmax)
            simtime = tmax        

    w_counts, binedges = reweight(weights_total, emissions_total)

    return state, emissions_total, w_counts, binedges

# Namegaurd -> only run the contents of if statement when python script
#               is run as the main script (not imported as module)
if __name__ == "__main__":
    # Half-life reference list, units of seconds
    # if None, then it is stable
    thalf = {
        'Iso1': 10, 
        'Iso2': None, 
        'Iso3': None,
    }

    # branching fractions of each parent into daughters
    # not normalized
    # multiple different edge cases you have to deal with here:
    #   -> not all isotopes listed in each daughter dict, so handle unlisted 
    #       daughters as having 0 as entry (preprocess bfrac inside of sim?)
    #   -> check for empty daghter dicts, handle as no possible decays
    #   -> branching fractions don't add up to 1, so they need to be normalized
    #       to be used as ratios
    #   -> check for sum of daughter fraction equal to 0, interpret that as no 
    #       possible decays for that parent
    # 

    # dict of possible transitions
    # this, combined with thalf information, is all that is needed to determine
    #   branching ratios or transition probabilities
    # Format: 'parent1': [('daughter1', bfrac1), ('daughter2', bfrac2)]
    #           where bfrac1 + bfrac2 = 1
    transitions = {
        'Iso1': [
            ('Iso2', 0.1), 
            ('Iso3', 0.9)
            ],
    }

    # # initial state
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


    Emax = 1.5 # MeV

    xbins = 10
    tmax = 1000
    final_state, original_emits, w_counts, binedges = sim(init, transitions, thalf, Qvec, Qvec_prime, tmax)
    # print(final_state)

    plt.hist(original_emits, bins=xbins)
    plt.xlim(0, Emax)
    plt.xlabel('Kinetic Energy (MeV)')
    plt.ylabel('dN/dE')
    plt.title('Original Energy distribution of Iso1 emissions')
    plt.show()

    binwidth = Emax / xbins
    binlower = binedges[:-1]

    plt.bar(binlower, w_counts, width=binwidth, align='edge')
    plt.xlabel('Kinetic Energy (MeV)')
    plt.ylabel('dN/dE')
    plt.title('Reconstructed Energy distribution of Iso1 emissions')
    plt.show()

