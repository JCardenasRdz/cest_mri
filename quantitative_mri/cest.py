import numpy as np
from scipy.optimize import least_squares
from numpy.random import uniform

def sim_offsets(initial,last,steps):
    offsets = np.linspace(initial,last,steps)
    offsets = np.squeeze(offsets)

def Lorentzian(pars,xdata_):
    """
    Lorentzian function for a single pool
    """
    fwhm = pars[1]**2/4
    L= (pars[0]*fwhm**2)/ ( (xdata_-pars[2])**2 + (fwhm)**2 )
    return L

def lorentzian_sim(xdata, Amp, Width, Center):
    """
    Estimates sum of Lorentzian functions where:
    Amp   = 1 X N lorentzian of amplitudes
    Width = 1 X N lorentzian of widths
    Center = 1 X N lorentzian of centers
    xdata = 1XN indepedent variable
    """
    # Convert to arrays (just in case):
    Amp = np.array(Amp)
    Width = np.array(Width)
    Center = np.array(Center)
    # Estimate number of pools
    Num_variables = Amp.__len__() + Width.__len__() + Center.__len__()

    # make sure it is divisible by 3.0
    assert (Num_variables % 3 == 0),"Please provide 3 variables per pool"

    # calculate final output
    num_pools = int(Num_variables/3)

    # Preallocate output
    Lsum = np.zeros( (xdata.shape[0]))

    for idx in range(num_pools):
        # assign each variable
        amp = Amp[idx]
        width = Width[idx]
        center = Center[idx]
        # estimate signal and sum
        Lsum += Lorentzian( [amp,width,center], xdata)

    return Lsum

def lorentzian_fit(x_data,experimental_data,initial_guess_offsets,repetitions = 10):
    """
    lorentzian_fit(x_data,experimental_data,initial_guess_offsets)

    returns:
    x_pred, Z_predicted, x_results
    """


    # correct xdata for B0 shift
    x_data = x_data - x_data[np.argmax(experimental_data)]
    initial_guess_offsets = np.array(initial_guess_offsets)

    # residual function
    L = lambda pars: lorentzian_sim(x_data,pars[0::3],pars[1::3], pars[2::3])
    def res_(p):
        return np.squeeze(experimental_data)- np.squeeze( L(p) )

    N_pools= initial_guess_offsets.__len__()
    Amplitudes_x0 = uniform(0,1,(N_pools,repetitions))
    Width_x0 = uniform(.1,5,(N_pools,repetitions))
    x0 = np.zeros((N_pools*3))
    x_results = np.zeros((N_pools*3,repetitions))

    up_lim, low_lim, init_offsets= build_limits(initial_guess_offsets)

    for i in range(10):
        x0[0::3] = Amplitudes_x0[:,i]
        x0[1::3] = Width_x0[:,i]
        x0[2::3] = init_offsets
        pars_predicted = least_squares(res_, x0, bounds=(low_lim,up_lim)).x
        x_results[:,i] = pars_predicted
    x_pred = np.mean(x_results,axis=1)
    Z_predicted = L(np.mean(x_results,axis=1))
    return x_pred, Z_predicted, x_results

def _build_limits(initial_guess_offsets):
    x0 = np.array(initial_guess_offsets)
    x0 = x0.astype(float)
    np.place(x0,x0==0,[0.1])

    N_pools= x0.__len__()
    lower_lim = np.zeros((N_pools*3))
    lower_lim[1::3]=0.1
    lower_lim[2::3]=x0 * 0.50

    upper_lim = np.ones((N_pools*3))
    upper_lim[1::3]=5
    upper_lim[2::3]=x0 * 1.50
    # swap if we have negative offsets
    true_lower = upper_lim[upper_lim<0]
    true_upper = lower_lim[lower_lim<0]

    upper_lim[upper_lim<0] = true_upper
    lower_lim[lower_lim<0] = true_lower

    return upper_lim, lower_lim, x0

def normalize(Zpectra,n):
    Zn = np.zeros_like(Zpectra.T)
    for idx, z in enumerate(Zpectra.T):
        Zn[idx,:] = z/z[n]
    return Zn.T
