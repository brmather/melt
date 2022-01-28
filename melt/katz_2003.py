import numpy as np
from scipy.optimize import root
import warnings

def T_solidus(P):
    return 1085.7 + 132.9*P - 5.1*P**2

def T_liquidus_lherz(P):
    return 1475.0 + 80.0*P - 3.2*P**2

def T_liquidus(P):
    return 1780.0 + 45.0*P - 2.0*P**2

def R_cpx(P):
    return 0.5 + 0.08*P

def X_sat(P):
    return 12.0*P**0.6 + 1.0*P

def X_H2O(X, F):
    return X/(0.01 + F*(1.0 - 0.01))

def delta_T(X):
    return 43.0*X**0.75


def F_dry(P,T,M):

    T_s   = T_solidus(P)
    T_l   = T_liquidus(P)
    T_lh  = T_liquidus_lherz(P)
    
    R = R_cpx(P)
    F_cpx_out = M/R
    T_cpx_out = F_cpx_out**(1.0/1.5) * (T_lh - T_s) + T_s
    
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        F_opx = F_cpx_out + (1.0 - F_cpx_out)*((T - T_cpx_out)/(T_l - T_cpx_out))**1.5
        F_cpx = ((T - T_s)/(T_lh - T_s))**1.5

    F = F_cpx.copy()
    F[F > F_cpx_out] = F_opx[F > F_cpx_out]
    return F

def F_wet(P,T,X):
    T_s   = T_solidus(P)
    T_l   = T_liquidus(P)
    T_lh  = T_liquidus_lherz(P)

    # Evaluate anhydrous melting
    R = R_cpx(P)
    M = 0.15
    F_cpx_out = M/R
    F_opx = F_dry(P, T, M)

    X_s = X_sat(P)
    delT_sat = delta_T(X_s)

    def F_iter(F):

        F[np.isnan(F)] = 0.0
        Xi = X_H2O(X, F)

        delT = delta_T(Xi)
        delT = np.minimum(delT, delT_sat)
        delT = np.clip(delT, 0, 1e99)

        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            F_new = ((T - (T_s - delT))/(T_lh - T_s))**1.5

        F_new[np.isnan(F_new)] = 0

        return F_new - F

    sol = root(F_iter, x0=F_opx, method='anderson')
    F = sol.x

    Xi = X_H2O(X, F)

    delT = delta_T(Xi)
    delT = np.minimum(delT, delT_sat)
    delT = np.clip(delT, 0, 1e99)

    T_s   = T_solidus(P) - delT
    T_l   = T_liquidus(P) - delT
    T_lh  = T_liquidus_lherz(P) - delT

    T_cpx_out = F_cpx_out**(1.0/1.5) * (T_lh - T_s) + T_s

    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        F_opx = F_cpx_out + (1.0 - F_cpx_out)*((T - T_cpx_out)/(T_l - T_cpx_out))**1.5

    F[F > F_cpx_out] = F_opx[F > F_cpx_out]
    return F