import numpy as np
from scipy.integrate import quad
from scipy.optimize import fsolve
from scipy.special import ker, kei

eps = 1e-16
pi = np.pi
theta = pi/2
square = np.square
sqrt = np.sqrt
cos = np.cos
sin = np.sin
power = np.power



def grant_madsen_1979(uc, vc, uw, vw, omega, kb=0.2, rho=1025, kappa=0.41):
         # compute velocity magnitudes due to waves and currents
    ua_mag = sqrt(square(uc)+square(vc))    
    ub_mag = sqrt(square(uw)+square(vw))
    
    # compute the angle
    phi_c = np.arccos(np.clip((uw*uc+vw*vc+eps)/((ua_mag+eps)*(ub_mag+eps)),-1,1)) #eqn

    # load functions for gx and gy as defined in equations 9 and 10
    gx = lambda _theta: sin(_theta)*ub_mag + ua_mag*cos(phi_c) #eqn 9 
    gy = lambda _theta: ua_mag*sin(phi_c)*(_theta*0 + 1) #eqn 10
    
#    # change of coordinates to wave direction
#    u = sin(theta)*ub_mag + ua_mag*cos(phi_c) 
#    v = ua_mag*sin(phi_c) 

    
    alpha = square(ub_mag) + square(ua_mag) + 2 * ua_mag* ub_mag*cos(phi_c) #eqn 20 alpha * ub^2
    Ab = ub_mag/omega
    if ua_mag/ub_mag > 1/cos(phi_c):
#        print('Warning: We gotta problem |ua|/|ub| > 1/cos φ. Loop: ', ub_mag)
        theta_star = pi/2
    else:
        theta_star = np.arcsin(ua_mag/(ub_mag+eps)*cos(phi_c)) #eqn 12
    
    # Integrals for equations 11, 14, and 17
    integrand_1 = lambda _theta: sqrt(power(gx(_theta),4) + square(gx(_theta)) * square(gy(_theta)))
    integrand_2 = lambda _theta: sqrt(power(gy(_theta),4) + square(gx(_theta)) * square(gy(_theta)))
    integral_a = quad(integrand_1, -1*theta_star, pi+theta_star)[0]
    integral_b = quad(integrand_1, pi+theta_star, 2*pi-theta_star)[0]
    integral_c = quad(integrand_2, 0, 2*pi)[0] 
    
    # Compute V2 and the simplfied V2 for large waves
    V2 = 1/ (2*pi) *sqrt(square(integral_a-integral_b) + square(integral_c)) #eqn 14
    V2_large_waves = 2/pi*(ua_mag*ub_mag)*sqrt(4-3*square(sin(phi_c))) #eqn. 16
    phi_bar = np.arctan((integral_c+eps)/(integral_a-integral_b + eps)) #eqn. 17    

    #This is the workhorse that solves for fcw by solving equation 54 
    def min_gm(_fcw):
#        if _fcw < 0:
#            return [np.inf]
        ucw_star_mag = sqrt(1/2*_fcw*alpha) #eqn 21
        l = kappa*ucw_star_mag/omega #eqn 29
        zeta_0 = kb/30/l #eqn 31

        ztemp = 2*sqrt(zeta_0) #eqn 55
        K = 1/ztemp*1/sqrt(square(ker(ztemp)) + square(kei(ztemp))) #eqn 55
    
        temp_54_a = 0.097*sqrt(kb/Ab)*K/power(_fcw,0.75) #eqn 54
        temp_54_b = V2/square(ub_mag)/2/power(alpha/square(ub_mag),0.25) #eqn 54
        
        rhs_54 = square(temp_54_a) + 2*temp_54_a*temp_54_b*cos(phi_bar) #eqn 54
        lhs_54 = power(alpha/square(ub_mag), 0.75)/4 - square(temp_54_b) #eqn 54
        return rhs_54 - lhs_54
    
    
    # Solve for fcw
    fcw = fsolve(min_gm,1e-6)[0]
  
    # Compute terms of interest
    tau_c = 1/2 * rho *fcw  *V2 #eqn 11
    uc_star_mag = sqrt(1/2*fcw*V2) #eqn 15
    tau_bmax = 1/2 * rho *fcw * alpha #eqn 19
    ucw_star_mag = sqrt(1/2*fcw*alpha) #eqn 21
    l = kappa*ucw_star_mag/omega #eqn 29
    zeta_0 = kb/30/l #eqn 31
    
    beta = 1 - uc_star_mag/ucw_star_mag #eqn 49
    kbc = kb*power(24*ucw_star_mag/omega/kb,beta) #eqn 49
    
    data = {
        'phi_c': phi_c,
        'phi_bar': phi_bar,
        'V2': V2/square(ub_mag),
        'alpha': alpha/square(ub_mag),
        'V2_large_waves': V2_large_waves/square(ub_mag),
        'ua': ua_mag,
        'ub': ub_mag,
        'fcw': fcw,
        'tau_c': tau_c,
        'uc_star': uc_star_mag,
        'ucw_star': ucw_star_mag,
        'tau_b': tau_bmax,
        'l': l,
        'zeta_0': zeta_0,
        'kbc': kbc,
        'kb': kb,
        'Ab': Ab
    }
    
    return data
    
   



if __name__ == "__main__":
    import matplotlib.pyplot as plt

    kb = 0.02
#    UW = np.array([1,3,5])
#    VW = np.array([1,4,12])
#
#    UC = np.array([.51,.41,1])
#    VC = np.array([0,0,0])
#    omega = np.array([3,2,1])
#    N = len(UC)
#    for i in range(0,N):
#       gm = grant_madsen_1979(UC[i], VC[i], UW[i], VW[i], omega[i])
#       print(gm['tau_b'])
#input params    
    N = 100
    kbAb = .002
    uaub = np.linspace(1e-3,1.4,N)
    theta = 0# np.linspace(0,pi/2-eps,N)
    
    UC = 5*np.ones(N)#*sqrt(2)/2
    VC = 5*np.ones(N)#*sqrt(2)/2
    UW = UC*cos(theta)/uaub - VC*sin(theta)/uaub
    VW = UC*sin(theta)/uaub + VC*cos(theta)/uaub
    ua = sqrt(square(UC)+square(VC))
    ub = sqrt(square(UW)+square(VW))
    Ab = kb/kbAb
    omega = ub/Ab
    ab = ub/omega
    

    fcw = np.zeros(N)
    phi_c = np.zeros(N)
    phi_bar = np.zeros(N)
    V2 = np.zeros(N)
    kbc = np.zeros(N)
    newuaub = np.zeros(N)

    for i in range(0,N):
       gm = grant_madsen_1979(UC[i], VC[i], UW[i], VW[i], omega[i], kb)
       V2[i] = gm['V2']
       fcw[i] = gm['fcw']
       phi_c[i] = gm['phi_c']
       phi_bar[i] = gm['phi_bar']
       kbc[i] = gm['kbc']
       newuaub[i] = gm['ua']/gm['ub']
       print(gm['phi_c'],gm['kb']/gm['Ab'],gm['fcw'])
#    phi_bar_est = np.arctan(1/2*np.tan(phi_c)) 

#%%
    plt.plot(newuaub,fcw)
#    
#    
#    
#    #%%
#    plt.plot(phi_c,phi_bar,phi_c,phi_bar_est)
    
