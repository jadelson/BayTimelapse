import numpy as np
from scipy.integrate import quad
from scipy.optimize import fsolve
from scipy.special import ker, kei

eps = 1e-16
pi = np.pi
square = np.square
sqrt = np.sqrt
cos = np.cos
sin = np.sin
power = np.power



def grant_madsen_1979(uc_star, vc_star, uw, vw, omega, kb=0.2, rho=1000, kappa=0.4):
    uc = uc_star
    vc = vc_star
    uc_star_mag = sqrt(square(uc)+square(vc))
    ub_mag = sqrt(square(uw)+square(vw))   
    
    # compute velocity magnitudes due to waves and currents
  
    class GM:
        def process(self, ua_mag):
            uc = uc_star/uc_star_mag*ua_mag
            vc = vc_star/uc_star_mag*ua_mag
            # compute the angle between vectors
            vec_dot = uw*uc+vw*vc
            vec_cross = uc*vw-vc*uw
            self.phi_c = np.arctan2(vec_cross, vec_dot)
        
            # load functions for gx and gy as defined in equations 9 and 10
            gx = lambda _theta: sin(_theta) + ua_mag/ub_mag*cos(self.phi_c) #eqn 9 
            gy = lambda _theta: ua_mag/ub_mag*sin(self.phi_c) #eqn 10
            
            self.alpha = 1 + square(ua_mag/ub_mag) + 2 * ua_mag/ub_mag*cos(self.phi_c) #eqn 20
            self.Ab = ub_mag/omega
        
            if ua_mag/ub_mag >= np.abs(1/cos(self.phi_c)):
        #        print('Warning: We gotta problem |ua|/|ub| > 1/cos φ. Loop: ', ub_mag)
                self.theta_star = pi/2
            else:
                self.theta_star = np.arcsin(ua_mag/ub_mag*cos(self.phi_c)) #eqn 12
            
            # Integrals for equations 11, 14, and 17
            integrand_1 = lambda _theta: sqrt(power(gx(_theta),4) + square(gx(_theta)) * square(gy(_theta)))
            integrand_2 = lambda _theta: sqrt(power(gy(_theta),4) + square(gx(_theta)) * square(gy(_theta)))
            integral_a = quad(integrand_1, -1*self.theta_star, pi+self.theta_star)[0]
            integral_b = quad(integrand_1, pi+self.theta_star, 2*pi-self.theta_star)[0]
            integral_c = quad(integrand_2, 0, 2*pi)[0] 
            
            # Compute V2 and the simplfied V2 for large waves
            self.V2 = 1/ (2*pi) *sqrt(square(integral_a-integral_b) + square(integral_c)) #eqn 14
            
            self.V2_large_waves = 2/pi*(ua_mag/ub_mag)*sqrt(4-3*square(sin(self.phi_c))) #eqn. 16
        
            temp_17_denom = integral_a - integral_b
            if temp_17_denom == 0:
                temp_17_denom = eps
                
            self.phi_bar = np.arctan((integral_c)/(temp_17_denom)) #eqn. 17    
            #This is the workhorse that solves for fcw by solving equation 54 
            def get_fcw(_fcw):
        #        if _fcw < 0:
        #            return [np.inf]
                ucw_star_mag = sqrt(1/2*_fcw*alpha)*ub_mag#eqn 50
                l = kappa*ucw_star_mag/omega #eqn 29
                zeta_0 = kb/30/l #eqn 31
        
                ztemp = 2*sqrt(zeta_0) #eqn 55
                K = 1/ztemp*1/sqrt(square(ker(ztemp)) + square(kei(ztemp)))#eqn 55
            
                temp_54_a = 0.097*sqrt(kb/Ab)*K/power(_fcw,3.0/4.0) #eqn 54
                temp_54_b = V2/2.0/power(alpha,1.0/4.0) #eqn 54
                
                rhs_54 = square(temp_54_a) + 2.0*temp_54_a*temp_54_b*cos(phi_c) #eqn 54
                lhs_54 = power(alpha, 3.0/4.0)/4.0 - square(temp_54_b) #eqn 54
                return rhs_54 - lhs_54
            
            
            # Solve for fcw
            self.fcw = fsolve(get_fcw,1e-3)[0]
            self.uc_star_mag_guess = sqrt(1/2*self.fcw*self.V2)*ub_mag #eqn 15
            
            #return stuff

        
        
    def min_gm(ua_mag):
        gm = GM()
        gm.process(ua_mag)
        return gm.uc_star_mag_guess - uc_star_mag
    
    ua_mag = fsolve(min_gm, uc_star_mag)
    
#    print(min_gm(fcw))
    # Compute terms of interest
    gm = GM()
    gm.process(ua_mag)
    tau_c = 1/2 * rho *gm.fcw  *gm.V2*square(ub_mag) #eqn 11
    uc_star_mag = sqrt(1/2*gm.fcw*gm.V2)*ub_mag #eqn 15
    tau_bmax = 1/2 * rho *gm.fcw * gm.alpha*square(ub_mag) #eqn 19
    ucw_star_mag = sqrt(1/2*gm.fcw*gm.alpha)*ub_mag #eqn 21
    l = kappa*ucw_star_mag/omega #eqn 29
    zeta_0 = kb/30/l #eqn 31

    beta = 1 - uc_star_mag/ucw_star_mag #eqn 49
    kbc = kb*power(24*ucw_star_mag/ub_mag*gm.Ab/kb,beta) #eqn 49
    
    data = {
        'phi_c': gm.phi_c,
        'phi_bar': gm.phi_bar,
        'V2': gm.V2,
        'alpha': gm.alpha,
        'V2_large_waves': gm.V2_large_waves,
        'ua': ua_mag,
        'ub': ub_mag,
        'fcw': gm.fcw,
        'tau_c': tau_c,
        'uc_star': uc_star_mag,
        'ucw_star': ucw_star_mag,
        'tau_b': tau_bmax,
        'l': l,
        'zeta_0': zeta_0,
        'kbc': kbc,
        'kb': kb,
        'Ab': gm.Ab,
        'theta_star': gm.theta_star
    }
    
    return data
    
   



if __name__ == "__main__":
    import matplotlib.pyplot as plt 
    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)


    KBAB = np.array([2e-4, 4e-4, 6e-4, 8e-4])
    for j in range(0,3):
        KBAB = np.append(KBAB, KBAB*10)
    KBAB = np.append(KBAB,1)

    #figure 1 setup
    UAUB = np.arange(.1,1.01,.1)
    KBAB = [1]
    THETA = np.arange(0,pi/2,.2)
    N = len(THETA)
    
#    #figure 4 setup
#    UAUB = np.arange(.001,1.61,.11)*10
#    KBAB = np.array([2e-4, 4e-4, 8e-4, 2e-3, 4e-3, 8e-3, 2e-2, 4e-2, 6e-2, 2e-1, 4e-1, 6e-1, 1])
#    THETA = [0]
#    N = len(UAUB)
    
#    #table 1 setup
#    KBAB = [0.2, 0.0002]
#    UAUB = [0.025, 0.1, 0.6, 1.0, 1.2]
#    THETA = [pi/2]
#    N = 1
    for kbAb in KBAB:
        
        
        UC = np.random.randn(N)#*sqrt(2)/2
        VC = np.random.randn(N)#*sqrt(2)/2
        fcw = np.zeros(N)
        phi_c = np.zeros(N)
        phi_bar = np.zeros(N)
        theta_star = np.zeros(N)
        V2 = np.zeros(N)
        alpha = np.zeros(N)
        kbc = np.zeros(N)
        newuaub = np.zeros(N)

        
        for uaub in UAUB: 
            i = 0
            for theta in THETA:
                H = 3
                UW = UC[i]*cos(theta)/uaub - VC*sin(theta)/uaub
                VW = UC[i]*sin(theta)/uaub + VC*cos(theta)/uaub
                ua = sqrt(square(UC)+square(VC))
                ub = sqrt(square(UW)+square(VW))
                
                kb = 0.2
                Ab = kb/kbAb
                omega = ub/Ab 
                
                ustar_c = ua#0.4*ua/(np.log(H/(kb/30)) + kb/30/H - 1)
                
                
                gm = grant_madsen_1979(UC[i]/ua[i]*
                ustar_c[i], VC[i]/ua[i]*ustar_c[i], UW[i], VW[i], omega[i], kb)
                V2[i] = gm['V2']
                fcw[i] = gm['fcw']
                phi_c[i] = gm['phi_c']
                phi_bar[i] = gm['phi_bar']
                alpha[i] = gm['alpha']
                theta_star[i] = gm['theta_star']
                kbc[i] = gm['kbc']
                newuaub[i] = gm['ua']/gm['ub']
#           prin/t(gm['phi_c'],gm['kb']/gm['Ab'],gm['fcw'])
#    phi_bar_est = np.arctan(1/2*np.tan(phi_c)) 
                
#                #table 1 printout
#                print(uaub, gm['kb']/gm['Ab'] ,gm['ua']/gm['ub'], gm['kbc']/gm['kb'])
                
                i = i+1
            #%% figure 1
            line = ax.plot(THETA, phi_c, lw=.5)
            ax.set_xlim([0,pi/2])      

#        ##figure 4
#        line = ax.plot(newuaub, fcw, lw=.5)
#        ax.set_yscale('log')
#        ax.set_xlim([0, 2])      
#%%  
  
#    line = ax.plot(phi_c, phi_c, 'k--')
    
#    
#    #%%
#    plt.plot(phi_c,phi_bar,phi_c,phi_bar_est)
    
