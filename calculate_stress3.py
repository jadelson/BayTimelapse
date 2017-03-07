#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 27 01:53:59 2017

@author: jadelson
"""

import pickle
import numpy as np
import SFBayRemoteSensing.conversions.calculate_fetch as fetch
import utm
from scipy.optimize import minimize

def calculate_stress(data):
    kappa = 0.41
    rho = 1023
    nu = 10**-6
    z0 = 10**-4
    g = 9.81
    
    
    
    m = data['a0'].shape
    
    data['H'] = np.empty(m)
    data['k'] = np.empty(m)
    data['A'] = np.empty(m)
    data['k_x'] = np.empty(m)
    data['k_y'] = np.empty(m)
    data['fw'] = np.empty(m)
    data['omega'] = np.empty(m)
    data['Rew'] = np.empty(m)
    data['Cd'] = np.empty(m)
    
    data['tau_m_x'] = np.empty(m)
    data['tau_w_x'] = np.empty(m)
    data['tau_m_y'] = np.empty(m)
    data['tau_w_y'] = np.empty(m)
    data['tau_x'] = np.empty(m)
    data['tau_y'] = np.empty(m)
    data['tau_w'] = np.empty(m)
    data['tau_m'] = np.empty(m)

    data['tau'] = np.empty(m)
    
    data['omega'] = 2*np.pi/data['T']

    depth = data['depth']
    surface = data['eta']
    u_vel_1 = data['u']
    v_vel_1 = data['v']
    
    #  utm.to_latlon(
    #  utm.from_latlon(
    for i in range(0,m[0]):
#        if np.mod(i,2500) == 0:
#            print(i*1./m[0])
        
        H = depth[i] + surface[i]
        Cd = kappa**2*(np.log(H/z0) + z0/H - 1)**-2
        u = u_vel_1[i]
        v = v_vel_1[i]
        mag_u = np.sqrt(u**2 + v**2)
        tau_m_x = rho*Cd*mag_u*u
        tau_m_y = rho*Cd*mag_u*v
            
        omega = data['omega'][i]
        
        a0 = data['a0'][i]
        angle = data['theta'][i]   
#        if np.isnan(Tw) or np.isnan(a0):
#    #        t.append(i)
#            continue
        

        def dispersion_relation(k):
            return np.abs(omega**2 - g*k*np.tanh(k*H))
        k = minimize(dispersion_relation,omega/(np.sqrt(g*H)),
                options={'disp': False})
        mag_k = np.abs(k.x)
        
        kx = mag_k*np.cos(angle)
        ky = -mag_k*np.sin(angle)
    
        A = a0/(np.sinh(mag_k*H))
        Rew = A**2*omega/nu
        fw = 2*Rew**(-0.5)
        mag_tau_w = 0.5*rho*fw*(A*omega)**2
        tau_w_x = mag_tau_w*kx/mag_k
        tau_w_y = mag_tau_w*ky/mag_k
    
        tau_x = tau_w_x + tau_m_x
        tau_y = tau_w_y + tau_m_y
    
        tau = np.sqrt(tau_x**2 + tau_y**2)    
    
        
        data['H'][i] = H
        data['k'][i] = mag_k
        data['k_x'][i] = kx
        data['k_y'][i] = ky
        data['A'][i] = A
        data['fw'][i] = fw
        
        data['Rew'][i] = Rew
        data['Cd'][i] = Cd
    
        data['tau_m_x'][i] = tau_m_x
        data['tau_w_x'][i] = tau_w_x
        data['tau_m_y'][i] = tau_m_y
        data['tau_w_y'][i] = tau_w_y
        data['tau_w'][i] = mag_tau_w
        data['tau_m'][i] = np.sqrt(tau_m_x*tau_m_x + tau_m_y*tau_m_y)
        data['tau_x'][i] = tau_x
        data['tau_y'][i] = tau_y
        data['tau'][i] = tau
    

    return data







