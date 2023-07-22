import numpy as np 
from utilities import *
from SpeciesData import *

def VTFlash(v,T,Z,Pres,Species,Species_names):
    P = Pres
    Ki = wilsonEquation(T,P,Species,Species_names)
    Sp_i = Species
    Sp_n = Species_names
    #theta = 0.7 #initial guess
    er = 100

    if(np.sum(Z*Ki) - 1 < 0):
        print("sub cooled liquid")
        theta = 0
        x = Z/(1+theta*(Ki-1))
        y = x*Ki 
        return x,y,theta
    elif ((1 - np.sum(Z/Ki) ) > 0):
        print("supercritical gas")
        theta = 1
        x = Z/(1+theta*(Ki-1))
        y = x*Ki 
        return x,y,theta
    else:
        while (er>1E-3):
            # Ki,Z,Sp_i,Sp_n = DescendingOrder(Ki,Z,Sp_i,Sp_n)
            print("Ki = "+str(Ki))
            print('Species indices = ' + str(Sp_i))
            # c = np.zeros(len(Ki))
            #d = np.zeros(len(Ki))
            # for i in range(len(Ki)):
            #     c[i] = 1/(1-Ki[i]) 
            # c = 1/(1-Ki)
            # d = (c[0] - c)/(c[len(Ki)-1] - c[0])
            #print("d = " +str(d))
            #print("c = " +str(c))

            #step 1 - solving Rachford Rice
            error_RR = 10
            theta = 0.3 #initial guess #
            # sigma = (theta - c[0])/(c[len(Ki)-1] - theta) 
            count = 0

            # while error_RR>1E-10:
            #     count += 1
            #     #print("Z = " +str(Z))
            #     #print("sigma = " +str(sigma))
            #     S = np.sum( Z/(d + sigma*(1+d)))
            #     #print("S = " +str(S))
            #     G = (1+sigma)*S
            #     H = -sigma*(1+sigma)*S

            #     Sd = np.sum( -Z*(1+d)/((d + sigma*(1+d))**2))
            #     Gd = S + (1+sigma)*Sd
            #     Hd = -(1+2*sigma)*S - sigma*(1+sigma)*Sd
                
            #     sigmaold = sigma
            #     if(G > 0):
            #         sigma = sigma - G/Gd
            #     elif(G < 0):
            #         sigma = sigma - H/Hd
            #     error_RR = np.abs(sigma - sigmaold)/sigma
            
            flag = -1
            if( np.sum(Z*(Ki - 1)/(1 + 0.5*(Ki - 1))) ):
                flag = 0
            else: 
                flag = 1
            if(flag = 0):
                while error_RR>1E-10:
                    count += 1
                    g = np.sum(Z*(Ki - 1)/(1 + theta*(Ki - 1)))
                    gd = -np.sum(Z*((Ki - 1)/(1 + theta*(Ki - 1)))**2)
                    theta = theta - g/gd
                    error_RR = np.abs(g/gd)
            elif(flag = 1):
                while error_RR>1E-10:
                    count += 1
                    g = np.sum(Z*(Ki - 1)/(theta_l + theta_l*(Ki - 1)))
                    gd = -np.sum(Z*((Ki - 1)/(1 + theta*(Ki - 1)))**2)
                    theta = theta - g/gd
                    error_RR = np.abs(g/gd)

            # print('sigma = '+str(sigma))
            print("count = "+str(count))
            # theta = (c[0] + sigma*c[len(Ki)-1]) / (1 + sigma)
            print('theta = '+str(theta))

            #step 2
            x = Z/(1+theta*(Ki-1))
            y = x*Ki
            print('x = '+str(x))
            print('y = '+str(y))
            a_L, b_L = Mixture_AB(x,T,P,Sp_i,Sp_n)
            a_V, b_V = Mixture_AB(y,T,P,Sp_i,Sp_n)

            #step 3 
            vl, vv = liq_vol(v,T,a_L,b_L,a_V,b_V,theta,Sp_i)
            print('liq vol = ' + str(vl))
            print('gas vol = ' + str(vv))

            #step 4
            fl = fugacity(x,T,P,vl,a_L,b_L,Sp_i,Sp_n)
            fv = fugacity(y,T,P,vv,a_V,b_V,Sp_i,Sp_n)
            print('fugacity liq = ' +str(fl))
            print('fugacity vap = ' +str(fv))
            #er = np.sum(np.abs(fl-fv))
            Kiold = Ki 
            Ki = fl/fv
            er = np.sum(np.abs(np.log(Ki/Kiold)))
            print("error = "+str(er))
            print("\n")
            #Ki = y/x
        
        # x,y = CorrectOrder(x,y,Sp_i)

        return x,y,theta

def wilsonEquation(T,P,Species,Species_names):
    Pc = np.zeros(len(Species))
    Tc = np.zeros(len(Species))
    omega = np.zeros(len(Species))

    for i in range(len(Species)):
        Pc[i]= SpeciesData[Species_names[i]]['Pc']
        Tc[i]= SpeciesData[Species_names[i]]['Tc']
        omega[i]= SpeciesData[Species_names[i]]['omega']

    Ki = (Pc/P) * np.exp( 5.373 * (1 + omega) * (1 - Tc/T ))
    return Ki

def Mixture_AB(x,T,P,Species,Species_names):
    ai = np.zeros(len(Species))
    bi = np.zeros(len(Species))
    beta = SpeciesData['beta']
    
    R = 8.314
    Pc = np.zeros(len(Species))
    Tc = np.zeros(len(Species))
    omega = np.zeros(len(Species))

    for i in range(len(Species)):
        Pc[i]= SpeciesData[Species_names[i]]['Pc']
        Tc[i]= SpeciesData[Species_names[i]]['Tc']
        omega[i]= SpeciesData[Species_names[i]]['omega']

        if(omega[i] < 0.5):
            c_omega = 0.37464 + 1.54226*omega[i] - 0.26992 * omega[i]**2
        elif(omega[i] >= 0.5):
            c_omega = 0.3796 + 1.485*omega[i] - 0.1644 * omega[i]**2 + 0.01667*omega[i]**3

        omega_a = 0.45724
        omega_b = 0.0778
        ai[i] = omega_a *((R*Tc[i])**2 / Pc[i]) * (1 + c_omega * (1 - np.sqrt(T/Tc[i]))) ** 2
        bi[i] = omega_b*R*Tc[i]/Pc[i]

    a = 0
    b = 0
    for i in range(len(Species)):
        for j in range(len(Species)):
            a += x[i]*x[j]*(1-beta[i][j])*np.sqrt(ai[i]*ai[j])
        b += x[i]*bi[i]

    return a,b

def liq_vol(v,T,a_l,b_l,a_v,b_v,theta,Sp_i):
    d1 = 1+np.sqrt(2)
    d2 = 1-np.sqrt(2)
    R = 8.314

    a1_L = b_l*(d1+d2)-a_l/(R*T)
    a2_L = b_l*(b_l*d1*d2 + a_l/(R*T))
    a3_L = b_l*(d1*d2 - 1)
    a4_L = (b_l**2)*(d1+d2 - d1*d2)
    a5_L = -(b_l**3)*d1*d2

    a1_V = b_v*(d1+d2)-a_v/(R*T)
    a2_V = b_v*(b_v*d1*d2 + a_v/(R*T))
    a3_V = b_v*(d1*d2 - 1)
    a4_V = (b_v**2)*(d1+d2 - d1*d2)
    a5_V = -(b_v**3)*d1*d2

    c0 = (a2_L*a5_V - a2_V*a5_L)*theta**3 - (a5_L*a1_V - a4_V*a2_L)*v*theta**2 + (a2_L*a3_V - a5_L)*theta*v**2 + a2_L*v**3
    c1 = (a1_L*a5_V - a5_L*a1_V + a4_V*a2_L - a2_V*a4_L)*theta**3 + (a1_L*a3_V + 3*a2_L - a4_L)*v**2*theta + (v*a1_L - 3*a2_L)*v**2 + 2*(a5_L-a2_L*a3_V)*v*theta + (a1_L*a4_V - a1_V*a4_L + 2*a2_L*a3_V - 2*a5_L)*v*theta**2 + (a5_L*a1_V - a4_V*a2_L)*theta**2
    c2 = (a1_L*a4_V - a1_V*a4_L + a2_L*a3_V - a2_V*a3_L - a5_L + a5_V)*theta**3 + (v**2 - 3*v*a1_L + 3*a2_L)*v + (2*a1_L*a3_V - a1_V*a3_L + 3*a2_L - 2*a4_L + a4_V)*v*theta**2 + (a1_V*a4_L - a1_L*a4_V - 2*a2_L*a3_V + 2*a5_L)*theta**2 + (3*a1_L - a3_L + a3_V)*v**2*theta + 2*(a4_L - a1_L*a3_V - 3*a2_L)*v*theta + (a2_L*a3_V - a5_L)*theta
    c3 = (a1_L*a3_V - a1_V*a3_L + a2_L - a2_V - a4_L + a4_V)*theta**3 + (3*a1_L - a1_V -2*a3_L + 2*a3_V)*v*theta**2 + (a1_V*a3_L - 2*a1_L*a3_V - a4_V - 3*a2_L + 2*a4_L)*theta**2 + (-6*a1_L + 2*a3_L - 2*a3_V)*v*theta + (2*v**2 + a1_L*a3_V + 3*a2_L - a4_L)*theta - 3*v**2 + 3*v*a1_L - a2_L
    c4 = ( (a1_L - a1_V - a3_L + a3_V)*theta**2 + (v - 2*a1_L + a3_L - a3_V)*theta - 3*v + a1_L) * (theta - 1)
    c5 = - (theta - 1)**2

    #using newton iteration to find vl
    vl = 0.99*v #initial guess
    error = 10
    # coeff = [c5,c4,c3,c2,c1,c0]
    # vol = np.roots(coeff)

    # for i in range(len(vol)):
    #     if(np.imag(vol[i]) == 0):
    #         vl = np.real(vol[i])
    #for i in range(2000):
    while error > 1E-5:    
        fx = c5*vl**5 + c4*vl**4 + c3*vl**3 + c2*vl**2 + c1*vl + c0
        fdx = 5*c5*vl**4 + 4*c4*vl**3 + 3*c3*vl**2 + 2*c2*vl + c1
        vl = vl - fx/fdx
        error = np.abs(fx/fdx)
    
    print(error)

    #update vv from vl  
    vv = (v - (1-theta)*vl)/theta
    return vl, vv

def fugacity(x,T,P,vl,a_L,b_L,Species,Species_names):
    R = 8.314
    Z = P*vl/(R*T)

    ai = np.zeros(len(Species))
    bi = np.zeros(len(Species))
    beta = SpeciesData['beta']

    Pc = np.zeros(len(Species))
    Tc = np.zeros(len(Species))
    omega = np.zeros(len(Species))

    for i in range(len(Species)):
        Pc[i]= SpeciesData[Species_names[i]]['Pc']
        Tc[i]= SpeciesData[Species_names[i]]['Tc']
        omega[i]= SpeciesData[Species_names[i]]['omega']

        if(omega[i] < 0.5):
            c_omega = 0.37464 + 1.54226*omega[i] - 0.26992 * omega[i]**2
        elif(omega[i] >= 0.5):
            c_omega = 0.3796 + 1.485*omega[i] - 0.1644 * omega[i]**2 + 0.01667*omega[i]**3

        omega_a = 0.45724
        omega_b = 0.0778
        ai[i] = omega_a *((R*Tc[i])**2 / Pc[i]) * (1 + c_omega * (1 - np.sqrt(T/Tc[i]))) ** 2
        bi[i] = omega_b*R*Tc[i]/Pc[i]

    Ai = ai*P/(R*T)**2
    Bi = bi*P/(R*T)

    Amix = 0
    Bmix = 0

    for i in range(len(Species)):
        for j in range(len(Species)):
            Amix += x[i]*x[j]*(1-beta[i][j])*np.sqrt(Ai[i]*Ai[j])
        Bmix += x[i]*Bi[i]
    
    #print(Z-Bmix)
    #print((Z+(1+np.sqrt(2))*Bmix)/(Z+(1-np.sqrt(2))*Bmix))
    
    f = P*x*np.exp( (Bi/Bmix)*(Z-1) - np.log(Z-Bmix) - ( Amix / (2*np.sqrt(2)*Bmix) ) * (2*np.sum(x*Ai)/Amix - Bi/Bmix ) * np.log( (Z+(1+np.sqrt(2))*Bmix)/(Z+(1-np.sqrt(2))*Bmix)) )

    return f
