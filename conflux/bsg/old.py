################
# Fermi theory
def nuclear_radius(A):
    result = 1.121*pow(A,1/3.0)+2.426*pow(A,-1/3.0)-6.614/A
    return result

WO = lambda energy: energy/(1.0*ELECTRON_MASS_MEV) + 1.0

p = lambda W: np.sqrt(W**2 - 1.)

y = lambda W, Z: 1.0*Z*W/p(W)*constants.fine_structure

gamma = lambda Z: np.sqrt(1.0 - (constants.fine_structure*Z*1.0)**2)

phasespace = lambda W, W0: p(W)*W*(W0-W)*(W0-W)

def neutrino(enu, Parameters):
    np.seterr(divide='ignore')
    buf = Parameters
    result = 0.
    R= 1.0*FERMI_to_W * nuclear_radius(buf.A)
    W0=WO(buf.e0)
    W=W0-enu/(ELECTRON_MASS_MEV*1.0)
    thresh=(W0-1.0-V0(buf.Z))*(ELECTRON_MASS_MEV*1.0)
    buf.thresh = thresh
    result = forbidden(W,W0,buf.WM,buf.ftype)*GN(W)*(L0(W, buf.Z, R, gamma(buf.Z))+L0b(W,buf.Z,R))*CC(R, buf.Z, W, W0)*phasespace(W, W0)*F(y(W,buf.Z), gamma(buf.Z), p(W), R)

    belowThresh = enu<thresh   # prepared for a list of comparison
    #print("E value: ",buf.e0, "calc: ", result, "stdev: ", np.std(result))
    return result*(S(buf.e0-enu, buf.Z)**belowThresh)

def electron(ebeta, Parameters):
    buf = Parameters
    result = 0.
    R = 1.0*FERMI_to_W * nuclear_radius(buf.A)
    W0=WO(buf.e0)
    W=WO(ebeta)
    thresh=V0(buf.Z)*(ELECTRON_MASS_MEV*1.0)
    buf.thresh = thresh
    result = forbidden(W,W0,buf.WM,buf.ftype)*G(W,W0)*(L0(W, buf.Z, R, gamma(buf.Z))+L0b(W,buf.Z,R))* CC(R, buf.Z, W, W0)*phasespace(W, W0)*F(y(W,buf.Z), gamma(buf.Z), p(W), R)

    aboveThresh = (ebeta>=thresh)   # prepared for a list of comparison
    return result*(S(ebeta,buf.Z)**aboveThresh)

def F(y, gamma, p, R):
    result = 0.
    res = special.loggamma(gamma+y*1j)
    absgamma2= np.exp(2*res.real)
    g2g = special.gamma(2*gamma+1)
    result = 2*(gamma+1)*absgamma2*np.exp(y*np.pi)/(pow(2*p*R,(2*(1-gamma)))*g2g**2)
    return result

###############
# Finite size

l0dat = [[0.115, -1.8123, 8.2498, -11.223, -14.854, 32.086],
		    [-0.00062, 0.007165, 0.01841, -0.53736, 1.2691, -1.5467],
		    [0.02482, -0.5975, 4.84199, -15.3374, 23.9774, -12.6534],
		    [-0.14038, 3.64953,-38.8143,172.1368,-346.708, 288.7873],
		    [0.008152,-1.15664,49.9663,-273.711,657.6292, -603.7033],
		    [1.2145, -23.9931,149.9718,-471.2985, 662.1909,-305.6804],
		    [-1.5632,33.4192,-255.1333,938.5297,-1641.2845,1095.358]]

def L0(W, Z, R, gamma):
    result = 0.
    alpha = constants.fine_structure
    result = 1. + (13/60.)*(alpha*Z)*(alpha*Z) - (W*R*alpha*Z*(41 - 26*gamma))/(15*(2*gamma - 1)) - (alpha*Z*R*gamma*(17 - 2*gamma))/(30*W*(2*gamma -1)) - 0.41*(R-0.0164)*pow(alpha*Z,4.5);
    return result

def af(j0, Z):
    j = j0+1
    result = 0.
    for i in range(0, 6):
        result += l0dat[j][i]*pow(Z*constants.fine_structure, i+1.0)
    return result

def L0b(W, Z, R):
    result = 0.
    for k in range(0, 6):
        result += af(k,Z)*pow(W*R,k)
    result=result+af(-1,Z)*R/W;
    return result

######################
# Screening Correction
#TODO: build a screening effect database to replace

Wb = lambda energy, V0: WO(energy)-V0

nn_energy =[1,8,13,16,23,27,29,49,84,92]
nn_mfp=[1,1.42,1.484,1.497,1.52,1.544,1.561,1.637,1.838,1.907]
nn_spline = interpolate.InterpolatedUnivariateSpline(nn_energy,nn_mfp)

def NN(Z):
    return nn_spline(Z)

def V0(Z):
    return NN(Z-1)*constants.fine_structure**2*pow(Z-1, 4/3.)

def S(energy, Z):
    result = 0.
    result = (Wb(energy,V0(Z))/WO(energy))*pow(p(Wb(energy,V0(Z)))/p(WO(energy)),2*gamma(Z)-1)*np.exp(-np.pi*(y(WO(energy),Z)-y(Wb(energy,V0(Z)),Z)))
    res = special.loggamma((gamma(Z) + 1j*y(Wb(energy,V0(Z)),Z)))
    absgamma2b= np.exp(2*res.real)
    res = special.loggamma((gamma(Z) + 1j*y(WO(energy),Z)))
    absgamma2= np.exp(2*res.real)
    result = result*absgamma2b/absgamma2
    return result

#####################
# Shape Corrections

def CC(R, Z, W, W0):
    result = 0.
    alpha = constants.fine_structure*1.0
    result += (-(233/630.0))*(alpha*Z)*(alpha*Z) - (1/5.0)*(W0*R)*(W0*R) + (2.0/35)*W0*R*alpha*Z
    result += W*((-(21/35.0))*R*alpha*Z + (4*W0*R*R)/9.0)
    result += -((4*R*R)/9.0)*W*W
    return result+1

###########################
# QED radiative corrections

def G(W, W0):
    result = 0.
    beta = p(W)/W
    result += 3*np.log(NUCLEON_MASS_W)-3/4.0+4*(np.arctanh(beta)/beta-1.0)*((W0-W)/(3*W)-3/2.0+np.log(2*(W0-W)))
    result += (4.0*special.spence(1-(2*beta)/(1+ beta)))/beta
    result += (np.arctanh(beta)*(2.0*(1.0+beta**2)+(W0-W)*(W0-W)/(6.0*W**2)-4.0*np.arctanh(beta)))/beta
    result = result*constants.fine_structure/(2.0*np.pi)+1.0
    return result

def GN(W):
    result=0.
    beta=np.sqrt(W*W-1.0)/W
    result += 3.0*np.log(PROTON_MASS_W)+23.0/4.0-8.0/beta * special.spence(1-2.0 * beta/(1+beta))
    result += 8.0*(np.arctanh(beta)/beta-1.0)*np.log(2.0*W*beta)
    result += 4.0*np.arctanh(beta)/beta*((7.0+3.0*beta*beta)/8.0-2.0*np.arctanh(beta))
    result = result*constants.fine_structure/(2.0*np.pi)+1.0
    return result

#######################
# Forbidden shapes

def forbidden(W, W0, WM, ftype):
    p1 = p(W)
    pn = W0-W
    result = 1.0
    enu = pn*ELECTRON_MASS_MEV*1.0
    ee = W*ELECTRON_MASS_MEV*1.0
    beta = p1/W
    shape = 0.
    wm = 0.
    pe = np.sqrt(ee**2 - ELECTRON_MASS_MEV**2*1.0)

    if ((ftype==1) or (ftype==-2)): # first unique, 2nd non-unique
        result = (pn*pn+p1*p1)*(1+W*ELECTRON_MASS_MEV*WM)
    if ((ftype==2) or (ftype==-3)): # 2nd unique, 3rd non-unique
        result=(pow(pn,4)+10.0/3*pn*pn*p1*p1+pow(p1,4))*(1+W*ELECTRON_MASS_MEV*WM)
    if ((ftype==3) or (ftype==-4)): # 3rd unique, 4th non-unique
        result=(pow(pn,6)+7.0*pow(pn,4)*p1*p1+7.0*pow(p1,4)*pn*pn+pow(p1,6))*(1+W*ELECTRON_MASS_MEV*WM)
    if (ftype == -10):  #   first non-unique 0-
        shape=(pe*pe+enu*enu+2*beta*beta*enu*ee)
        wm=0
        wm=wm*WM
        result=(1+wm)*shape
    if (ftype==-11):    #   first non-unique 1-
        shape= pe*pe + enu*enu - 4.0/3.0*beta*beta*enu*ee
        wm=((pe*pe + enu*enu)*(beta*beta*ee - enu) + 2.0*beta*beta*ee*enu*(enu - ee)/3.0)/((pe*pe + enu*enu - 4.0*beta*beta*enu*ee/3.0))
        wm=wm*WM
        result=(1+wm)*shape
    if (ftype == 10):   #   first unique, 2-
        shape=pe**2 + enu**2
        wm=3.0/5.0*((pe*pe + enu*enu)*(beta*beta*ee - enu) + 2.0*beta*beta*ee*enu*(enu - ee)/3.0)/((pe*pe + enu*enu))
        wm=wm*WM
        result=(1+wm)*shape

    return result
