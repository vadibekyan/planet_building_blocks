#!/usr/bin/python
# Imports
import numpy as np
import os
import pandas as pd

def abundances_solar_to_absolute(c, o, mg, si, fe):
    """
    Converts the relative-to-Sun abundances to absolute abundances.
    
    Parameters:
    c (float): Relative-to-Sun abundance of Carbon.
    o (float): Relative-to-Sun abundance of Oxygen.
    mg (float): Relative-to-Sun abundance of Magnesium.
    si (float): Relative-to-Sun abundance of Silicon.
    fe (float): Relative-to-Sun abundance of Iron.

    Returns:
    tuple: Absolute abundances of Carbon, Oxygen, Magnesium, Silicon, and Iron.
    
    Note:
    Solar reference abundances are from Asplund et al. (2021).
    """

    # Convert relative-to-Sun abundances to absolute abundances using solar reference values
    C_abs = 10**(c + 8.46)
    O_abs = 10**(o + 8.69)
    Mg_abs = 10**(mg + 7.55)
    Si_abs = 10**(si + 7.51)
    Fe_abs = 10**(fe + 7.46)

    return C_abs, O_abs, Mg_abs, Si_abs, Fe_abs


def planet_composition(c,o,mg,si,fe):
    """
    Calculate the planetary composition based on the observed abundances of elements.

    Parameters:
    c (float): Relative-to-Sun abundance of Carbon.
    o (float): Relative-to-Sun abundance of Oxygen.
    mg (float): Relative-to-Sun abundance of Magnesium.
    si (float): Relative-to-Sun abundance of Silicon.
    fe (float): Relative-to-Sun abundance of Iron.

    Returns:
    tuple: Mass fractions of refractories, silicates, total metal, iron, and water.
    """


    #Atomic masses
    muH= 1.00794
    muHe=4.002602
    muC=12.0107
    muN=14.0067
    muO=15.9994
    muMg=24.305
    muSi=28.0855
    muFe=55.845

    #Primordial abundances
    Xprim=0.76
    Yprim=0.24

    #Molecular weights
    muH2O=muO+2.*muH
    muCO=muC+muO
    muCO2=muC+2.*muO
    muCH4=muC+4.*muH
    muNH3=muN+3.*muH
    muMgSiO3=muMg+muSi+3.*muO
    muMg2SiO4=2.*muMg+muSi+4.*muO
    muSiO2=muSi+2.*muO

    #Absolute abundances
    H_abs = 1*10**12
    He_abs = 7.93*10**10


    C_abs, O_abs, Mg_abs, Si_abs, Fe_abs = abundances_solar_to_absolute(c,o,mg,si,fe)

    #number of atoms
    nH = H_abs
    nHe = He_abs
    nC = C_abs
    nO = O_abs
    nMg = Mg_abs
    nSi = Si_abs
    nFe = Fe_abs

    #total number of atoms
    NNTOT=nH+nHe+nC+nO+nMg+nSi+nFe
    MMTOT=nH*muH+nHe*muHe+nC*muC+nO*muO+nMg*muMg+nSi*muSi+nFe*muFe

    #number fractions of molecules
    nnH=nH/NNTOT
    nnHe=nHe/NNTOT
    nnC=nC/NNTOT
    nnO=nO/NNTOT
    nnMg=nMg/NNTOT
    nnSi=nSi/NNTOT
    nnFe=nFe/NNTOT

    #total mass
    MTOT=nnH*muH+nnHe*muHe+nnC*muC+nnO*muO+nnMg*muMg+nnSi*muSi+nnFe*muFe

    #mass fraction
    maH=nnH*muH/MTOT
    maHe=nnHe*muHe/MTOT
    mC=nnC*muC/MTOT
    mO=nnO*muO/MTOT
    mMg=nnMg*muMg/MTOT
    mSi=nnSi*muSi/MTOT
    maFe=nnFe*muFe/MTOT

    if(nMg>nSi):
        nMg2SiO4=nMg-nSi
        nMgSiO3=2.*nSi-nMg
        nH2O=nO-(3.*nMgSiO3+4.*nMg2SiO4)
        nCH4=nC
        nSiO2=0.

    else:
        nMg2SiO4=0.
        nMgSiO3=nMg
        nSiO2=nSi-nMgSiO3
        nH2O=nO-(3.*nMgSiO3+2.*nSiO2)
        nCH4=nC


    nH=nH-2.*nH2O-4.*nCH4

    mMg2SiO4=muMg2SiO4*nMg2SiO4
    mMgSiO3=muMgSiO3*nMgSiO3
    mH2O=muH2O*nH2O
    mCH4=muCH4*nCH4
    mSiO2=muSiO2*nSiO2
    mFe=muFe*nFe
    mH=muH*nH
    mHe=muHe*nHe

    mmMg2SiO4=mMg2SiO4/MMTOT
    mmMgSiO3=mMgSiO3/MMTOT
    mmH2O=mH2O/MMTOT
    mmCH4=mCH4/MMTOT
    mmSiO2=mSiO2/MMTOT
    mmH=mH/MMTOT
    mmFe=mFe/MMTOT
    mmHe=mHe/MMTOT

    refractories = (mmFe+mmMgSiO3+mmMg2SiO4+mmSiO2)*100.
    silicates = (mmMgSiO3+mmMg2SiO4+mmSiO2)*100.
    z = (mmCH4+mmH2O+mmFe+mmMgSiO3+mmMg2SiO4+mmSiO2)*100.
    iron_frac = mmFe/(mmFe+mmMgSiO3+mmMg2SiO4+mmSiO2)*100.
    water_frac = mmH2O/(mmFe+mmMgSiO3+mmMg2SiO4+mmSiO2+mmH2O)*100.

    return refractories, silicates, z, iron_frac, water_frac


if __name__ == '__main__':
    """
    Example: Calculate the planetary composition for solar abundances.
    """
    refractories, silicates, z, iron_frac, water_frac = planet_composition(0,0,0,0,0)
    print (f'Iron mass fraction = {iron_frac:.2f}')