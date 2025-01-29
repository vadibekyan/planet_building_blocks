#!/usr/bin/python
#Imports

import numpy as np
import os
import pandas as pd
import argparse


def planet_composition(c,o,mg,si,fe):
    """
    Input parameters are the absolute abundances of the elements
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


    C_abs, O_abs, Mg_abs, Si_abs, Fe_abs = abundances_solarTOabsolute(c,o,mg,si,fe)

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

    # in case nMg>nSi

    if(nMg>nSi):

        """
        In this case we must have: H2O, MgSiO3, Mg2SiO4

        We have the following master equations:
        1) nO=nH2O+3D0*nMgSiO3+4D0*nMg2SiO4
        2) nMg=nMgSiO3+2D0*nMg2SiO4
        3) nSi=nMgSiO3+nMg2SiO4
        Thus
        1) nMgSiO3=nSi-nMg2SiO4
        2) nMg=nSi-nMg2SiO4+2D0*nMg2SiO4=nSi+nMg2SiO4 => nMg2SiO4=nMg-nSi
        so nMgSiO3=nSi-(nMg-nSi)=2*nSi-nMg
        3) nO=nH2O+3D0*(2D0*nSi-nMg)+4D0*(nMg-nSi)
        so nH2O=nO-(3D0*(2D0*nSi-nMg)+4D0*(nMg-nSi))
        """

        nMg2SiO4=nMg-nSi
        nMgSiO3=2.*nSi-nMg
        nH2O=nO-(3.*nMgSiO3+4.*nMg2SiO4)
        nCH4=nC
        nSiO2=0.

    else:
        """
        In this case we must have: H2O, MgSiO3, SiO2

        Master equations
        1) nO=nH2O+3D0*nMgSiO3+2D0*SiO2
        2) nMg=nMgSiO3
        3) nSi=nMgSiO3+nSiO2
        thus
        nMgSiO3=nMg
        nSiO2=nSi-nMgSiO3
        nH2O=nO-(3D0*nMgSiO3+2D0*SiO2)    (nO-(3D0*nMgSiO3+2D0*(nSi-nMg))=nO-2D0*nMg-2D0*nSi
        """
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

def abundances_solarTOabsolute(c, o, mg, si, fe):
    """ Convert solar abundances to absolute abundances """
    C_abs = 10**(c + 8.46)  # Asplund et al. (2021)
    O_abs = 10**(o + 8.69)  # Asplund et al. (2021)
    Mg_abs = 10**(mg + 7.55)  # Asplund et al. (2021)
    Si_abs = 10**(si + 7.51)  # Asplund et al. (2021)
    Fe_abs = 10**(fe + 7.46)  # Asplund et al. (2021)
    return C_abs, O_abs, Mg_abs, Si_abs, Fe_abs

def bootstrap_composition(input_table_path):
    """
    Compute planet compositions for given input table.
    """
    input_sample = pd.read_csv(input_table_path)
    #Number of bootstraps
    N_bootstrap = 10000

    output_table = pd.DataFrame(input_sample.star, columns=['star'])
    output_table[['Ref', 'Ref_err', 'Silc', 'Silc_err', 'Z', 'Z_err', 'f_iron', 'f_iron_err', 'wf', 'wf_err']] = 0.0

    for n, star in enumerate(input_sample.star):
        output_composition = np.zeros((N_bootstrap, 5))
        for i in range(N_bootstrap):
            rand_norm_c = float(np.random.normal(input_sample.C[n], input_sample.C_err[n], 1)[0])
            rand_norm_o = float(np.random.normal(input_sample.O[n], input_sample.O_err[n], 1)[0])
            rand_norm_mg = float(np.random.normal(input_sample.Mg[n], input_sample.Mg_err[n], 1)[0])
            rand_norm_si = float(np.random.normal(input_sample.Si[n], input_sample.Si_err[n], 1)[0])
            rand_norm_fe = float(np.random.normal(input_sample.feh[n], input_sample.feh_err[n], 1)[0])
            
            output_composition[i] = planet_composition(rand_norm_c, rand_norm_o, rand_norm_mg, rand_norm_si, rand_norm_fe)
            
        means = np.mean(output_composition, axis=0)
        stds = np.std(output_composition, axis=0)
        output_table.loc[n, ['Ref', 'Silc', 'Z', 'f_iron', 'wf']] = np.round(means, 2)
        output_table.loc[n, ['Ref_err', 'Silc_err', 'Z_err', 'f_iron_err', 'wf_err']] = np.round(stds, 2)
        print (star)
        
    output_file = os.path.join(os.path.dirname(input_table_path), 'planet_compositions.csv')
    output_table.to_csv(output_file, index=False)
    print(f"Saved compositions to {output_file}")
            

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_table", type=str, help="Path to input CSV file")
    parser.add_argument("--starname", type=str, default="all", help="Specific star name or 'all'")
    args = parser.parse_args()
    bootstrap_composition(args.input_table)