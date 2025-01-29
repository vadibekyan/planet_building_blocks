#!/usr/bin/python
# Imports
import numpy as np
import os
import pandas as pd
import argparse

def planet_composition(c, o, mg, si, fe):
    """
    Input parameters are the absolute abundances of the elements
    """
    # Atomic masses
    muH = 1.00794
    muHe = 4.002602
    muC = 12.0107
    muN = 14.0067
    muO = 15.9994
    muMg = 24.305
    muSi = 28.0855
    muFe = 55.845

    # Primordial abundances
    Xprim = 0.76
    Yprim = 0.24

    # Absolute abundances
    H_abs = 1 * 10**12
    He_abs = 7.93 * 10**10
    
    C_abs, O_abs, Mg_abs, Si_abs, Fe_abs = abundances_solarTOabsolute(c, o, mg, si, fe)
    
    # Number of atoms
    nH = H_abs
    nHe = He_abs
    nC = C_abs
    nO = O_abs
    nMg = Mg_abs
    nSi = Si_abs
    nFe = Fe_abs
    
    # Total number of atoms
    NNTOT = nH + nHe + nC + nO + nMg + nSi + nFe
    MMTOT = nH * muH + nHe * muHe + nC * muC + nO * muO + nMg * muMg + nSi * muSi + nFe * muFe
    
    # Mass fraction
    maH = nH * muH / MMTOT
    maHe = nHe * muHe / MMTOT
    mC = nC * muC / MMTOT
    mO = nO * muO / MMTOT
    mMg = nMg * muMg / MMTOT
    mSi = nSi * muSi / MMTOT
    maFe = nFe * muFe / MMTOT
    
    # Compute planet composition fractions
    refractories = (mMg + mSi + maFe) * 100.
    silicates = (mMg + mSi) * 100.
    z = (mC + mO + mMg + mSi + maFe) * 100.
    iron_frac = maFe / (mMg + mSi + maFe) * 100.
    water_frac = mO / (mC + mO + mMg + mSi + maFe) * 100.
    
    return refractories, silicates, z, iron_frac, water_frac

def abundances_solarTOabsolute(c, o, mg, si, fe):
    """ Convert solar abundances to absolute abundances """
    C_abs = 10**(c + 8.46)  # Asplund et al. (2021)
    O_abs = 10**(o + 8.69)  # Asplund et al. (2021)
    Mg_abs = 10**(mg + 7.55)  # Asplund et al. (2021)
    Si_abs = 10**(si + 7.51)  # Asplund et al. (2021)
    Fe_abs = 10**(fe + 7.46)  # Asplund et al. (2021)
    return C_abs, O_abs, Mg_abs, Si_abs, Fe_abs

def bootstrap_composition(input_table_path, starname="all"):
    """
    Compute planet compositions for given input table.
    """
    input_sample = pd.read_csv(input_table_path)
    option = '%s' % starname
    N_bootstrap = 10000
    
    if option == 'all':
        output_table = pd.DataFrame(input_sample.star, columns=['star'])
        output_table['Ref'] = 0.0
        output_table['Silc'] = 0.0
        output_table['Z'] = 0.0
        output_table['f_iron'] = 0.0
        output_table['wf'] = 0.0
        
        for n, star in enumerate(input_sample.star):
            output_composition = np.zeros((N_bootstrap, 5))
            for i in range(N_bootstrap):
                rand_norm_c = float(np.random.normal(input_sample.C[n], input_sample.C_err[n], 1)[0])
                rand_norm_o = float(np.random.normal(input_sample.O[n], input_sample.O_err[n], 1)[0])
                rand_norm_mg = float(np.random.normal(input_sample.Mg[n], input_sample.Mg_err[n], 1)[0])
                rand_norm_si = float(np.random.normal(input_sample.Si[n], input_sample.Si_err[n], 1)[0])
                rand_norm_fe = float(np.random.normal(input_sample.feh[n], input_sample.feh_err[n], 1)[0])
                
                refractories, silicates, z, iron_frac, water_frac = planet_composition(rand_norm_c, rand_norm_o, rand_norm_mg, rand_norm_si, rand_norm_fe)
                output_composition[i] = [refractories, silicates, z, iron_frac, water_frac]
                
            output_table.loc[n, 'Ref'] = np.mean(output_composition[:, 0])
            output_table.loc[n, 'Silc'] = np.mean(output_composition[:, 1])
            output_table.loc[n, 'Z'] = np.mean(output_composition[:, 2])
            output_table.loc[n, 'f_iron'] = np.mean(output_composition[:, 3])
            output_table.loc[n, 'wf'] = np.mean(output_composition[:, 4])
        
        output_file = os.path.join(os.path.dirname(input_table_path), 'compositions_all.csv')
        output_table.to_csv(output_file, sep=',', index=False)
        print(f"Saved compositions to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute planet compositions from input table")
    parser.add_argument("input_table", type=str, help="Path to the input CSV file with abundances")
    parser.add_argument("--starname", type=str, default="all", help="Specific star name or 'all'")
    args = parser.parse_args()
    
    bootstrap_composition(args.input_table, args.starname)
