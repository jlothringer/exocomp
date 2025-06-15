#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  6 14:04:20 2025

@author: jlothringer
"""

import easychem.easychem as ec
from scipy.optimize import curve_fit,fmin
import numpy as np
import matplotlib.pyplot as plt
import random
#import emcee

#Tasks
#   1) Convert molecular measurements to bulk abundances [mostly done, make version with refractories]
#        e.g., take log10(H2O)_vmr and convert
#   2) Convert VMRs to MMRs and vice versa
#   3) Convert elemental ratios relative to solar to different definitions of solar or stellar [mostly done, add more solar defs]
#       a) Add stellar abund option
#   4) Given [M/H], C/O, elemental abund, P, and T, give expected molecular abundances
#       4.5) And vice-versa - given molecular abundances, P, T, and solar definition, give best-fitting [M/H] & C/O [done, make version with refractories]
#           -A bit better than just adding up molecules... because you may be missing some.
#           -So given H2O and CO2, you can get [M/H] and C/O!

#todo: 1) add refractory to volatile fits [done]
#2) add stellar [done]
#3) do some real-world examples [done]
#4) write up 
#5) make PEP compliant 
#6) look at obs and find pop means 
#7) make table publishable [and add uncertainty to mass and temp] 
#8) Make mcmc version 
#9) make posterior version [done]
#9) Make version where neither C or O are individually adjusted for C/O, but rather both! preserving metallicity! [done]
#10) add goodness of fit info... will want to choose between ways to handle C/O and whether adding ref/vol is justified [done]
#11) enhanced errorbars - upper/lower limits and asymm bars
#12) create test block/routine
#13) look at posteriors of W-178b
#14) The presence of gaseous H2O sets a limit of ref/vol for cloudy atmospheres
#   -lol if only ref-O enriched a planet atmosphere, we'd only see O in refractory phase (more the most part)

##Test
## M/H= 1.0 and C/O = 0.5 at T=1000 and P = 1
#{'H2O': -2.193862599311742, 'CO': -2.9951664932143895,'CO2': -4.957449831706032} for Asplund 09
#{'H2O': -2.196454926078874, 'CO': -2.95559683423187, 'CO2': -4.9118293333092335} for Lodders 25

#Depends on how you change C/O! Esp. when there are other elements setting M/H
#g = Abund(species_abunds={'H2O':-2.96,'CH4':-9.05,'CO2':-5.81,'K':-8.07,'Na':-9.31},species_errs={'H2O':0.31,'CH4':2.1,'CO2':1.31,'K':0.58,'Na':1.77}) #W17b
#g.convert_species_abunds(1271, 5e-3)
#versus 
#g.convert_species_abunds(1271, 5e-3,adjust_C=False)


class Abund:
    """
    Top-level class for abundance calculations.
    
    Parameters
    ----------
    
    species_abunds : dict 
        (Optional) Input/output for bulk_abundance calculations. Of form {species:VMR} 
        (e.g., {'H2O':1e-3})
        
    species_errs : dict 
        (Optional) Input/output for bulk_abundance calculations. Of form {species:VMR_err} for symmetric errors (e.g., {'H2O':1e-4}),  
        {species:[VMR_low,VMR_high]} for asymetric error bars, or {species:[VMR_0,VMR_1,...]} for posteriors.
        
    solar : str
        (Optional) Which set of solar abundances to use. Default is Lodders25. Possible options listed in
        self.possible_solars.
    
    retrieval : str
        (Optional) Which retrieval code was used to produce metallicity and C/O constraints. This is used to
        choose the correct conversion calculation. Defines self.mh_type, self.co_type, and self.solar.
        Not all retrieval codes are included and you may need a different set up, so you can just define
        self.mh_type, self.co_type, and self.solar yourself.
        
        
    Returns
    -------
    
    bulk_abunds : dict 
        Where output of fits to species abundances or conversion to different solar compositions is placed.
        Will be of the format {'O':8.0} for an O/H ratio of 1e8 parts per trillion (i.e., H = 1e12).
        
    bulk_errs : dict 
        (Optional) Where output of fits to species abundances or conversion to different solar compositions is placed. 
        Of form {species:bulk_errs} for symmetric errors (e.g., {'O/H':0.5}),  
        {species:[abund_low,abund_high]} for asymetric error bars, or {species:[abund_0,abund_1,...]} for posteriors.
     
    
    
    Attributes
    ----------
    possible_solars() : list possible solar compositions
    define_stellar() : define host star composition to use instead of solar
    VMR_to_MMR() : convert a dict of volume mixing ratios to mass fraction
    MMR_to_VMR() : convert a dict of mass fractions to volume mixing ratios
    convert_bulk_abundance() : given bulk abundance properties, calculate elemental abundances
    convert_solar()  : convert self.bulk_abundances to another definition of solar
    convert_species_abunds() : given species abundances, calculate bulk abundances
    predict_abund() : given bulk abundance properties, calculate species abundances
    init_ec() : utility function to initialize easyCHEM with proper solar abundances
    import_masses() : read in molecular masses for MMW calculation 
        (e.g., conversion between VMR and MMR)
    
    
    species_type : VMR or MMR
    solar_abundances : elemental abundances for given self.solar
    solar_abundances_orig : solar abundances first defined in class definition
    possible_solars : list of possible solar compositions
    best_fit_mh : best-fit metallicity from convert_species_abunds
    best_fit_co : best-fit C/O from convert_species_abunds
    best_fit_refvol : best-fit ref/vol from convert_species_abunds
    best_fit_abund : best-fit species abundances from convert_species_abunds
    Asplund09 : elemental abundance definitions from Asplund, M., et al. (2009). ARA&A, 47, 481
    Asplund21 : elemental abundance definitions from Asplund, M., et al. (2021). A&A, 653, A141  
    Asplund05 : elemental abundance definitions from 
    Lodders10 : elemental abundance definitions from 
    Lodders20 : elemental abundance definitions from Lodders, K. (2020). Space Sci. Rev., 216, 44
    Caffau11 : elemental abundance definitions from 
    
    Examples
    --------
    Setting up for a bulk abundance calculation with symmetric errors
    g = Abund(species_abunds={'H2O':-2.96,'CH4':-9.05,'CO2':-5.81,'K':-8.07,'Na':-9.31},
              species_errs={'H2O':0.31,'CH4':2.1,'CO2':1.31,'K':0.58,'Na':1.77})
    
    Setting up for a bulk abundance calculation with posteriors
    g = Abund(species_abunds={'H2O':-2.96,'CH4':-9.05,'CO2':-5.81,'K':-8.07,'Na':-9.31},
              species_errs={'H2O':np.random.normal(-2.96,0.31,100),'CH4':np.random.normal(-9.05,2.1,100),
                            'CO2':np.random.normal(-5.81,1.31,100),
                            'K':np.random.normal(-8.07,0.58,100),'Na':np.random.normal(-9.31,1.77,100)})
    """
    
    def __init__(self, bulk_abunds=None, species_abunds=None, 
                 bulk_errs=None, species_errs=None, 
                 species_type='VMR', 
                 solar='Lodders25',retrieval=None):
        
        self.__version__ = "0.1.0"

        self.bulk_abunds = bulk_abunds
        self.species_abunds = species_abunds
        self.bulk_errs = bulk_errs
        self.species_errs = species_errs
        self.species_type = species_type
        self.solar = solar
        self.retrieval = retrieval
        self.solar_abundances = None
        
        self.possible_solars = ['Lodders25','Asplund09',
                                'Lodders20','Asplund21',
                                'Lodders10','Caffau11',
                                'Asplund05']
        
        match retrieval:
            case 'pRT': #documentation
                self.mh_type = 'C/H'
                self.co_type = 'O/H'
                self.solar = 'Asplund09'
            case 'PICASO': #Mukherjee et al. 2023
                self.mh_type = '(O+C)/H'
                self.co_type = 'MH_Preserve' 
                self.solar = 'Lodders10'
            case 'PETRA': #Lothringer & Barman 2020
                self.mh_type = 'O/H'
                self.co_type = 'C/H'
                self.solar = 'Asplund05'
            case 'POSEIDON': #Kirk et al. 2025
                self.mh_type = 'O/H'
                self.co_type = 'C/H'
                self.solar = 'Asplund09'
                print('Sometimes POSEIDON varies O/H for C/O (Meech et al. 2025)')
            case 'CHIMERA': #assuming same as sCHIMERA
                self.mh_type = '(O+C)/H'
                self.co_type = 'MH_Preserve' 
                self.solar = 'Asplund09'
            case 'sCHIMERA': #Mansfield et al. 2024 
                self.mh_type = '(O+C)/H'
                self.co_type = 'MH_Preserve' 
                self.solar = 'Asplund09'
            case 'SCARLET':  #Benneke 2015, now Pelletier et al. 2025 -> fastchem
                self.mh_type = 'O/H' 
                self.co_type = 'C/H' 
                self.solar = 'Asplund09' #unknown      
            case 'ATMO':
                self.mh_type = '(O+C)/H' 
                self.co_type = 'MH_Preserve'              
                self.solar = 'Asplund09'
            case 'PLATON': ##Kirk et al. 2025
                #uses ggchem but hard to tell how
                #built on fastchem which uses asplund09
                self.mh_type = 'O/H'
                self.co_type = 'C/H'            
                self.solar = 'Asplund09'
            case 'Gibson': #Ramkumar et al. 2025
                self.mh_type = '(O+C)/H' 
                self.co_type = 'MH_Preserve'     
                self.solar = 'Asplund09'   
            #case 'Aurora': #Welbanks & Madhusudhan 2021
                #self.mh_type = 'O/H'
                #self.co_type = 'C/H'
                #self.solar = 'Asplund09'
           # case 'NEMESIS':
            
           # case 'Hydra':
                
                   
        match solar:
            case 'Lodders25':
                self.solar_abundances = self.Lodders25()
                self.solar_abundances_orig = self.solar_abundances
            case 'Asplund09':
                self.solar_abundances = self.Asplund09()
                self.solar_abundances_orig = self.solar_abundances
            case 'Asplund05':
                self.solar_abundances = self.Asplund09()
                self.solar_abundances_orig = self.solar_abundances
            case 'Asplund21':
                self.solar_abundances = self.Asplund21()
                self.solar_abundances_orig = self.solar_abundances
            case 'Lodders20':
                self.solar_abundances = self.Lodders20()
                self.solar_abundances_orig = self.solar_abundances
            case 'Lodders10':
                self.solar_abundances = self.Lodders10()
                self.solar_abundances_orig = self.solar_abundances
            case 'Caffau11':
                print('Using abundances from Caffau11, filling in \
                      elements not included with Asplund09')
                self.solar_abundances = self.Caffau11()
                self.solar_abundances_orig = self.solar_abundances
            case _:
                print('Please define Abund.solar_abundances')
                
    def define_stellar(self, abund_dict):
        #will use self.solar as backup

        for key in abund_dict.keys():
            solar_abundance[key] = abund_dict[key]
        return
    
    def VMR_to_MMR(self, vmr_dict, mmw=2.3, replace_abunds=False):
        self.import_masses()
        mmr_dict = {}
        for key in vmr_dict.keys():
            if key not in self.molecular_mass.keys():
                print(f'{key} not in molecular mass dictionary. Skipping. Please add.')
                continue
            mmr_dict[key] = np.log10(10**vmr_dict[key] * self.molecular_mass[key] / mmw)
        if replace_abunds:
            self.species_abunds = mmr_dict
        return mmr_dict
    
    def MMR_to_VMR(self,mmr_dict, mmw=2.3, replace_abunds=False):
        self.import_masses()
        vmr_dict = {}
        for key in mmr_dict.keys():
            if key not in self.molecular_mass.keys():
                print(f'{key} not in molecular mass dictionary. Skipping. Please add.')
                continue
            vmr_dict[key] = np.log10(10**mmr_dict[key] * mmw / self.molecular_mass[key])
        if replace_abunds:
            self.species_abunds = vmr_dict
        return vmr_dict
    
    def convert_bulk_abundance(self,mh,co,rv=0.0,mh_err=0.0,co_err=0.0,rv_err=0.0,solar = None):
        """Convert bulk abundance parameters to elemental abundances with error propagation.
        
        This method transforms fundamental stellar/planetary parameters (metallicity, C/O ratio, 
        and refractory enhancement) into individual elemental abundances relative to hydrogen.
        It supports both symmetric/asymmetric error propagation and posterior sampling for 
        uncertainty quantification.
        
        
        Parameters
        ----------
        mh : float or array-like
            Metallicity [M/H] in dex. Can be a single value or array for posterior sampling.
        co : float or array-like
            Carbon-to-oxygen ratio. Interpretation depends on `co_type` attribute:
            
            * 'O/H': C/O ratio, oxygen abundance adjusted
            * 'C/H': C/O ratio, carbon abundance adjusted  
            * 'MH_Preserve': C/O ratio, preserves total C+O abundance
            
        rv : float or array-like, optional
            Refractory enhancement factor in dex for elements heavier than oxygen.
            Default is 0.0 (no enhancement).
        mh_err : float, list, or array-like, optional
            Uncertainty in metallicity. Can be:
            
            * Single float: symmetric error
            * 2-element list: [lower, upper] asymmetric errors
            * Array: posterior samples
            
            Default is 0.0.
        co_err : float, list, or array-like, optional
            Uncertainty in C/O ratio, same format as `mh_err`. Default is 0.0.
        rv_err : float, list, or array-like, optional
            Uncertainty in refractory enhancement, same format as `mh_err`. Default is 0.0.
        solar : str, optional
            Solar abundance reference to use. If provided, overrides current solar abundance.
            Supported values: 'Lodders25', 'Asplund09', 'Asplund05', 'Asplund21', 
            'Lodders20', 'Lodders10', 'Caffau11'. Default uses current `self.solar`.
        
        Returns
        -------
        bulk_dict : dict
            Dictionary of elemental abundances in the format log10(X/H) + 12, where
            keys are element symbols and values are abundances.
        err_dict : dict
            Dictionary of abundance uncertainties with same keys as `bulk_dict`.
            
            * For symmetric/asymmetric errors: contains propagated uncertainties
            * For posterior sampling: contains arrays of posterior abundance samples
        
        
        Notes
        -----
        The method requires the following class attributes to be set:
        
        * `mh_type`: Metallicity parameterization type
        * `co_type`: C/O ratio treatment ('O/H', 'C/H', or 'MH_Preserve')
        * `solar`: Solar abundance reference (if `solar` parameter not provided)
        
        The method automatically detects whether inputs represent:
        
        1. **Single values with errors**: Standard error propagation using Gaussian approximation
        2. **Posterior samples**: Full posterior propagation through the abundance calculation
        
        For elements heavier than oxygen (Z > 8), the refractory enhancement `rv` is applied 
        as an additional multiplicative factor in dex space.
        
        Examples
        --------
        >>> # Single value with symmetric errors
        >>> bulk_dict, err_dict = obj.convert_bulk_abundance(
        ...     mh=0.1, co=0.55, rv=0.2, mh_err=0.05, co_err=0.03, rv_err=0.1
        ... )
        
        >>> # Asymmetric errors
        >>> bulk_dict, err_dict = obj.convert_bulk_abundance(
        ...     mh=0.1, co=0.55, mh_err=[-0.03, 0.07], co_err=[-0.02, 0.04]
        ... )
        
        >>> # Posterior sampling
        >>> import numpy as np
        >>> mh_samples = np.random.normal(0.1, 0.05, 1000)
        >>> co_samples = np.random.normal(0.55, 0.03, 1000)
        >>> bulk_dict, err_dict = obj.convert_bulk_abundance(
        ...     mh=0.1, co=0.55, mh_err=mh_samples, co_err=co_samples
        ... )
        
        See Also
        --------
        init_ec : Initialize elemental composition object
        
        References
        ----------
        Solar abundance compilations used depend on the `solar` parameter choice.
        """
        
        if self.mh_type is None:
            print('Please define retrieval-type or mh_type, co_type, and solar_type')
            return
        match self.co_type:
            case 'O/H':
                adjust_C = False
            case 'C/H':
                adjust_C = True
            case 'MH_Preserve':
                adjust_C = 'Neither'
        
        # Re-define solar abundances       
        if solar is not None:            
            self.solar = solar
        match self.solar:
            case 'Lodders25':
                self.solar_abundances = self.Lodders25()
                self.solar_abundances_orig = self.solar_abundances
            case 'Asplund09':
                self.solar_abundances = self.Asplund09()
                self.solar_abundances_orig = self.solar_abundances
            case 'Asplund05':
                self.solar_abundances = self.Asplund09()
                self.solar_abundances_orig = self.solar_abundances
            case 'Asplund21':
                self.solar_abundances = self.Asplund21()
                self.solar_abundances_orig = self.solar_abundances
            case 'Lodders20':
                self.solar_abundances = self.Lodders20()
                self.solar_abundances_orig = self.solar_abundances
            case 'Lodders10':
                self.solar_abundances = self.Lodders10()
                self.solar_abundances_orig = self.solar_abundances
            case 'Caffau11':
                print('Using abundances from Caffau11, filling in \
                      elements not included with Asplund09')
                self.solar_abundances = self.Caffau11()
                self.solar_abundances_orig = self.solar_abundances
            case _:
                print('Please define Abund.solar_abundances')
        print(f"Using {self.solar} solar abundances")
        
        exo = self.init_ec()
        
        def get_comp(mh,co,rv=0.0):
        
            exo.metallicity = mh
            #Adjust C/O either by changing C/H or O/H:
            if adjust_C is True:
                exo.atomAbunds[2] = exo.atomAbunds[4] * co
            elif adjust_C == 'Neither':
                tot = exo.atomAbunds[2] + exo.atomAbunds[4]
                #co = exo.atomAbunds[2] / exo.atomAbunds[4]
                exo.atomAbunds[2] = co * tot / (co+1)
                exo.atomAbunds[4] = tot / (co+1)
            else:
                exo.co = co #easychem changes O/H by default
                
            #Enhance elements above O by rv
            for i in range(5,len(exo.atoms)):
                exo.atomAbunds[i] = 10**(np.log10(exo.atomAbunds[i])+rv)
            
            return
       
        get_comp(mh,co,rv)
        
        bulk_dict = {}
        err_dict = {}
        for i,key in enumerate(exo.atoms):
            bulk_dict[key] = np.log10((exo.atomAbunds[i] / exo.atomAbunds[0])*1e12)
            err_dict[key] = mh_err
            
        try:
            s = len(mh_err)
            if s == 2: # asymmetric errors
                s = None # don't treat as posterior
        except TypeError:
            s = None

        if s is None: # If we have errors (symmetric or asymmetric)     
            if self.co_type == 'O/H':
                #Then error in O is set also by error in C/O
                try:
                    tmp = len(mh_err)
                    err_dict['O'] = list(np.sqrt(np.array(mh_err)**2 + np.array(co_err)**2))
                except:
                    err_dict['O'] = np.sqrt(mh_err**2+co_err**2)
            if self.co_type == 'C/H':
                #Error in C is set also by error in C/O
                try:
                    tmp = len(mh_err)
                    err_dict['C'] = list(np.sqrt(np.array(mh_err)**2 + np.array(co_err)**2))
                except:
                    err_dict['C'] = np.sqrt(mh_err**2+co_err**2)
            if self.co_type == 'Neither':
                #Unsure here, conservatively, let's place C/O error on both...
                #I think one could split the C/O error but need to think about that
                try:
                    tmp = len(mh_err)
                    err_dict['C'] = list(np.sqrt(np.array(mh_err)**2 + np.array(co_err)**2))
                    err_dict['O'] = list(np.sqrt(np.array(mh_err)**2 + np.array(co_err)**2))
                except:
                    err_dict['C'] = np.sqrt(mh_err**2+co_err**2)
                    err_dict['O'] = np.sqrt(mh_err**2+co_err**2)
            if np.array(rv_err == 0.0).all():
                #Again, will be correlated usually, but without posteriors, assume Gauss
                for i in range(5,len(exo.atoms)):
                    try:
                        tmp = len(mh_err)
                        err_dict[exo.atoms[i]] = list(np.sqrt(np.array(mh_err)**2 + np.array(rv_err)**2))
                    except:
                        err_dict[exo.atoms[i]] = np.sqrt(mh_err**2+rv_err**2)
    
            self.bulk_abunds = bulk_dict
            self.bulk_errs = err_dict
        else: # If we have posteriors
            err_dict = {}
            for j in exo.atoms:
                err_dict[j] = []
            for i in range(len(mh_err)):
                try:
                    get_comp(mh_err[i],co_err[i],rv_err[i])
                except TypeError:
                    get_comp(mh_err[i],co_err[i],0.0)
                for j in range(len(exo.atoms)):
                    err_dict[exo.atoms[j]].append(np.log10((exo.atomAbunds[j] / exo.atomAbunds[0])*1e12))
            self.bulk_errs = err_dict
        return bulk_dict,err_dict

    def convert_solar(self, convert_to):
        """Convert existing bulk abundances to a different solar abundance reference.

        This method transforms previously calculated elemental abundances from one solar 
        abundance scale to another while preserving the relative abundances. The conversion
        accounts for differences between solar abundance compilations by applying appropriate
        offset corrections to maintain consistency in the abundance ratios.
        
        Parameters
        ----------
        convert_to : str
            Target solar abundance reference scale. Supported values:
            
            * 'Lodders25': Lodders (2025) solar photospheric abundances
            * 'Asplund09': Asplund et al. (2009) solar photospheric abundances  
            * 'Asplund21': Asplund et al. (2021) solar photospheric abundances
            * 'Lodders20': Lodders (2020) solar photospheric abundances
            * 'Lodders10': Lodders (2010) solar photospheric abundances
        
        Returns
        -------
        None
            Method modifies class attributes (self.bulk_abunds and sel.fsolar_abundances) 
            in-place and prints status messages.
    
        
        Notes
        -----
        This method modifies the following class attributes:
        
        * `solar_abundances`: Updated to the new solar abundance reference
        * `bulk_abunds`: Element abundances converted to the new solar scale
        * `solar_abundances_orig`: Changes at the end
        
        The conversion preserves the physical meaning of abundance ratios by applying
        the transformation:
        
        .. math::
            
            [X/H]_{new} = [X/H]_{old} - ([X/H]_{solar,new} - [X/H]_{solar,old})
        
        where the solar abundance difference accounts for systematic offsets between
        different solar abundance compilations.
        
        **Behavior when no bulk abundances exist:**
        
        If `bulk_abunds` is None (no previous abundance calculation), the method only
        updates the solar abundance reference without performing any conversions.
        
        **Workflow Integration:**
        
        This method is typically used after `convert_bulk_abundance()` when you want
        to compare results using different solar abundance scales or when literature
        values use a different solar reference than your initial calculation.
        
        Examples
        --------
        >>> # After calculating abundances with one solar reference
        >>> obj.convert_bulk_abundance(mh=0.1, co=0.55, solar='Asplund09')
        >>> 
        >>> # Convert to a different solar abundance scale
        >>> obj.convert_solar('Lodders25')
        Changed solar abunds to Lodders25
        
        >>> # If no bulk abundances calculated yet
        >>> obj.convert_solar('Asplund21') 
        Changed solar abunds to Asplund21
        No bulk abundances to convert...
        
        >>> # Invalid solar reference
        >>> obj.convert_solar('InvalidRef')
        Input solar not recognized
        
        See Also
        --------
        convert_bulk_abundance : Calculate initial elemental abundances
        
        References
        ----------
        Solar abundance scales and their systematic differences are discussed in:
        
        * Asplund, M., et al. (2009). ARA&A, 47, 481
        * Asplund, M., et al. (2021). A&A, 653, A141  
        * Lodders, K. (2010). Astrophys. Space Sci., 334, 1
        * Lodders, K. (2020). Space Sci. Rev., 216, 44
        """
        match convert_to:
            case 'Lodders25':
                solar_new = self.Lodders25()
                print('Changed solar abunds to Lodders25')
            case 'Asplund09':   
                solar_new = self.Asplund09()
                print('Changed solar abunds to Asplund09')
            case 'Asplund21':
                solar_new = self.Asplund21()
                print('Changed solar abunds to Asplund21')
            case 'Lodders20':
                solar_new = self.Lodders20()
                print('Changed solar abunds to Lodders20')
            case 'Lodders10':
                solar_new = self.Lodders10()
                print('Changed solar abunds to Lodders10')
            case _:
                print('Input solar not recognized')
                
        if self.bulk_abunds == None:
            self.solar_abundances = solar_new
            self.solar_abundances_orig = solar_new
            print('No bulk abundances to convert...')
            return
        
        for element in self.bulk_abunds:
            #Conversion factor
            factor = solar_new[element][0] - self.solar_abundances_orig[element][0] 
            #print(factor,solar_new[element][0],self.solar_abundances_orig[element][0])
            
            #Convert
            self.bulk_abunds[element] = self.bulk_abunds[element] - factor #minus b/c increase/decrease in planet *relative to solar/star*
            
        self.solar_abundances = solar_new
        self.solar_abundances_orig = solar_new
        return
    

    def init_ec(self):
        '''
        Function to initialize easyCHEM with the elemental abundances from the
        present definition of solar.
        
        Returns
        -------
        exo : obj
            easyCHEM object
        
        '''
        #Init easychem
        exo = ec.ExoAtmos()
        
        reactants = exo.reactants.copy()
        reactants = np.append(reactants,'Fe+')
        #reactants = np.append(reactants,'Mg+')
        #reactants = np.append(reactants,'C2H2')    
        #reactants = np.append(reactants,'C2H4')                                                                                                                                                                                                        
        exo.updateReactants(reactants)
        
        #Update easychem abunds
        d = [10**self.solar_abundances[element][0] for element in exo.atoms]
        for i,element in enumerate(exo.atoms):
            exo.atomAbunds[i] = d[i]/np.sum(d)
        exo.updateAtomAbunds(exo.atomAbunds)
        return exo
    
    
    def predict_abund(self,mh,co, T, P, species=['H2O','CO','CO2'],verbose=True,
                      mode='VMR', adjust_C=True):
        """Predict molecular abundances from bulk parameters and atmospheric conditions.
        
        Calculates equilibrium molecular abundances for specified species given 
        metallicity, C/O ratio, temperature, and pressure.
        Uses the definition of solar in self.solar!
        
        Parameters
        ----------
        mh : float
            Metallicity [M/H] in dex.
        co : float
            Carbon-to-oxygen ratio.
        T : float
            Temperature in Kelvin.
        P : float
            Pressure in bar.
        species : list of str, optional
            Molecular species to calculate. Default is ['H2O','CO','CO2'].
        verbose : bool, optional
            Print abundances to console. Default is True.
        mode : {'VMR', 'MMR'}, optional
            Return volume mixing ratios ('VMR') or mass mixing ratios ('MMR'). 
            Default is 'VMR'.
        adjust_C : bool, optional
            If True, adjust carbon abundance to match C/O ratio. If False, 
            adjust oxygen abundance. Default is True.
        
        Returns
        -------
        dict
            Dictionary with species names as keys and log10 abundances as values.
        
        Examples
        --------
        >>> # Predict H2O, CO, CO2 abundances
        >>> abunds = obj.predict_abund(mh=0.0, co=0.55, T=1500, P=1.0)
        H2O: -3.2
        CO: -3.8
        CO2: -7.1
        
        >>> # Custom species list, silent mode
        >>> abunds = obj.predict_abund(mh=0.5, co=0.8, T=2000, P=10.0, 
        ...                          species=['CH4', 'NH3'], verbose=False)
        """
        
        exo = self.init_ec()
        
        exo.metallicity = mh
        #Adjust C/O either by changing C/H or O/H:
        if adjust_C is True:
            exo.atomAbunds[2] = exo.atomAbunds[4] * co
        elif adjust_C == 'Neither':
            tot = exo.atomAbunds[2] + exo.atomAbunds[4]
            #co = exo.atomAbunds[2] / exo.atomAbunds[4]
            exo.atomAbunds[2] = co * tot / (co+1)
            exo.atomAbunds[4] = tot / (co+1)
        else:
            exo.co = co #easychem changes O/H by default
        exo.solve(P,T)
        
        if mode == 'VMR':
            result = exo.result_mol()
        elif mode == 'MMR':
            result = exo.result_mass()
        ret_dict = {}
        for s in species:
            if verbose:
                print(f"{s}: {np.log10(result[s])}")
            ret_dict[s] = np.log10(result[s])
        return ret_dict
    
    def convert_species_abunds(self, T, P, adjust_C=True, plot_it=True,
                               posterior_samples=100, fit_refvol = False):
        """Fit bulk atmospheric parameters to observed molecular abundances.
        
        Performs inverse retrieval to determine metallicity [M/H], C/O ratio, and 
        optionally refractory enhancement from observed molecular species abundances
        using thermochemical equilibrium modeling. Supports both error propagation
        and full posterior sampling.
        
        IMPORTANT: This assumes chemical equilibrium, of course! Results can be very
        sensitive to the assumed T and P. And atmospheres are not necessarily in 
        equilibrium. Interpret these results with care.
        
        Parameters
        ----------
        T : float
            Atmospheric temperature in Kelvin for which the species abundances are representative.
        P : float
            Atmospheric pressure in bar for which the species abundances are representative.
        adjust_C : bool, optional
            Carbon adjustment method for C/O ratio. If True, adjust carbon abundance;
            if False, adjust oxygen abundance. If string 'Neither' conserve [M/H].
            Default is True.
        plot_it : bool, optional
            Generate diagnostic plots showing fit results and parameter distributions.
            Default is True.
        posterior_samples : int, optional
            Number of posterior samples to use when `species_errs` contains posteriors.
            Default is 100.
        fit_refvol : bool, optional
            Include refractory enhancement (R/V) as a free parameter in the fit.
            Default is False.
        
        Returns
        -------
        popt : list
            Best-fit parameters [M/H, C/O] or [M/H, C/O, R/V] if `fit_refvol=True`.
        pcov : array-like or list
            Parameter covariance matrix (for error-based fitting) or list of 
            posterior samples (for posterior-based fitting).
        best_fit : array
            Model predictions for molecular abundances at best-fit parameters.
        
        Notes
        -----
        This method requires the following class attributes to be set:
        
        * `species_abunds`: Dictionary of observed molecular abundances
        * `species_errs`: Uncertainties (symmetric, asymmetric, or posteriors)
        
        **Fitting Methods:**
        
        1. **Symmetric errors**: Uses `scipy.optimize.curve_fit` with Gaussian likelihood
        2. **Asymmetric errors**: Uses `scipy.optimize.fmin` with custom chi-squared
        3. **Posterior samples**: Fits each sample individually, returns percentiles
        
        **Parameter Bounds:**
        
        * Metallicity: -2.0 < [M/H] < 3.0
        * C/O ratio: 0.0 < C/O < 2.0
        * Model returns -1000 for abundances outside these bounds
        
        The method updates several class attributes:
        
        * `best_fit_mh`, `best_fit_co`: Best-fit bulk parameters
        * `best_fit_refvol`: Best-fit refractory enhancement (if fitted)
        * `best_fit_abund`: Dictionary of best-fit molecular abundances  
        * `bulk_abunds`: Updated elemental abundances from best-fit parameters
        * `fits`: Posterior samples (for posterior-based fitting)
        
        Examples
        --------
        >>> # Basic fitting with symmetric errors
        >>> obj.species_abunds = {'H2O': -3.2, 'CO': -3.8, 'CO2': -7.1}
        >>> obj.species_errs = {'H2O': 0.1, 'CO': 0.15, 'CO2': 0.3}
        >>> popt, pcov, fit = obj.convert_species_abunds(T=1500, P=1.0)
        
        >>> # Include refractory enhancement
        >>> popt, pcov, fit = obj.convert_species_abunds(T=1500, P=1.0, fit_refvol=True)
        
        >>> # Posterior sampling (when species_errs contains arrays)
        >>> popt, pcov, fit = obj.convert_species_abunds(T=1500, P=1.0, 
        ...                                             posterior_samples=500)
        
        Raises
        ------
        AttributeError
            If `species_abunds` is not defined in the class instance.
        ValueError
            If optimization fails to converge or input data is inconsistent.
        
        See Also
        --------
        predict_abund : Forward modeling of molecular abundances
        convert_bulk_abundance : Convert bulk parameters to elemental abundances
        """
        
        if self.species_errs == None:
            print("Consider adding errors/uncertainty to the provided measurements!")
        
        #Init easychem
        exo = self.init_ec()
        
        def fit_species(x,mh,co,rv=0):
            #Put some bounds in so easychem doesn't flip out
            if (mh < -2.0) or (mh > 3.0):
                return np.array([-1000.0 for el in self.species_abunds.keys()])
            if (co < 0) or (co > 2.0):
                return np.array([-1000.0 for el in self.species_abunds.keys()])
            exo.metallicity = mh #resets atomAbunds
            #Adjust C/O either by changing C/H or O/H:
            if adjust_C is True:
                exo.atomAbunds[2] = exo.atomAbunds[4] * co
            elif adjust_C == 'Neither':
                tot = exo.atomAbunds[2] + exo.atomAbunds[4]
                #co = exo.atomAbunds[2] / exo.atomAbunds[4]
                exo.atomAbunds[2] = co * tot / (co+1)
                exo.atomAbunds[4] = tot / (co+1)
            else:
                exo.co = co #easychem changes O/H by default
            #Enhance elements above O by rv
            for i in range(5,len(exo.atoms)):
                exo.atomAbunds[i] = 10**(np.log10(exo.atomAbunds[i])+rv)
            exo.solve(P,T)
            result = exo.result_mol()
            return np.array([np.log10(result[el]) for el in self.species_abunds.keys()])
        
        def fit_species_loss(params, x, y, errs):

            mh = params[0]
            co = params[1]
            if len(params) == 3:
                rv = params[2]
            else:
                rv = 0.0
            results = fit_species(x,mh,co,rv)
            
            chi2 = 0
            for i, key in enumerate(errs.keys()):
                # Asymmetric error handling
                if results[i] > y[key]:
                    chi2 += ((results[i] - y[key]) / errs[key][1]) ** 2
                elif errs[key][0] == 0.0:
                    #if only upper limit, don't penalize if below...
                    continue
                else:
                    chi2 += ((results[i] -y[key]) / errs[key][0]) ** 2
                    
            return chi2
        
        # See if we have errors or posteriors
        # Dumb logic finds if we have 2 elements in errors, which means we were given asymmetric error bars
        try:
            s = len(self.species_errs[sorted(self.species_errs)[0]])
            if s == 2:
                s = None
        except TypeError:
            s = None
                
        # If just given errors
        if s is None:
            #
            #if asymmetric errorbars given, use fmin
            if all(isinstance(value, list) for value in self.species_errs.values()):
                if fit_refvol:
                    popt = fmin(fit_species_loss, [0.5,0.5,-0.1], args = (None,self.species_abunds,self.species_errs))
                    best_fit = fit_species(range(len(self.species_abunds)),popt[0],popt[1],popt[2])       
                    errs = np.zeros(len(best_fit))
                    print('---')
                    print(f'Best fit [M/H]: {popt[0]:.3f} ')
                    print(f'Best fit C/O: {popt[1]:.3f}')
                    print(f'Best fit R/V: {popt[2]:.3f}')
                else:
                    popt = fmin(fit_species_loss, [0.5,0.5], args = (None,self.species_abunds,self.species_errs))
                    best_fit = fit_species(range(len(self.species_abunds)),popt[0],popt[1])       
                    errs = np.zeros(len(best_fit))
                    print('---')
                    print(f'Best fit [M/H]: {popt[0]:.3f} ')
                    print(f'Best fit C/O: {popt[1]:.3f}')
                chi = fit_species_loss(best_fit,None,self.species_abunds,self.species_errs)
                print(f'Chi^2 = {chi}, Reduced Chi^2 = {chi/(len(best_fit)-1)}')
                #set sigma and pcov to None for below since they are not applicable
                sigma = None
                pcov = None

            # If errors are not lists, we must have been given symmetric uncertainties
            else:                
                if self.species_errs == None:
                    sigma = None
                else:
                    sigma = [vals for key, vals in self.species_errs.items()]
                if fit_refvol:
                    popt,pcov = curve_fit(fit_species,range(len(self.species_abunds)),
                              [vals for key, vals in self.species_abunds.items()],
                              sigma=[vals for key, vals in self.species_errs.items()],
                              absolute_sigma=True,
                              p0=[0.5,0.5,-0.1])
                    best_fit = fit_species(range(len(self.species_abunds)),popt[0],popt[1],popt[2])
                    errs = np.diag(pcov)**2
                    
     
                    print('---')
                    print(f'Best fit [M/H]: {popt[0]:.3f}  +/- {errs[0]:.3f}')
                    print(f'Best fit C/O: {popt[1]:.3f}  +/- {errs[1]:.3f}')
                    print(f'Best fit R/V: {popt[2]:.3f}  +/- {errs[2]:.3f}')
                else:
                    # same thing but without refvol
                    popt,pcov = curve_fit(fit_species,range(len(self.species_abunds)),
                              [vals for key, vals in self.species_abunds.items()],
                              sigma=sigma,
                              absolute_sigma=True,
                              p0=[0.5,0.5])
                    best_fit = fit_species(range(len(self.species_abunds)),popt[0],popt[1])
                    errs = np.diag(pcov)**2
                    
                    print('---')
                    print(f'Best fit [M/H]: {popt[0]:.3f}  +/- {errs[0]:.3f}')
                    print(f'Best fit C/O: {popt[1]:.3f}  +/- {errs[1]:.3f}')

            
            print('---')
            print('Best Fit Abundances')
            print('---')
            for i, key in enumerate(self.species_abunds.keys()):
                print(f'{key}: {best_fit[i]}')
                
            if sigma is not None:
                chi = 0.0
                for i,key in enumerate(self.species_abunds.keys()):
                    chi += (best_fit[i] - self.species_abunds[key])**2 / sigma[i]**2
                print(f'Chi^2 = {chi}, Reduced Chi^2 = {chi/(len(best_fit)-1)}')
            
            if plot_it:
                fig,ax = plt.subplots()
                if self.species_errs is None:
                    ax.plot(self.species_abunds.keys(),self.species_abunds.values(),'ko',label='Data')
                elif all(isinstance(value, list) for value in self.species_errs.values()):
                    lower_err = [self.species_errs[key][0] for key in self.species_errs.keys()]
                    upper_err = [self.species_errs[key][1] for key in self.species_errs.keys()]
                    asymmetric_errors = [lower_err, upper_err]
                    ax.errorbar(self.species_abunds.keys(),self.species_abunds.values(),yerr=asymmetric_errors,fmt='o',color='k',label='Data')
                else:
                    ax.errorbar(self.species_abunds.keys(),self.species_abunds.values(),yerr=list(self.species_errs.values()),fmt='o',color='k',label='Data')
                    
                
                if fit_refvol:
                    highco = fit_species(range(len(self.species_abunds)),popt[0],popt[1]+errs[1],popt[2])
                    lowco = fit_species(range(len(self.species_abunds)),popt[0],popt[1]-errs[1],popt[2])
                    highmh = fit_species(range(len(self.species_abunds)),popt[0]+errs[0],popt[1],popt[2])
                    lowmh = fit_species(range(len(self.species_abunds)),popt[0]-errs[0],popt[1],popt[2])
                    highrv = fit_species(range(len(self.species_abunds)),popt[0],popt[1],popt[2]+errs[2])
                    lowrv = fit_species(range(len(self.species_abunds)),popt[0],popt[1],popt[2]-errs[2])
                    
                   #print((best_fit-np.min([lowco,highco,lowmh,highmh,lowrv,highrv],axis=0),
                   #         np.max([lowco,highco,lowmh,highmh,lowrv,highrv],axis=0)-best_fit))
                    try:
                        ax.errorbar(self.species_abunds.keys(),best_fit,fmt='X',color='b',
                                    yerr = (best_fit-np.min([lowco,highco,lowmh,highmh,lowrv,highrv],axis=0),
                                            np.max([lowco,highco,lowmh,highmh,lowrv,highrv],axis=0)-best_fit),
                                    label=f'Best-Fit +/- 1-sig Chem. Eq. \n([M/H]={popt[0]:.2f}, C/O={popt[1]:.2f}, R/V={popt[2]:.2f})\n at {T}K and {P} bars')        
                    except ValueError:
                        print('Unable to plot errorbars from M/H, C/O, R/V fit')
                        ax.plot(self.species_abunds.keys(),best_fit,'X',color='b',
                                    label=f'Best-Fit Chem. Eq. \n([M/H]={popt[0]:.2f}, C/O={popt[1]:.2f}, R/V={popt[2]:.2f})\n at {T}K and {P} bars')     
          
                else:
                    highco = fit_species(range(len(self.species_abunds)),popt[0],popt[1]+errs[1])
                    lowco = fit_species(range(len(self.species_abunds)),popt[0],popt[1]-errs[1])
                    highmh = fit_species(range(len(self.species_abunds)),popt[0]+errs[0],popt[1])
                    lowmh = fit_species(range(len(self.species_abunds)),popt[0]-errs[0],popt[1])

                    try:
                        ax.errorbar(self.species_abunds.keys(),best_fit,fmt='X',color='b',
                                    yerr = (best_fit-np.min([lowco,highco,lowmh,highmh],axis=0),np.max([lowco,highco,lowmh,highmh],axis=0)-best_fit),
                                    label=f'Best-Fit +/- 1-sig Chem. Eq. \n([M/H]={popt[0]:.2f}, C/O={popt[1]:.2f})\n at {T}K and {P} bars')
                    except ValueError:
                        print('Unable to plot errorbars from M/H, C/O, R/V fit')
                        ax.plot(self.species_abunds.keys(),best_fit,'X',color='b',
                                    label=f'Best-Fit Chem. Eq. \n([M/H]={popt[0]:.2f}, C/O={popt[1]:.2f})\n at {T}K and {P} bars')                          
        
                ax.set_xlabel('Species',fontsize=14)
                ax.set_ylabel('Abundance (VMR)',fontsize=14)
                ax.tick_params(labelsize=12)
                ax.legend(fontsize=12) 
            
        # If given posteriors
        else:
            inds = np.arange(s)
            samp = random.sample(list(inds),posterior_samples)
            fits = []
            abunds = [[vals[i] for key, vals in self.species_errs.items()] for i in samp]
            #print(np.shape(abunds))
            #print(abunds[0])
            print(f'Running on {s} samples... May take a moment if >10')
            for i, n in enumerate(samp):
                a = abunds[i]
                if fit_refvol:
                    popt,pcov = curve_fit(fit_species,range(len(self.species_abunds)),
                              a,
                              p0=[0.5,0.5,-0.1])
                else:
                    popt,pcov = curve_fit(fit_species,range(len(self.species_abunds)),
                              a,
                              p0=[0.5,0.5])
                fits.append(list(popt))
            self.fits = fits
            mhs = [fits[i][0] for i in range(len(fits))]
            cos = [fits[i][1] for i in range(len(fits))]
            
            mhs_low = np.percentile(mhs,15.8655)
            mhs_high = np.percentile(mhs,100-15.8655)
            mhs_med = np.median(mhs)

            cos_low = np.percentile(cos,15.8655)
            cos_high = np.percentile(cos,100-15.8655)    
            cos_med = np.median(cos)
            
            print('---')
            print(f"Metallicity - Median: {mhs_med:.4f}, Lower bound: {mhs_low:.4f}, Upper bound: {mhs_high:.4f}")
            print(f"C/O - Median: {cos_med:.4f}, Lower bound: {cos_low:.4f}, Upper bound: {cos_high:.4f}")
            
            if fit_refvol:
                refvol = [fits[i][2] for i in range(len(fits))]
                refvol_low = np.percentile(refvol,15.8655)
                refvol_high = np.percentile(refvol,100-15.8655)    
                refvol_med = np.median(refvol)   
                
                popt = [mhs_med,cos_med,refvol_med]
                pcov = fits
                best_fit = fit_species(range(len(self.species_abunds)),mhs_med,cos_med)
                
                print(f"Ref/Vol - Median: {refvol_med:.4f}, Lower bound: {refvol_low:.4f}, Upper bound: {refvol_high:.4f}")
            else:
           
                popt = [mhs_med,cos_med]
                pcov = fits
                best_fit = fit_species(range(len(self.species_abunds)),mhs_med,cos_med)
                

            print('---')
            print('Median Fit Abundances')
            print('---')
            for i, key in enumerate(self.species_abunds.keys()):
                print(f'{key}: {best_fit[i]}')
           
                
            if plot_it:
                if not fit_refvol:
                    fig,ax = plt.subplots(1,2,figsize=(15,6))
                else:
                    fig,ax = plt.subplots(1,3,figsize=(12,6))
                    p = ax[2].hist(refvol)
                    ax[2].vlines(refvol_med,0,1.2*np.max(p[0]),'k')
                    ax[2].vlines(refvol_low,0,1.2*np.max(p[0]),'k',linestyle = '--')
                    ax[2].vlines(refvol_high,0,1.2*np.max(p[0]),'k',linestyle = '--')
                    ax[2].set_xlabel('Ref/Vol')


                    
                p = ax[0].hist(mhs)
                ax[0].vlines(mhs_med,0,1.2*np.max(p[0]),'k')
                ax[0].vlines(mhs_low,0,1.2*np.max(p[0]),'k',linestyle = '--')
                ax[0].vlines(mhs_high,0,1.2*np.max(p[0]),'k',linestyle = '--')

                p = ax[1].hist(cos)
                ax[1].vlines(cos_med,0,1.2*np.max(p[0]),'k',linestyle = '--')
                ax[1].vlines(cos_low,0,1.2*np.max(p[0]),'k',linestyle = '--')
                ax[1].vlines(cos_high,0,1.2*np.max(p[0]),'k',linestyle = '--')
                                    

                ax[0].set_xlabel('[M/H]')
                ax[1].set_xlabel('C/O')
                ax[0].set_ylabel('Probability')      
                
            fig,ax = plt.subplots()
            ax.plot(self.species_abunds.keys(),self.species_errs.values(),'ko',alpha=0.025)
    

            if fit_refvol:
                highco = fit_species(range(len(self.species_abunds)),mhs_med,cos_med+cos_high,refvol_med)
                lowco = fit_species(range(len(self.species_abunds)),mhs_med,cos_med-cos_low,refvol_med)
                highmh = fit_species(range(len(self.species_abunds)),mhs_med+mhs_high,cos_med,refvol_med)
                lowmh = fit_species(range(len(self.species_abunds)),mhs_med-mhs_low,cos_med,refvol_med)
                highrv = fit_species(range(len(self.species_abunds)),mhs_med,cos_med,refvol_med+refvol_high)
                lowrv = fit_species(range(len(self.species_abunds)),mhs_med,cos_med,refvol_med-refvol_low)
            
                ax.errorbar(self.species_abunds.keys(),best_fit,fmt='X',color='b',
                            yerr = (best_fit-np.min([lowco,highco,lowmh,highmh,lowrv,highrv],axis=0),
                                    np.max([lowco,highco,lowmh,highmh,lowrv,highrv],axis=0)-best_fit),
                            label=f'Best-Fit +/- 1-sig Chem. Eq. \n([M/H]={popt[0]:.2f}, C/O={popt[1]:.2f}, R/V={popt[2]:.2f})\n at {T}K and {P} bars')  
            else:
                highco = fit_species(range(len(self.species_abunds)),mhs_med,cos_med+cos_high)
                lowco = fit_species(range(len(self.species_abunds)),mhs_med,cos_med-cos_low)
                highmh = fit_species(range(len(self.species_abunds)),mhs_med+mhs_high,cos_med)
                lowmh = fit_species(range(len(self.species_abunds)),mhs_med-mhs_low,cos_med)
            
                ax.errorbar(self.species_abunds.keys(),best_fit,fmt='X',color='b',
                            yerr = (best_fit-np.min([lowco,highco,lowmh,highmh],axis=0),
                                    np.max([lowco,highco,lowmh,highmh],axis=0)-best_fit),
                            label=f'Best-Fit +/- 1-sig Chem. Eq. \n([M/H]={popt[0]:.2f}, C/O={popt[1]:.2f})\n at {T}K and {P} bars')  
                
            ax.set_xlabel('Species',fontsize=14)
            ax.set_ylabel('Abundance (VMR)',fontsize=14)
            ax.tick_params(labelsize=12)
            ax.legend(fontsize=12) 
            
                     
        #Save best fit M/H and C/O to object
        self.best_fit_mh = popt[0]
        self.best_fit_co = popt[1]
        if fit_refvol:
            self.best_fit_refvol = popt[2]

        
        #Save best fit abundances to object
        best_fit_dict = {}
        for i, key in enumerate(self.species_abunds.keys()):
            best_fit_dict[key] = best_fit[i]
        self.best_fit_abund = best_fit_dict

        #Place this new inferred bulk density into bulk abund attribute
        print('Updating bulk abundances with mean/median fit')
        
        #Init easychem
        exo = self.init_ec()

        exo.metallicity = popt[0]
        #Adjust C/O either by changing C/H or O/H:
        if adjust_C is True:
            exo.atomAbunds[2] = exo.atomAbunds[4] * popt[1]
        elif adjust_C == 'Neither':
            tot = exo.atomAbunds[2] + exo.atomAbunds[4]
            #co = exo.atomAbunds[2] / exo.atomAbunds[4]
            exo.atomAbunds[2] = popt[1] * tot / (popt[1]+1)
            exo.atomAbunds[4] = tot / (popt[1]+1)
        else:
            exo.co = popt[1] #easychem changes O/H by default
            
        if fit_refvol:
            #Enhance elements above O by rv
            for i in range(5,len(exo.atoms)):
                exo.atomAbunds[i] = 10**(np.log10(exo.atomAbunds[i])+popt[2])
        
        bulk_dict = {}
        for i,key in enumerate(exo.atoms):
            bulk_dict[key] = np.log10((exo.atomAbunds[i] / exo.atomAbunds[0])*1e12)
        self.bulk_abunds = bulk_dict

        return popt,pcov,best_fit

    def Lodders25(self):
        solar_abundances = {
            "H": (12.00, 0.00),
            "He": (10.984, 0.02),
            "Li": (3.35, 0.06),
            "Be": (1.48, 0.08),
            "B": (2.85, 0.04),
            "C": (8.46, 0.04),
            "N": (7.90, 0.11),
            "O": (8.76, 0.05),
            "F": (4.53, 0.06),
            "Ne": (7.95, 0.10),
            "Na": (6.37, 0.03),
            "Mg": (7.62, 0.02),
            "Al": (6.54, 0.02),
            "Si": (7.61, 0.02),
            "P": (5.54, 0.04),
            "S": (7.26, 0.04),
            "Cl": (5.33, 0.06),
            "Ar": (6.62, 0.08),
            "K": (5.18, 0.05),
            "K*": (5.18, 0.05),
            "Ca": (6.41, 0.03),
            "Sc": (3.15, 0.04),
            "Ti": (5.00, 0.03),
            "V": (4.07, 0.03),
            "Cr": (5.72, 0.05),
            "Mn": (5.58, 0.03),
            "Fe": (7.54, 0.03),
            "Co": (4.98, 0.03),
            "Ni": (6.29, 0.03),
            "Cu": (4.34, 0.06),
            "Zn": (4.70, 0.04),
            "Ga": (3.17, 0.06)
        }
        return solar_abundances

    def Asplund09(self):
        solar_abundances = {
            "H": (12.00, 0.00),
            "He": (10.93, 0.01),
            "Li": (1.05, 0.10),
            "Be": (1.38, 0.09),
            "B": (2.70, 0.20),
            "C": (8.43, 0.05),
            "N": (7.83, 0.05),
            "O": (8.69, 0.05),
            "F": (4.56, 0.30),
            "Ne": (7.93, 0.10),
            "Na": (6.24, 0.04),
            "Mg": (7.60, 0.04),
            "Al": (6.45, 0.03),
            "Si": (7.51, 0.03),
            "P": (5.41, 0.03),
            "S": (7.12, 0.03),
            "Cl": (5.50, 0.30),
            "Ar": (6.40, 0.13),
            "K": (5.03, 0.09),
            "Ca": (6.34, 0.04),
            "Sc": (3.15, 0.04),
            "Ti": (4.95, 0.05),
            "V": (3.93, 0.08),
            "Cr": (5.64, 0.04),
            "Mn": (5.43, 0.04),
            "Fe": (7.50, 0.04),
            "Co": (4.99, 0.07),
            "Ni": (6.22, 0.04),
            "Cu": (4.19, 0.04),
            "Zn": (4.56, 0.05),
            "Ga": (3.04, 0.09),
            "Ge": (3.65, 0.10)
        }
        return solar_abundances

    def Lodders20(self):
        solar_abundances = {
            "H": 12.00,
            "He": (10.994, 0.02),
            "Li": (3.35, 0.03),
            "Be": (1.40, 0.04),
            "B": (2.85, 0.03),
            "C": (8.56, 0.06),
            "N": (7.94, 0.12),
            "O": (8.82, 0.07),
            "F": (4.79, 0.09),
            "Ne": (8.24, 0.10),
            "Na": (6.36, 0.03),
            "Mg": (7.61, 0.02),
            "Al": (6.51, 0.03),
            "Si": (7.60, 0.01),
            "P": (5.52, 0.03),
            "S": (7.24, 0.03),
            "Cl": (5.32, 0.06),
            "Ar": (6.59, 0.10),
            "K": (5.16, 0.02),
            "Ca": (6.36, 0.03),
            "Sc": (3.13, 0.03),
            "Ti": (4.99, 0.03),
            "V": (4.04, 0.03),
            "Cr": (5.72, 0.02),
            "Mn": (5.56, 0.03),
            "Fe": (7.54, 0.02),
            "Co": (4.95, 0.02),
            "Ni": (6.29, 0.03),
            "Cu": (4.33, 0.04),
            "Zn": (4.70, 0.60),
            "Ga": (3.16, 0.02),
            "Ge": (3.68, 0.03)  # Presumed from periodic order
        }
            
        return solar_abundances
    
    def Asplund21(self):
        solar_abundances = {
            "H": (12.00, 0.00),
            "He": (10.914, 0.013),
            "Li": (0.96, 0.06),
            "Be": (1.38, 0.09),
            "B": (2.70, 0.20),
            "C": (8.46, 0.04),
            "N": (7.83, 0.07),
            "O": (8.69, 0.04),
            "F": (4.40, 0.25),
            "Ne": (8.06, 0.05),
            "Na": (6.22, 0.03),
            "Mg": (7.55, 0.03),
            "Al": (6.43, 0.03),
            "Si": (7.51, 0.03),
            "P": (5.41, 0.03),
            "S": (7.12, 0.03),
            "Cl": (5.31, 0.20),
            "Ar": (6.38, 0.10),
            "K": (5.07, 0.03),
            "Ca": (6.30, 0.03),
            "Sc": (3.14, 0.04),
            "Ti": (4.97, 0.05),
            "V": (3.90, 0.08),
            "Cr": (5.62, 0.04),
            "Mn": (5.42, 0.06),
            "Fe": (7.46, 0.04),
            "Co": (4.94, 0.05),
            "Ni": (6.20, 0.04),
            "Cu": (4.18, 0.05),
            "Zn": (4.56, 0.05),
            "Ga": (3.02, 0.05),
            "Ge": (3.62, 0.10)
        }   
        return solar_abundances
    
    def Lodders10(self):
        solar_abundances = {
            "H": (12.00,0.0),
            "He": (10.925, 0.02),
            "Li": (3.28, 0.05),
            "Be": (1.32, 0.03),
            "B": (2.81, 0.04),
            "C": (8.39, 0.04),
            "N": (7.86, 0.12),
            "O": (8.73, 0.07),
            "F": (4.44, 0.06),
            "Ne": (8.05, 0.10),
            "Na": (6.29, 0.04),
            "Mg": (7.54, 0.06),
            "Al": (6.46, 0.07),
            "Si": (7.53, 0.06),
            "P": (5.45, 0.05),
            "S": (7.16, 0.02),
            "Cl": (5.25, 0.06),
            "Ar": (6.50, 0.10),
            "K": (5.11, 0.04),
            "Ca": (6.31, 0.02),
            "Sc": (3.07, 0.02),
            "Ti": (4.93, 0.03),
            "V": (3.99, 0.03),
            "Cr": (5.65, 0.02),
            "Mn": (5.50, 0.01),
            "Fe": (7.46, 0.08),
            "Co": (4.90, 0.08),
            "Ni": (6.22, 0.04),
            "Cu": (4.27, 0.04),
            "Zn": (4.65, 0.04),
            "Ga": (3.10, 0.02),
            "Ge": (3.59, 0.06)}
        return solar_abundances
    
    def Lodders10(self):
        solar_abundances = {
            "H": (12.00, 0.00),
            "He": (10.93, 0.01),
            "Li": (1.05, 0.10),
            "Be": (1.38, 0.09),
            "B": (2.70, 0.20),
            "C": (8.39, 0.05),
            "N": (7.78, 0.06),
            "O": (8.66, 0.05),
            "F": (4.56, 0.30),
            "Ne": (7.84, 0.06),
            "Na": (6.17, 0.04),
            "Mg": (7.53, 0.09),
            "Al": (6.37, 0.06),
            "Si": (7.51, 0.04),
            "P": (5.36, 0.04),
            "S": (7.14, 0.05),
            "Cl": (5.50, 0.30),
            "Ar": (6.18, 0.08),
            "K": (5.08, 0.07),
            "Ca": (6.31, 0.04),
            "Sc": (3.05, 0.08),
            "Ti": (4.90, 0.06),
            "V": (4.00, 0.02),
            "Cr": (5.64, 0.10),
            "Mn": (5.39, 0.03),
            "Fe": (7.45, 0.05),
            "Co": (4.92, 0.08),
            "Ni": (6.23, 0.04),
            "Cu": (4.21, 0.04),
            "Zn": (4.60, 0.03),
            "Ga": (2.88, 0.10),
            "Ge": (3.58, 0.05)
        }
        return solar_abundances
    
    def Caffau11(self):
        #Asplund09 fill in the rest
        solar_abundances = {
            "H": (12.00, 0.00),
            "He": (10.93, 0.01),
            "Li": (1.03, 0.03),
            "Be": (1.38, 0.09),
            "B": (2.70, 0.20),
            "C": (8.50, 0.06),
            "N": (7.86, 0.12),
            "O": (8.76, 0.07),
            "F": (4.56, 0.30),
            "Ne": (7.93, 0.10),
            "Na": (6.24, 0.04),
            "Mg": (7.60, 0.04),
            "Al": (6.45, 0.03),
            "Si": (7.51, 0.03),
            "P": (5.46, 0.04),
            "S": (7.16, 0.05),
            "Cl": (5.50, 0.30),
            "Ar": (6.40, 0.13),
            "K": (5.11, 0.09),
            "Ca": (6.34, 0.04),
            "Sc": (3.15, 0.04),
            "Ti": (4.95, 0.05),
            "V": (3.93, 0.08),
            "Cr": (5.64, 0.04),
            "Mn": (5.43, 0.04),
            "Fe": (7.52, 0.06),
            "Co": (4.99, 0.07),
            "Ni": (6.22, 0.04),
            "Cu": (4.19, 0.04),
            "Zn": (4.56, 0.05),
            "Ga": (3.04, 0.09),
            "Ge": (3.65, 0.10)
        }
        return solar_abundances
    
    def import_masses(self):
        self.molecular_mass = {
            # Atoms
            'H': 1.00784,
            'He': 4.002602,
            'C': 12.0107,
            'N': 14.0067,
            'O': 15.999,
            'F': 18.998403163,
            'Ne': 20.1797,
            'Na': 22.98976928,
            'Mg': 24.305,
            'Al': 26.9815385,
            'Si': 28.085,
            'P': 30.973761998,
            'S': 32.06,
            'Cl': 35.45,
            'Ar': 39.948,
            'K': 39.0983,
            'Ca': 40.078,
            'Fe': 55.845,
            'Ti': 47.867,
            'V': 50.9415,
        
            # Molecules
            'H2': 2 * 1.00784,
            'H2O': 2 * 1.00784 + 15.999,
            'OH': 1.00784 + 15.999,
            'O2': 2 * 15.999,
            'O3': 3 * 15.999,
            'CO': 12.0107 + 15.999,
            'CO2': 12.0107 + 2 * 15.999,
            'CH4': 12.0107 + 4 * 1.00784,
            'C2H2': 2 * 12.0107 + 2 * 1.00784,
            'C2H4': 2 * 12.0107 + 4 * 1.00784,
            'C2H6': 2 * 12.0107 + 6 * 1.00784,
            'NH3': 14.0067 + 3 * 1.00784,
            'N2': 2 * 14.0067,
            'NO': 14.0067 + 15.999,
            'NO2': 14.0067 + 2 * 15.999,
            'HCN': 1.00784 + 12.0107 + 14.0067,
            'SiO': 28.085 + 15.999,
            'SO': 32.06 + 15.999,
            'SO2': 32.06 + 2 * 15.999,
            'H2S': 2 * 1.00784 + 32.06,
            'PH3': 30.973761998 + 3 * 1.00784,
            'CS': 12.0107 + 32.06,
            'HCl': 1.00784 + 35.45,
            'FeH': 1.00784 + 55.845,
            'TiO': 47.867 + 15.999,
            'VO': 50.9415 + 15.999
        }
        return

