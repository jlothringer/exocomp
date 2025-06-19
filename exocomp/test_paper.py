#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  7 15:21:49 2025

@author: jlothringer
"""

##Test
## M/H= 1.0 and C/O = 0.5 at T=1000 and P = 1
#{'H2O': -2.193862599311742, 'CO': -2.9951664932143895,'CO2': -4.957449831706032} for Asplund 09
#{'H2O': -2.196454926078874, 'CO': -2.95559683423187, 'CO2': -4.9118293333092335} for Lodders 25
g = Abund(species_abunds={'H2O': -2.196454926078874, 'CO': -2.95559683423187, 'CO2': -4.9118293333092335})
g.convert_species_abunds(1000, 1, plot_it=True) # Need to add errors or you get weird errors!

g = Abund(species_abunds={'H2O': -2.196454926078874, 'CO': -2.95559683423187, 'CO2': -4.9118293333092335},
          species_errs={'H2O': 0.1, 'CO': 0.1, 'CO2': 0.1})
g.convert_species_abunds(1000, 1, plot_it=True) # Doesn't matter if C/O adjusts C or O here
g.convert_species_abunds(1000, 1, plot_it=True, adjust_C=False) 
g.convert_species_abunds(1000, 1, plot_it=True, adjust_C='Neither') 

# Or you can give it posteriors
# A little slow with posteriors and many samples
g = Abund(species_abunds={'H2O': -2.196454926078874, 'CO': -2.95559683423187, 'CO2': -4.9118293333092335},
          species_errs={'H2O':np.random.normal(-2.19645,0.1,100),'CO':np.random.normal(-2.9556,0.1,100),'CO2':np.random.normal(-4.9118,0.1,100)})
g.convert_species_abunds(1000, 1, plot_it=True, posterior_samples=100)
g.convert_species_abunds(1000, 1, plot_it=True, adjust_C=False, posterior_samples=100) 
g.convert_species_abunds(1000, 1, plot_it=True, adjust_C='Neither', posterior_samples=100) 



# Or you can try fitting the equivalent in Asplund 09 abundances
g = Abund(species_abunds={'H2O': -2.196454926078874, 'CO': -2.95559683423187, 'CO2': -4.9118293333092335},
          species_errs={'H2O': 0.1, 'CO': 0.1, 'CO2': 0.1}, solar='Asplund09')
g.convert_species_abunds(1000, 1, plot_it=True) # Doesn't matter if C/O adjusts C or O here
g.convert_species_abunds(1000, 1, plot_it=True, adjust_C=False) 

# But let's try what actually corresponds to M/H=1 and C/O=0.5 for Asplund09:
g = Abund(species_abunds={'H2O': -2.193862599311742, 'CO': -2.9951664932143895,'CO2': -4.957449831706032},
          species_errs={'H2O': 0.1, 'CO': 0.1, 'CO2': 0.1}, solar='Asplund09')
g.convert_species_abunds(1000, 1, plot_it=True) # Doesn't matter if C/O adjusts C or O here
g.convert_species_abunds(1000, 1, plot_it=True, adjust_C=False) 


#WASP-17b
g = Abund(species_abunds={'H2O':-2.96,'CH4':-9.05,'CO2':-5.81},
          species_errs={'H2O':0.31,'CH4':2.1,'CO2':1.31})
g.convert_species_abunds(1271, 5e-3,plot_it=True,fit_refvol=False)
g.convert_species_abunds(1271, 5e-3,plot_it=True,fit_refvol=False,adjust_C=False)
g.convert_species_abunds(1271, 5e-3,plot_it=True,fit_refvol=False,adjust_C='Neither')



g = Abund(species_abunds={'H2O':-2.96,'CH4':-9.05,'CO2':-5.81,'K':-8.07,'Na':-9.31},
          species_errs={'H2O':0.31,'CH4':2.1,'CO2':1.31,'K':0.58,'Na':1.77})
g.convert_species_abunds(1271, 5e-3,plot_it=True,fit_refvol=False)
g.convert_species_abunds(1271, 5e-3,plot_it=True,fit_refvol=False,adjust_C=False)
g.convert_species_abunds(1271, 5e-3,plot_it=True,fit_refvol=False,adjust_C='Neither')

g.convert_species_abunds(1271, 5e-3,plot_it=True,fit_refvol=True)
g.convert_species_abunds(1271, 5e-3,plot_it=True,fit_refvol=True,adjust_C=False)
g.convert_species_abunds(1271, 5e-3,plot_it=True,fit_refvol=True,adjust_C='Neither')



#posterior test:
g = Abund(species_abunds={'H2O':-2.96,'CH4':-9.05,'CO2':-5.81,'K':-8.07,'Na':-9.31},
          species_errs={'H2O':np.random.normal(-2.96,0.31,100),'CH4':np.random.normal(-9.05,2.1,100),'CO2':np.random.normal(-5.81,1.31,100),
                        'K':np.random.normal(-8.07,0.58,100),'Na':np.random.normal(-9.31,1.77,100)})

g.convert_species_abunds(1271, 5e-3,posterior_samples=100,plot_it=True,fit_refvol=False)

# W39b test
# Aurora free
#Must remove C2H4 and SO2
#g = Abund(species_abunds={'H2O':-2.96,'CO':-2.6,'CH4':-8.5,'CO2':-3.7,'C2H2':-8.7,'HCN':-7.9,
#                          'NH3':-9.5,'H2S':-6.5,'SO2':-5.5,'K':-8.07,'Na':-9.31},
#          species_errs={'H2O':0.2,'CO':0.5,'CH4':0.5,'CO2':0.3,'C2H2':2.1,'HCN':-2.5,
#                        'NH3':1.8,'H2S':3.4,'SO2':0.2,'K':0.6,'Na':0.5})
g = Abund(species_abunds={'H2O':-2.96,'CO':-2.6,'CH4':-8.5,'CO2':-3.7,'HCN':-7.9,
                          'NH3':-9.5,'H2S':-6.5,'K':-8.07,'Na':-9.31},
          species_errs={'H2O':0.2,'CO':0.5,'CH4':2.2,'CO2':0.3,'HCN':2.5,
                        'NH3':1.8,'H2S':3.4,'K':0.6,'Na':0.5})

g.convert_species_abunds(1200, 5e-3,plot_it=True,fit_refvol=False)
g.convert_species_abunds(1200, 5e-3,plot_it=True,fit_refvol=True)
g.convert_species_abunds(1200, 5e-3,plot_it=True,fit_refvol=False,adjust_C=False)
g.convert_species_abunds(1200, 5e-3,plot_it=True,fit_refvol=True,adjust_C=False)
g.convert_species_abunds(1200, 5e-3,plot_it=True,fit_refvol=False,adjust_C='Neither')
g.convert_species_abunds(1200, 5e-3,plot_it=True,fit_refvol=True,adjust_C='Neither')


# WASP-121b - example with dissociation and RV
# PETRA PC
g = Abund(species_abunds={'H2O':-5.74,'CO':-2.34,'SiO':-4.60},species_errs={'H2O':0.12,'CO':0.34,'SiO':0.17})
g.convert_species_abunds(2000, 1e-4, plot_it=True,fit_refvol=True)

# PETRA TR
g = Abund(species_abunds={'H2O':-4.95,'CO':-2.9,'SiO':-3.87},species_errs={'H2O':0.58,'CO':0.67,'SiO':0.62})
g.convert_species_abunds(2000, 1e-4, plot_it=True,fit_refvol=True)

# NEMESIS PC
g = Abund(species_abunds={'H2O':-8.3,'CO':-1.08,'SiO':-3.48},species_errs={'H2O':1.0,'CO':0.08,'SiO':0.37})
g.convert_species_abunds(3500, 1e-3, plot_it=True,fit_refvol=True)


# NEMESIS TR
g = Abund(species_abunds={'H2O':-7.5,'CO':-1.27,'SiO':-3.51},species_errs={'H2O':1.0,'CO':0.13,'SiO':0.47})
g.convert_species_abunds(3000, 1e-3, plot_it=True,fit_refvol=True)

# W121b Dayside - PETRA
g = Abund(species_abunds={'H2O':-3.92,'CO':-2.03,'SiO':-3.15},species_errs={'H2O':0.29,'CO':0.36,'SiO':0.28})
g.convert_species_abunds(3000, 1e-3, plot_it=True,fit_refvol=True) #>3000 breaks b/c of dissociation
g.convert_species_abunds(3000, 1e-3, plot_it=True,fit_refvol=True,adjust_C=False)


# W121b Dayside - NEMESIS
g = Abund(species_abunds={'H2O':-3.31,'CO':-1.29,'SiO':-2.66},species_errs={'H2O':0.26,'CO':0.19,'SiO':0.21})
g.convert_species_abunds(3000, 1e-3, plot_it=True,fit_refvol=True) #>3000 breaks b/c of dissociation <-- interesting! quantify T vs. M/H relation/bias
g.convert_species_abunds(3000, 1e-3, plot_it=True,fit_refvol=True,adjust_C=False)


#WASP-178b
g = Abund(species_abunds={'H2O':-6.26,'CO':-5.44,'SiO':-5.08},species_errs={'H2O':0.14,'CO':0.24,'SiO':0.18})
g.convert_species_abunds(3300, 1e-3, plot_it=True,fit_refvol=True)
g.convert_species_abunds(3000, 1e-3, plot_it=True,fit_refvol=True) #Difference between 3000 and 3300K is crazy! UHJs are hard!

#all atoms
#no Mg + in easychem
# g = Abund(species_abunds={'H2O':-6.26,'CO':-5.44,'SiO':-5.08,'CO2':-9.58,
#                           'Fe':-6.56,'Fe+':-10.06,'FeH':-6.4,'Mg':-6.0,
#                           'Mg+':-3.35,'TiO':-10.73,'VO':-10.29},
#           species_errs={'H2O':0.14,'CO':0.24,'SiO':0.18,'CO2':1.09,
#                         'Fe':3.37,'Fe+':0.93,'FeH':3.4,'Mg':3.35,
#                         'Mg+':0.17,'TiO':0.58,'VO':0.83})
g = Abund(species_abunds={'H2O':-6.26,'CO':-5.44,'SiO':-5.08,'CO2':-9.58,
                          'Fe':-6.56,'Fe+':-10.06,'FeH':-6.4,'Mg':-6.0,
                          'TiO':-10.73,'VO':-10.29},
          species_errs={'H2O':0.14,'CO':0.24,'SiO':0.18,'CO2':1.09,
                        'Fe':3.37,'Fe+':0.93,'FeH':3.4,'Mg':3.35,
                        'TiO':0.58,'VO':0.83})
g.convert_species_abunds(3300, 1e-3, plot_it=True,fit_refvol=True) #again, very different if T=3000 <- small changes in scale height can result in big changes to chem eq. interpretation (even if molecular abundances are the same)
g.convert_species_abunds(3300, 1e-3, plot_it=True,fit_refvol=True,adjust_C=False) #again, very different if T=3000 <- small changes in scale height can result in big changes to chem eq. interpretation (even if molecular abundances are the same)
g.convert_species_abunds(3300, 1e-3, plot_it=True,fit_refvol=True,adjust_C='Neither') #again, very different if T=3000 <- small changes in scale height can result in big changes to chem eq. interpretation (even if molecular abundances are the same)


##
## Convert chem eq results to atomic abundnaces
##

#W178b
g = Abund(retrieval = 'pRT')
g.mh_type = 'O/H'
g.co_type = 'C/H'
bulk = g.convert_bulk_abundance(1.47,0.01,-1.56) #without err
bulk = g.convert_bulk_abundance(1.47,0.01,-1.56, mh_err = 1.1, co_err=0.01, rv_err=0.97) #sym err
bulk = g.convert_bulk_abundance(1.47,0.01,-1.56, mh_err = [1.1,0.28], co_err=[0.01,0.01], rv_err=[0.3,0.97]) #asym err
bulk = g.convert_bulk_abundance(1.47,0.01,-1.56, mh_err = np.random.normal(1.47,0.5,100), co_err=np.random.normal(0.4,0.1,100), rv_err=np.random.normal(-0.5,0.5,100)) #posterior

#JWST-only Chem eq (Fe/H, C/O)
samples_eq = np.genfromtxt('/Users/jlothringer/Research/extreme_irrad/wasp178/W178_eqchem_JWST_retrieval_10bin_1000LP_iso_post_equal_weights.dat')
mh = samples_eq[:,3]
co = samples_eq[:,4]

bulk = g.convert_bulk_abundance(np.median(mh), np.median(co), mh_err = mh, co_err=co) #posterior

fig,ax = subplots(2,2,figsize=(12,6))
ax[0,0].hist(mh)
ax[0,1].hist(co)
ax[1,0].hist(bulk[1]['O'])
ax[1,1].hist(bulk[1]['C'])

#full spectrum
samples_free = np.genfromtxt('/Users/jlothringer/Research/extreme_irrad/wasp178/evaluate_W178_eqchem_JWST_retrieval_10bin_1000LP_iso_new_offs_hot/W178_eqchem_JWST_retrieval_10bin_1000LP_iso_new_offs_hot_post_equal_weights.dat')


##
## Posterior test
##

h2o = np.random.random(10000)*-11 -1
co = np.random.random(10000)*-11 -1
co2 = np.random.random(10000)*-11 -1
ch4 = np.random.random(10000)*-11 -1
na = np.random.random(10000)*-11 -1
k = np.random.random(10000)*-11 -1

co_ratio = (10**co*1 + 10**co2*1 + 10**ch4*1) / (10**h2o*1+10**co*1+10**co2*2)
refvol = (10**na+10**k) / (10**h2o*1+10**co*2+10**co2*3+10**ch4*1)
refvol[np.where(refvol == 0.0)[0]] = 1e-14

figure()
hist(co_ratio,bins=50,range=(0,1.5))

figure()
hist(np.log10(refvol),bins=50)

