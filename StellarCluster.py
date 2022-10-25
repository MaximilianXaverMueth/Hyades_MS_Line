#Finds Distance estimates for a specified cluster
#Error propigation via taylor approximation (assuming symettric gaussian errors)
    #TODO: Correct Spectroscopic parallax distance estimate

    #TODO: add cepheid PL(C) distance estimate
    #TODO: get virial mass
    #TODO: add pipeline for queries (gaia, simbad)

    #TODO: get better metallicity source
    #TODO: add ZAMS MS fit
    #TODO: add isochrone fit, also get age
    #TODO: add better cluster membership filtering with algorithms


#Imports
from ssl import PEM_cert_to_DER_cert
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy import spatial
from scipy import stats as stats_sci
import statistics as stat
import sys
from uncertainties import ufloat
from uncertainties.umath import *

#Stopping printing of warnings
pd.set_option("mode.chained_assignment", None)

#Main
def main():
#----------SET USER PARAMETERS HERE---------
    #Set z-score for pm exclusion
    z=1
    #Choose GAIA CSV input file name
    GAIAname = "Pleiades-result.csv"
    #SIMBADname = "pleiadesSim.csv"
    Cepheidname= "Cepheid-table.odt"
    GaiaCepheid = "Cepheid-result.csv"
#----------SET USER PARAMETERS HERE---------
    gaia = pd.read_csv(GAIAname)
    cepheidT1 = pd.read_csv(Cepheidname, sep=" ")
    cepheid= pd.read_csv(GaiaCepheid)
    #Filtering gaia data for nonvalues and similar apparant magnitudes
    maxMag = 10
    gaia = gaia.dropna(subset=["phot_g_mean_mag","phot_bp_mean_mag","phot_rp_mean_mag","parallax"])
    gaia = gaia[(gaia["phot_g_mean_mag"] < maxMag) & (gaia["phot_bp_mean_mag"] < maxMag) & (gaia["phot_rp_mean_mag"] < maxMag) & (gaia["parallax"] > 0)]

    #Get cluster members by removing pm and parallax z-std outliers
    gaia = getMems(gaia,z)

    #Get Spectral Class test case for B
    SpecB= gaia[gaia["spectraltype_esphs"]=="B"]

    #Apply Variable Extinction Method; correct for extinction and reddening
    R, R_err, dist_ext, err_dist_ext = varExt(gaia)
    gaia = ext_correct(gaia,R,R_err)

    #Make line of best fit
    a ,b = np.polyfit(gaia["(B-V)_intr"],gaia["V_rel", 1])

    #Plot corrected Hyades CMD
    plt.figure()
    plt.gca().invert_yaxis()
    plt.title("Corrected CMD")
    plt.xlabel("(B-V)_intr [Mag]")
    plt.ylabel("V_rel [Mag]")
    plt.scatter(gaia["(B-V)_intr"],gaia["V_rel"], color="red", label = "Cluster Stars")
    plt.plot(x, a*x+b)
    plt.legend()
    plt.show()

   
    #cephid distnce
    #outside of the table cepheid apparent magnitude can be found many different ways and the code should be updated accordingly
    AbsMag=[-1*(2.76*(np.log10(cepheid["PF"]))-1.0)-4.16]
    apparentMag = [cepheid["INT_AVERAGE_G"]]
    arrayAbsMag = np.array(AbsMag)
    arrayapparentMag = np.array(apparentMag)
    dist = [10**((arrayapparentMag - arrayAbsMag +5)/5)]
    AbsMagError = 2.76*(np.log10(cepheid["PF_ERROR"]))
    ApparentMagError= cepheid["INT_AVERAGE_G_ERROR"]
    arrayAbsMagError = np.array(AbsMagError)
    arrayApparentMagError = np.array(ApparentMagError)
    DError= [10**((arrayApparentMagError- arrayAbsMagError)/5)]

    #make the cephid distance table I have no clue how to express the error in the same table as the distance value I don't even know if this code will result in a tbale I'm just kinda guessing
    data= dist
    collumns = ["Distance in pc",]
    rows = [cepheid["SOURCE_ID"]]
    the_table= plt.table(cellText=data, rowLabels=rows, colLabels=collumns)
    plt.show    
   




def getMems(stars,z):
    stars = stars[(np.abs(stats_sci.zscore(stars["pmra"]))<z) & (np.abs(stats_sci.zscore(stars["pmdec"]))<z) & (np.abs(stats_sci.zscore(stars["dist_trig_parallax"]))<z)]
    return stars

#Applies Variable Extinction Method and returns R value and distance estimate
def varExt(stars):
        #Regresson to get R
        #errors estimated by taking sigma=sqrt(err_dMod_app**2 + (3.1)**2 * err_colorExcess**2)
    popt,pcov = curve_fit(ext_fit, stars["colorExcess"], stars["dMod_app"],sigma=np.sqrt(stars["err_dMod_app"]**2 + (3.1)**2 * stars["err_colorExcess"]**2))
        #get fit parameters
    perr = np.sqrt(np.diag(pcov))
    r = popt[0]
    r_err = perr[0]
    dMod_real = popt[1]
    dMod_real_err = perr[1]
        #Calculate distance estimate with error
    dMod = ufloat(dMod_real,dMod_real_err)
    dist = 10**((dMod + 5) / 5)
    dist_err = dist.std_dev
    dist = dist.nominal_value

    #Print Results
    print("R = " + str(sciRound(r,r_err)))
        #Use R=3.1 if regression R values doesn't make sense
    if r<0 or r>5:
        r,r_err = 3.1,.1
        print("Assuming R=3.1 +/- .1")

    #Create variable Extinciton graph
    plt.figure()
    plt.title("Variable Extinction (Non-Reddened Removed from R calculation)")
    plt.xlabel("E(B-V) [Mag]")
    plt.ylabel("dMod_app [Mag]")
    plt.scatter(stars["colorExcess"],stars["dMod_app"], color="red", label = "Cluster Stars")
    es = np.linspace(np.min(stars["colorExcess"]),np.max(stars["colorExcess"]),10000)
    plt.plot(es,(dMod_real + es*r), color="blue", label = "Variable Extinction Trend")
    plt.legend()
    plt.show()

    #Return R and extinction distance estimate
    return r, r_err, dist, dist_err

#Returns pandas dataframe to correct for extinction and reddening
def ext_correct(stars,r,r_err):
    #Use intrinsic color from now on to account for reddening
    stars["V_rel"] = stars["V_app"] - r*stars["colorExcess"]
    stars["err_V_rel"] = np.sqrt(stars["err_V_app"]**2 + r**2 * stars["err_colorExcess"]**2)
    return stars
    
#Equation for variable extinction regression
def ext_fit(e, r, dMod_real):
    dMod_app = dMod_real + r*e
    return dMod_app

#round num to 2 sig figs
def figRound(num):
    return float("%.2g" % num)

#get number decimals
def numDec(num):
    if num == 0:
        return 2
    else:
        return -int(np.floor(np.log10(abs(num)))) + 1

#do value and error rounding (to 2 sig figs in error)
def sciRound(num,err):
    err = figRound(err)
    num = round(num, numDec(err))
    return [num,err]
