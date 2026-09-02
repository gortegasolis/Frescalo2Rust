NOTES ON RUNNING FRESCALO - Mark Hill March 2011

FRESCALO and supporting programs SAMPDIST and NEIGHSIM work from data files to estimated time factors.  If there is only a single time, which may be blank, time factors are still produced but are irrelevant.  It is suggested that you try the programs on the trial datasets to see how the calculations proceed.  The two supporting programs allow you to calculate neighbourhood weightings for use in FRESCALO.

*** CAUTION ABOUT LOCATION AND SPECIES NAMES ***
The programs assume that location names and species names have no blanks, and that they are not more than 10 characters long.  If you have longer names than this, or names with blanks, then you must use an abbreviation or a code.  Numbers are acceptable as codes but are not recommended, as it is easier to work with recognizable names.

0.  DATA
In the example, the data are as follows (for eastern England);  

Test.txt – this dataset is a summary by quinquennium of the dataset Test_B, referenced in the supplementary data to the 2011 local frequency paper. There are 10 date classes.  A date such as 1997 refers to the quinquennium 1995-1999, of which 1997 is the middle year.  The data are in the form location, species, dateclass

samples.txt - positons of locations, specified here as kilometres east and north of the GB Ordnance Survey SW corner.  Each location is a 10-km square or 'hectad'.

Training_vasc.txt - vascular plant data for hectads in eastern England (used only as a training set to define neighbourhoods).

1.  CALCULATION OF PHYSICAL DISTANCES FROM POSITIONS OF LOCATIONS
This calculation is done by SAMPDIST.  The output lists the nearest hectads to the target hectad.
Input data - samples.txt
User-supplied parameter - Number of neighbours = 200
Output - dist.txt

2.  CALCULATION OF NEIGHBOURHOOD WEIGHTS FOR LOCATIONS
This is done by NEIGHSIM, and requires distances to be input from previous stage.  It also requires a list of species (here vascular plants) for each location, to judge their biological similarity as well as their distance.
Input data - dist.txt (from previous stage), Training_vasc.txt
User-supplied parameter - Number of neighbours = 100 (should be lower than number in SAMPDIST)
Output - sim.txt, weights.txt
sim.txt – floristic similarity between samples and neighbours; not used in subsequent calculations, but may be useful for reference
weights.txt - neighbourhood weighting for each hectad;  it defines the neighbourhoods and is used in FRESCALO

3.  CALCULATION OF SPECIES PROBABILITIES AND TIME FACTORS
This is done by FRESCALO; the input parameters are all listed in the log file.
Input data - Test.txt (bryophyte data in date classes), weights.txt (from previous stage)
User-supplied parameters - both mean local frequency and benchmark limit are set to default values
Output – 
      log.txt - the log file showing parameters and file names fed to FRESCALO
      samples.txt - statistics for each location
      frequencies.txt - rescaled frequencies for each location and species
      trends.txt - time factors for all species
Note 1.  It is possible to recalculate time factors, marking particular species as unsuitable to be used as benchmarks.  A list is supplied here (but not used in the example) consisting of those species for which the mean time factors in the early period (1960-1984) differed from those in the later period (1985-2009) by more than a factor of 2.  Omitting these ‘dynamic’ species has almost no effect on calculated time factors.
Note 2.  Where, as here, a locality is present in the species data (Test.txt) but not in the supplied weights (weights.txt), it cannot be used.  Such localities are listed in the log file.  Here there is one such hectad, TA25, which has perhaps 200 square metres of land above sea level; it is simply thrown out of the calculation.

EXPLANATION OF HEADERS IN OUTPUT FROM FRESCALO

(a) Log file (log.txt)
This should be self-explanatory

(b) Location report (samples.txt)
Location - Name of location; in this case locations are hectads of the GB National Grid
Loc_no - Numbering (added) of locations in alphanumeric order
No_spp - Number of species at that location; the actual number which may be zero
Phi_in - Initial value of phi, the frequency-weighted mean frequency
Alpha - Sampling effort multiplier (to achieve standard value of phi)
Wgt_n2 – ‘effective number’ N2 for the neighbourhood weights; this is small if there are few floristically similar hectads close to the target hectad.  It is (sum weights)^2 / (sum weights^2)
Phi_out - Value of phi after rescaling; constant, if the algorithm has converged
Spnum_in - Sum of neighbourhood frequencies before rescaling
Spnum_out - Estimated species richness, i.e. sum of neighbourhood frequencies after rescaling
Iter - Number of iterations for algorithm to converge

(c) Listing of rescaled species frequencies (frequencies.out)
Location - Name of location
Species - Name of species (a number can be used as a name if that is convenient)
Pres - Record of species in location (1 = recorded, 0 = not recorded)
Freq - Frequency of species in neighbourhood of location
Freq_1 - Estimated probabilty of occurrence, i.e. frequency of species after rescaling
SD_Frq1 – Standard error of Freq_1, calculated on the assumption that Freq is a binomial variate with standard error sqrt(Freq*(1-Freq)/ Wgt_n2), where Wgt_n2 is as defined for samples.txt in section (b)
Rank - Rank of frequency in neighbourhood of location
Rank_1 - Rescaled rank, defined as Rank/Estimated species richness

(d) Listing of time factors for species
Species - Name of species
Time - Time period, specified as a class (e.g. 1970); times need not be numeric and are indexed as character strings
TFactor - Time factor, the estimated relative frequency of species at the time
St_Dev - Standard deviation of the time factor, given that spt (defined below) is a weighted sum of binomial variates
Count - Number of occurrences of species at the time period
spt - Number of occurrences, given reduced weight of locations having very low sampling effort
est - Estimated number of occurrences; this should be equal to spt if the algorithm has converged
N>0.00 - Number of locations with non-zero probability of the species occurring
N>0.98 - Number of locations for which the probability of occurrence was estimated as greater than 0.98
Note 3. In FRESCALO, data from times at a location in which the proportion of benchmark species (recording effort) is less than 0.1 are given low weight, ranging linearly from a minimum of 0.05 to 1.0 if sampling effort is 0.095.  The logic of this is that where no systematic sample is taken, observations are isolated and unsupported, and are rather unsuitable for calculating time factors.


APPENDIX - technical notes on FRESCALO

There are two key subroutines
- fresca calculates sampling effort and probabilities of species occurrence based on frequencies in the neighbourhood
- tfcalc calculates time factors based on sampling effort and probability of occurrence 

subroutine fresca(i,n,itot,jocc,f,ff,jrank,samp,sp,phibig,blim,fmax,fmin,wn2,spnum,tol,irepmx)
Parameters are:-
i - serial number of location
n - total number of species
itot - number of species at location i
jocc(n) - recorded or not recorded (1 or 0)
f(n) - input species frequencies
ff(n) - workspace, mainly for rescaled frequencies
jrank(n) - integer for sorting to calculate rank order
samp - name of Location
sp(n) - names of species
phibig - target value of phi, the frequency-weighted mean frequency
blim - limit for accepting benchmark species; this is 0.2703 by default
fmax - maximum value of frequency, set as parameter to fmax=0.99999
fmin - minimum value of frequency, set as parameter to fmin=1.0E-10
wn2 - effective species number, the reciprocal of Simpson's index
spnum – probability sum of species at location i
tol - tolerance for convergence of phi to phibig, set as parameter to 0.0003
irepmx - maximum number of iterations, set to 100

The algorithm works by seeking a value of alpha, the Sampling Effort Multiplier, such that when frequencies are rescaled as
   (1)   ff(j)=1-exp(-ff(j)*alpha)
Then the value of phi, the frequency-weighted mean frequency, achieves its traget value phibig. For iterations 1 to 19, alpha is adjusted by the empirical formula
   (2)   alpha=alpha*exp(1.86*(log(1-phi)-log(1-phibig)))
This greatly increases the speed of convergence if the frequencies correspond to a typically-shaped rank-frequency curve.  If the data are of a less tractable shape, convergence may fail, so for iterations 20 to 100 alpha is adjusted by the formula
   (3)   alpha=alpha*phibig/phi
This normally has much slower convergence, but should work for all data, including cases where there is an unusual shape of the rank-frequency curve.  If the reported number of iterations is 100, then convergence has failed, and the degree to which it has failed can be judged from the location report.

subroutine tfcalc(tf,sd,sptot,jtot,esttot,iocc,smpint,fff,m,ic1,ic2)
Parameters, specified for a species j at time t are:-
tf - Time factor
sd - Standard deviation of time factor
sptot - Number of occurrences, given reduced weight of locations having very low sampling effort
jtot - Actual number of occurrences of species j at time t
esttot - Estimated number of occurrences of species j at time t, given fitted model
iocc(i) - is 1 if species j is found at location i at time t, 0 otherwise
smpint(i) - is the sampling intensity at location i and time t
fff(i) - smoothed time-independent frequency of species j at location i
m - number of locations
ic1 - Number of locations with non-zero probability of the species occurring
ic2 - Number of locations for which the probability of occurrence was estimated as greater than 0.98

In this subroutine, the time factor is adjusted so that the estimated number of occurrences of the species is equal to the actual number.
