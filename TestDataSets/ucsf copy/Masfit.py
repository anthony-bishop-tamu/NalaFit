#MASfit v2.0.0
#Goal of this script is to improve the generality of the fitting software and work with python 3.0+
import matplotlib.pyplot as plt
import numpy
from matplotlib.ticker import ScalarFormatter
import matplotlib.mlab as mlab
import scipy.optimize
from scipy.stats import norm, chi2
import warnings
import scipy.optimize as opt
import pandas as pd
import re
import math
import numpy as np
import os
from concurrent.futures import ThreadPoolExecutor
import multiprocessing
import sys
import cProfile

VERSION_NUMBER="v2.2.0"
FIRSTDATACOLUMN=1
def denormalizeDataFrame(dataFrame,normalizationFactors):
    result = normalizationFactors[:,np.newaxis]* np.array(dataFrame.iloc[:,FIRSTDATACOLUMN:].values)
    #result = dataFrame.iloc[:,3:].values
    dataFrame.iloc[:,FIRSTDATACOLUMN:] = pd.DataFrame(result)
    return dataFrame, normalizationFactors
#
def normalizeDataFrame(dataFrame: pd.DataFrame):
    normalizationFactors = np.max(dataFrame.iloc[:,FIRSTDATACOLUMN:].values,axis=1)
    result = dataFrame.iloc[:, FIRSTDATACOLUMN:].values/normalizationFactors[:,np.newaxis]
    print(dataFrame)
    print(result)
    dataFrame.iloc[:, FIRSTDATACOLUMN:] = pd.DataFrame(result)
    return normalizationFactors
#

def calculateNoise(timeValues,dataFrame: pd.DataFrame,outputFile):
    index_dict = {}
    #determine duplicates
    print(dataFrame)
    for i, value in enumerate(timeValues):
        if value not in index_dict:
            index_dict[value] = [i+FIRSTDATACOLUMN]
        else:
            index_dict[value].append(i+FIRSTDATACOLUMN)


    #calculate rms
    count = 0;
    rms = 0;
    errorDict = {}
    multiIndexes = np.array([ len(value) for value in index_dict.values() ])
   # if len(multiIndexes[multiIndexes > 1]) < 3:
  #      raise Exception(f"Error was expecting at least 3 duplicate points only found: {len(multiIndexes[multiIndexes > 1])}")
    for indexList in index_dict.items():
        indexList = indexList[1]
        if len(indexList) < 2:
            continue
        elif len(indexList) > 2:
            raise Exception("Cannot have more than two instances of the same time point")
        else:
            time = timeValues[indexList[0]-1]
            difference = dataFrame.iloc[:,indexList[0]] - dataFrame.iloc[:,indexList[1]]
            difference = np.array(difference)
            errorDict[time]=difference
            difference = difference*difference
            count += difference.shape[0]
            rms += np.sum(difference)
    #
    plotErrorHistogram(errorDict,outputFile)
    rms = math.sqrt(rms/count)/math.sqrt(2)
    #rms = 1.1E5/math.sqrt(2)
    #The /sqrt(2) comes from the fact that the rms is a gaussian variable with a variance that is the sum
    # of the variances of the gaussians that sampled by the original points. Since were are assuming, that the variance
    # is the same across all points we need to divide by the sqrt(2)
    rms = rms*np.ones(dataFrame.iloc[:,FIRSTDATACOLUMN:].values.shape)
    normalizedValues = normalizeDataFrame(dataFrame)
    rms /= normalizedValues[:,np.newaxis]
    return rms
#
def plotErrorHistogram(errorsDict: dict, outputFile: str):

    figText = []
    fig = plt.Figure()
    ax1 = fig.add_subplot(111)
    fig.suptitle(f"Analysis of duplicate error")
    for key in errorsDict:
        mean = np.mean(errorsDict[key])
        std = np.std(errorsDict[key],ddof=1)
        SEM = std/np.sqrt(errorsDict[key].shape[0])
        figText.append(f"TimePoint: {key}, Mean= {mean: .2E}, Std= {std: 0.2E}, SEM: {SEM: 0.2E}")
        binnum = int(2*(errorsDict[key].shape[0]**(1.0/3.0)))
        n,bins,patches = ax1.hist(errorsDict[key],bins=binnum,density=True, alpha=0.5,label=f"time: {key}")
        pmfxvals = np.linspace(min(bins),max(bins),500)
        histfit = norm.pdf(pmfxvals,loc=mean,scale=std)
        ax1.plot(pmfxvals,histfit,'r--',linewidth=2.0)

    #
    total = np.concatenate(tuple(np.array(errorsDict[key]) for key in errorsDict))
    totalMean = np.mean(total)
    totalSTD = np.std(total,ddof=1)
    totalSEM = totalSTD/np.sqrt(total.shape[0])
    figText.append(f"Total, Mean= {totalMean: .2E}, Std= {totalSTD: 0.2E}, SEM: {totalSEM: 0.2E}")
    binnum = int(2 * (total.shape[0] ** (1.0 / 3.0)))
    n, bins, patches = ax1.hist(total, bins=binnum, density=True, alpha=0.5, label=f"total")
    pmfxvals = np.linspace(min(bins), max(bins), 500)
    histfit = norm.pdf(pmfxvals, loc=totalMean, scale=totalSTD)
    ax1.plot(pmfxvals, histfit, 'r--', linewidth=2.0)
    ax1.set_xlabel('Error')
    ax1.set_ylabel('Probability Density')
    count = 0
    ax1.legend(loc='lower right')
    for text in figText:
        fig.text(0.3,0.85-count*0.05,text,fontsize=8)
        count += 1

    fig.savefig(outputFile)
#
class relaxationModel:
    def __init__(self,modelName, parameterNames,bounds,function):
        self.parameterNames = parameterNames
        self.function = function
        self.name = modelName
        self.bounds = bounds
    #
    def __call__(self,x,p):
        if len(p) != len(self.parameterNames):
            raise Exception("The number of parameters used when calling: " + self.name + " "
                            "does not equal the expected number: " + len(self.parameterNames))
        return self.function(x,p)
    #
    def numParameters(self):
        return len(self.parameterNames)
#
def fit2(x,p):
    assert(len(p) == 2)
    with warnings.catch_warnings():
        warnings.filterwarnings("error")
        try:
            return (p[0]*(np.exp((x*-1.0/p[1]))))
        except Warning as e:
            raise OverflowError(f"Overflow error in fit2: {p}, {x}")
        #


def fit3(x,p):
    assert len(p) == 3
    with warnings.catch_warnings():
        warnings.filterwarnings("error")
        try:
            return ((p[0]*np.exp(x*-1.0/p[1]))+p[2])
        except Warning as e:
            raise OverflowError(f"Overflow error in fit3: {p}, {x}")
        #

def residuals(p,y,x,u,model):
    err=((y-model(x,p))/u)
    return err

def sumSquared(p,y,x,u,model):
    return np.sum(np.square(residuals(p,y,x,u,model)))

def fitFunction(model,timeValues,ydata,rms,initialParams=None):
    if initialParams is None:
        ranges = [ (np.max(ydata),2*np.max(ydata)), (np.min(timeValues),5*np.max(timeValues))]
        if model.numParameters() == 3:
            ranges.append((-1,1))
        #
        initialParams = scipy.optimize.brute(sumSquared,ranges,args=(ydata,timeValues,rms,model),Ns=5,finish=None)
    params_opt, cov_x, infodict, msg, ier = scipy.optimize.leastsq(residuals,initialParams,args=(ydata,timeValues,rms,model),full_output=True)
    return params_opt
#
def monteCarloFit(model,initialParams,timeValues,ydata,rms,numIterations):
    parametersForFit = None
    i = 0
    while i < numIterations:
        randomizedYdata = np.random.normal(ydata,rms)
        try:
            result = fitFunction(model,timeValues,randomizedYdata,rms,initialParams)
        except OverflowError as e:
            continue
        #

        if(parametersForFit is None):
            parametersForFit = np.zeros((numIterations,result.shape[0]))
        #
        parametersForFit[i,:] = result
        i=i+1
    #
    return parametersForFit
#
def generateFitPlot(model,timeValues,ydata,rmsError,fit,siteID,destination):
    fitParameters = fit
    fitParameters = tuple(fitParameters)
    chi2 = np.sum(np.square((ydata - model(timeValues,fitParameters))/rmsError))
    normChi2 = chi2/(ydata.shape[0] - model.numParameters())

    xrange = np.linspace(0,np.max(timeValues),1000)
    yrange = model(xrange,fitParameters)
    residuals = ydata - model(timeValues,fitParameters)

    fittedParamText = []
    fittedParamText.append(f"Chi2: {chi2: .2E}")
    fittedParamText.append(f"Norm Chi2: {normChi2: .2E}")
    for i in range(model.numParameters()):
        fittedParamText.append(f"Fitted {model.parameterNames[i]}: {fitParameters[i]:.2E}")
    #
    fig = plt.Figure()
    formatter = ScalarFormatter(useOffset=False, useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((0,0))
    fig.suptitle(f" {model.name} for {str(siteID)}")
    ax1 = fig.add_axes([0.18, 0.82, 0.72, 0.08]) #ax1 is the residual plot
    ax1.plot(timeValues, residuals, 'ok', xrange, np.zeros(len(xrange)), '--r')
    ax1.errorbar(timeValues, residuals, xerr=None, yerr=rmsError, color='red', barsabove=True, fmt=',',
                 linewidth=1.5)
    ax1.yaxis.set_major_formatter(formatter)

    ax2 = fig.add_axes([0.18, 0.09, 0.72, 0.67])
    #ax2.plot(xrange, yrange, '-k', timeValues, ydata, 'ok', rx, oy, '--r', linewidth=1.5)
    ax2.plot(xrange, yrange, '-k', timeValues, ydata, 'ok')
    ax2.errorbar(timeValues, ydata, xerr=None, yerr=rmsError, color='red', barsabove=True, fmt=',',
                 linewidth=1.5)
    formatter = ScalarFormatter(useOffset=False,useMathText=True)
    formatter.set_scientific(True)
    ax2.yaxis.set_major_formatter(formatter)
    ax2.set_xlabel('Delay Time (s)')
    ax2.set_ylabel('Intensity (a.u.)')
    ax2.set_ylim(0, 1.1*max(np.max(ydata),np.max(yrange)))
    count = 0
    for text in fittedParamText:
        fig.text(0.48,0.70-count*0.05,text,fontsize=12)
        count +=1
    #"""
    os.makedirs(destination+"/"+model.name,exist_ok=True)
    fig.savefig(destination+"/"+model.name+"/"+str(siteID)+"_"+model.name+"_fitPlot.png")
#
def generateMonteCarloHistogram(parameterName, parameterValues, modelName, siteID,outDir):
    mean = np.mean(parameterValues)
    std = np.std(parameterValues)
    figText =[]
    figText.append(f"Mean= {mean: .2E}")
    figText.append(f"Std= {std: 0.2E}")
    binnum = int(2.0*(parameterValues.shape[0]**(1.0/3.0)))

    fig = plt.Figure()
    ax1 = fig.add_subplot(111)
    n,bins,patches = ax1.hist(parameterValues,bins=binnum,density=True, alpha=0.8)
    pmfxvals =  np.linspace(min(bins),max(bins),500)
    histfit = norm.pdf(pmfxvals,loc=mean,scale=std)
    ax1.plot(pmfxvals,histfit,'r--',linewidth=2.0)
    fig.suptitle(f"Monte Carlo results for parameter {parameterName}")
    ax1.set_xlabel('Fitted Value')
    ax1.set_ylabel('Probability Density')
    count = 0
    for text in figText:
        fig.text(0.68,0.85-count*0.05,text,fontsize=16)
        count += 1

    outputDirectory = f"{outDir}/{modelName}/{siteID}"
    os.makedirs(outputDirectory,exist_ok=True)
    fig.savefig(f"{outputDirectory}/{siteID}_{parameterName}_{modelName}_monteCarlo.png")


#
def plotChi2Distribution(chi2Array, dof, outputDirectory):
    chi2Array = np.array(chi2Array)
    fig = plt.Figure()
    ax1 = fig.add_subplot(111)
    lower = 0
    upper = chi2.ppf(0.999,dof)


    lowerRange = chi2Array[chi2Array <= upper]
    nbins = int(2.0*(lowerRange.shape[0]**(1.0/3.0)))

    # Compute the histogram using numpy.histogram
    counts, bin_edges = np.histogram(lowerRange, bins=nbins)
    bin_edges = numpy.append(bin_edges,np.max(chi2Array)+1)
    # Plot the histogram
    ax1.hist(chi2Array,bins=bin_edges,density=True)
    pmfxvals = np.linspace(lower, upper, 500)
    histfit = chi2.pdf(pmfxvals, dof)
    ax1.plot(pmfxvals, histfit, 'r--', linewidth=2.0)

    ax1.set_xlim((lower,upper))
    ax1.set_xlabel('Chi2 Value')
    ax1.set_ylabel('Probability Density')

    os.makedirs(outputDirectory,exist_ok=True)
    fig.savefig(f"{outputDirectory}/fit_Chi2_distribution.png")
#
def processResidueModel(model,site,timeValues,ydata,rmsError,outDir):
    fit = fitFunction(model, timeValues, ydata, rmsError)
    generateFitPlot(model, timeValues, ydata, rmsError, fit, site, outDir)
    parameterArray = monteCarloFit(model, fit, timeValues, dataDictionary[site], rmsError, 500)
    idx = 0
    for parameter in model.parameterNames:
        generateMonteCarloHistogram(parameter,parameterArray[:,idx], model.name, site, outDir)
        idx+=1
    #
    chi2 = sumSquared(fit,ydata,timeValues,rmsError,model)
    BIC = chi2 + model.numParameters()*math.log(ydata.shape[0])
    return BIC, chi2, np.mean(parameterArray,axis=0), np.std(parameterArray,axis=0)


twoParameterFit = relaxationModel("TwoParameterFit",
                        ["Intensity", "Tau"],[(1E-9,1E100),(1E-9,1E9)],fit2)
threeParameterFit = relaxationModel("ThreeParameterFit",
                                  ["Intensity", "Tau", "Baseline"], [(1E-9,1E100),(1E-9,1E9),(-1E100,1E100)] ,fit3)
def processResidue(site,timeValues,ydata,rmsError,outDir,ThreeParamFitFlag=False):
    twoParamResult = processResidueModel(twoParameterFit,site,timeValues,ydata,rmsError,outDir)
    best = twoParamResult
    if ThreeParamFitFlag:
        threeParamResult = processResidueModel(threeParameterFit,site,timeValues,ydata,rmsError,outDir)
        print(f"Site {site}: twoParameterFit BIC {twoParamResult[0]}, threeParameterFit {threeParamResult[0]}")

        if(twoParamResult[0] < threeParamResult[0]):
            print("Chose two parameter")
        else:
            best = threeParamResult
            print("chose three parameter")
        #
    #
    return best
#

splitCharacter = r'[ \t]+|,'

print(f"version number {VERSION_NUMBER}")
if(len(sys.argv) != 4):
    print("Error: syntax: \n python Masfit.py inputFileName outputDirectory ThreeParamFit?(True or False)")
    exit(1)

inputFileName = sys.argv[1]
outputDirectory = sys.argv[2]
threeParameterFitFlag = sys.argv[3]
if(threeParameterFitFlag != "True" and threeParameterFitFlag != "False"):
    raise ("Specify True or False if you want to try Three parameter fits")
if(threeParameterFitFlag=="True"):
    threeParameterFitFlag=True
else:
    threeParameterFitFlag=False

if os.path.exists(outputDirectory):
    print(f" Error: Directory {outputDirectory} exists, delete or choose new directory")
    exit(1)

os.makedirs(outputDirectory)
with open(f"{outputDirectory}/runParameters.txt", 'w') as parameterFile:
    print(VERSION_NUMBER,sys.argv,file=parameterFile)

with open(inputFileName, 'r') as inputFile:
    df = pd.read_csv(inputFile, sep=splitCharacter, engine='python', skiprows=1, header=None)

with open(inputFileName, 'r') as inputFile:
    header = re.split(splitCharacter,inputFile.readline().strip())

timeValues = np.array([ float(value) for value in header[FIRSTDATACOLUMN:] ])
print(timeValues)
df.iloc[:,FIRSTDATACOLUMN:] = df.iloc[:,FIRSTDATACOLUMN:].astype(float)
print(df)

rmsError = calculateNoise(timeValues,df,f"{outputDirectory}/errorHistogram.png")

keys = df.iloc[:,0]
values = df.iloc[:,FIRSTDATACOLUMN:].apply(np.array,axis = 1)


dataDictionary = dict(zip(keys,values))



profiler = cProfile.Profile()
profiler.enable()
count=0
with open(f"{outputDirectory}/fittedTaus.txt",'w') as outputFile:
    print("Site\tTau\tTau_err",file=outputFile)
    chi2List=[]
    for site in dataDictionary:
        fittedParameters = processResidue(site,timeValues,dataDictionary[site],rmsError[count],outputDirectory,threeParameterFitFlag)
        print(f"{site}\t{fittedParameters[2][1]: 0.2E}\t{fittedParameters[3][1]:0.2E}",file=outputFile)
        count+=1
        print(f" completed {count} out of {len(dataDictionary.keys())}")
        chi2List.append(fittedParameters[1])

    #
    if not threeParameterFitFlag:
        plotChi2Distribution(chi2List, timeValues.shape[0] - twoParameterFit.numParameters(), outputDirectory)
#
profiler.disable()
profiler.dump_stats("output.prof")



#








