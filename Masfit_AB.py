#MASfit v2.0.0
#Goal of this script is to improve the generality of the fitting software and work with python 3.0+
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import scipy.optimize
import scipy.optimize as opt
import pandas as pd
import re
import math
import numpy as np

def denormalizeDataFrame(dataFrame):
    result = dataFrame.iloc[:,[2]].values*dataFrame.iloc[:,3:].values
    dataFrame.iloc[:,3:] = pd.DataFrame(result)
    dataFrame = dataFrame.drop(columns=dataFrame.columns[2])
    return dataFrame
#

def calculateNoise(timeValues,dataFrame,firstCDataColumnIndex):
    index_dict = {}
    #determine duplicates
    for i, value in enumerate(timeValues):
        if value not in index_dict:
            index_dict[value] = [i+firstCDataColumnIndex]
        else:
            index_dict[value].append(i+firstCDataColumnIndex)


    #calculate rms
    count = 0;
    rms = 0;
    for indexList in index_dict.items():
        indexList = indexList[1]
        if len(indexList) < 2:
            continue
        elif len(indexList) > 2:
            raise Exception("Cannot have more than two instances of the same time point")
        else:
            difference = dataFrame.iloc[:,indexList[0]] - dataFrame.iloc[:,indexList[1]]
            difference = np.array(difference)
            difference = difference*difference
            count += difference.shape[0]
            rms += np.sum(difference)
    #
    return math.sqrt(rms/count)
#
class relaxationModel:
    def __init__(self,modelName, parameterNames,function):
        self.parameterNames = parameterNames
        self.function = function
        self.name = modelName
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
    return (p[0]*(np.exp((x*-1.0/p[1]))))

def fit3(x,p):
    assert len(p) == 3
    return ((p[0]*np.exp(x*-1.0/p[1]))+p[2])

def residuals(p,y,x,u,model):
    err=((y-model(x,p))/u)**2
    return np.sum(err)

def fitFunction(model,timeValues,ydata,rms,initialParams):

    initialParams[0] = np.max(ydata)
    initialParams[1] = np.mean(timeValues)
    result = scipy.optimize.minimize(residuals,initialParams,args=(ydata,timeValues,rms,model))
    return result
#
def monteCarloFit(model,timeValues,ydata,rms,numIterations):
    parametersForFit = None
    for i in range(numIterations):
        result = fitFunction(model,timeValues,ydata,rms)
        if(parametersForFit is None):
            parametersForFit = np.zeros((numIterations,result.x.shape[0]))
        #
        parametersForFit[i,:] = result.x
    #
    return parametersForFit;
#
def generateFitPlot(model,timeValues,ydata,rmsError,fit,siteID,destination):
    fitParameters = fit.x
    chi2 = fit.fun
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
    fig.savefig(destination+"/"+str(siteID)+"_"+model.name+"_fitPlot.png")

#
twoParameterFit = relaxationModel("TwoParameterFit",
                                  ["Intensity", "Time constant"], fit2)
threeParameterFit = relaxationModel("ThreeParameterFit",
                                  ["Intensity", "Time constant", "Baseline"], fit3);

splitCharacter = r'[ \t]+|,'
firstDataColumnIndex = 2
inputFileName = 'T2600_4fitting.txt'
with open(inputFileName, 'r') as inputFile:
    df = pd.read_csv(inputFile, sep=splitCharacter, engine='python', skiprows=1, header=None)

with open('T2600_4fitting.txt', 'r') as inputFile:
    header = re.split(splitCharacter,inputFile.readline().strip())

timeValues = np.array([ float(value) for value in header[3:] ])
print(timeValues)
df.iloc[:,firstDataColumnIndex:] = df.iloc[:,firstDataColumnIndex:].astype(float)
df = denormalizeDataFrame(df)

keys = df.iloc[:,1]
values = df.iloc[:,firstDataColumnIndex:].apply(np.array,axis = 1)

dataDictionary = dict(zip(keys,values))
rmsError = calculateNoise(timeValues,df,firstDataColumnIndex)
rmsError = np.ones(len(timeValues))*rmsError

for site in dataDictionary:
    initalParamValues = [ np.max(dataDictionary[site]) , np.max(timeValues) ]
    fit = fitFunction(twoParameterFit,timeValues,dataDictionary[site],rmsError,initalParamValues)
    generateFitPlot(twoParameterFit,timeValues,dataDictionary[site],rmsError,fit,site,".")
#







