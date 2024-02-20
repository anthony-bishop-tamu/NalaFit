#MASfit v2.0.0
#Goal of this script is to improve the generality of the fitting software and work with python 3.0+
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
def fit2(x,p):
    return (p[0]*(scipy.exp((scipy.array(x)*-1.0/p[1]))))

def fit3(x,p):
    return ((p[0]*scipy.exp(scipy.array(x)*-1.0/p[1]))+p[2])

def residuals(p,y,x,u,model):
    err=((y-model(x,p))/u)
    return err

def fitFunction(model,timeValues,ydata,rms):
    initialParams = np.zeros(2)
    if model == fit2:
        initialParams = np.zeros(2)
    elif model == fit3:
        initalParams = np.zeros(3)
    else:
        assert False

    initialParams[0] = np.max(ydata)
    initialParams[1] = np.mean(timeValues)
    result = scipy.optimize.minimize(residuals,initalParams,args=(ydata,timeValues,rms,model))
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
#

#
splitCharacter = r'[ \t]+|,'
firstDataColumnIndex = 2

with open('T2600_4fitting.txt', 'r') as inputFile:
    df = pd.read_csv(inputFile, sep=splitCharacter, engine='python', skiprows=1, header=None)

with open('T2600_4fitting.txt', 'r') as inputFile:
    header = re.split(splitCharacter,inputFile.readline().strip())

timeValues = [ float(value) for value in header[3:] ]
print(timeValues)
df.iloc[:,firstDataColumnIndex:] = df.iloc[:,firstDataColumnIndex:].astype(float)
df = denormalizeDataFrame(df)
print(calculateNoise(timeValues,df,firstDataColumnIndex))





