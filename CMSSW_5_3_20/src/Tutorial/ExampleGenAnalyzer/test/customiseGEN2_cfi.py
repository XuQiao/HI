import FWCore.ParameterSet.Config as cms

def getPy6ProcParameters(genType):
    para = {
        # NSD (beyon pt>5~10 GeV/c, it's safe to say NSD = MB, but below, it's not)
        'NSD':['MSEL=1'], 
        'NSD_no_ISR':['MSEL=1','PARP(67)=0'], 
        # MB
        'MB':['MSEL=0','MSUB(11)=1','MSUB(12)=1','MSUB(13)=1','MSUB(28)=1','MSUB(53)=1','MSUB(68)=1','MSUB(92)=1','MSUB(93)=1','MSUB(94)=1','MSUB(95)=1'],
        'HMB':['MSEL=0','MSUB(11)=1','MSUB(12)=1','MSUB(13)=1','MSUB(28)=1','MSUB(53)=1','MSUB(68)=1','MSUB(92)=1','MSUB(93)=1','MSUB(94)=1'],
        }
    print "PYTHIA process parameters: ",para[genType]
    return para[genType]

def getPy6PthatLow(pthatLow):
    print "Minimum pt-hat: "+str(pthatLow)
    return ['CKIN(3)='+str(pthatLow)] 

def getPy6PthatHigh(pthatHigh):
    print "Minimum pt-hat: "+str(pthatHigh)
    return ['CKIN(4)='+str(pthatHigh)]


def updatePy6ProcParameters(gen,genType,low,high,sqrtS):
    print "Center-of-mass energy: ",sqrtS
    gen.PythiaParameters.processParameters = getPy6ProcParameters(genType)
    gen.PythiaParameters.processParameters += getPy6PthatLow(low)
    gen.PythiaParameters.processParameters += getPy6PthatHigh(high)
    gen.comEnergy = sqrtS
