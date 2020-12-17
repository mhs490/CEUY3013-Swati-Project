#!/usr/bin/python
import pandas as pd
import numpy as np
import argparse
from os import path
import matplotlib.pyplot as plt
import random
from argparse import RawTextHelpFormatter
from pandas.core.frame import DataFrame


def disinfectionSys(data, virus, rate, flow, decayRate, order, count): #defining the inactivation or disinfection system
    """
HRT, required volume for both CMFR and PFR for both orders are calculated and plotted in the output image here. 
the values of  HRT and required volume are calculated once for zero order since they are the same. 
while for the first order reaction they are different since calculated for PFR and CMFR from two different formual.
    """
    plt.text(-8.3, -1, 'CMFR: ', weight='bold') #writing CMFR heading on the output image.
    plt.text(0, -1, 'PFR: ', weight='bold') #writing PFR heading on the output image.

    virusData = data[virus] #getting the orginal virus data that are been below in the code. 
    time = 0
    if order == 'zeroOrder': 
        C0 = virusData.iloc[0] #getting the first concentration value from  orgignal data. 
        kInv = (1 / decayRate) #defining the 1/K which is the inverse of decay rate used below in calculations.
        rate = (100 - rate) / 100 # converting desired percentage inactivation rate to a constant value.

        time = kInv * ((1 - rate) * C0) #calculating the hydraulic residence time (HTR) for the virus.
        plt.text(0, -2, 'Values same as CMFR') #for zero order the value of HTR for CMFR and PFR are the same.
        
        plt.text(-10, -10, 'Any of the reactor can be chosen since the required vol and HTR is the same for both.')

    elif order == 'firstOrder': #if not zero then first order.
        ctbyc0 = 1 / ((100 - rate)/100) #calculatiing Ct/Co
        kInv = 1 / decayRate #defining the 1/K which is the inverse of decay rate used below in calculations.
        ctbyc0 = ctbyc0 - 1 

        time = int(kInv * ctbyc0) #calculating the HRT for first order reactions. 
        plt.text(-8.3, -1.8-((count+1.5)*count), 
                 'Hydralic Residence Time of ' + str(virus)+':') #wrting the HRT heading for CMFR
        plt.text(-4.2, -1.8-((count+1.5)*count),
                 str(time) + ' hours', weight='bold') #writing the HRT calculated value for CMFR

        timePFR = round(-(kInv) * (np.log(((100 - rate)/100))), 2) #calculating the HRT for PFR when first order reaction.
        volPFR = flow * (10 ** 6) * timePFR * 0.0001577088 #calculating the required volume for PFR. 

        plt.text(0, -2.2-((count+1.5)*count),
                 'Reactor Volume of ' + str(virus) + ':') #writing the heading for required volume for PFR.
        plt.text(3, -2.2-((count+1.5)*count),
                 "{:e}".format(volPFR)+' m^3', weight='bold') #writing the calculated volume required for PFR

        plt.text(0, -3.5-((count+1.5)*count),
                 'Hydralic Residence Time of ' + str(virus)+':') #wrting the HRT heading for PFR.
        plt.text(4, -3.5-((count+1.5)*count),
                 str(timePFR) + ' hours', weight='bold')  #writing the HRT calculated value for PFR.
        
        plt.text(-5, 16-count, str(virus) + ': The chosen reactor should be the one with lower HTR and required volume')

    volume = flow * (10 ** 6) * time * 0.0001577088  #calculating the required volume CMFR. 

    plt.text(-8.3, -2.6-((count+2)*count),
             'Reactor Volume of ' + str(virus)+':') #writing the heading for required volume CMFR.
    plt.text(-4.2, -2.6-((count+2)*count),
             "{:e}".format(volume)+' m^3', weight='bold') #writing the calculated volume required for CMFR.
    


def firstOrderSlope(data, rate, flow, count): #calculation done in order to find if the data is first order
    """
    here the data of the LN(concentraton) and time is checked against the condition of fist order reaction.
    the system will check for first order only if the zero order requirments are not met.
    """
    newData = data.copy() #copying the data from orignal data in the input file so that the orignal data remains untouched.
    newData.iloc[:, 1:] = np.log(newData.iloc[:, 1:]) #finding the natural log of virus concentration from the copied data. 

    diff = newData.diff(periods=1) # calculating the difference i.e y2-y1 and x2-x1 of graph of ln(virus concentration) vs time. 
    diff.dropna(inplace=True) #droping empyty rows.

    viruses = newData.iloc[:, 1:] #extracting the virus concentration from the input copy.
   
    for virus in viruses: #going thorugh the copy ln(virus conc) data
    
        slope = DataFrame() #defining new object that represents slope.
        slope['slope ' + str(virus)] = diff[virus].div(diff.iloc[:, 0])  # dividing the difference calculated above to find the slope of ln(conc vs time). 

        slopeDev = round(slope.std()[0], 2) #calculating the standard deviation of the slopes calculated above. 

        slopeMean = round(slope.mean()[0], 2) #calculating the mean of the slopes calculated above.

        lowerCutoff = slopeMean - slopeDev #figuring out the outliers of the slope
        upperCutOff = slopeMean + slopeDev #figuring out the outliers of the slope

        currOrder = False  #checking if the virus has first order reaction 
        if slope.iloc[2][0] > lowerCutoff and slope.iloc[2][0] < upperCutOff: #difining the condition for first order reaction.
            currOrder = True

        plt.text(-8.3, 2-((count+1)*count),
                 'Reaction order of ' + str(virus)+':')  #creating heading for reaction order. 
        
        if currOrder: #for first order
            decayRate = round(abs(slope['slope '+str(virus)].iloc[0]), 2) #calculatiing decay order
            
            plt.text(-4.5, 2-((count+1)*count),
                     'First-order decay', weight='bold') # writing the calculated reaction order
            plt.text(0, 2-((count+1)*count), 'K ' +
                     str(virus)+': ' + str(decayRate)) # writing the calculated decay rate
            
            disinfectionSys(data, virus, rate, flow,
                            decayRate, 'firstOrder', count) #calculating the inactivation or disinfection rate for the system
        else:
            print("We do not handle second order") 
        count += 1


def calcOrder(data, rate, flow):#defining the function which will be responsible for checking the reaction order of the viruses.
    """
here the data of the concentraton and time is checked against the condition of zero order reaction.
if the data fulfil the zero order reaction condition it will not go for first order conditions.
if the data is not zero order reaction it will go to first order conditions.
    """
    newData = data.copy() #copying the data from orignal data in the input file so that the orignle data remains untouch.
    diff = newData.diff(periods=1)# calculating the difference i.e y2-y1 and x2-x1 of graph of virus concentration vs time. 
    diff.dropna(inplace=True)# dropping an empty row. 

    viruses = newData.iloc[:, 1:] #extracting the virus concentration from the input copy.
    count = 0 
    for virus in viruses: #going thorugh the copy virus data
        slope = DataFrame()  #defining new object that represents slope. 
        slope['slope ' + str(virus)] = diff[virus].div(diff.iloc[:, 0]) # dividing the difference calculated above to fine the slope. 

        slopeDev = round(slope.std()[0], 2) #calculating the standard deviation of the slopes calculated above. 

        slopeMean = round(slope.mean()[0], 2)#calculating the mean of the slopes calculated above.

        lowerCutoff = slopeMean - slopeDev #figuring out the outliers of the slope
        upperCutOff = slopeMean + slopeDev #figuring out the outliers of the slope

        currOrder = False #checking if the virus has zero order reaction 
        if slope.iloc[0][0] > lowerCutoff and slope.iloc[0][0] < upperCutOff: #difining the condition for zero order reaction.
            currOrder = True
        plt.text(-8.3, 2-((count+1)*count),
                 'Reaction order of ' + str(virus)+':')  #creating heading for reaction order. 
        plt.text(0, 2-((count+1)*count), 'K ' + str(virus)+':') #creating heading for decay rate

        if currOrder: #for zero order 
            decayRate = round(abs(slope['slope '+str(virus)].iloc[0]), 2) #calculatiing decay order
        
            plt.text(-4.5, 2-((count+1)*count), 
                     'Zero-order decay', weight='bold') # writing the calculated reaction order
            plt.text(0, 2-((count+1)*count), 'K ' +
                     str(virus)+': ' + str(decayRate)) # writing the calculated decay rate 

            disinfectionSys(data, virus, rate, flow,
                            decayRate, 'zeroOrder', count) #calculating the inactivation rate for the system
        else:
            firstOrderSlope(data, rate, flow, count) #if the data is not zero order the system will move to calculation for first order reaction.
            break
        count += 1


def plotGraphs(data):# starting the plotting of the data from input.
    """
the data from input file is extracted in steps, first time then virus concentration. 
after that the data is plotted in two different graph one for concentratioin vs time and the other for ln(conc) vs time. 
    """
    timePeriod = data.iloc[:, 0] #extracting time period columns from the inpiut file.
    viruses = data.iloc[:, 1:] #extracting the virus columns from the input file 
    markers = ['^', '.', '+', 'x', 'D', 's'] #deciding the markers for the graph ploted from the columns above.

    plt.figure(figsize=(10, 8)) #deciding the size of the figure
    plt.subplot(221) #ploting virus concentration vs time graph.
    for virus in viruses:#ploting the graphs of all given viruses in the input folder. 
        plt.plot(timePeriod, viruses[virus],#ploting concentration vs time graphs
                 random.choice(markers), label=str(virus))

    plt.xlabel('Time (hours)') #labeling the x axis
    plt.ylabel('Virus Concentration') #labeling the y axis
    plt.legend(loc='best') #placing the legend

    # natual log plot
    plt.subplot(222) #plotting the ln(virus concentration) vs time graph
    for virus in viruses:
        plt.plot(timePeriod, np.log(viruses[virus]),
                 random.choice(markers)+'-', label=str(virus), linewidth=1.0)#deciding the marker and labeling the virus on graph.

    plt.xlabel('Time (hours)')#labeling the x axis for ln graph 
    plt.ylabel('ln (Virus Concentration)')#labeling the y axis for ln graph
    plt.legend(loc='best')#placing the legend


def main(filePath, rate, flow): 
    """ 
    checking if the desired input file in suggested is present in the input folder.
    if thats not the case it will give error to the user.
    calculating the decay rate, HTR and flowrate which are required to get to inactivation rate of the virus.
    after that the answers are saved in output folder in .png format

    """
    
    FinalFilePath = 'input/' + filePath + '.xlsx' # all input files are assumed to be in the input folder.
    
    if not path.exists(FinalFilePath): # checking if the file is present.if not it will print error 
        print('File does not exist')
        exit()

    
    data = pd.read_excel(FinalFilePath)  # extracting data from excel to data function.
    plotGraphs(data) # ploting the data from data function. and checking if the given data is zero order or firdt order.
    calcOrder(data, rate, flow) # calculating the decay rate, HTR and flowrate.
    

    plt.savefig('output/' + filePath + '.png')# saving the output result in image.png format in the output folder.

    print('Output written to: ', 'output/' + filePath + '.png')# printing the location of output image. 

    plt.show()# opening the output image file


if __name__ == "__main__": #Calling function with argument from command line
    parser = argparse.ArgumentParser(description=
'''
Virus file should be an excel file in the following format:
        Time (hours) | Virus 1 
        0            | 1 x 10^7 
        1            | 1 x 10^6 
        2            | 1 x 10^4 
        3            | 1 x 10^2 
        4            | 1 x 10^2 
''',
        formatter_class=RawTextHelpFormatter
        )

    parser.add_argument('-i', '--input', type=str, required=True,
                        help='(str) Add the path to the data file in the input directory without extension.'  )

    parser.add_argument('-r', '--rate', type=float, required=True,
                        help='(float) Provide inactivation rate percentage upto 3 precison, e.g: 99.9')

    parser.add_argument('-f', '--flow', type=int, required=True,
                        help='(int) Provide flowrate in MGD (million gallons per day) e.g: 1 will be taken as 1 (MGD)')

    args = parser.parse_args() # Resolving the above arguments
    main(args.input, args.rate, args.flow) # Calling main function
