#CARO: run as python outputfolder datasetname datatype &
#python3 Masfit_pipe_v1.py ca08 dz 7
#for Taylor python Masfit.py
#masfit v11
#Corrected single exponential relaxation fitting script
#Corrected Rate to Time Constant
#by Matt Stetz
#
###Usage
#This script is used to fit single or bi-exponetial relaxation data derived 
#from my nmrPipe processing scripts (NUS or cartesian)
#If you want to use Felix processed data, run my script felix2pipe.py first
#
#I'm sorry for the use of global variables but it is easier for users
#who may be unfamiliar with python to edit.
#
###References
#For determining peak height uncertainties and errors in estimators:
#Palmer AG 3rd, Rance M, and Wright PE JACS 113 (12) 1991 
#Skelton NJ, Palmer AG 3rd, Akke M, Kordel J, Rance M, and Chazin WJ JMR B 102 1993
#Fushman D, Cahill S, and Cowburn D JMB 266(1) 1997 
#Keith G. Inman (Weber lab) Ph.D. Dissertation 2003
#
#For best practices with respect to model selection:
#Viles JH, Wright Pe et al. JBNMR 21 2001
#
#Bootstrapping method for obtaining errors in fitted parameters:
#Kamath U and Shrived JW JBC 264(10) 1989 (purely a resampling method; no replacement)
#Farrow NA, Kay LE et al. Biochemistry 33 1994
#Press WH et al. 'Numerical Recipes 3rd Edition: The Art of Scientific Computing'
#Rule GS and Hitchens TK 'Fundamentals of Protein NMR Spectroscopy'
#
#Basic statistics:
#Taylor JR 'An Introduction to Error Analysis The Study of Uncertainties in Physical Measurements'

import numpy
import numpy.random

import pylab as plt
import matplotlib.mlab as mlab 

import scipy
import scipy.optimize

import os
import errno

import sys

###Users must edit the information below###

#Specify path to input data (should be .txt file) below
#here=os.getcwd()
#handle=sys.argv[1]+sys.argv[3]+sys.argv[2]
datapath='/home/nmrbox/tcole/nmrdata/FLParkin_backbonedynamics/FLParkin_backboneT1T2ratecalc/T2600_4fitting.txt'

#excludepath=here+'/'+handle+'/'+handle+'_forfitting.txt.excluded'
#if os.path.isfile(excludepath):
#	excludeflag=1
#else:
#	excludeflag=0
#	print excludepath
#	print 'NO EXCLUDED RESIDUES LIST EXISTS! \n'

#Specify path to output directory below
#here=os.getcwd()
#try:
#	os.makedirs(here+'/'+handle)
#	#avoid race condition
#except OSError:
#	pass
	#if not os.path.isdir(here):
	#	raise
#return 0
#os.makedirs(here+'/'+sys.argv[2])
#outpath=here+'/'+handle+'/'
outpath='/home/nmrbox/tcole/nmrdata/FLParkin_backbonedynamics/FLParkin_backboneT1T2ratecalc/T2_600/'

#Define tau values as list
taulist=0


#Do you want to do a 2-parameter fit? (2p_choice='y' for yes, ='n' for no)
twop_choice='y'

#Specify folder name for 2 parameter fit plots
twopfolder='2_param_fit_plots'

#Specify file name for 2 parameter fit
ext2='2param'

#Do you want to do bootstrapping for the 2-parameter fits (twoboot='y' for yes, ='n' for no)?
twoboot='y'

#How many simulations for bootstrapping 2-parameter fits? (500-1000 recommended)
twobootit=1000

#Specify folder name for 2 parameter bootstrap plots
twobfolder='2_param_boot_plots'

#Specify file name for 2 parameter fit bootstrapping results
ext4='bootstrap_2param'

#Do you want to do a 3-parameter fit? (3p_choice='y' for yes, ='n' for no)
threep_choice='n'

#Specify folder name for 3 parameter fit plots
threepfolder='3_param_fit_plots'

#Specify file name for 3 parameter fit
ext3='3param'

#Do you want to do bootstrapping for the 3-parameter fits (threeboot='y' for yes, ='n' for no)?
threeboot='n'

#How many simulations for bootstrapping 3-parameter fits? (500-1000 recommended)
threebootit=500

#Specify folder name for 3 parameter bootstrap plots
threebfolder='3_param_boot_plots'

#Specify file name for 3 parameter fit bootstrapping results
ext5='bootstrap_3param'

#Would you like to do a bi-exponential fit?
bi_choice='n'

#Specify folder name for bi-exponential fit plots
bifolder='biexp'

#Specify file name for bi-exponential fit
extb='biexp'

#Use normalized data heights? (normalized='y' for yes, ='n' for no)
normalized='y'

#Delay times used in experiment. List in ascending order and whatever units you want
tau=[0.008,0.012,0.016,0.020,0.024,0.028,0.032]
#taulen=len(tau)+1

#tau=[float(i) for i in tau]
#tau.sort()                         #SORTS TAU. IF DELAYS AREN'T KOSHER, THIS WILL BUG OUT
########### 
#if not 9 delay points, delete the delay here. also omit stderror on line 358 to match array size.
tau=tau[:-1]
###########
taulen=len(tau)+1
print tau

#tau=[0.016, 0.048, 0.08, 0.112, 0.08, 0.145, 0.177, 0.225, 0.257]
#0.016064,0.032128,0.048192,0.064256,0.080320,0.096384,0.112448,0.128512,0.144576]

#Which data POINTS were duplicated (first delay=1, 8th delay=8, etc)? List in ascending order.
taudup=[1,4,6]

#Specify maximum value of x-axis for plot figures (Only used to make plots nice looking, ~2*max delay time looks good)
maxx=tau[-1]*2

#Specify scaling factor for sigma (uncertainties in peak heights)
#Start with 1.0 but change if you get an error message
sigmascale=1

#Specfiy name of output text file containing a summary of all fitted rates, errors, and stats
txtname='FLParkin_t2600'

#### single exponential fitting below ####

#Initial guess for intensity (=1.0 if using normalized data heights)
io=1.0

#Initial guess for rate
r=(tau[4])
#r=1.2
#Initial guess for offset (used in three parameter fits only)
asym=0.01

#### bi exponential fitting guess below for aqueous proton t1rho ####

##Initial guess for intensity 1 for bi-expfit (=1.0 if using normalized data heights)
io1=1.0

##Initial guess for intensity 2 for bi-expfit (=1.0 if using normalized data heights)
io2=1.0

#Initial guess for rate 1 (sec) for bi-expfit
r1=25

#Initial guess for rate 2 (sec) for bi-expfit
r2=2500.0

##!!!No more info from user needed!!!###

###Define Functions to Organize Input###

#function to read formatted input data
def read_file(datafile):
        data=file(datafile,'r')
        contents=data.readlines()
        data.close()
        return contents

#build list of excluded residues due to bad fits
#if excludeflag==1:
#	excludedata=read_file(excludepath)
#	excludelist=[]
#	for i in excludedata:
#		if len(i)>0:
#			i=i[0:-1]
#			excludelist.append(i)
#else:
#	excludelist=[0]


#function to extract relevant data from input
def parse_data(data):
        final=[]
        delays=[]
        #check which data height values to use
        if normalized == 'y' or normalized == 'Y':
            i=3
        elif normalized == 'n' or normalized == 'N':
            i=2
            print '\nWARNING: You have opted to use absolute data heights.'
            print 'nmrPipe data will experience a loss in numerical precision!'
        else:
            print 'ERROR: Please specify to use either normalized or absolute data heights.'
            return 0
        #go through data and parse
        for x in range(len(data)):
                split=data[x].split()
#		if split[1] in excludelist:
#			print split[1], 'excluded'
#			pass
#		else:
		delays.append(str(split[1]))
                while i < len(split):
                	delays.append(float(split[i]))
                	i+=1
                	#reset loop counter (probably a more elegant way to do this)
                	
		if normalized == 'y' or normalized == 'Y':
			i=3
                elif normalized == 'n' or normalized == 'N':
                	i=2
                else:
                	pass
		Ampl=float(split[2])
		Ampl='%.3e'%Ampl
		delays.append(Ampl)
		final.append(delays)
		delays=[]
	print '\n*** *** *** *** *** ***\n'
	return final

#function to get absolute data heights--is not called if using normalized data heights
#Note: since absolute data heights are back calculated from the normalized output of nmrPipe
#there is a loss in numerical precision
def get_dataheight(data):
        temp=[]
        heights=[]
        for x in range(len(data)):
                i=2
                temp.append(data[x][0])
                while i < len(data[x]):
                        temp.append(data[x][i]*data[x][1])
                        i+=1
                heights.append(temp)
                temp=[]
        return heights


##produce a list of outliers
def printoutliers(data):
        output=open(outpath+'uncertainty_outliers.txt','w')
        for x in range(len(data)):
                for c in range(len(data[x])):
                        output.write("%-20s"%str(data[x][c]))
                        #output.write('\t')
                output.write('\n')
        print '\nuncertainty_outliers summary file saved.\n'
        output.close()
        return 0

#function to get uncertainty values from duplicate spectra
#this is done as described in the Skelnar paper (see references above)
#this function makes no assumptions about the total number of experimental delay times collected
#but it does assume that three duplicate times were collected
def get_sigma(data,c):
        cutoff=0.05
        temp1,temp2,temp3,temp4=[],[],[],[]
        outliers=[]
        final_outliers=[]
        total=len(taudup)+len(tau)
        #calculate differences between reference and duplicate peak heights
        #then take standard deviation of all the differences and scale by root 2
        for x in range(len(data)):
		diff1=data[x][taudup[0]]-data[x][total-2]
                if abs(diff1) >= cutoff:
                        outliers.append(data[x][0])
                        outliers.append(diff1)
                        outliers.append('duplicate 1')
                        final_outliers.append(outliers)
                        outliers=[]
                diff2=data[x][taudup[1]]-data[x][total-1]
                if abs(diff2) >= cutoff:
                        outliers.append(data[x][0])
                        outliers.append(diff2)
                        outliers.append('duplicate 2')
                        final_outliers.append(outliers)
                        outliers=[]
                diff3=data[x][taudup[2]]-data[x][total]
                if abs(diff3) >= cutoff:
                        outliers.append(data[x][0])
                        outliers.append(diff3)
			outliers.append('duplicate 3')
                        final_outliers.append(outliers)
                        outliers=[]
                temp1.append(diff1)
                temp2.append(diff2)
                temp3.append(diff3)
        temp4.append(temp1)
        temp4.append(temp2)
        temp4.append(temp3)
        printoutliers(final_outliers)
        #save histograms of differences between ref and duplicate data
        i=0
        while i < len(taudup):
               print i
               meanplot=numpy.mean(temp4[i])
               stdvplot=numpy.std(temp4[i],ddof=1)
               binnum=int(2.0*(len(temp4[i])**(1.0/3.0)))
               n,bins,patches=plt.hist(temp4[i],binnum,normed=1, alpha=0.8)
               pdfsmooth=numpy.linspace(min(bins),max(bins),500)
               hisfit=mlab.normpdf(pdfsmooth,meanplot,stdvplot)
               plt.plot(pdfsmooth,hisfit,'r--',linewidth=2.0)
               plt.title('Distribution of Differences in Peak Heights for Duplicate '+str(i+1))
               plt.xlabel('Difference in Peak Heights (Reference - Duplicate)')
               plt.ylabel('Probability Density')
               plt.savefig(outpath+'/'+'uncertainty_'+str(i+1)+'.png')
               plt.clf()
	       #CARO
	       plt.close('all')
               i+=1               
        std1=(numpy.std(temp1, ddof=1))/numpy.sqrt(2.0)*c
        std2=(numpy.std(temp2, ddof=1))/numpy.sqrt(2.0)*c
        std3=(numpy.std(temp3, ddof=1))/numpy.sqrt(2.0)*c
        #if you don't want to use interpolated uncertainties, 'ulist' below will contain
        #promoted uncertainties (using the same uncertainty for different delays)
        temp=[]
        ulist=[]
        for x in range(len(data)):
                temp.append(data[x][0])
                temp.append(std1)
                temp.append(std1)
                temp.append(std1)
                temp.append(std2)
                temp.append(std2)
                temp.append(std2)
#                temp.append(std3)
#                temp.append(std3)
#                temp.append(std3)     ###MUTE THIS IF DELAY #9 IS ABSENT
                ulist.append(temp)
                temp=[]
        return std1,std2,std3,ulist





###Fitting Function Definitions###

#Two-parameter fit, keep as global
p2=[io,r]

#two parameter single exponential decay function with parameters in "p2" vector
def fit2(x,p):
    return (p[0]*(scipy.exp((scipy.array(x)*-1.0/p[1]))))

#residuals scaled by uncertainties for minimization (will be squared)
def residuals2(p,y,x,u):
    err=((y-fit2(x,p))/u)
    return err

#Three-parameter fit, keep as global
p3=[io,r,asym]

#three parameter single exponential decay function with parameters in "p3" vector
def fit3(x,p):
    return ((p[0]*scipy.exp(scipy.array(x)*-1.0/p[1]))+p[2])

#residuals scaled by uncertainties for minimization (will be squared)
def residuals3(p,y,x,u):
    err=((y-fit3(x,p))/u)
    return err

####

#bi-exponential 2-parameter fit, keep as global
pb=[io1, io2, r1, r2]

#bi-exponential 2-parameter fit with parameters in "pb"
def fitbi(x,p):
        return (p[0]*scipy.exp(scipy.array(x)*-p[1]))+(p[2]*scipy.exp(scipy.array(x)*-p[3]))

#residuals scaled by uncertainties for minimization (will be squared)
def residualsbi(p,y,x,u):
        err=((y-fitbi(x,p))/u)
        return err

####

###Define Functions to Generate Output###

#generic function to generate a text file which summarizes everything
def printout(data,ext):
        output=open(outpath+txtname+'_'+ext+'.txt','w')
        for x in range(len(data)):
                for c in range(len(data[x])):
                        output.write("%-20s"%str(data[x][c]))
                        #output.write('\t')
                output.write('\n')
        print '\n'+ext+' summary file saved.\n'
        output.close()
        return 0

#generic function to make new directories
def check_path(path,ext):
        try:
		os.makedirs(path+ext+'/old')
	except OSError:
		pass
	try:
                os.makedirs(path+ext+'/')
        #avoid race condition
        except OSError:
	##REMOVES ANY PREVIOUS PLOTS#########
		os.system('mv -f '+path+ext+'/*.png '+path+ext+'/old')
		if not os.path.isdir(path):
                        raise
        return 0

###Author's Note###
#i've written separate functions for generating output plots. there is a lot of reduncdancy among the functions
#and they can certainly be combined/separated/recursively called. i've opted against doing that here because it
#is easier for me to troubleshoot if i'm not using recursive functions but feel free to combine them if you want.

#function to generate 2-param fit plots
def gen_2param_plot(data,uncert,ext):
        #these definitions are for plotting residuals vs. fitted parameter
        rx=numpy.linspace(0.0,maxx,100)
        ry=numpy.zeros(len(rx))
        #this will make a smooth decay curve
        smoothtau=numpy.linspace(0.0,maxx,1000)        
        #define buffer vector
        temp=[]
        chi=[]
        #define output vector
        rates=[]
        print '\nInitiating 2-parameter fitting'
        temp.append('ID')
        #temp.append('Int(0)')
        temp.append('TIME CONSTANT')
        temp.append('TIME ERROR')
        temp.append('CHISQ')
        temp.append('RedCHISQ')
	#CARO
	temp.append('RESID')
        rates.append(temp)
        temp=[]
	coun=0
        dof=float((len(tau)-len(p2)))
	print 'degrees of freedom = '+str(dof)
        for x in range(len(data)):
                print '\n2-parameter Fitting '+str(data[x][0])
		#do the fit, calculate chi-squared
                plsq2=scipy.optimize.leastsq(residuals2,p2,args=(data[x][1:taulen],tau,uncert[x][1:]),full_output=1)
                chi2=((plsq2[2]['fvec'])**2).sum()
		chi.append(chi2)
                #define text to put on plots
		
		print len(data),coun,x, data[x]
		coun+=1
		text='Fitted TC = '+str(numpy.round(plsq2[0][1],3))
                text2='Error in Fitted Time ='+str(numpy.round((numpy.sqrt((chi2/dof)*plsq2[1][1,1])),5))
                text3='Reduced Chi-Squared ='+str(numpy.round((chi2/dof),3))
                text4='Initial Intensity ='+str(data[x][-1])
		#save rates, errors, and stats in buffer
                temp.append(data[x][0])
                #temp.append(str(numpy.round(plsq2[0][0],3)))
                temp.append(str(numpy.round(plsq2[0][1],3)))
                temp.append(str(numpy.round((numpy.sqrt((chi2/dof)*plsq2[1][1,1])),5)))
                temp.append(str(numpy.round(chi2,3)))
                temp.append(str(numpy.round(chi2/dof,3)))
		#CARO RESID
		temp.append(str(data[x][0][-3:]))
                rates.append(temp)
                #generate plot figures
                fig=plt.figure()
                #plot of residuals vs. fitted parameter
                ax1=fig.add_axes([0.18, 0.82, 0.72, 0.08])
                ax1.plot(tau, plsq2[2]['fvec'],'ok',rx,ry,'--r')
                plt.title('2-Parameter Fit for Res: '+str(data[x][0]))
                ax1.set_ylabel('Residuals')
                #plot decay and fit with error bars
		ax2=fig.add_axes([0.18, 0.09, 0.72, 0.67])
                ax2.plot(smoothtau,fit2(smoothtau,plsq2[0]),'-k',tau,data[x][1:taulen],'ok',linewidth=1.5)
                ax2.errorbar(tau,data[x][1:taulen],xerr=None,yerr=uncert[x][1:],color='red',barsabove=True,fmt=',',linewidth=1.5)
                ax2.set_xlabel('Delay Time (s)')
                ax2.set_ylabel('Intensity (a.u.)')
		#plot text
                plt.figtext(0.48,0.70,text,fontsize=16)
                plt.figtext(0.48,0.65,text2,fontsize=16)
                plt.figtext(0.48,0.60,text3,fontsize=16)
		plt.figtext(0.48,0.55,text4,fontsize=16)
                #CARO
		#set y-axis limits
		plt.ylim(0,1.2)
		#write plot to disk
                plt.savefig(outpath+twopfolder+'/'+str(data[x][0][0:])+'_plot_'+ext+'.png')
                plt.clf()
		#CARO
		plt.close('all')
                temp=[]
        print '\n2-parameter fitting is finished.'
        chimean=numpy.mean(chi)
        return rates, chimean

#function to generate 3-param fit plots
def gen_3param_plot(data,uncert,ext):
        #these definitions are for plotting residuals vs. fitted parameter
        rx=numpy.linspace(0.0,maxx,100)
        ry=numpy.zeros(len(rx))
        #this will make a smooth decay curve
        smoothtau=numpy.linspace(0.0,maxx,1000)        
        #define buffer vector
        temp=[]
        chi=[]
        #define output vector
        rates=[]
        temp.append('ID')
       # temp.append('Int(0)')
        temp.append('TIME CONSTANT')
        temp.append('TIME ERROR')
        temp.append('OFFSET')
        temp.append('OFFSET ERROR')
        temp.append('CHISQ')
        temp.append('RedCHISQ')
        temp.append('RESID')
	rates.append(temp)
        temp=[]
        dof=float((len(tau)-len(p3)))
        print '\nInitiating 3-parameter fitting'
        for x in range(len(data)):
                print '\n3-parameter Fitting '+str(data[x][0])
                plsq3=scipy.optimize.leastsq(residuals3,p3,args=(data[x][1:taulen],tau,uncert[x][1:]),full_output=1)
                chi2=((plsq3[2]['fvec'])**2).sum()
                chi.append(chi2)
                #define text to put on plots
		
		print len(data),x, data[x]
		text='Fitted TC = '+str(numpy.round(plsq3[0][1],3))
                text2='Error in Fitted Time ='+str(numpy.round((numpy.sqrt((chi2/dof)*plsq3[1][1,1])),3))
                text3='Offset = '+str(numpy.round(plsq3[0][2],3))
                text4='Error in Offset ='+str(numpy.round((numpy.sqrt((chi2/dof)*plsq3[1][2,2])),3))
                text5='Reduced Chi-Squared ='+str(numpy.round((chi2/dof),3))
                text6='Initial Intensity ='+str(data[x][-1])
		#save rates, errors, and stats in buffer
                temp.append(data[x][0])
                #temp.append(str(numpy.round(plsq3[0][0],3)))
                temp.append(str(numpy.round(plsq3[0][1],3)))
                temp.append(str(numpy.round((numpy.sqrt((chi2/dof)*plsq3[1][1,1])),3)))
                temp.append(str(numpy.round(plsq3[0][2],3)))
                temp.append(str(numpy.round((numpy.sqrt((chi2/dof)*plsq3[1][2,2])),3)))
                temp.append(str(numpy.round(chi2,3)))
                temp.append(str(numpy.round(chi2/dof,3)))
                temp.append(str(data[x][0][-3:]))
		rates.append(temp)
                #for plotting the offset
                oy=numpy.ones(len(rx))*plsq3[0][2]
                #generate plot figures
                fig=plt.figure()
                #plot of residuals vs. fitted parameter
                ax1=fig.add_axes([0.18, 0.82, 0.72, 0.08])
                ax1.plot(tau, plsq3[2]['fvec'],'ok',rx,ry,'--r')
                plt.title('3-Parameter Fit for Res: '+str(data[x][0]))
                ax1.set_ylabel('Residuals')
                #plot decay and fit with error bars
                ax2=fig.add_axes([0.18, 0.09, 0.72, 0.67])
                ax2.plot(smoothtau,fit3(smoothtau,plsq3[0]),'-k',tau,data[x][1:taulen],'ok',rx,oy,'--r',linewidth=1.5)
                ax2.errorbar(tau,data[x][1:taulen],xerr=None,yerr=uncert[x][1:],color='red',barsabove=True,fmt=',',linewidth=1.5)
                ax2.set_xlabel('Delay Time (s)')
                ax2.set_ylabel('Intensity (a.u.)')
                #plot text
                plt.figtext(0.48,0.70,text,fontsize=16)
                plt.figtext(0.48,0.65,text2,fontsize=16)
                plt.figtext(0.48,0.60,text3,fontsize=16)
                plt.figtext(0.48,0.55,text4,fontsize=16)
                plt.figtext(0.48,0.50,text5,fontsize=16)
                plt.figtext(0.48,0.45,text6,fontsize=16)
		#CARO
		#set y-axis limits
		plt.ylim(0,1.2)
		#write plot to disk
                plt.savefig(outpath+threepfolder+'/'+str(data[x][0][0:])+'_plot_'+'_'+ext+'.png')
                plt.clf()
		#CARO
		plt.close('all')
                temp=[]
        print '\n3-parameter fitting is finished'
        chimean=numpy.mean(chi)
        return rates, chimean

#generate bi-exponential fit
def gen_biexp_plot(data,uncert,ext):
        #these definitions are for plotting residuals vs. fitted parameter
        rx=numpy.linspace(0.0,maxx,100)
        ry=numpy.zeros(len(rx))
        #this will make a smooth decay curve
        smoothtau=numpy.linspace(0.0,maxx,1000)        
        #define buffer vector
        temp=[]
        chi=[]
        #define output vector
        rates=[]
        temp.append('ID')
        temp.append('Int 1(0)')
        temp.append('RATE 1')
        temp.append('RATE 1 ERROR')
        temp.append('Int 2(0)')
        temp.append('RATE 2')
        temp.append('RATE 2 ERROR')
        temp.append('CHISQ')
        temp.append('RedCHISQ')
        rates.append(temp)
        temp=[]
        dof=float((len(tau)-len(pb)))
        for x in range(len(data)):
                print '\nBi-exp-2-parameter Fitting '+str(data[x][0])
                plsqb=scipy.optimize.leastsq(residualsbi,pb,args=(data[x][1:taulen],tau,uncert[x][1:]),full_output=1)
                chi2=((plsqb[2]['fvec'])**2).sum()
                chi.append(chi2)
                #define text to put on plots
                print plsqb[0]
                text='Fitted Rate 1 = '+str(numpy.round(plsqb[0][2],3))
                text1='Error in Fitted Rate 1 =' +str(numpy.round((numpy.sqrt((chi2/dof)*plsqb[1][2,2])),3))
                text2='Fitted Rate 2 = '+str(numpy.round(plsqb[0][3],3))
                text3='Error in Fitted Rate 2 = '+str(numpy.round((numpy.sqrt((chi2/dof)*plsqb[1][3,3])),3))
                text4='Reduced Chi-Squared = '+str(numpy.round((chi2/dof),3))
                #save rates, errors, and stats in buffer
                temp.append(data[x][0])
                temp.append(str(numpy.round(plsqb[0][0],3)))
                temp.append(str(numpy.round(plsqb[0][2],3)))
                temp.append(str(numpy.round((numpy.sqrt((chi2/dof)*plsqb[1][2,2])),3)))
                temp.append(str(numpy.round(plsqb[0][1],3)))
                temp.append(str(numpy.round(plsqb[0][3],3)))
                temp.append(str(numpy.round((numpy.sqrt((chi2/dof)*plsqb[1][3,3])),3)))
                temp.append(str(numpy.round(chi2,3)))
                temp.append(str(numpy.round(chi2/dof,3)))
                rates.append(temp)
                #generate plot figures
                fig=plt.figure()
                #plot of residuals vs. fitted parameter
                ax1=fig.add_axes([0.18, 0.82, 0.72, 0.08])
                ax1.plot(tau, plsqb[2]['fvec'],'ok',rx,ry,'--r')
                plt.title('Bi-exp-2-parameter for Res: '+str(data[x][0]))
                ax1.set_ylabel('Residuals')
                #plot decay and fit with error bars
                ax2=fig.add_axes([0.18, 0.09, 0.72, 0.67])
                ax2.plot(smoothtau,fitbi(smoothtau,plsqb[0]),'-k',tau,data[x][1:taulen],'ok')
                ax2.errorbar(tau,data[x][1:taulen],xerr=None,yerr=uncert[x][1:],color='red',barsabove=True,fmt=',',linewidth=1.5)
                ax2.set_xlabel('Delay Time (s)')
                ax2.set_ylabel('Intensity (a.u.)')
                #plot text
                plt.figtext(0.48,0.70,text,fontsize=16)
                plt.figtext(0.48,0.65,text1,fontsize=16)
                plt.figtext(0.48,0.60,text2,fontsize=16)
                plt.figtext(0.48,0.55,text3,fontsize=16)
                plt.figtext(0.48,0.50,text4,fontsize=16)
                #write plot to disk
                plt.savefig(outpath+bifolder+'/'+str(data[x][0][1:])+'_plot_'+str(data[x][0])+'_'+ext+'.png')
                plt.clf()
                #CARO
		plt.close('all')
		temp=[]
        print '\nBi-exp-2-parameter fitting is finished'
        chimean=numpy.mean(chi)
        return rates, chimean

###Define Bootstrapping Functions###
#2-parameter fit rate bootstrapping
def boot_2param(data,uncert):
        i=1
        simrate=[]
        outboot=[]
        simdata=[]
        finalboot=[]
        finalboot.append('ID')
        finalboot.append('# MC')
        finalboot.append('MEAN TIME')
        finalboot.append('STDV TIME')
        outboot.append(finalboot)
        finalboot=[]
        for x in range(len(data)):
                print '\nStarting monte carlo for 2-parameter fit of '+str(data[x][0])
                for b in range(twobootit):
                        simdata.append(data[x][0])
                        while i < (len(data[x][1:taulen])+1):
                            offset=numpy.random.normal(data[x][i],uncert[x][i])
                            simdata.append(offset)
                            i+=1
                        i=1
                        plsq2=scipy.optimize.leastsq(residuals2,p2,args=(simdata[1:],tau,uncert[x][1:]),full_output=1)
                        chi2=((plsq2[2]['fvec'])**2).sum()
                        #save estimated rate for simulated data set
                        simrate.append(plsq2[0][1])
                        simdata=[]
                bootmean=numpy.mean(simrate)
                bootstdv=numpy.std(simrate,ddof=1)
                finalboot.append(data[x][0])
                finalboot.append(twobootit)
                finalboot.append(numpy.round(bootmean,3))
                finalboot.append(numpy.round(bootstdv,3))
                outboot.append(finalboot)
                #txt to put on plots
                mctxt1='Mean = '+str(numpy.round(bootmean,3))
                mctxt2='STDV = '+str(numpy.round(bootstdv,3))
                mctxt3='-- Fit to PDF'
                #generate images for bootstrap
                binnum=int(2.0*(len(simrate)**(1.0/3.0)))
                n,bins,patches=plt.hist(simrate,binnum,normed=1, alpha=0.8)
                pdfsmooth=numpy.linspace(min(bins),max(bins),500)
                hisfit=mlab.normpdf(pdfsmooth,bootmean,bootstdv)
                plt.plot(pdfsmooth,hisfit,'r--',linewidth=2.0)
                plt.title('Monte Carlo Results for 2-parameter fit of Res: '+str(data[x][0]))
                plt.xlabel('Fitted Time Constant of Simulated Data')
                plt.ylabel('Probability Density')
                plt.figtext(0.68,0.85,mctxt1,fontsize=16)
                plt.figtext(0.68,0.80,mctxt2,fontsize=16)
                plt.figtext(0.68,0.75,mctxt3,fontsize=16)
                plt.savefig(outpath+twobfolder+'/'+'mc_hist_'+str(data[x][0])+'_'+ext4+'.png')
                plt.clf()
                #CARO
		plt.close('all')
		finalboot=[]
                simrate=[]
#                print '\n'+str(twobootit)+' simulations for '+str(data[x][0])+' completed\n'
        return outboot

#3-parameter fit rate bootstrapping
def boot_3param(data,uncert):
        i=1
        simrate=[]
        simasym=[]
        outboot=[]
        simdata=[]
        finalboot=[]
        finalboot.append('ID')
        finalboot.append('# MC')
        finalboot.append('MEAN TIME')
        finalboot.append('STDV TIME')
        finalboot.append('MEAN OFFSET')
        finalboot.append('STDV OFFSET')
        outboot.append(finalboot)
        finalboot=[]
        for x in range(len(data)):
#                print '\nStarting monte carlo for 3-parameter fit of '+str(data[x][0])
                for b in range(threebootit):
                        simdata.append(data[x][0])
                        while i < (len(data[x][1:taulen])+1):
                            offset=numpy.random.normal(data[x][i],uncert[x][i])
                            simdata.append(offset)
                            i+=1
                        i=1
                        plsq3=scipy.optimize.leastsq(residuals3,p3,args=(simdata[1:],tau,uncert[x][1:]),full_output=1)
                        chi2=((plsq3[2]['fvec'])**2).sum()
                        #save estimated rate for simulated data set
                        simrate.append(plsq3[0][1])
                        simasym.append(plsq3[0][2])
                        simdata=[]
                bootmean=numpy.mean(simrate)
                bootstdv=numpy.std(simrate,ddof=1)
                bootasymmean=numpy.mean(simasym)
                bootasymstdv=numpy.std(simasym,ddof=1)            
                finalboot.append(data[x][0])
                finalboot.append(threebootit)
                finalboot.append(numpy.round(bootmean,3))
                finalboot.append(numpy.round(bootstdv,3))
                finalboot.append(numpy.round(bootasymmean,3))
                finalboot.append(numpy.round(bootasymstdv,3))
                outboot.append(finalboot)
                #txt to put on plots
                mctxt1='Mean = '+str(numpy.round(bootmean,3))
                mctxt2='STDV = '+str(numpy.round(bootstdv,3))
                mctxt3='-- Fit to PDF'
                #generate images for bootstrap
                binnum=int(2.0*(len(simrate)**(1.0/3.0)))
                n,bins,patches=plt.hist(simrate,binnum,normed=1, alpha=0.8)
                pdfsmooth=numpy.linspace(min(bins),max(bins),500)
                hisfit=mlab.normpdf(pdfsmooth,bootmean,bootstdv)
                plt.plot(pdfsmooth,hisfit,'r--',linewidth=2.0)
                plt.title('Monte Carlo Results for 3-parameter fit of Res: '+str(data[x][0]))
                plt.xlabel('Fitted Time Constant of Simulated Data')
                plt.ylabel('Probability Density')
                plt.figtext(0.68,0.85,mctxt1,fontsize=16)
                plt.figtext(0.68,0.80,mctxt2,fontsize=16)
                plt.figtext(0.68,0.75,mctxt3,fontsize=16)
                plt.savefig(outpath+threebfolder+'/'+'mc_hist_'+str(data[x][0])+'_'+ext5+'.png')
                plt.clf()
                #CARO
		plt.close('all')
		finalboot=[]
                simrate=[]
                simasym=[]
#                print '\n'+str(threebootit)+' simulations for '+str(data[x][0])+' completed\n'
        return outboot

###General Function Calls###

print 'MASFIT v.9 by Matt Stetz\n'

#create directories

if twop_choice == 'y' or twop_choice == 'Y':
        print 'Creating '+twopfolder+' directory\n'
        check_path(outpath,twopfolder)
else:
	pass

if threep_choice == 'y' or threep_choice == 'Y':
        print 'Creating '+threepfolder+' directory\n'
        check_path(outpath,threepfolder)
else:
        pass

if twoboot=='y' or twoboot=='Y':
        print 'Creating '+twobfolder+' directory\n'
        check_path(outpath,twobfolder)
else:
        pass

if threeboot=='y' or threeboot=='Y':
        print 'Creating '+threebfolder+' directory\n'
        check_path(outpath,threebfolder)
else:
        pass

if bi_choice=='y' or bi_choice=='Y':
        print 'Creating '+ bifolder+' directory\n'
        check_path(outpath,bifolder)
else:
        pass

print '\n*** *** *** *** *** ***\n'

#the following two functions get data from the txt and clean it up
pdata=read_file(datapath)
pdata=parse_data(pdata)

#only call the get_dataheight function if using absolute data heights
if normalized == 'n' or normalized == 'N':
     pdata=get_dataheight(pdata)
else:
     pass

#get the uncertainties
s1,s2,s3,ulist=get_sigma(pdata,sigmascale)

###Fitting Function Calls###

#two-parameter fit
if twop_choice == 'y' or twop_choice == 'Y':
        twoparamrates,twoparamchimean=gen_2param_plot(pdata,ulist,ext2)
        printout(twoparamrates,ext2)
        print '\n*** *** *** *** *** ***\n'
elif twop_choice == 'n' or twop_choice == 'N':
        print 'ATTN: You have chosen not to do a 2-parameter fit.'
        print '\n*** *** *** *** *** ***\n'
else:
        print 'ERROR: You have entered an invalid choice for the 2-parameter fit flag!'
        print '\n*** *** *** *** *** ***\n'

#three-parameter fit
if threep_choice == 'y' or threep_choice == 'Y':
        threeparamrates,threeparamchimean=gen_3param_plot(pdata,ulist,ext3)
        printout(threeparamrates,ext3)
        print '\n*** *** *** *** *** ***\n'
elif threep_choice == 'n' or threep_choice == 'N':
        print 'ATTN: You have chosen not to do a 3-parameter fit.'
        print '\n*** *** *** *** *** ***\n'
else:
        print 'ERROR: You have entered an invalid choice for the 3-parameter fit flag!'
        print '\n*** *** *** *** *** ***\n'

#two-parameter fit bootstrap
if twoboot == 'y' or twoboot == 'Y':
        twoparambootresults=boot_2param(pdata,ulist)
        printout(twoparambootresults,ext4)
        print '\n*** *** *** *** *** ***\n'
elif twoboot == 'n' or twoboot == 'N':
        print 'ATTN: You have chosen not to do bootstrapping for 2-parameter fits.'
        print '\n*** *** *** *** *** ***\n'
else:
        print 'ERROR: You have entered an invalid choice for the 2-parameter fit bootstrap flag!'
        print '\n*** *** *** *** *** ***\n'

#three-parameter fit bootstrap        
if threeboot == 'y' or threeboot == 'Y':
        threeparambootresults=boot_3param(pdata,ulist)
        printout(threeparambootresults,ext5)
        print '\n*** *** *** *** *** ***\n'
elif threeboot == 'n' or threeboot == 'N':
        print 'ATTN: You have chosen not to do bootstrapping for 3-parameter fits.'
        print '\n*** *** *** *** *** ***\n'
else:
        print 'ERROR: You have entered an invalid choice for the 3-parameter fit bootstrap flag!'
        print '\n*** *** *** *** *** ***\n'


#bi-exp-2-parameter fit
if bi_choice == 'y' or bi_choice == 'Y':
        biexpparamrates,biexpchimean=gen_biexp_plot(pdata,ulist,extb)
        printout(biexpparamrates,extb)
        print '\n*** *** *** *** *** ***\n'
elif bi_choice == 'n' or bi_choice == 'N':
        print 'ATTN: You have chosen not to do a bi-exp-2-parameter fit.'
        print '\n*** *** *** *** *** ***\n'
else:
        print 'ERROR: You have entered an invalid choice for the 3-parameter fit flag!'
        print '\n*** *** *** *** *** ***\n'

###I have opted to put the error messages at the very end so they won't get lost in the log
###Error Messages and Warnings###
print '\nChecking for errors and warnings...'

if twop_choice == 'y' or threep_choice == 'Y':
        if twoparamchimean > (len(tau)+1):
                print '\nWARNING: Chi-squared values for 2-parameter fits are too high!'
                print '\nCheck data and/or consider re-scaling your uncertainties!'
                print '\nYour average chi-squared value = '+str(twoparamchimean)
                print '\nYour current uncertainty scaling factor is = '+str(sigmascale)
        else:
                print '\n2-parameter fitting finished without any warnings.'
else:
        print '\n2-parameter fitting was not performed.'
if threep_choice == 'y' or threep_choice == 'Y':
        if threeparamchimean > (len(tau)+1):
                print '\nWARNING: Chi-squared values for 3-parameter fits are too high!'
                print '\nCheck data and/or consider re-scaling your uncertainties!'
                print '\nYour average chi-squared value = '+str(threeparamchimean)
                print '\nYour current uncertainty scaling factor is = '+str(sigmascale)
        else:
                print '\n3-parameter fitting finished without any warnings.'
else:
        print '\n3-parameter fitting was not performed.'

if bi_choice == 'y' or bi_choice == 'Y':
        if biexpchimean > (len(tau)+1):
                print '\nWARNING: Chi-squared values for bi-exp-2-parameter are too high!'
                print '\nCheck data and/or consider re-scaling your uncertainties!'
                print '\nYour average chi-squared value = '+str(biexpchimean)
                print '\nYour current uncertainty scaling factor is = '+str(sigmascale)
        else:
                print '\nbi-exp-2-parameter fitting finished without any warnings.'
else:
        print '\nbi-exp-2-parameter fitting was not performed.'

