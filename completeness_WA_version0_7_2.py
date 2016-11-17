#! /usr/bin python
# -*- coding: utf-8 -*-

"""
Name: Completeness_WA
Version 0.7.1 (first release, beta, point-like source version)
	Edited by clsj at CEFCA.
	0.7.1 : Including the realistic galaxy models.
Description: Provide completeness percentages, over a range of magnitudes
and types of objects, for a given image
Authors: Walter Santos & Alberto Molino
Created on: May 14th 2015
Last Updated: Nov 23th 2015
Instructions: Run the code with a -h/--help option to list all the arguments,
necessary and optional, on the command line.
Latest Changes:
-
to change that to a galaxy/extended library later
Requirements: Numpy, SciPy, Astropy, PyRAF, Sextractor, pymangle
TODO:
- Implement more safeguards, exception handlings, ...
- Define a ___main___ so people can import the code elsewhere
- Follow the defined conventions for a python code, to be compatible with 
the rest of the J-PAS pipeline codes
"""

"""
IMPORTS
"""
import argparse
import os,sys
import commands
import numpy as N
#from pyraf.iraf import psfmatch
from scipy import interpolate
from astropy.io import fits
from astropy.stats import sigma_clip
from astropy import wcs
from scipy.signal import convolve2d
import numpy.ma as ma
from scipy.ndimage.filters import gaussian_filter
#import useful as U
import random
# import matplotlib.pyplot as plt
import pymangle
###

SExroot = '/usr/share/sextractor/'

"""
CONSTANTS
TODO: Probably move them from script constants to input parameters in the
main.
"""
MODEL_OBJECTS_PIXSCALE = 0.065 #arcsec/pixel
WEIGHT_MAP_CUT = 0.5
PSF_CUT = 0.01
GAIN = 13.42

"""
FUNCTIONS
"""
def generateSexInputModels(modelsTable, firstLine='# X Y magID modelID mag', filename='input_model.coo'):
    """
    Generate Sextractor input txt file that contains the list of objects to be found in the image.
    """
    f = open(filename, 'w')
    f.write(firstLine+'\n')
    
    for model in modelsTable:
        for column in model:
            f.write(str(column)+' ')
        f.write('\n')
    f.close()
    

def generateSexParams(filename='assoc.param'):
    """
    Generate Sextractor input param file: the list of info it should output.
    """
    f = open(filename, 'w')
    
    outParams = ['NUMBER',
              'XWIN_IMAGE',
              'YWIN_IMAGE',
              'VECTOR_ASSOC(1)', #assoc input x
              'VECTOR_ASSOC(2)', #assoc input y
              'VECTOR_ASSOC(3)',  #assoc input id
              'VECTOR_ASSOC(4)', #assoc mag
              'VECTOR_ASSOC(5)',
              'ISOAREA_IMAGE',
              'MAG_APER[3]',
              'MAG_AUTO',
              'FWHM_WORLD[1]'
             ]

    for param in outParams:
        f.write(param+'\n')
    f.close()
    
    
def generateSexConfig(outFile, configFile,
                      paramsFile, inputModelsfile, weightImage, sexFile,
                      matchRadius):
    """
    Generate Sextractor input config file, with information for the run
    """                      
    f = open(configFile, 'w')

    #First generate a standard config options...
    """
    sexConfig = dict(CATALOG_NAME='matched.cat', 
                    CATALOG_TYPE='ASCII_HEAD',  
                    PARAMETERS_NAME='assoc.param', 
                    DETECT_TYPE='CCD',
                    DETECT_MINAREA='8.78733031674',
                    DETECT_THRESH='1.25',
                    ANALYSIS_THRESH='1.25',
                    FILTER='Y',
                    FILTER_NAME=SExroot+'tophat_3.0_3x3.conv',
                    DEBLEND_NTHRESH='64',
                    DEBLEND_MINCONT='0.0002',
                    CLEAN='Y',
                    CLEAN_PARAM='1.',
                    MASK_TYPE='CORRECT',
                    PHOT_APERTURES='14.0',
                    PHOT_AUTOPARAMS='2.5,3.5',
                    SATUR_LEVEL='50000.0',
                    MAG_ZEROPOINT='31.2671969775',
                    MAG_GAMMA='4.0',
                    GAIN='13.4213208113',
                    PIXEL_SCALE='0.221',
                    SEEING_FWHM='1.060202',
                    STARNNW_NAME=SExroot+'default.nnw',
                    BACK_SIZE='64',
                    BACK_FILTERSIZE='3',
                    BACKPHOTO_TYPE='LOCAL',
                    BACKPHOTO_THICK='102',
                    MEMORY_OBJSTACK='15000',
                    MEMORY_PIXSTACK='2600000',
                    MEMORY_BUFSIZE='4600',
                    ASSOC_NAME='input_model.coo',
                    ASSOC_PARAMS='1,2',
                    ASSOC_RADIUS='5.0',
                    ASSOCSELEC_TYPE='MATCHED',
                    ASSOC_TYPE='NEAREST',
                    ASSOC_DATA='1,2,3,4,5',
                    VERBOSE_TYPE='NORMAL',
                    #WEIGHT_GAIN='Y',
                    WEIGHT_IMAGE='/home/walter/Dropbox/Completeness/f02p01_F814W_2.swp.weight.fits',
                    WEIGHT_TYPE='MAP_WEIGHT'
                    )
    """

    sexdata = open(sexFile, "r")
    
    sexdict = {l.split()[0] : l.split()[1] for l in sexdata.readlines() if l[0] != " "}
    print sexdict
    sexdata.close()
    
    sexConfig = sexdict

    #... then change the information with some input, image-specific, info
    sexConfig['CATALOG_NAME']=outFile
    sexConfig['PARAMETERS_NAME']=paramsFile
    #sexConfig['PIXEL_SCALE']=sexdict['PIXEL_SCALE']
    #sexConfig['SEEING_FWHM']=sexdict['SEEING_FWHM']
    sexConfig['ASSOC_RADIUS']=str(matchRadius)
    sexConfig['ASSOC_NAME']=inputModelsfile
    sexConfig['GAIN']=sexdict['GAIN']
    #sexConfig['DETECT_MINAREA']=sexdict['DETECT_MINAREA']
    #sexConfig['DETECT_THRESH']=sexdict['DETECT_THRESH']
    #sexConfig['ANALYSIS_THRESH']=sexdict['ANALYSIS_THRESH']
    #sexConfig['MAG_ZEROPOINT']=sexdict['MAG_ZEROPOINT']
    aperture =  float(sexConfig['SEEING_FWHM'])/float(sexConfig['PIXEL_SCALE'])
    sexConfig['PHOT_APERTURES']=str(aperture)+','+str(2*aperture)+','+str(3*aperture)
    sexConfig['WEIGHT_IMAGE']=str(weightImage)
    sexConfig['FILTER_NAME']=SExroot+sexConfig['FILTER_NAME']
    sexConfig['STARNNW_NAME']=SExroot+sexConfig['STARNNW_NAME']
    sexConfig['ASSOC_PARAMS']='1,2'
    sexConfig['ASSOCSELEC_TYPE']='MATCHED'
    sexConfig['ASSOC_TYPE']='NEAREST'
    sexConfig['ASSOC_DATA']='1,2,3,4,5'
    
    #Now write them out in a config file
    for key, value in sexConfig.iteritems():
        f.write(key+' '+value+'\n')
    f.close()

    
def changepixelscale(imagein,finalpix,imageinHDU=0,oripix=MODEL_OBJECTS_PIXSCALE):
    """
    Change the pixelscale ('origpix') from a model image ('imagein') to be 
    inserted to the pixelscale ('finalpix') of the final image.
    """
    interpol = 'linear' # 'cubic','quintic'
    fill_value = 0.
    imageinHDU = fits.open(imagein)[imageinHDU]
    header=imageinHDU.h-eader
    data = imageinHDU.data
    
    ny = N.shape(data)[0]
    nx = N.shape(data)[1]
    xx = N.linspace(0,oripix*nx,nx)
    yy = N.linspace(0,oripix*ny,ny)
    
    interpmat = interpolate.interp2d(yy,xx,data,interpol,fill_value)
    
    nx2=(oripix/finalpix)*nx
    ny2=(oripix/finalpix)*ny
    xnew = N.linspace(0,finalpix*nx2,nx2)
    ynew = N.linspace(0,finalpix*ny2,ny2)
    
    imagein_path, _ = os.path.splitext(imagein)
    imageout = imagein_path+'.repixeld.fits'
    newmat = interpmat(ynew, xnew)
    newheader = header
    fits.writeto(imageout,newmat,newheader,clobber=True)
    
    return imageout
    
    
def get_nickname(image):
    """
    It returns the nickname for an input file.
    """
    extension = len(image.split('/')[-1:][0].split('.')[-1:][0])
    return image.split('/')[-1:][0][:-extension-1]
    

def get_nickname_mine(image):
    return os.path.splitext(os.path.basename(image))[0]

def get_filepath_mine(image):
    return (os.path.dirname(image)+'/')
    
    
def get_filepath(image):
    """
    It returns a path file.
    """
    extension = (len(image.split('/')[-1:][0]))
    return image[:-extension]


def psfconvolve(image,psfimage,psfnew,newimage):
     """
     It convolves a sample of images with a sample PSF-models.
     """
     # Creating the kernel to degrade a PSF to another
     path = get_filepath(image)
     nick_psfimage = get_nickname(psfimage)
     nick_psfnew   = get_nickname(psfnew)
     kernel_psf = '%skernelPSF_%s_to_%s.fits'%(path,nick_psfimage,nick_psfnew)
     
     #Overwrite the kernel if it already exists
     if os.path.exists(kernel_psf):
        os.remove(kernel_psf)
                
     psfmatch(input=psfimage,reference=psfnew,
            psfdata=psfimage,kernel=kernel_psf,
            convolution="psf",background='none',verbose='yes')
     
     #Overwrite the new image if it already exists
     if os.path.exists(newimage):
        os.remove(newimage)

     # Degrading image to match new PSF condition.
     psfmatch(input=image,reference=psfnew,psfdata=psfimage,
            kernel=kernel_psf,output=newimage,
            convolution="kernel", verbose='yes')
    
     return newimage
     

def get_psffwhm(image,pixelscale):
    """
    It estimates the FWHM from a PSF-model
    ---USAGE:
    image = 'alhambra.psf.fits'
    pixelscale = 0.221 #"/pix
    seeing = get_psffwhm(image,pixelscale)
    
    """
    data = fits.open(image)[0].data # Open image
    pixels = N.shape(data)[0] # Matrix size
    base = N.arange(pixels) # base for the radial distrib.
    dim = pixels/2   # Matrix dimension
    vect_x = data[dim,:] # Horizontal signal
    vect_y = data[:,dim] # Vertical signal
    signal = (vect_x+vect_y)/2. # Averaged signal
    signal /= signal.max() # Normaliz. signal
    pos_peak = N.where(signal==max(signal))[0][0] # Peak position
    # Position where distribution drops 50% its signal (FWHM)
    pos_50 = N.where(abs(signal-0.5)==min(abs(signal-0.5)))[0][0]
    fwhm_pix = abs(base[pos_50]-base[pos_peak])*2. # Both sides 
    # Converting to "
    return fwhm_pix, fwhm_pix * pixelscale


def isInsideImage(x, y, nx, ny, imageNx, imageNy):
    """
    Simple check o see if a model image fits in a field image, based on their
    size
    """
    return ( ((x+nx) < imageNx) and ((y+ny) < imageNy) )


def getWeightAverage(x, y, nx, ny, weight):
    """
    Get average values in a weight map
    """
    weightData = fits.open(weight)[0].data
    return N.mean(weightData[y:y+ny,x:x+nx])


def sumModelIntoImage(x, y, nx, ny, modelData, imageData):
    """
    Insert (sum) model image into the field image
    """
    imageData[y:y+ny,x:x+nx] += modelData


def isModelDetected(modelID, sexOutData):
    """
    Based on sextractor's output matched catalog, check whether a specific
    inserted model was detected, based on its inserted id
    """
    detectedIDList = [row[5] for row in sexOutData]
    
    for detectedID in detectedIDList:
        if modelID == int(float(detectedID)):
            return True 
            

def generateMaskedImage(imagein):
    """
    Generate a new masked image with a sigma clipping cut and with only positive
    values, from an original 'imagein'.
    """
    orig_subdata = fits.open(imagein)[0].data
    imageinHDR = fits.open(imagein)[0].header
    
    orig_filtered_data = sigma_clip(orig_subdata, sig=2.9, iters=5)
    sigma_mask = ma.getmask(orig_filtered_data)
   
    detections_pos_data = ma.masked_greater_equal(orig_filtered_data, 0.0)
    pos_mask = ma.getmask(detections_pos_data)
    final_mask = sigma_mask & pos_mask
    
    outdata = N.zeros((5001,5001))
    dataout = ma.masked_array(outdata, final_mask).filled(1)
    
    imagein_path, _ = os.path.splitext(imagein)
    imageout = imagein_path+'masked29_gauss.fits'
    fits.writeto(imageout,dataout,imageinHDR,clobber=True)
    
    
def getEffectiveNumberPixels(weightMapData):
    """
    From a weight map image, count the number of effective pixels,
    based on a predefined weight cut: 0.5 as a constant, for now
    """
    count = 0
    for value in N.nditer(weightMapData):
        if value >= WEIGHT_MAP_CUT:
            count += 1
    return count
        

def findMaxModelsToInsert(imageData, weightMapData, psfFWHM_pix, signalDensity, maxpercent=0.03, times=4):
    """
    Come up with a maximum number of model objects to insert into a field image,
    so that the resulting image doesn't have the background signal changed
    by more than 'maxpercent=0.03'.
    """    
    psfArea = N.pi*((psfFWHM_pix)/2)**2
    galAvgArea = times*psfArea
    galDensity = maxpercent*signalDensity
    
    totalEffectivePixels = getEffectiveNumberPixels(weightMapData)
    galTotalPixels = totalEffectivePixels*galDensity
    
    return galTotalPixels/galAvgArea
    
        
def getPixelSignalDensity(imageData, weightMapData=None, size=1000):
    """
    Calculate a signal density (signals over a certain background, calculated
    iteratively). The analysis is done considering only a sub-sample of image
    around its center.
    """    
    #get a sub-square around the center
    ny = N.shape(imageData)[0]
    nx = N.shape(imageData)[1]
    origSubdata = imageData[(nx/2-size/2):(nx/2+size/2),(ny/2-size/2):(ny/2+size/2)]
    
    #invert the image, so that we can also have the "background signals"    
    invSubdata = (-1)*origSubdata
    hasWeightMap = False    
    
    if weightMapData is not None:
        hasWeightMap = True
    
    #constants and initial conditions
    SIGMA_INITIAL = 3.0
    MIN_RATIO = 0.03
    MAX_RATIO = 0.10
    ratio = 0.0
    SIGMA_FILTER = 1.7 #for the gaussian filter to eliminate 1-2 pixels background
    
    # consider the weight map, if given: make a mask out of it, based on the weight cut constant
    if hasWeightMap:
        weightMapSubdata = weightMapData[(nx/2-size/2):(nx/2+size/2),(ny/2-size/2):(ny/2+size/2)]
        weightMapNormal = N.copy(weightMapSubdata)
        for value in N.nditer(weightMapNormal, op_flags=['readwrite']):
            if value >= WEIGHT_MAP_CUT:
                value[...] = 0
            else:
                value[...] = 1
        weightMapMask = ma.make_mask(weightMapNormal)
        origSubdata = ma.array(origSubdata, mask=weightMapMask)
        
    #eliminate the 'too small' sized signals with a sigma filter
    origSubdata = gaussian_filter(origSubdata, SIGMA_FILTER)    
    
    #number of total pixels, of the sub-image, considering the weight map mask
    totalPixels = ma.count(origSubdata)
    
    #iteratively sigma clipping to derive ratios of "detections" and
    #"background" signals
    sigma = SIGMA_INITIAL+0.1
    while (ratio >= MAX_RATIO or ratio < MIN_RATIO):
        if (ratio >= MAX_RATIO):
            sigma += 0.1
        if (ratio < MIN_RATIO):
            sigma -= 0.1
    
        origSignals = sigma_clip(origSubdata, sig=sigma, iters=5)
        signals = origSignals.data[origSignals.mask]
        origPosSignals_count = 0
        for value in signals:
            if value > 0:
                origPosSignals_count += 1
        
        invSignals = sigma_clip(invSubdata, sig=sigma, iters=5)
        signals = invSignals.data[invSignals.mask]
        invPosSignals_count = 0
        for value in signals:
            if value > 0:
                invPosSignals_count += 1
    
        ratio = float(invPosSignals_count)/float(origPosSignals_count)
             
    
    #signal density (ratio of signals per pixel), to be used to define 
    #maximum number of inserted galaxies
    return (float(origPosSignals_count)/float(totalPixels))
    

def buildModelsList(modelsDir, modelTypes='All'):
    """
    build models list with:
    path,modelID,modelType,modelIndex,modelFlux
    TODO: Revise this function later once we have libraries of galaxy models
    """
    modelsList = []
    modelPaths = []
    try:
        modelPaths = [os.path.join(modelsDir,fn) for fn in next(os.walk(modelsDir))[2]] 
    except StopIteration:
        pass
    
    for path in modelPaths:
        if (os.path.splitext(os.path.basename(path))[1].lower() == '.fits'): #if extension is fits
            filename = os.path.splitext(os.path.basename(path))[0] #filename without extension
            filenameSplit = filename.split('.')
            if (filenameSplit[0][:3].lower() == 'gal'): #if model filename is galaxy
                modelType = 'galaxy'                
                modelID = filenameSplit[0][3:]
                modelIndex = 3 #change to the Sersic index later
                modelFlux = 0.0
                if (filenameSplit[1][:4].lower() == 'flux'):
                    modelFlux = float(filenameSplit[1][4:])
                modelsList.append([path, modelID, modelType, modelIndex, modelFlux])
            elif (filenameSplit[0][:4].lower() == 'star'): #if model filename is star
                modelType = 'star'                
                modelID = filenameSplit[0][4:]
                modelIndex = 1 #change to the Sersic index later
                modelFlux = 0.0
                if (filenameSplit[1][:4].lower() == 'flux'):
                    modelFlux = float(filenameSplit[1][4:])
                modelsList.append([path, modelID, modelType, modelIndex, modelFlux])
    return modelsList


def guessHDU(filename):
    """
    Try to guess which hdu contains the image, but the default overall is always
    the primary (hdu=0)
    """
    hdulist = fits.open(filename)
    if(len(hdulist)==1):
        hdulist.close()
        return 0
    else:
        for n,hdu in enumerate(hdulist):
            if (hdu.name == 'SCIENCE'):
                hdulist.close()
                return n
        for n,hdu in enumerate(hdulist):
            if (hdu.name == 'PRIMARY'):
                hdulist.close()
                return n
    hdulist.close()    
    return 0


def outputPlot(xpoints, ypoints, xline, yline, outfile='completeness.png', magType='Instrumental'):
    """
    Output the completeness results as a plot figure.
    """
    plt.figure(1, figsize=(9,8),dpi=80, facecolor='w', edgecolor='k') #New window for the plot
    plt.clf() #Cleaning the figure.
    plt.plot(xpoints,ypoints,'ko',alpha=0.5,ms=8) # This plots the individual points
    plt.plot(xline,yline,'r-',lw=5,alpha=0.5) # This plots the averaged line
    plt.grid() #This adds a grid
    plt.ylabel('Completeness Factor',size=30) #Ylabel
    plt.xlabel(magType+' Magnitude',size=30) #Xlabel
    plt.xlim(min(xline)-0.5,max(xline)+0.5) #X-ranges
    plt.ylim(0.,1.) #Y-ranges
    plt.xticks(fontsize=20) #Thickness of X-axis
    plt.yticks(fontsize=20) #Thickness of Y-axis
    plt.savefig(outfile,dpi=80) # Saving figure
    
    
def outputTable(xpoints, ypoints, ypointse, outfile='completeness.txt', magType='Instrumental'):
    """
    Output the completeness results (the 'mean' values in magnitude bins) in text file
    """
    fout = open(outfile, 'w')
    fout.write('# '+magType+'F814W fc fce\n')

    for i,_ in enumerate(xpoints):
        fout.write(str(xpoints[i])+' '+"%1.3f" %ypoints[i]+' '+"%1.3f" %ypointse[i]+'\n')
    fout.close()


def outputResults(mags, comp, compe, filename, magTypeInstorObs='Instrumental', magBin=0.2):
    """
    Output the results, both to a plot and a table file, of the completeness
    ratios data points and a 'mean' line averaged in magnitude bins    
    """
    #dm = magBin
    
    #m1 = min(mags) #inferior limit
    #m2 = max(mags) #superior limit
    #basem = N.arange(m1,m2+dm,dm) #base
    #line = U.bin_stats(mags,comp,basem,'mean_robust') #the averaged numbers.
    
    filename_path, _ = os.path.splitext(filename)
    #outfilePlot = filename_path+'.completeness.png'
    outfileTable = filename_path+'.completeness.txt'
    
    #outputPlot(mags, comp, basem, line, outfile=outfilePlot, magType=magTypeInstorObs)
    outputTable(mags,comp,compe,outfile=outfileTable,magType=magTypeInstorObs)
    
    
def buildFluxScales(fluxOrig, magBin=0.1, magMin=20, magMax=26, magZero=30.64):
    magScale = N.arange(magMin-magZero,magMax-magZero,magBin)
    
    fluxScales = []
    for mag in magScale:
        flux = 10**(mag/(-2.5))
        fluxScales.append(flux/fluxOrig)
   
    return fluxScales, magScale

def readWCS(header):
    """
    Module to load the WCS parameters from the header
    """
    # Creating the WCS instance
    w = wcs.WCS(header = header)
    # Filling the WCS
    w.wcs.crval = N.array([float(header["CRVAL1"]), float(header["CRVAL2"])])
    w.wcs.crpix = N.array([float(header["CRPIX1"]), float(header["CRPIX2"])])
    w.wcs.cdelt = N.array([float(header["CD1_1"]) , float(header["CD2_2"])])
    # Returning the filled instance
    return w

"""
MAIN
"""

parser = argparse.ArgumentParser(description='Image photometric completeness')
parser.add_argument('-i', '--input_image', required=True, help='required \
                    FITS file image', metavar='INPUT_IMAGE.fits')
parser.add_argument('--hdu_number', '--hdu', default=0, type=int, 
                    help='HDU number for the FITS images')
parser.add_argument('--mag_range_min', default=21, type=int, choices=range(-10, 30), 
                    help='Minimun value for the magnitude range')
parser.add_argument('--mag_range_max', default=30, type=int, choices=range(-10, 30), 
                    help='Maximun value for the magnitude range')
parser.add_argument('--model_type', default='STAR', help='Type of \
                    model objects, this version accepts STAR or GALAXY')
parser.add_argument('-w','--weight_map_file', default='', help='Weight map \
                     FITS to be applied in the image', metavar='WEIGHT_MAP.fits')
parser.add_argument('-seg','--segmentation_map_file', default='', help='Segmentation map \
                     Map to avoid object areas', metavar = '')
parser.add_argument('--mask', default="/home/CEFCA/clsj/Trabajo/Proyectos/ALHAMBRA/catalogs/Mascaras/allsky.pol", help='Masked area \ MANGLE file with the high queality area')
parser.add_argument('-sex','--sex_file', default='', help='SExtractor file\
                     Original SExtractor parameters', metavar = '')
parser.add_argument('--psf_models', default='/home/walterjr/Documents/Completeness/HST.PSF.0.065.23x23.fits', 
                    help='PSF models FITS that was applied to the models', metavar='PSF_MODELS.fits')
parser.add_argument('--psf_image', default='', help='PSF models FITS \
                    to be applied to the models', metavar='PSF_IMAGE.fits')
parser.add_argument('--models_path', default='/home/walterjr/Documents/Completeness/Models/Stars/', 
                    help='path to where all the models images are')
parser.add_argument('--mag_zero_point', default=None) #TODO: do a if for the instrumental or observed mag later
parser.add_argument('--models_per_mag', default=100, type = int,
		    help = 'Number of models per magnitude bin.')
parser.add_argument('--iter', default=2, type = int,
		    help = 'Number of iterations.')

#check the pixscale from image header, if not make the user provide it with a error message

#command line example:
"""
command_line = '-i /home/walter/Dropbox/Completeness/f02p01_F814W_2.swp.fits \
                -w /home/walter/Dropbox/Completeness/f02p01_F814W_2.swp.weight.fits \
                --psf_image /home/walter/Dropbox/Completeness/f02p01_F814W_2.swp.psf.fits \
                --mag_zero_point 30.64'

python completeness_WA_version0_6_1.py -i /Users/albertomolino/doctorado/photo/imagenes/f08p01_F814W_1.swp.fits -w /Users/albertomolino/doctorado/photo/imagenes/f08p01_F814W_1.swp.weight.fits --psf_image /Users/albertomolino/doctorado/photo/imagenes/f08p01_F814W_1.swp.psf.fits --mag_zero_point 30.64

python completeness_WA_version0_7_1.py -i ../Images/F814W/alhambra.f02p01_F814W_1.swp.fits --psf_image PSF/f02p01_F814W_1.swp.psf.fits --mask ../catalogs/Mascaras/ALHAMBRA_v5_mask/mask_f02p01c01.pol --w ../Images/F814W/f02p01_F814W_1.swp.weight.fits -s SegMaps/f02p01_F814W_1.swp.seg.unit.fits --mag_zero_point 30.64 --mag_range_min 21 --mag_range_max 27 --models_per_mag 10 --iter 10 --models_path modelos/

python completeness_WA_version0_7_2.py -i A2589.riz.sci.fits --psf_image A2589.PSF.fits --mask allsky.pol --w A2589.riz.wht.fits --segmentation_map_file A2589.riz.segm.fits --mag_zero_point 30.64 --mag_range_min 21 --mag_range_max 27 --models_per_mag 10 --iter 10 --model_type STAR --sex_file doublecluster.sex

"""
args = parser.parse_args()


"""
Parse the input parameters, define and declare variables
"""
#args = parser.parse_args()
settings = vars(args)

nHDU = settings['hdu_number']
imageinHDUList = fits.open(settings['input_image'])
imageinHDU = imageinHDUList[nHDU]
imageinHDR = imageinHDU.header
w = readWCS(imageinHDR)
imageinData = imageinHDU.data
imageinPixscale = float(imageinHDR['HIERARCH OAJ INS PIXSCALE'])
imageinPSF = settings['psf_image']
imageinWeight = settings['weight_map_file']
dataSeg = fits.open(settings['segmentation_map_file'])[0].data
dataSex = settings['sex_file']
modelsFolder = settings['models_path']
if modelsFolder[-1] == "/": modelsFolder = modelsFolder[:-1]
modelsOrigPSF = settings['psf_models']
magZeroPoint = float(settings['mag_zero_point'])
MODELS_PER_MAG = int(settings['models_per_mag'])
iterations = int(settings['iter'])
modelsType = settings['model_type'] # right now, 'STAR' (from psf)
                                    # or 'GALAXY' chef galaxies
####
# I read the mask to avoid bad regions. By default it reads a full sky mask.
maskfile = settings['mask']
#print maskfile
mng = pymangle.Mangle(maskfile)
####

modelsSexTable = []
modelsInsertedTable = []
modelsOutTable = []
imageoutData = imageinData.copy()
modelInsertedID = 0

imageinNy = N.shape(imageinData)[0]
imageinNx = N.shape(imageinData)[1]

"""
Models handling.
NOTE: This version is for CHEFS galaxy models models.
"""

if(modelsType=='GALAXY'):
    if not(os.access(modelsFolder+"/summary.model",os.F_OK)):
        print "Creating the summary file of the models..."
        imagelist = [image for image in commands.getoutput("ls "+modelsFolder).split() if ".fits" in image]
        s = "# Model Flux\n"
        for image in imagelist:
        	Data = fits.open(modelsFolder+"/"+image)[0].data
        	flux = N.sum(Data)
        	s += image.split("_")[0]+"\t"+str(flux)+"\n"
        f = open(modelsFolder+"/summary.model", "a")
        f.write(s)
        f.close()	
    # The fluxes of the models in the folder!
    datamodel = N.loadtxt(modelsFolder+"/summary.model")
    magsmodel = N.where(datamodel[:,1] >= 0, -2.5*N.log10(datamodel[:,1]) + 25.5, 0)

dm = 0.25
magScale = N.arange(float(settings['mag_range_min'])-magZeroPoint,float(settings['mag_range_max'])-magZeroPoint,dm)

dataPSF = fits.open(imageinPSF)[0].data
fluxPSF = N.sum(dataPSF)
#print fluxPSF

if(modelsType=='STAR'):
    modelDataOrig = fits.open(imageinPSF)[0].data
    
    modelNy = N.shape(modelDataOrig)[0]
    modelNx = N.shape(modelDataOrig)[1]
    
    maxFlux = N.max(modelDataOrig)
    countsOrig = N.sum(modelDataOrig) # total flux of the model.
    
    modelDataOrig = modelDataOrig/maxFlux #normalizing with the max. Then PSF_CUt makes sense as a relative value.
    for value in N.nditer(modelDataOrig, op_flags=['readwrite']): #setting zero low values (0.01) (why????)
        if value < PSF_CUT:
            value[...] = 0
    
    modelDataExtended = N.zeros((modelNy+10,modelNx+10))
    sumModelIntoImage(5, 5, modelNx, modelNy, modelDataOrig, modelDataExtended)
    
    modelNy = N.shape(modelDataExtended)[0]
    modelNx = N.shape(modelDataExtended)[1]
    
    modelDataExtended = modelDataExtended*maxFlux #de-normalize
    counts = N.sum(modelDataExtended)
    
   
    
    #build flux scales for a model based on its counts and the desired magnitude range
    #fluxScales = buildFluxScales(counts, magBin=0.05, magMin=20, magMax=26, magZero=30.64)
    fluxScales, mags = buildFluxScales(counts, magBin=0.25, magMin=float(settings['mag_range_min']), magMax=float(settings['mag_range_max']), magZero=magZeroPoint)
    print 'number of fluxscales: ',len(fluxScales)
    ids = range(len(fluxScales))
    datamodel = N.column_stack((ids,fluxScales))
    magsmodel = N.column_stack((ids,mags+magZeroPoint))
    print datamodel, magsmodel
    #pausa = raw_input('paused')


print 'Entering Models insertion...'

"""
Models insertion into the field image, keeping track of them by IDs.
"""

# Matrix to store the results
fcomp = N.zeros((iterations, len(magScale)))
print N.shape(fcomp)
fbias = N.zeros((iterations, len(magScale)))
radius = 12   # Object exclusion radius

for i in range(iterations):
    print "Starting iteration ", i+1
    modelsSexTable = []
    modelsInsertedTable = []
    modelsOutTable = []
    modelMagID = 0
    # I reset the image and the segmenation map in each interation
    imageoutData = imageinData.copy()
    dataSegout = dataSeg.copy()
    #print magScale
    #print magsmodel
    #print magZeroPoint

    for magnitude in magScale:
        # Targeted flux
        fluxMag = 10**(-0.4*magnitude)
        #print magnitude, magnitude - dm/2. + magZeroPoint, magnitude + dm/2. + magZeroPoint
        # The models that are in the desired magnitude range
        modelTargets = N.where((magsmodel >= magnitude - dm/2. + magZeroPoint) & (magsmodel < magnitude + dm/2. + magZeroPoint))[0]
        
        
        print "\t", magnitude + magZeroPoint, " ..."
        print "\t", N.where((magsmodel >= magnitude - dm/2. + magZeroPoint) & (magsmodel < magnitude + dm/2. + magZeroPoint))
        
        #if(len(modelTargets))>0:
        for j in range(MODELS_PER_MAG):
        # I read the model and scale to the targeted flux
            ModelIndex = random.choice(modelTargets)
            if(modelsType=='STAR'):
                modelDataOrig = modelDataExtended
            else:
                modelDataOrig = fits.open(modelsFolder+"/"+str(int(datamodel[ModelIndex][0]))+"_model.fits")[0].data
            modelData = modelDataOrig*datamodel[ModelIndex][1]
            #print N.sum(modelDataExtended)
            flux = N.sum(modelData)
            #print N.sum(modelDataOrig), flux
            
    
            #print -2.5*N.log10(flux) + magZeroPoint
            
            # PSF convolution and final shape of the model
            modelDataPSF = convolve2d(modelData, dataPSF)
            #modelNy = N.shape(modelDataPSF)[0]
            #modelNx = N.shape(modelDataPSF)[1]
            
            modelNy = N.shape(modelData)[0]
            modelNx = N.shape(modelData)[1]
            	
            # This includes the Pisson noise in the model.
            # I checked that several "poisson" curves are very close to the "no poisson" one. So I skip this step.
            #modelDataP = N.random.poisson(modelData*GAIN)/GAIN

            # random position of the source, leaving a strip near the boundaries
            posx = N.random.randint(radius,imageinNx - radius)    
            posy = N.random.randint(radius,imageinNy - radius)
            # RA and DEC of the source. Needed to apply the mangle mask
            ra, dec = w.wcs_pix2world(posx+1+modelNx/2, posy+1+modelNy/2, 1)
            
            while(not(isInsideImage(posx,posy,modelNx,modelNy,imageinNx,imageinNy) & (mng.polyid([ra], [dec])[0] >= 0) & (N.sum(dataSeg[posy+modelNy/2-radius:posy+modelNy/2+radius, posx+modelNx/2-radius:posx+modelNx/2+radius]) == 0) )):
                # If the position is not OK, new position :D
                posx = N.random.randint(radius,imageinNx - radius)    
                posy = N.random.randint(radius,imageinNy - radius)
                ra, dec = w.wcs_pix2world(posx+1+modelNx/2, posy+1+modelNy/2, 1)
                		
            # Add the model to the segmentation map with value 2
            dataSegout[posy+modelNy/2-radius:posy+modelNy/2+radius, posx+modelNx/2-radius:posx+modelNx/2+radius] = 2
            sumModelIntoImage(posx, posy, modelNx, modelNy, modelData, imageoutData)
            modelInsertedID += 1

            modelsSexTable.append([int(posx+1+modelNx/2), int(posy+1+modelNy/2), modelMagID, modelInsertedID, magnitude])
            modelsInsertedTable.append([modelInsertedID, 0, modelMagID, magnitude])
            
        modelMagID += 1
            #modelsOutTable.append([0, modelMagID, magnitude, 0.0, 0.0, magnitude])

    #print dataSegout
    fits.writeto("prueba_seg.fits",dataSegout, fits.open(settings['segmentation_map_file'])[0].header,clobber=True)
	
    """
    Handling and running sextractor on the field image+model objects
    """
    input_image_path, _ = os.path.splitext(settings['input_image'])
    imageout = input_image_path+'toSext2.fits'
    imageoutHeader = imageinHDR
    fits.writeto(imageout,imageoutData,imageoutHeader,clobber=True)
    
    #sys.exit(0)
    
    seeingFWHM_pix, seeingFWHM_arcsec = get_psffwhm(imageinPSF,imageinPixscale)
    psfArea = N.pi*((seeingFWHM_pix)/2)**2
    
    generateSexInputModels(modelsSexTable)
    generateSexParams()
    generateSexConfig('matched.cat', 'assoc.sex',
    'assoc.param', 'input_model.coo', imageinWeight, dataSex,
    7.)
    
    sexConfigFile = 'assoc.sex'
    sexOutputFile = 'matched.cat'
    sexCMD = 'sex '+imageout+' -c '+sexConfigFile+' -MAG_ZEROPOINT '+str(magZeroPoint)
    print 'sexCMD', sexCMD
    os.system(sexCMD)
    
    #sys.exit(0)

    """
    From the sextractor output file, count the retrieved models from their id,
    and calculate the completeness ratios.
    """
    sexOutputData = N.loadtxt(sexOutputFile)
    print modelMagID
    for magID in range(modelMagID):
        okmagid = N.where(sexOutputData[:,5] == magID)
        Ndetect = len(set(sexOutputData[:,6][okmagid]))
        fdetect = float(Ndetect)/float(MODELS_PER_MAG)
        print magID, Ndetect, fdetect
        if fdetect <= 1.:
            fcomp[i][magID] = fdetect
        else:
            fcomp[i][magID] = 1.
            
        fbias[i][magID] = N.median(sexOutputData[:,7][okmagid] + magZeroPoint - sexOutputData[:,-1][okmagid])


"""
Handling output
"""
magType = 'Instrumental'
if magZeroPoint is not None:
   magType = 'Observed'
   magScale += magZeroPoint
  
comp = N.mean(fcomp, axis = 0)
compe = N.std(fcomp, axis = 0)
compe = N.where(compe < 0.005, 0.005, compe)  # To avoid "0" errors
bias = N.mean(fbias, axis = 0)

print 'mags: ', list(magScale)
print 'comp: ', list(comp)
print "compe: ", list(compe)
print "bias: ", list(bias)

f = open("biasfile.txt", "w")
f.write(str(N.median(mags))+" "+str(N.mean(compe))+"\n")
f.close()

imageinHDUList.close()
