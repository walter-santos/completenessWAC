#! /usr/bin python
# -*- coding: utf-8 -*-

"""
Name: completeness_WAC
Version 0.8 (first release for J-PLUS, beta, point-like source version)
Description: Provide completeness percentages, over a range of magnitudes
and types of objects, for a given image
Authors: Walter Santos, Alberto Molino, Carlos San Juan
Created on: May 14th 2015
Last Updated: Nov 18th 2016
Instructions: Run the code with a -h/--help option to list all the arguments,
necessary and optional, on the command line.
Latest Changes:
- Adapted to run on J-PLUS images, either with stars (already tested) or galaxies
- Cleanup legacy old code
Requirements: Numpy, SciPy, Astropy, PyRAF, Sextractor, pymangle
TODO:
- Revisit galaxy models
- Implement more safeguards, exception handlings, ...
- Define a ___main___ so people can import the code elsewhere
- Follow the defined conventions for a python code, to be compatible with 
the rest of the J-PAS pipeline codes
"""

"""
IMPORTS
"""
import argparse
import os
import commands
import numpy as N
from astropy.io import fits
from astropy import wcs
from scipy.signal import convolve2d
import random
#import matplotlib.pyplot as plt
import pymangle

SExroot = '/usr/share/sextractor/'

"""
CONSTANTS
TODO: Probably move them from script constants to input parameters in the
main.
"""
#MODEL_OBJECTS_PIXSCALE = 0.065 #arcsec/pixel
WEIGHT_MAP_CUT = 0.5
PSF_CUT = 0.01
#GAIN = 13.42
MODELS_ZERO_POINT = 25.5

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
    sexdata.close()
    
    sexConfig = sexdict

    #... then change the information with some input, image-specific, info
    sexConfig['CATALOG_NAME']=outFile
    sexConfig['PARAMETERS_NAME']=paramsFile
    sexConfig['ASSOC_RADIUS']=str(matchRadius)
    sexConfig['ASSOC_NAME']=inputModelsfile
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

def isInsideImage(x, y, nx, ny, imageNx, imageNy):
    """
    Simple check o see if a model image fits in a field image, based on their
    size
    """
    return ( ((x+nx) < imageNx) and ((y+ny) < imageNy) )

def sumModelIntoImage(x, y, nx, ny, modelData, imageData):
    """
    Insert (sum) model image into the field image
    """
    imageData[y:y+ny,x:x+nx] += modelData
    
def outputTable(xpoints, ypoints, ypointse, outfile='completeness.txt', magType='Instrumental'):
    """
    Output the completeness results (the 'mean' values in magnitude bins) in text file
    """
    fout = open(outfile, 'w')
    fout.write('# '+magType+'F814W fc fce\n')

    for i,_ in enumerate(xpoints):
        fout.write(str(xpoints[i])+' '+"%1.3f" %ypoints[i]+' '+"%1.3f" %ypointse[i]+'\n')
    fout.close()
    
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
parser.add_argument('--mag_range_min', default=17, type=int, choices=range(10, 30), 
                    help='Minimun value for the magnitude range')
parser.add_argument('--mag_range_max', default=22, type=int, choices=range(10, 30), 
                    help='Maximun value for the magnitude range')
parser.add_argument('--model_type', default='STAR', help='Type of \
                    model objects, this version accepts STAR or GALAXY')
parser.add_argument('-w','--weight_map_file', default='', help='Weight map \
                     FITS to be applied in the image', metavar='WEIGHT_MAP.fits')
parser.add_argument('-seg','--segmentation_map_file', default='', help='Segmentation map \
                     Map to avoid object areas', metavar = '')
parser.add_argument('--mask', default="allsky.pol", help='Masked area \ MANGLE file with the high queality area')
parser.add_argument('-sex','--sex_file', default='', help='SExtractor file\
                     Original SExtractor parameters', metavar = '')
parser.add_argument('--psf_image', default='', help='PSF models FITS \
                    to be applied to the models', metavar='PSF_IMAGE.fits')
parser.add_argument('--models_path', default='models/', 
                    help='path to where all the models images are')
parser.add_argument('--mag_zero_point', default=None) #TODO: do a if for the instrumental or observed mag later
parser.add_argument('--models_per_mag', default=10, type = int,
		    help = 'Number of models per magnitude bin.')
parser.add_argument('--iter', default=2, type = int,
		    help = 'Number of iterations.')

#command line example:
"""
python completeness_WAC.py -i Alberto/1500029-A2593_rSDSS_swp.fits --mag_range_min 17 --mag_range_max 22 --model_type STAR --weight_map_file Alberto/1500029-A2593_rSDSS_swpweight.normed.fits --segmentation_map_file Alberto/A2593.riz.segm.fits --mask allsky.pol --sex_file Alberto/doublecluster.sex --psf_image Alberto/A2589.riz.sci.psf.fits --mag_zero_point 23.133 --models_per_mag 10 --iter 2
"""
args = parser.parse_args()


"""
Parse the input parameters, define and declare variables
"""
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
magZeroPoint = float(settings['mag_zero_point'])
MODELS_PER_MAG = int(settings['models_per_mag'])
iterations = int(settings['iter'])
modelsType = settings['model_type'] # right now, 'STAR' (from psf)
                                    # or 'GALAXY' chef galaxies
####
# I read the mask to avoid bad regions. By default it reads a full sky mask.
maskfile = settings['mask']
mng = pymangle.Mangle(maskfile)
####

modelsSexTable = []
modelsInsertedTable = []
imageoutData = imageinData.copy()
modelInsertedID = 0

imageinNy = N.shape(imageinData)[0]
imageinNx = N.shape(imageinData)[1]

"""
Models handling.
NOTE: This version is for CHEFS galaxy models models.
"""
i=0
if(modelsType=='GALAXY'):
    if not(os.access(modelsFolder+"/summary.model",os.F_OK)):
        print "Creating the summary file of the models..."
        imagelist = [image for image in commands.getoutput("ls "+modelsFolder).split() if ".fits" in image]
        s = "# Model Flux\n"
        for image in imagelist:
            Data = fits.open(modelsFolder+"/"+image)[0].data
            flux = N.sum(Data)
            s += str(i)+"\t"+str(flux)+"\t"+image+"\n"
            i=i+1
        f = open(modelsFolder+"/summary.model", "a")
        f.write(s)
        f.close()	
    # The fluxes of the models in the folder!
    datamodel = N.genfromtxt(modelsFolder+"/summary.model", dtype=str)
    magsmodel = N.where(datamodel[:,1].astype(float) >= 0, -2.5*N.log10(datamodel[:,1].astype(float)) + 25.5, 0)

dm = 0.25
magScale = N.arange(float(settings['mag_range_min'])-magZeroPoint,float(settings['mag_range_max'])-magZeroPoint,dm)

dataPSF = fits.open(imageinPSF)[0].data
fluxPSF = N.sum(dataPSF)

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
    fluxScales, mags = buildFluxScales(counts, magBin=0.25, magMin=float(settings['mag_range_min']), magMax=float(settings['mag_range_max']), magZero=magZeroPoint)
    ids = range(len(fluxScales))
    datamodel = N.column_stack((ids,fluxScales))
    magsmodel = N.column_stack((ids,mags+magZeroPoint))

print 'Entering Models insertion...'

"""
Models insertion into the field image, keeping track of them by IDs.
"""

# Matrix to store the results
fcomp = N.zeros((iterations, len(magScale)))
fbias = N.zeros((iterations, len(magScale)))
radius = 12   # Object exclusion radius

for i in range(iterations):
    print "Starting iteration ", i+1
    modelsSexTable = []
    modelsInsertedTable = []
    modelMagID = 0
    # I reset the image and the segmenation map in each interation
    imageoutData = imageinData.copy()
    dataSegout = dataSeg.copy()

    for magnitude in magScale:
        # Targeted flux
        fluxMag = 10**(-0.4*magnitude)
        # The models that are in the desired magnitude range
        modelTargets = N.where((magsmodel[:,1] >= magnitude - dm/2. + magZeroPoint) & (magsmodel[:,1] < magnitude + dm/2. + magZeroPoint))[0]
        
        print "\t", magnitude + magZeroPoint, " ..."
        
        if(len(modelTargets))>0:
            for j in range(MODELS_PER_MAG):
            # I read the model and scale to the targeted flux
                ModelIndex = random.choice(modelTargets)
                if(modelsType=='STAR'):
                    modelData = modelDataExtended*N.float(datamodel[ModelIndex][1])
                else:
                    modelData = fits.open(modelsFolder+"/"+datamodel[ModelIndex][2])[0].data
                
                flux = N.sum(modelData)
                        
                # PSF convolution and final shape of the model
                modelDataPSF = convolve2d(modelData, dataPSF)
                
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

    print 'Number of inserted objects: ', len(modelsInsertedTable)
    fits.writeto("prueba_seg.fits",dataSegout, fits.open(settings['segmentation_map_file'])[0].header,clobber=True)
	
    #Handling and running sextractor on the field image+model objects
    input_image_path, _ = os.path.splitext(settings['input_image'])
    imageout = input_image_path+'toSext2.fits'
    imageoutHeader = imageinHDR
    fits.writeto(imageout,imageoutData,imageoutHeader,clobber=True)
    
    generateSexInputModels(modelsSexTable)
    generateSexParams()
    generateSexConfig('matched.cat', 'assoc.sex',
    'assoc.param', 'input_model.coo', imageinWeight, dataSex,
    7.)
    
    sexConfigFile = 'assoc.sex'
    sexOutputFile = 'matched.cat'
    sexCMD = 'sex '+imageout+' -c '+sexConfigFile+' -MAG_ZEROPOINT '+str(magZeroPoint)
    os.system(sexCMD)

    
    #From the sextractor output file, count the retrieved models from their id,
    #and calculate the completeness ratios.
    sexOutputData = N.loadtxt(sexOutputFile)
    for magID in range(modelMagID):
        okmagid = N.where(sexOutputData[:,5] == magID)
        Ndetect = len(set(sexOutputData[:,6][okmagid]))
        fdetect = float(Ndetect)/float(MODELS_PER_MAG)
        if fdetect <= 1.:
            fcomp[i][magID] = fdetect
        else:
            fcomp[i][magID] = 1.
            
        fbias[i][magID] = N.median(sexOutputData[:,7][okmagid] + magZeroPoint - sexOutputData[:,-1][okmagid])


"""
Handling output
"""
comp = N.mean(fcomp, axis = 0)
compe = N.std(fcomp, axis = 0)
compe = N.where(compe < 0.005, 0.005, compe)  # To avoid "0" errors
bias = N.mean(fbias, axis = 0)

print 'mags: ', list(magScale)
print 'comp: ', list(comp)
print "compe: ", list(compe)
print "bias: ", list(bias)

f = open("compfile.txt", "w")
for i,_ in enumerate(comp):
        f.write(str(magScale[i]+magZeroPoint)+' '+str(comp[i])+' '+str(compe[i])+'\n')
f.close()

f = open("biasfile.txt", "w")
for i,_ in enumerate(bias):
        f.write(str(magScale[i]+magZeroPoint)+' '+str(bias[i])+'\n')
f.close()

imageinHDUList.close()