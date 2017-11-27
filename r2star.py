#!/usr/bin/python
import sys,  os,  getopt
import timeit
from timeit import default_timer as timer 
from scipy.optimize import curve_fit
import dicom
import numpy as np
import dicom.UID
import datetime
from lmfit import *
# import tqdm
import matplotlib.pyplot as plt

class DataWrapper:
    def __init__(self, inputFolder, outputFolder):
        self.inputDir = inputFolder
        self.outputDir = outputFolder

    def run(self):
        self.parseFiles()
        self.dataCheck01()
        self.parseHeaders()
        self.dataCheck02()
        self.parseData()

    def parseFiles(self):
        self.lstFilesDCM = []
        for dirName,  subdirList,  fileList in os.walk(self.inputDir):
            for filename in fileList:
                if ".dcm" in filename.lower(): 
                    self.lstFilesDCM.append(os.path.join(dirName, filename))

    def dataCheck01(self):
        print('')
        print("MESSAGE: Performing data check ...")
        if len(self.lstFilesDCM) != 75:
            print("ERROR: Sanity check failed! Unexpected number of DICOM images in the folder specified!")
            print("ERROR: Expected 75 (15 images x 5 sub-folers), got " + str(len(self.lstFilesDCM)))
            print('')
            sys.exit(os.EX_DATAERR)
        print("MESSAGE: number of files seems to be a-ok!")
        print('')

    def parseHeaders(self):
        self.TEs = []
        self.Slices = []
        for file in self.lstFilesDCM:
            curD = dicom.read_file(file, stop_before_pixels=True)
            TE = curD[0x18, 0x81].value
            Slice = curD[0x20, 0x32].value[2]
            if TE not in self.TEs:
                self.TEs.append(TE)
            if Slice not in self.Slices:
                self.Slices.append(Slice)
        self.TEs.sort()
        self.Slices.sort()

    def dataCheck02(self):
        print("MESSAGE: Performing data check ...")
        if (len(self.TEs) != 15):
            print("ERROR: Sanity check failed! Number of TEs do not match!")
            print("ERROR: Expected 15, got " + len(self.TEs))
            print('')
            sys.exit(os.EX_DATAERR)
        if (len(self.Slices) != 5 ):
            print("ERROR: Sanity check failed! Number of Slices do not match!")
            print("ERROR: Expected 5, got " + len(self.Slices))
            print('')
            sys.exit(os.EX_DATAERR)
        print("MESSAGE: ... data seems to be a-ok!")
        print('')

        print("DEBUG: TEs and Slice Locations:")
        print(self.TEs)
        print(self.Slices)
        print('')

    def parseData(self):
        self.RefDs = dicom.read_file(self.lstFilesDCM[0])
        # Load dimension based on the number of rows, columns, and slices
        self.ConstPixelDims = (int(self.RefDs.Rows),  int(self.RefDs.Columns), len(self.Slices), len(self.TEs))
        # The array is sized based on 'ConstPixelDims'
        self.ArrayDicom = np.zeros(self.ConstPixelDims, dtype=self.RefDs.pixel_array.dtype)
        # loop through all the DICOM files
        for filenameDCM in self.lstFilesDCM:
            # read the file
            ds = dicom.read_file(filenameDCM)
            TE = ds[0x18, 0x81].value
            Slice = ds[0x20, 0x32].value[2]
            idx_TE = self.TEs.index(TE)
            idx_Slice = self.Slices.index(Slice)
            if idx_TE is 0:
                if idx_Slice is 0:
                    self.RefSlice1 = dicom.read_file(filenameDCM)
                elif idx_Slice is 1:
                    self.RefSlice2 = dicom.read_file(filenameDCM)
                elif idx_Slice is 2:
                    self.RefSlice3 = dicom.read_file(filenameDCM)
                elif idx_Slice is 3:
                    self.RefSlice4 = dicom.read_file(filenameDCM)
                elif idx_Slice is 4:
                    self.RefSlice5 = dicom.read_file(filenameDCM)
            # store the raw image data
            self.ArrayDicom[:, :, idx_Slice,  idx_TE] = ds.pixel_array 
        self.RefSlices = [self.RefSlice1, self.RefSlice2, self.RefSlice3, self.RefSlice4, self.RefSlice5] 

    def getDims(self, idx):
        return int(self.ConstPixelDims[idx])

    def getOutput(self, pArrayOutput):
        self.ArrayOutput = pArrayOutput

    def writeDicomFromTemplate(self, pixel_array, filename, refd):
        ds = refd
        #ds.ContentDate = str(datetime.date.today()).replace('-','')
        ds.SecondaryCaptureDeviceManufctur = 'python'
        ds.SeriesNumber = '99'
        tmpNr = ds.AccessionNumber
        ds.StudyID = tmpNr[-2:]
        ds.AccessionNumber = tmpNr[:-2]
        ds.ProtocolName = 'R2Star_Map_minDelta'
        
        n = len(ds.SOPInstanceUID)
        if n != 52:
            print('WARNING: Is this really a dataset from a Siemens Aera?')
        uuid = dicom.UID.generate_uid()
        ds.SOPInstanceUID = ds.SOPInstanceUID[0:35] + uuid[ (len(uuid)-18):(len(uuid)-1) ]

        ds.SamplesPerPixel = 1
        ds.PhotometricInterpretation = "MONOCHROME2"
        ds.PixelRepresentation = 0
        ds.HighBit = 15
        ds.BitsStored = 16
        ds.BitsAllocated = 16
        ds.SmallestImagePixelValue = np.min(pixel_array)
        ds.LargestImagePixelValue = np.max(pixel_array)
        if pixel_array.dtype != np.uint16:
            pixel_array = np.round(pixel_array, 0)
            pixel_array = pixel_array.astype(np.uint16)
        ds.PixelData = pixel_array.tostring()
        ds.save_as(filename)

    def writeOutput(self):
        print('MESSAGE: Exporting output maps ...')
        if not os.path.isdir(self.outputDir):
            os.makedirs(self.outputDir)
        for i in range(0, 5):
            filename = '\\slice' + str(i+1) + '.dcm'
            self.writeDicomFromTemplate(self.ArrayOutput[:,:,i,1], self.outputDir+filename, self.RefSlices[i])
        print('MESSAGE: ... export completed!')
        print('')

def fit_fct(params,x,y):   #cost function to minimize
    a=params['a']
    r2s=params['r2s']
    N=params['N']
    y_est=a*np.exp(-r2s*x)+N
    return y_est - y      

def fit_minimizer(x,y,seeds,lbounds,hbounds):   #lmfit optimiser
    params=Parameters() 
    params.add('a',value=seeds[0],min=lbounds[0],max=hbounds[0])   
    params.add('r2s',value=seeds[1],min=lbounds[1],max=hbounds[1])   
    params.add('N',value=seeds[2],min=lbounds[2],max=hbounds[2])   
    minner=Minimizer(fit_fct,params,fcn_args=(x,y))
    result = minner.minimize()
    return result.params['a'].value,result.params['r2s'].value,result.params['N'].value

class R2StarFitting:
    def __init__(self, i_nRows, i_nColumns, i_nSlices, i_TEs, i_ArrayDicom):
        self.nRows = int(i_nRows)
        self.nColumns = int(i_nColumns)
        self.nSlices = int(i_nSlices)
        self.nvxl=self.nRows*self.nColumns*self.nSlices
        self.TEs = i_TEs
        self.ArrayDicom = i_ArrayDicom

    def run(self):    #debug: self=r2sf   ; i=15000
        print('MESSAGE: Performing R2* fitting using non-parametric least squares ...')
        print('MESSAGE: Monoexponential with plateau')
        x = np.array([float(ii) for ii in self.TEs])
        y=self.ArrayDicom.reshape(self.nvxl,len(self.TEs))
        self.ArrayOutput=np.zeros((self.nvxl,3))
        # for i in tqdm.tqdm(range(self.nvxl)):
        for i in range(self.nvxl):
            if np.mean(y[i,0:5]) < 10.0:
                self.ArrayOutput[i] = 0
            else:
                seeds=[y[i,0],1,0]
                lbounds=[0,0,-10]
                hbounds=[100*np.mean(y[i,0:5]),5,10]  
                a,r2s,N=fit_minimizer(x,y[i],seeds,lbounds,hbounds)
                self.ArrayOutput[i]=a,r2s*1000,N
        self.ArrayOutput.shape=(self.nRows,self.nColumns,self.nSlices,3)
        print('MESSAGE R2* fitting completed!')
        print('')

def main(argv):
    inputFolder = r'C:\Users\User\Desktop\relaxometry\test_data\Mag_t1_fl2d_15p5o_graphcut_LEVER_minDelta_tra_bh_1'
    outputFolder = r'C:\Users\User\Desktop\relaxometry\output'
    try:
        opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
    except getopt.GetoptError:
        print( 'r2star.py -i <input directory> -o <output directory>' )
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('r2star.py -i <input directory> -o <output director>' )
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputFolder = arg
        elif opt in ("-o", "--ofile"):
            outputFolder = arg

    print( 'Input folder is "', inputFolder )
    print( 'Output folder is "', outputFolder )
    
    dw = DataWrapper(inputFolder, outputFolder)
    dw.run()
    r2sf = R2StarFitting(dw.getDims(0), dw.getDims(1), dw.getDims(2), dw.TEs, dw.ArrayDicom)
    r2sf.run()
    dw.getOutput(r2sf.ArrayOutput)
    dw.writeOutput()

if __name__ == "__main__":
    print("")
    print("")
    print("Initializing R2STAR analysis for Iron quantification in LIVER")
    print("Scripted for Dept. Central Radiology, MRI Unit")
    print("Authored by Yanlu Wang, MR Physics, Karolinska University Hospital, Solna, Sweden")
    print("02-2017")
    print("")
    main(sys.argv[1:])
    print("")
    print("Script finished sucessfully, exiting normally.")
    print("")
    # sys.exit()

