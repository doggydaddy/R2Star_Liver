#!/usr/bin/python
import sys,  os,  getopt
import timeit
from timeit import default_timer as timer 
from scipy.optimize import curve_fit
import dicom
import numpy
import dicom.UID
import datetime

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
        self.ArrayDicom = numpy.zeros(self.ConstPixelDims, dtype=self.RefDs.pixel_array.dtype)
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
        ds.SmallestImagePixelValue = numpy.min(pixel_array)
        ds.LargestImagePixelValue = numpy.max(pixel_array)
        if pixel_array.dtype != numpy.uint16:
            pixel_array = numpy.round(pixel_array, 0)
            pixel_array = pixel_array.astype(numpy.uint16)
        ds.PixelData = pixel_array.tostring()
        ds.save_as(filename)

    def writeOutput(self):
        print('MESSAGE: Exporting output maps ...')
        if not os.path.isdir(self.outputDir):
            os.makedirs(self.outputDir)
        for i in range(0, 5):
            filename = 'slice' + str(i+1) + '.dcm'
            self.writeDicomFromTemplate(self.ArrayOutput[:,:,i,1], self.outputDir+filename, self.RefSlices[i])
        print('MESSAGE: ... export completed!')
        print('')

class R2StarFitting:
    def __init__(self, i_nRows, i_nColumns, i_nSlices, i_TEs, i_ArrayDicom):
        self.nRows = int(i_nRows)
        self.nColumns = int(i_nColumns)
        self.nSlices = int(i_nSlices)
        self.TEs = i_TEs
        self.ArrayDicom = i_ArrayDicom
        OutputPixelDims = (self.nRows, self.nColumns, self.nSlices, 3)
        self.ArrayOutput = numpy.zeros(OutputPixelDims, dtype=float)

    def func(self, x, a, b, c):
        return a * numpy.exp(-x * b) + c

    def run(self):
        print('MESSAGE: Performing R2* fitting using non-parametric least squares ...')
        print('MESSAGE: Monoexponential with plateau')
        x = [float(ii) for ii in self.TEs]
        for s in range(0, self.ArrayDicom.shape[2]):
            start = timer()
            for i in range(0, self.ArrayDicom.shape[0]):
                for j in range(0,  self.ArrayDicom.shape[1]):
                    yn = self.ArrayDicom[i, j, s, :]
                    if numpy.mean(yn[0:3]) < 10.0:
                        self.ArrayOutput[i, j, s, :] = 0
                    else:
                        A_bound =  100 * numpy.mean(yn[0:2])
                        R2Star_bound = 5.0
                        C_bound = 10.0
                        bounds = [A_bound,  R2Star_bound,  C_bound]
                        try:
                            popt, pcov = curve_fit(self.func, x, yn,  bounds=(0.0, bounds))
                        except RuntimeError:
                            popt = [0., 0., 0.] 
                        self.ArrayOutput[i, j, s, 0] = popt[0]
                        self.ArrayOutput[i, j, s, 1] = popt[1]*1000
                        self.ArrayOutput[i, j, s, 2] = popt[2]
            end = timer()
            print('MESSAGE: Slice ' + str(s) + ' done! Time taken: ' + str(end - start) + ' seconds')
        print('MESSAGE R2* fitting completed!')
        print('')

def main(argv):
    inputFolder = ''
    outputFolder = ''
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
    sys.exit()

