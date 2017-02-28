#/!/usr/bin/python

from timeit import default_timer as timer
import sys,  os
import dicom
import numpy
from matplotlib import pyplot
from matplotlib.widgets import Slider
from pylab import * 

print("")
print("")
print("Initializing R2STAR for LIVER data visualization tool")
print("Scripted for Dept. Central Radiology, MRI Unit")
print("Authored by Yanlu Wang, MR Physics, Karolinska University Hospital, Solna, Sweden")
print("02-2017")
print("")

start = timer()

DataPath = "/home/doggydaddy/Documents/pydicom_test/Data_From_Aera/"
lstFilesDCM = [] #create empty list
for dirName,  subdirList,  fileList in os.walk(DataPath):
    for filename in fileList:
        if ".dcm" in filename.lower(): 
            lstFilesDCM.append(os.path.join(dirName, filename))

print("MESSAGE: Performing data check ...")
if len(lstFilesDCM) != 75:
     print("ERROR: Sanity check failed! Unexpected number of DICOM images in the folder specified!")
     print("ERROR: Expected 75 (15 images x 5 sub-folers), got " + str(len(lstFilesDCM)))
     sys.exit(os.EX_DATAERR)
print("MESSAGE: number of files seems to be a-ok!")

TEs = []
Slices = []
for file in lstFilesDCM:
    curD = dicom.read_file(file,  stop_before_pixels=True)
    TE = curD[0x18, 0x81].value
    Slice = curD[0x20, 0x32].value[2]
    if TE not in TEs:
        TEs.append(TE)
    if Slice not in Slices:
        Slices.append(Slice)

TEs.sort()
Slices.sort()
print("MESSAGE: Performing data check ...")
if (len(TEs) != 15):
    print("ERROR: Sanity check failed! Number of TEs do not match!")
    print("ERROR: Expected 15, got " + len(TEs))
    sys.exit(os.EX_DATAERR)
if (len(Slices) != 5 ):
    print("ERROR: Sanity check failed! Number of Slices do not match!")
    print("ERROR: Expected 5, got " + len(Slices))
    sys.exit(os.EX_DATAERR)
print("MESSAGE: ... data seems to be a-ok!")

print(TEs)
print(Slices)

RefDs = dicom.read_file(lstFilesDCM[0])
# Load dimension based on the number of rows, columns, and slices
ConstPixelDims = (int(RefDs.Rows),  int(RefDs.Columns), 5,  15)
# The array is sized based on 'ConstPixelDims'
ArrayDicom = numpy.zeros(ConstPixelDims, dtype=RefDs.pixel_array.dtype)
# loop through all the DICOM files
for filenameDCM in lstFilesDCM:
    # read the file
    ds = dicom.read_file(filenameDCM)
    TE = ds[0x18, 0x81].value
    Slice = ds[0x20, 0x32].value[2]
    idx_TE = TEs.index(TE)
    idx_Slice = Slices.index(Slice)
    if idx_TE is 0:
        if idx_Slice is 0:
            RefSlice1 = dicom.read_file(filenameDCM)
        elif idx_Slice is 1:
            RefSlice2 = dicom.read_file(filenameDCM)
        elif idx_Slice is 2:
            RefSlice3 = dicom.read_file(filenameDCM)
        elif idx_Slice is 3:
            RefSlice4 = dicom.read_file(filenameDCM)
        elif idx_Slice is 4:
            RefSlice5 = dicom.read_file(filenameDCM)
    # store the raw image data
    ArrayDicom[:, :, idx_Slice,  idx_TE] = ds.pixel_array 

def signal(data,  slice_number,  te_number):
    if slice_number >= len(Slices):
        slice_number = len(Slices) - 1
    if te_number >= len(TEs):
        te_number = len(TEs) - 1
    return data[:, :, slice_number,  te_number]

def signal_plotdata(dat, id_x, id_y, id_z):
    return dat[id_y, id_x, id_z, :] # note x and y is switched

fig = pyplot.figure()
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)
fig.subplots_adjust(bottom=0.25)

# Draw the plot
slice_no = 0
te_no = 0
pyplot.set_cmap(pyplot.gray())
img = ax1.imshow(ArrayDicom[:, :, slice_no, te_no])

coords = [30,  42,  0] # set up init values
x = [float(i) for i in TEs]
data_x = numpy.asarray(x)
data_y = signal_plotdata(ArrayDicom,  coords[0],  coords[1],  coords[2])
#line = pyplot.matplotlib.lines.Line2D( data_x,  data_y )
#ax2.add_line( line )
scat = ax2.scatter( data_x,  data_y )
ax2.set_autoscale_on(True)
fig.show()

# Add two sliders for tweaking the parameters
slice_slider_ax  = fig.add_axes([0.25, 0.15, 0.65, 0.03])
slice_slider = Slider(slice_slider_ax, 'Slice', 0.0, 5.0, valinit=slice_no,  valfmt='%0.0f')
te_slider_ax = fig.add_axes([0.25, 0.1, 0.65, 0.03])
te_slider = Slider(te_slider_ax, 'TE', 0.0, 15.0, valinit=te_no,  valfmt='%0.0f')

def sliders_on_changed(val):
    img.set_data(signal(ArrayDicom,  int(slice_slider.val), int(te_slider.val)))
    fig.canvas.draw()

def onclick(event):
    event.xdata, event.ydata
    coords = [int(event.xdata), int(event.ydata),  int(slice_slider.val)]
    #line.set_ydata(signal_plotdata(ArrayDicom, coords[0], coords[1], coords[2]))
    data_y = signal_plotdata(ArrayDicom, coords[0], coords[1], coords[2])
    new_offsets = list(zip(data_x,  data_y))
    scat.set_offsets( new_offsets )

    corners = (min(data_x), min(data_y)), (max(data_x), max(data_y))
    ax2.update_datalim(corners)
#    ax2.relim()    
    ax2.autoscale_view()
    fig.canvas.draw()

slice_slider.on_changed(sliders_on_changed)
te_slider.on_changed(sliders_on_changed)
cid = fig.canvas.mpl_connect('button_press_event', onclick)

pyplot.show()
