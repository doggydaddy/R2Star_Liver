import numpy
import os
import dicom
from scipy.optimize import curve_fit
import matplotlib.pyplot as pyplot
from matplotlib.widgets import Slider
from pylab import *
from timeit import default_timer as timer

print("")
print("")
print("Initializing R2STAR for LIVER curve-fitting visualization tool")
print("Scripted for Dept. Central Radiology, MRI Unit")
print("Authored by Yanlu Wang, MR Physics, Karolinska University Hospital, Solna, Sweden")
print("02-2017")
print("")
print('')
print('fittingInspector.py')
print('not program proper, intended for debug use')
print('')
print('')

DataPath = '/home/doggydaddy/Documents/pydicom_test/Data_From_Aera/'
print('parsing folder: ' + DataPath)
lstFilesDCM = [] #create empty list
for dirName,  subdirList,  fileList in os.walk(DataPath):
    for filename in fileList:
        if ".dcm" in filename.lower(): 
            lstFilesDCM.append(os.path.join(dirName, filename))

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

print("DEBUG: TEs and Slice Locations:")
print(TEs)
print(Slices)
print('')

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
    # store the raw image data
    ArrayDicom[:, :, idx_Slice,  idx_TE] = ds.pixel_array 

# -----------------------------------------
# Performing the calculations
# -----------------------------------------
def func(x, a, b, c):
    return a * numpy.exp(-x * b) + c

x = [float(i) for i in TEs]
OutputPixelDims = (int(RefDs.Rows),  int(RefDs.Columns), 5,  3)
ArrayOutput = numpy.zeros(OutputPixelDims, dtype=float)
start = timer() # starting timer
for s in range(0, ArrayDicom.shape[2]):
    for i in range(0, ArrayDicom.shape[0]):
        for j in range(0,  ArrayDicom.shape[1]):
            yn = ArrayDicom[i, j, s, :]
            if numpy.mean(yn) < 10.0:
                ArrayOutput[i, j, s, :] = 0
            else:
                A_bound =  100 * numpy.mean(yn[0:2])
                R2Star_bound = 5.0
                C_bound = 10.0
                bounds = [A_bound,  R2Star_bound,  C_bound]
                err = 1e-02
                try:
                    popt, pcov = curve_fit(func, x, yn,  bounds=(0.0, bounds))
                except RuntimeError:
                    popt = [0., 0., 0.]
                #if any( popt >= [i - err for i in bounds] ):
                #    popt =[0., 0., 0.]
                ArrayOutput[i, j, s, 0] = popt[0]
                ArrayOutput[i, j, s, 1] = popt[1]*1000
                ArrayOutput[i, j, s, 2] = popt[2]
    # -----------------------------------------
end = timer()
print("DEBUG: Calculations done, time taken is:")
print(end-start)
print('')

# -----------------------------
# Plotting the results
# -----------------------------
def signal(data,  slice_number,  te_number):
    if slice_number >= data.shape[2]:
        slice_number = data.shape[2]- 1
    return data[:, :, slice_number,  te_number]

def signal_plotdata(dat, id_x, id_y, id_z):
    return dat[id_y, id_x, id_z, :]

def signal_plotmodel(dat, id_x, id_y, id_z,  x_range):
    model_params = dat[id_y, id_x, id_z, :]
    model_y = [func(i,model_params[0], model_params[1]/1000., model_params[2]) for i in x_range]
    model_y = numpy.asarray(model_y)
    return model_y

fig = pyplot.figure()
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)
fig.subplots_adjust(bottom=0.25)

slice_no = 0
pyplot.set_cmap(pyplot.gray())
img = ax1.imshow(signal(ArrayOutput, slice_no, 1),  vmax=1000.)

coords = [30,  42,  0] # set up init values
data_x = numpy.asarray(x)
data_y = signal_plotdata(ArrayDicom,  coords[0],  coords[1],  coords[2])
scat = ax2.scatter( data_x,  data_y )
model_x = numpy.linspace(0.1, 4.99)
model_y = signal_plotmodel(ArrayOutput,  coords[0],  coords[1],  coords[2],  model_x)
line = matplotlib.lines.Line2D( model_x,  model_y )
ax2.add_line( line )

# Add slider for slices and maps
slice_slider_ax  = fig.add_axes([0.25, 0.15, 0.65, 0.03])
slice_slider = Slider(slice_slider_ax, 'Slice', 0.0, 4.0, valinit=slice_no,  valfmt='%0.0f')
var_slider_ax = fig.add_axes([0.25, 0.1, 0.65, 0.03])
var_slider = Slider(var_slider_ax, 'Variable', 0.0, 2.0, valinit=1,  valfmt='%0.0f')

def sliders_on_changed(val):
    img.set_data(signal(ArrayOutput,  int(slice_slider.val),  int(var_slider.val)))
    if int(var_slider.val) == 1:
        img.set_clim(vmax=1000.)
    else:
        img.set_clim( vmax=numpy.max( ArrayOutput[:, :, int(slice_slider.val), int(var_slider.val)]) ) 
    fig.canvas.draw()

def onclick(event):
    event.xdata, event.ydata
    coords = [int(event.xdata), int(event.ydata),  int(slice_slider.val)]
    line.set_ydata(signal_plotmodel(ArrayOutput, coords[0], coords[1], coords[2], model_x))
    data_y = signal_plotdata(ArrayDicom, coords[0], coords[1], coords[2])
    new_offsets = list(zip(data_x,  data_y))
    scat.set_offsets( new_offsets )
    #corners = (min(data_x), min(data_y)), (max(data_x), max(data_y))
    #ax2.update_datalim(corners)
    ax2.relim()    
    ax2.autoscale_view()
    fig.canvas.draw()

slice_slider.on_changed(sliders_on_changed)
var_slider.on_changed(sliders_on_changed)
cid = fig.canvas.mpl_connect('button_press_event', onclick)

pyplot.show()
