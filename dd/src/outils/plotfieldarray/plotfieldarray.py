####################################
### ------------- Libraries --------
####################################
import numpy             as np
import matplotlib        as mpl
import matplotlib.cm     as cm
import matplotlib.pyplot as plt
import sys

from matplotlib.offsetbox import AnchoredText

####################################
### ------------- Input ------------
####################################

header      = True # Is there a header (1 line) ?
headerlines = 1    # Number of lines for the heqder
col_start   = 2    # Starting Column (python numbering)
col_end     = 7    # Ending Column (python numbering)
z_min       = -40  # Field rang min value
z_max       = 40   # Field rang max value
name        = ['A','B']   # Name of field components
nbcolplot    = 3   # Number of columns for plot array)

rot90  = False     # Rotation of 90
fliplr = False     # Flip left / right
flipud = False     # Flip up / down

# Label (/; for spaces)
ptitle = 'Field \; array \; plot'
xlab = 'Common \; label \; \sigma'
ylab = 'Common \; label'
cbarlab = 'Cbar \; lab'

# If file does not contain i,j indices
noindices = False
indi = 100 # number of points on i
indj = 100 # number of points on j

zoombox = False          # if true, set indices boundaries for plotting
box =[[70,215],[88,228]] # box = [[imin,imax],[jmin,jmax]]

####################################
### -------- Some functions --------
####################################

def reverse_colourmap(cmap, name = 'my_cmap_r'):
    """
    In:
    cmap, name
    Out:
    my_cmap_r

    Explanation:
    t[0] goes from 0 to 1
    row i:   x  y0  y1 -> t[0] t[1] t[2]
                   /
                  /
    row i+1: x  y0  y1 -> t[n] t[1] t[2]

    so the inverse should do the same:
    row i+1: x  y1  y0 -> 1-t[0] t[2] t[1]
                   /
                  /
    row i:   x  y1  y0 -> 1-t[n] t[2] t[1]
    """
    reverse = []
    k = []

    for key in cmap._segmentdata:
        k.append(key)
        channel = cmap._segmentdata[key]
        data = []

        for t in channel:
            data.append((1-t[0],t[2],t[1]))
        reverse.append(sorted(data))

    LinearL = dict(zip(k,reverse))
    my_cmap_r = mpl.colors.LinearSegmentedColormap(name, LinearL)
    return my_cmap_r

####################################
### ------------- Code -------------
####################################

#-------------------------
# Some input verifications
#--
if (len(sys.argv) == 1 ): sys.exit("!> You must set a filemname as first argument : python plotfieldarray 'filename'")
if (len(sys.argv)  > 2 ): sys.exit("!> Only one argument")

if (col_start < 2 and not noindices):
  print ">>>"
  print ">>> WARNING : If you do not have the 2 indices columns,  activate flag 'noindices' and fill indi and indj values <<< "
  print ">>>"
  quit()

filename = sys.argv[1]

#-----------
# Load file
#--
if (not header):
  headerlines = 0

field = np.loadtxt(filename,skiprows=headerlines)

#---------------------
# Set up array indices
#--
if (noindices): # user input indices

 if (indi*indj != np.shape(field)[0]):

   raise ValueError('!> indi * indj and file column length are not equal: indi*indj: %i , column length: %i '% (indi*indj,np.shape(field)[0]))

 else :

   imin = 0
   imax = indi
   jmin = 0
   jmax = indj

else : # use indices of the first 2 columns of the input file

  imin = np.int(np.min(field[:,0]))
  jmin = np.int(np.min(field[:,1]))

  if (imin < 0):
    imax = np.int(np.max(field[:,0])) - imin + 1 # max i indice
  else:
    imax = np.int(np.max(field[:,0])) # max i indice

  if (jmin < 0):
    jmax = np.int(np.max(field[:,1])) -jmin + 1 # max j indice
  else :
    jmax = np.int(np.max(field[:,1])) # max j indice


#---------------------
# Field reshape
#--
nbfieldcomp = col_end+1-col_start   # number of field component

# Field reshape
field = np.reshape(field[:,col_start:col_end+1],(imax,jmax,nbfieldcomp))

# Zoom box
if (zoombox):
  field = field[box[0][0]:box[0][1],box[1][0]:box[1][1]]
  imax = box[0][1]-box[0][0]
  jmax = box[1][1]-box[1][0]


#---------------------
# Mesh Grid
#--
x,y = np.meshgrid(np.linspace(0,imax, imax),np.linspace(0,jmax, jmax))

#---------------------
# Mesh rotations
#--
# 90deg rotation
if (rot90):
  field = np.rot90(field)
  x = np.rot90(x)
  y = np.rot90(y)

#Flip left/right
if (fliplr):
  field = np.fliplr(field)
  x = np.fliplr(x)
  y = np.fliplr(y)

#Flip up/down
if (flipud):
  field = np.flipud(field)
  x = np.flipud(x)
  y = np.flipud(y)

#---------------------
# Corner labels
#--
if (name == []):
  name=['F'+ np.str(i) for i in range(nbfieldcomp)] #defaut names
elif (len(name) <nbfieldcomp) :
  for i in range(len(name),nbfieldcomp):
   name.append('F'+ np.str(i)) #defaut names


#---------------------
# Subplots grid
#--
if (nbfieldcomp == 1) : nbcolplot = 1
nblineplot = np.int(np.ceil(nbfieldcomp/np.float(nbcolplot))) # Calculate number of lines for subplot array


#---------------------------------
# Generate subplots and attributes
#--
# SubplotsPlot
f, axes = plt.subplots(nblineplot,nbcolplot,figsize=(5,5),dpi=150,sharex='col', sharey='row')

num=0

#In case we only have 1 fied component, we need to transform axes into an array to be iterative
if (nbfieldcomp==1) : axes = np.array([axes])

for ax in axes.flat:

  ax.set_xlim([0,imax])
  ax.set_ylim([0,jmax])

  ax.set(adjustable='box-forced', aspect='equal')

  ax.tick_params(
    axis='both',          # changes apply to the axis
    which='both',         # both major and minor ticks are affected
    bottom='on',          # ticks along the bottom edge are on/off
    left='on',            # ticks along the bottom edge are on/off
    right='on',           # ticks along the bottom edge are on/off
    top='on',             # ticks along the top edge are on/off
    labelbottom='on',     # labels along the bottom edge are on/off
    labeltop='off',       # labels along the top edge are on/off
    labelleft='on',       # labels along the left edge are on/off
    labelright='off',     # labels along the right edge are on/off
    )

  if (num < nbfieldcomp) :
    anchored_text = AnchoredText(r"$\mathbf{"+name[num]+"}$", loc=1)
    ax.add_artist(anchored_text)
    fieldtest=field[:,:,num]
    num=num+1
    im = ax.pcolormesh(x.T, y.T, fieldtest, cmap= reverse_colourmap(cm.RdBu), vmin=z_min, vmax=z_max)
  else :
    f.delaxes(ax)


#---------------------
# Plot and attributes
#--
plt.subplots_adjust(wspace=0, hspace=0)
f.subplots_adjust(right=0.8)
cbar_ax = f.add_axes([0.85, 0.15, 0.05, 0.7])
cbar_ax.set_xlabel(r"$\mathbf{"+cbarlab+"}$",labelpad=10) #labelpad = label offset
cbar_ax.xaxis.set_label_position('top')

f.colorbar(im, cax=cbar_ax)
f.suptitle(r"$\mathbf{"+ptitle+"}$", size = 26)

# Set common labels
f.text(0.5, 0.04, r"$\mathbf{"+xlab+"}$", ha='center', va='center')
f.text(0.06, 0.5, r"$\mathbf{"+ylab+"}$", ha='center', va='center', rotation='vertical')

#---------------------
# Show it
#--
plt.show()
