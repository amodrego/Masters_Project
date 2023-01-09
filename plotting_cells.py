
# Importing libraries
from datetime import time
#from os import times, listdir, startfile, remove
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.art3d as art3d
import pandas as pd
import seaborn as sns
from scipy import io as sio
from pathlib import Path
from scipy.spatial import distance as dist
from sklearn.cluster import DBSCAN
from matplotlib.patches import Ellipse
import sys
import time

start = time.time()

#%matplotlib inline
sns.set_context('notebook')
sns.set_palette(['#0071BC', '#d95319', '#efb320'])



##### METHODS USED HERE #####
def build_path_name(timestep, folder_name, data_type):
    """Returns a Path object with the adequate PhysiCell output name format.
    
    Uses the standard PhysiCell output filename structure,
    (output{data_type}_{time_id}.mat)
    
    Parameters
    ----------
    timestep : int
        The time point at which the output was recorded
    folder_path: Path
        The path to the folder where the output (.mat, .xml) files are stored
    data_type: string
        The type of data to be retrieved (cells for cell-based data and
        microenvironment for the continuum variables)
    """

    # All possible file types written by PhysiCell
    data_type_name = {
        'cells': 'cells_physicell'
    }
    # Variable definition
    file_name = data_type_name[data_type]

    time_str = str(timestep).zfill(8)

    file_name = 'output{}_{}.mat'.format(time_str, file_name)



    # Build path name
    path_name = folder_name / file_name

    return path_name


def get_cell_data(timestep, folder_name, variables='all'):
    """Returns a dictionary with the cell output data for the given variables.
    
    Parameters
    ----------
    timestep : int
        The time point at which the output was recorded
    folder_path: Path
        The path to the folder where the output (.mat, .xml) files are stored
    variables : list
        The variables to be extracted from the output files. If variables
        are not defined, all the available outputs will be saved.
    """

    # All possible output variables written by PhysiCell
    data_labels = [
        'ID',
        'position_x', 'position_y', 'position_z',
        'total_volume',
        'cell_type',
        'cycle_model', 'current_phase', 'elapsed_time_in_phase',
        'nuclear_volume', 'cytoplasmic_volume',
        'fluid_fraction', 'calcified_fraction',
        'orientation_x', 'orientation_y', 'orientation_z',
        'polarity',
        'migration_speed',
        'motility_vector_x', 'motility_vector_y', 'motility_vector_z',
        'migration_bias',
        'motility_bias_direction_x', 'motility_bias_direction_y', 'motility_bias_direction_z',
        'persistence_time',
        'motility_reserved'
    ]

    # Variable definition
    cells = {}
    file_type = 'cells'

    if variables == 'all':
        variables = data_labels

    # Build path name
    path_name = build_path_name(timestep, folder_name, file_type)


    # Read output file
    cell_data = sio.loadmat(path_name)['cells']

    # Select and save the variables of interest
    variables_indexes = [data_labels.index(var) for var in variables]

    for index, var, in zip(variables_indexes, variables):
        cells[var] = cell_data[index, :]

    return cells


def cells_in_z_area_of_interest(cells, height_of_interest):
    """Selects cells in a defined area and returns them in a DataFrame.
    
    
    Parameters
    ---------
    cells : DataFrame
        The ax object to be stylized
    height_of_interest : int
        The threshold that defines which cells to select. Assumed to be
        equal for positive and negative z coordinates
    """
    
    # Define conditions 
    cells_above_neg_height_of_interest = cells['position_z'] > -height_of_interest
    cells_below_height_of_interest = cells['position_z'] < height_of_interest
    area_of_interest = cells_above_neg_height_of_interest & cells_below_height_of_interest
    
    # Select cells
    cells_of_interest = cells[area_of_interest]

    return cells_of_interest


def classify_cells_into_clusters(cells, DBSCAN_radius=18, DBSCAN_min_cells=3):
    """Classifies cells into clusters using the DBSCAN algorithm.
    
    Uses the cells' spatial information (x and y coordinates) through sklearn's
    DBSCAN algorithm. More information can be found at:
    https://scikit-learn.org/stable/modules/generated/sklearn.cluster.DBSCAN.html
    Returns the input DataFrame, with an additional column ('cluster') corrsponding to
    the ID of the cluster the cell belongs to. Outlier cells are given an ID of -1.
    
    Parameters
    ----------
    cells : DataFrame
        The DataFrame with cell data to be classified
    DBSCAN_radius : int
        The radius to be considered by the DBSCAN algorithm
        (default is 18 [microns])
    DBSCAN_min_cells : int
        The minimum number of cells to be used by the DBSCAN algorithm
        (default is 3 cells)
    """
    
    # Classify cells based on spatial information
    cells_positions = cells[['position_x', 'position_y']]
    dbscan_clusters = DBSCAN(eps=DBSCAN_radius, min_samples=DBSCAN_min_cells).fit(cells_positions)
    
    # Update the cells DataFrame with the corresponding cluster labels
    cells['cluster'] = dbscan_clusters.labels_

    return cells

def classify_cells_into_death(cells):
    cell_cycles = cells['cycle_model']


    return cells


def draw_cells_as_spheres(ax, cells, y_variable='position_y', z_variable='position_z', clusterized=False, cell_radius=8):
    """Draws cells as circles of a defined radius.
    
    Parameters
    ----------
    ax : ax
        The ax object to draw in.
    cells : DataFrame
        The DataFrame with the cells to be drawn.
    y_variable : string
        The coordinate to be represented in the y axis (y or z, default
        is the y coordinate)
    clusterized : Boolean
        Defines if cells should be colorized based on the cluster they
        belong to (default is False, all cells get the same color)
    cell_radius : int
        The radius of the cells (default is 8 [microns])
    """
    
    # When cells are classified into clusters, cells from the same cluster are plotted with the same color.
    # Independent cells are plotted as black.
    if clusterized:
        num_colors = len(cells['cluster'].unique())

        my_palette = sns.color_palette('Spectral', num_colors)
        my_palette.insert(len(my_palette), (0,0,0))

        for x_pos, y_pos, z_pos, color, model in zip(cells['position_x'], cells[y_variable], cells[z_variable], cells['cluster'], cells['cycle_model']):
            #cell_contour = plt.Circle((x_pos, y_pos),
            #                          cell_radius,
            #                          facecolor=my_palette[color], edgecolor='black')

            #ax.add_patch(cell_contour)
            #print(model)
            if model == 101.0 :
                ax.scatter(x_pos, y_pos, z_pos, s=cell_radius, color='white', edgecolors='black')
            else:
                ax.scatter(x_pos, y_pos, z_pos, s=cell_radius, color=my_palette[color], edgecolors='black')

    # When cells are not classified into clusters, all cells are plotted with the same color.
    else:
        for x_pos, y_pos, z_pos, model in zip(cells['position_x'], cells[y_variable], cells[z_variable], cells['cycle_model']):
            #cell_contour = plt.Circle((x_pos, y_pos),
            #                          cell_radius,
            #                          facecolor='C0', edgecolor='black')


            #ax.add_patch(cell_contour)
            if model == 101.0 :
                ax.scatter(x_pos, y_pos, z_pos, s=cell_radius, color='white', edgecolors='black')
            else:
                ax.scatter(x_pos, y_pos, z_pos, s=cell_radius, color='black', edgecolors='black')


def set_cell_view_style_axes(ax):
    """Sets axes style to predefined settings.
    
    Includes a light gray dashed grid and all spines set as black, 
    with ticks pointing inwards.
    
    Parameters
    ---------
    ax : ax
        The ax object to be stylized.
    """

    # Grid
    ax.yaxis.grid(color='lightgray', linestyle='--', linewidth=0.5, zorder=-1)
    ax.xaxis.grid(color='lightgray', linestyle='--', linewidth=0.5, zorder=-1)
    ax.zaxis.grid(color='lightgray', linestyle='--', linewidth=0.5, zorder=-1)

    # Ticks
    ax.tick_params(axis="y", direction="in", right=True)
    ax.tick_params(axis="x", direction="in", top=True)
    ax.tick_params(axis="z", direction="in", top=True)

    # Spines
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(1)
        spine.set_edgecolor('black')

#############################

##### INITIAL VARIABLES DECLARATION HERE #####

variables = ['ID', 'position_x', 'position_y', 'position_z', 'cycle_model']
files = []

#if (timestep == ''):
timestep = int(sys.argv[1])

folder = str(sys.argv[2]) 

if(str(timestep) == ""):
	timestep = 100

if(folder == ""):
    folder = "output"

base_folder_name = Path('./' + str(folder) + '/')  # Cluster Path
#base_folder_name = Path('D:/TFM/output/output_o2dens_100/')    # Local Path
cells_folder = base_folder_name / 'cells_plot/'
#############################


# Create fig object
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
fig.set_figwidth(10)
fig.set_figheight(5)

# Bucle para recorrer todos los timesteps y obtener desarrollo visual 

#for file in listdir(base_folder_name):
#    if file.endswith('.xml'):
#         files.append(file)
#counter = counter - 2


#while timestep < counter:
#for timestep, xml_name in enumerate(files):
    # Extract cell data from output folders
    #while timestep < counter :
    #    if timestep%10 == 0:
    #        print('Timestep number {}'.format(timestep))
    #    get_save_fig(timestep, SUBSTANCE, BASE_FILE, ax)
    #    timestep +=1
    

cells = get_cell_data(timestep, base_folder_name, variables)
    
# Se almacenan los datos seleccionados de las celulas de todo el volumen
cells_df = pd.DataFrame(cells, columns=variables)

# Se eliminan las celulas que no esten en el plano pasado al metodo
#cells_df = cells_in_z_area_of_interest(cells_df, z_height)

# Clusterization
cells_df = classify_cells_into_clusters(cells_df, DBSCAN_min_cells=3)
    
# Define ax to draw in
#ax = axes[timestep]
    
# Plot cells with an accurate representation of the assumed geometry
draw_cells_as_spheres(ax, cells_df, y_variable='position_y', z_variable='position_z', clusterized=True)
    
# Figure aesthetics
ax.set_ylim(-300,300)
ax.set_xlim(-300,300)
ax.set_zlim(-300,300)

# Use the standard axes style
set_cell_view_style_axes(ax)
    
#ax.set_title("Collagen concentration: {} [mg/mL]".format(concentrations[index]), y=1.05, fontsize=15)
    
ax.set_xlabel("X Position [$\mu$m]", labelpad=15, fontsize=15)
ax.set_ylabel("Y Position [$\mu$m]", labelpad=10, fontsize=15)
ax.set_zlabel("Z Position [$\mu$m]", labelpad=10, fontsize=15)
ax.set_title("Step number " + str(timestep), loc='center')

end = time.time()
plt.show()

print('Elapsed time (s): '+ str(round(end-start, 1)))

