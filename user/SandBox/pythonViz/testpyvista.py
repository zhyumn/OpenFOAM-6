import pyvista as pv

# read the data
grid = pv.read('/home/zhy/new_drive/OpenFOAM-6/user/tutorials/sprayA1DIG/VTK/sprayA1DIG_0.vtk')

# plot the data with an automatically created Plotter
grid.plot(show_scalar_bar=False, show_axes=False)