

def vtk_show(renderer, size=(400,300)):
    """
    Takes vtkRenderer instance and returns an IPython Image with the rendering.
    from https://pyscience.wordpress.com/2014/09/03/ipython-notebook-vtk/
    """

    from IPython.display import Image
    import vtk

    renderWindow = vtk.vtkRenderWindow()
    renderWindow.SetOffScreenRendering(1)
    renderWindow.AddRenderer(renderer)
    renderWindow.SetSize(*size)
    renderWindow.Render()

    windowToImageFilter = vtk.vtkWindowToImageFilter()
    windowToImageFilter.SetInput(renderWindow)
    windowToImageFilter.Update()

    writer = vtk.vtkPNGWriter()
    writer.SetWriteToMemory(1)
    writer.SetInputConnection(windowToImageFilter.GetOutputPort())
    writer.Write()
    data = str(buffer(writer.GetResult()))

    return Image(data)

def create_vtk_volume(data, color = None, alpha = None):
    import vtk
    import numpy as np

    if not np.isrealobj(data):
        raise ValueError("Can only plot real valued data")

    data_matrix = data * 255./np.max(data)
    data_matrix = data_matrix.astype(np.uint8)
    data_matrix = np.fliplr(np.transpose(data_matrix,axes=(1,0,2)))

    # We begin by creating the data we want to render.
    # For this tutorial, we create a 3D-image containing three overlaping cubes.
    # This data can of course easily be replaced by data from a medical CT-scan or anything else three dimensional.
    # The only limit is that the data must be reduced to unsigned 8 bit or 16 bit integers.
    #data_matrix = zeros([75, 75, 75], dtype=uint8)
    #data_matrix[0:35, 0:35, 0:35] = 50
    #data_matrix[25:55, 25:55, 25:55] = 100
    #data_matrix[45:74, 45:74, 45:74] = 150

    # For VTK to be able to use the data, it must be stored as a VTK-image. This can be done by the vtkImageImport-class which
    # imports raw data and stores it.
    dataImporter = vtk.vtkImageImport()
    # The preaviusly created array is converted to a string of chars and imported.
    data_string = data_matrix.tostring(order='C')
    dataImporter.CopyImportVoidPointer(data_string, len(data_string))
    # The type of the newly imported data is set to unsigned char (uint8)
    dataImporter.SetDataScalarTypeToUnsignedChar()
    # Because the data that is imported only contains an intensity value (it isnt RGB-coded or someting similar), the importer
    # must be told this is the case.
    dataImporter.SetNumberOfScalarComponents(1)
    # The following two functions describe how the data is stored and the dimensions of the array it is stored in. For this
    # simple case, all axes are of length 75 and begins with the first element. For other data, this is probably not the case.
    # I have to admit however, that I honestly dont know the difference between SetDataExtent() and SetWholeExtent() although
    # VTK complains if not both are used.
    dataImporter.SetDataExtent(0, data_matrix.shape[2]-1, 0, data_matrix.shape[1]-1, 0, data_matrix.shape[0]-1)
    dataImporter.SetWholeExtent(0, data_matrix.shape[2]-1, 0, data_matrix.shape[1]-1, 0, data_matrix.shape[0]-1)

    # The following class is used to store transparencyv-values for later retrival. In our case, we want the value 0 to be
    # completly opaque whereas the three different cubes are given different transperancy-values to show how it works.
    alphaChannelFunc = vtk.vtkPiecewiseFunction()
    alphaChannelFunc.AddPoint(0, 0.0)

    if color is None:
        alphaChannelFunc.AddPoint(50, 0.05)
        alphaChannelFunc.AddPoint(100, 0.1)
        alphaChannelFunc.AddPoint(255, 0.4)

    else:
        alphaChannelFunc.AddPoint(255, alpha)

    # This class stores color data and can create color tables from a few color points. For this demo, we want the three cubes
    # to be of the colors red green and blue.
    colorFunc = vtk.vtkColorTransferFunction()
    if color is None:
        colorFunc.AddRGBPoint(0,   0.0, 0.0, 1.0)
        colorFunc.AddRGBPoint(50,  0.0, 1.0, 0.0)
        colorFunc.AddRGBPoint(100, 1.0, 0.0, 0.0)
        colorFunc.AddRGBPoint(255, 1.0, 1.0, 1.0)
    else:
        colorFunc.AddRGBPoint(0,*color)
        colorFunc.AddRGBPoint(255,*color)

    # The preavius two classes stored properties. Because we want to apply these properties to the volume we want to render,
    # we have to store them in a class that stores volume prpoperties.
    volumeProperty = vtk.vtkVolumeProperty()
    volumeProperty.SetColor(colorFunc)
    volumeProperty.SetScalarOpacity(alphaChannelFunc)

    # This class describes how the volume is rendered (through ray tracing).
    compositeFunction = vtk.vtkVolumeRayCastCompositeFunction()
    # We can finally create our volume. We also have to specify the data for it, as well as how the data will be rendered.
    volumeMapper = vtk.vtkVolumeRayCastMapper()
    #volumeMapper = vtk.vtkSmartVolumeMapper()
    #volumeMapper = vtk.vtkGPUVolumeRayCastMapper()

    volumeMapper.SetVolumeRayCastFunction(compositeFunction)
    volumeMapper.SetInputConnection(dataImporter.GetOutputPort())

    # The class vtkVolume is used to pair the preaviusly declared volume as well as the properties to be used when rendering that volume.
    volume = vtk.vtkVolume()
    volume.SetMapper(volumeMapper)
    volume.SetProperty(volumeProperty)

    return volume


def interactive_volumetric_plot(data,background_color = (1,1,1)):
    import vtk

    volume = create_vtk_volume(data)

    # With almost everything else ready, its time to initialize the renderer and window, as well as creating a method for exiting the application
    renderer = vtk.vtkRenderer()
    renderWin = vtk.vtkRenderWindow()
    renderWin.AddRenderer(renderer)
    renderInteractor = vtk.vtkRenderWindowInteractor()
    renderInteractor.SetRenderWindow(renderWin)

    # We add the volume to the renderer ...
    renderer.AddVolume(volume)
    # ... set background color to white ...
    renderer.SetBackground(*background_color)
    # ... and set window size.
    renderWin.SetSize(400, 400)

    # A simple function to be called when the user decides to quit the application.
    def exitCheck(obj, event):
        if obj.GetEventPending() != 0:
            obj.SetAbortRender(1)

    # Tell the application to use the function as an exit check.
    renderWin.AddObserver("AbortCheckEvent", exitCheck)

    renderInteractor.Initialize()
    # Because nothing will be rendered without any input, we order the first render manually before control is handed over to the main-loop.
    renderWin.Render()
    renderInteractor.Start()


