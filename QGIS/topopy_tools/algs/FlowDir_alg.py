from qgis.PyQt.QtCore import QCoreApplication
from qgis.core import QgsProcessingAlgorithm, QgsProcessingParameterRasterLayer, QgsProcessingParameterRasterDestination, QgsProcessingParameterBoolean
from topopy import DEM, Flow
from qgis import processing

class FlowDir(QgsProcessingAlgorithm):
    # Constants used to refer to parameters and outputs They will be
    # used when calling the algorithm from another algorithm, or when
    # calling from the QGIS console.

    INPUT_DEM = 'INPUT_DEM'
    OUTPUT_FD = 'OUTPUT_FD'
    VERBOSE = 'VERBOSE'
 
    def __init__(self):
        super().__init__()

    def createInstance(self):
        return type(self)()
 
    def name(self):
        """
        Rerturns the algorithm name, used to identify the algorithm.
        Must be unique within each provider and should contain lowercase alphanumeric characters only.
        """
        return "flowdir"
     
    def displayName(self):
        """
        Returns the translated algorithm name, which should be used for any
        user-visible display of the algorithm name.
        """
        return self.tr("Flow Direction") 
    
    def groupId(self):
        """
        Returns the unique ID of the group this algorithm belongs to.
        """
        return "drainage_net_processing"

    def group(self):
        """
        Returns the name of the group this algoritm belongs to.
        """
        return self.tr("Drainage Network Processing")

    def shortHelpString(self):
        """
        Returns a localised short helper string for the algorithm. 
        """
        texto = """
                    This script creates a flow direction raster. This is not an usual raster, but a raster with a special format used by topopy
                    Input DEM : Input pit-filled Digital Elevation Model (DEM)
                    Output Flow: Output flow direction raster
                    Show Messages: Show progress messages (useful for big rasters)
                    """
        return texto
 
    def tr(self, string):
        return QCoreApplication.translate('Processing', string)

    def helpUrl(self):
        return "https://qgis.org"
         

    def initAlgorithm(self, config=None):
        """
        Here we define the inputs and output of the algorithm, along
        with some other properties.
        """
        self.addParameter(QgsProcessingParameterRasterLayer(self.INPUT_DEM,  self.tr("Input DEM")))
        self.addParameter(QgsProcessingParameterRasterDestination(self.OUTPUT_FD, "Output filled DEM", None, False))
        self.addParameter(QgsProcessingParameterBoolean(self.VERBOSE, "Show Messages", False))

 
    def processAlgorithm(self, parameters, context, feedback):
        """
        Here is where the processing itself takes place.
        """
        input_dem = self.parameterAsRasterLayer(parameters, self.INPUT_DEM, context)
        output_fd = self.parameterAsOutputLayer(parameters, self.OUTPUT_FD, context)
        verbose = self.parameterAsBool(parameters, self.VERBOSE, context)
        
        dem = DEM(input_dem.source())
        fd = Flow(dem, filled =True, verbose=verbose, verb_func=feedback.setProgressText)
        fd.save(output_fd)
        
        results = {self.OUTPUT_FD : output_fd}
        return results