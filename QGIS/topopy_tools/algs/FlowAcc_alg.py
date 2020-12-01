from qgis.PyQt.QtCore import QCoreApplication
from qgis.core import QgsProcessingAlgorithm, QgsProcessingParameterRasterLayer, QgsProcessingParameterRasterDestination, QgsProcessingParameterBoolean
from topopy import DEM, Flow, Grid
from qgis import processing

class FlowAccumulation(QgsProcessingAlgorithm):
    # Constants used to refer to parameters and outputs They will be
    # used when calling the algorithm from another algorithm, or when
    # calling from the QGIS console.

    INPUT_FD = 'INPUT_FD'
    #INPUT_WG = 'INPUT_WG'
    OUTPUT_FAC = 'OUTPUT_FAC'
 
    def __init__(self):
        super().__init__()

    def createInstance(self):
        return type(self)()
 
    def name(self):
        """
        Rerturns the algorithm name, used to identify the algorithm.
        Must be unique within each provider and should contain lowercase alphanumeric characters only.
        """
        return "flowacc"
     
    def displayName(self):
        """
        Returns the translated algorithm name, which should be used for any
        user-visible display of the algorithm name.
        """
        return self.tr("Flow Accumulation") 
    
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
                    This script creates a flow accumulation raster. 
                    Flow direction : Input flow direction raster
                    Flow accumulation: Output raster that shows the accumulated flow to each cell.
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
        self.addParameter(QgsProcessingParameterRasterLayer(self.INPUT_FD,  self.tr("Flow direction")))
        #self.addParameter(QgsProcessingParameterRasterLayer(self.INPUT_WG,  self.tr("Weight raster"), optional=True))
        self.addParameter(QgsProcessingParameterRasterDestination(self.OUTPUT_FAC, self.tr("Flow accumulation"), None, False))
 
    def processAlgorithm(self, parameters, context, feedback):
        """
        Here is where the processing itself takes place.
        """
        input_fd = self.parameterAsRasterLayer(parameters, self.INPUT_FD, context)
        #input_wg = self.parameterAsRasterLayer(parameters, self.INPUT_WG, context)
        output_fac = self.parameterAsOutputLayer(parameters, self.OUTPUT_FAC, context)

        #if input_wg is None:
        #    wg = None
        #else:
        #    wg = Grid(input_wg.source())
        
        fd = Flow(input_fd.source())
        fac = fd.get_flow_accumulation()
        fac.save(output_fac)
        
        results = {self.OUTPUT_FAC : output_fac, }
        return results