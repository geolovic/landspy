from qgis.PyQt.QtCore import QCoreApplication
from qgis.core import QgsProcessingAlgorithm, QgsProcessingParameterRasterLayer, QgsProcessingParameterRasterDestination
from qgis.core import QgsProcessingParameterNumber, QgsProcessingParameterString
from topopy import Flow, Network
from qgis import processing

class Streams(QgsProcessingAlgorithm):
    # Constants used to refer to parameters and outputs They will be
    # used when calling the algorithm from another algorithm, or when
    # calling from the QGIS console.

    INPUT_FD = 'INPUT_FD'
    THRESHOLD = 'THRESHOLD'
    UNITS = 'UNITS'
    TYPE = 'TYPE'
    NETWORK = 'NETWORK'

 
    def __init__(self):
        super().__init__()

    def createInstance(self):
        return type(self)()
 
    def name(self):
        """
        Rerturns the algorithm name, used to identify the algorithm.
        Must be unique within each provider and should contain lowercase alphanumeric characters only.
        """
        return "streams"
     
    def displayName(self):
        """
        Returns the translated algorithm name, which should be used for any
        user-visible display of the algorithm name.
        """
        return self.tr("Get Streams") 
    
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
                    This script creates a Network raster that represents a drainage network.
                    Flow Direction: Input flow direction raster
                    Threshold: Threshold to initiate a channel (can be entered in cells o map unit)
                    Units: Threshold units, can be "CELLS" or "MAP"
                    Raster type: Type of raster, can be "str" (normal stream raster), "links" (segments links will be coded with a id), "strahler" (segments will be coded with their strahler order)
                    Streams Raster: Output network raster
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
        self.addParameter(QgsProcessingParameterRasterLayer(self.INPUT_FD, self.tr("Flow Direction")))
        self.addParameter(QgsProcessingParameterNumber(self.THRESHOLD, self.tr("Threshold")))
        self.addParameter(QgsProcessingParameterString(self.UNITS, self.tr("Units"), defaultValue="CELLS"))
        self.addParameter(QgsProcessingParameterString(self.TYPE, self.tr("Type"), defaultValue="str"))
        self.addParameter(QgsProcessingParameterRasterDestination(self.NETWORK, self.tr("Streams raster"), None, False))

 
    def processAlgorithm(self, parameters, context, feedback):
        """
        Here is where the processing itself takes place.
        """
        input_fd = self.parameterAsRasterLayer(parameters, self.INPUT_FD, context)
        threshold = self.parameterAsInt(parameters, self.THRESHOLD, context)
        units = self.parameterAsString(parameters, self.UNITS, context)
        ras_type = self.parameterAsString(parameters, self.TYPE, context)
        output_str = self.parameterAsOutputLayer(parameters, self.NETWORK, context)
        
        fd = Flow(input_fd.source())
        
        if units == "MAP":
            cellsize = -fd.get_cellsize()[1] * fd.get_cellsize()[0]
            threshold = int(threshold / cellsize)
        
        net = Network(fd, threshold, gradients=False)
        
        if ras_type == "links":
            streams = net.get_stream_segments()
        elif ras_type == "strahler":
            streams = net.get_stream_orders()
        else:
            streams = net.get_streams()
            
        streams.save(output_str)
        
        results = {self.NETWORK : output_str}
        return results