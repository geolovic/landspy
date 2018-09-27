from qgis.PyQt.QtCore import QCoreApplication, QVariant
from qgis.core import (QgsField, QgsFeature, QgsFeatureSink, QgsFeatureRequest, QgsProcessing, QgsProcessingAlgorithm, QgsProcessingParameterFeatureSource, QgsProcessingParameterFeatureSink, QgsProcessingParameterRasterLayer, QgsProcessingParameterString, QgsProcessingParameterRasterDestination,QgsProcessingParameterBoolean)

from topopy import DEM, Flow

class ProcessFlow(QgsProcessingAlgorithm):
    INPUT_DEM = 'INPUT_DEM' # Input DEM
    OUTPUT_FAC = 'OUTPUT_FAC'
    OUTPUT_FD = 'OUTPUT_FD'
    VERBOSE = 'VERBOSE'
    AUXTOPO = 'AUXTOPO'
    FILLED = 'FILLED'
 
    def __init__(self):
        super().__init__()
 
    def name(self):
        return "processflow"
     
    def tr(self, text):
        return QCoreApplication.translate("processflow", text)
         
    def displayName(self):
        return self.tr("Process Flow")
 
    def group(self):
        return self.tr("topopy")
 
    def groupId(self):
        return "topopy"
 
    def shortHelpString(self):
        texto = """
                    This script process the flow from a Digital Elevation Model:
                    Input DEM : Elevation raster (DEM) to process
                    Use auxiliar topography : Uses auxliar topography to route flow in flats areas (slow)
                    Filled DEM : True if the input DEM is already pit-filled
                    Show messages : Print progress messages (usefull for large grids)
                    Output Flow Direction Raster : Raster with flow directions (not an usual raster)
                    Output Flow Accumulation Raster : Raster with flow accumulation values
                    """
        return texto
 
    def helpUrl(self):
        return "https://qgis.org"
         
    def createInstance(self):
        return type(self)()
   
    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterRasterLayer(self.INPUT_DEM,  self.tr("Input DEM")))
        self.addParameter(QgsProcessingParameterRasterDestination(self.OUTPUT_FD, "Output Flow Direction raster", None, False))
        self.addParameter(QgsProcessingParameterRasterDestination(self.OUTPUT_FAC, "Output Flow Accumulation raster", None, False))
        self.addParameter(QgsProcessingParameterBoolean(self.AUXTOPO, "Use auxiliar topography", False))
        self.addParameter(QgsProcessingParameterBoolean(self.FILLED, "Filled DEM", False))
        self.addParameter(QgsProcessingParameterBoolean(self.VERBOSE, "Show messages", False))
 
    def processAlgorithm(self, parameters, context, feedback):
        dem_raster = self.parameterAsRasterLayer(parameters, self.INPUT_DEM, context)
        fac_raster = self.parameterAsOutputLayer(parameters, self.OUTPUT_FAC, context)
        fd_raster = self.parameterAsOutputLayer(parameters, self.OUTPUT_FD, context)
        verbose = self.parameterAsBool(parameters, self.VERBOSE, context)
        filled = self.parameterAsBool(parameters, self.FILLED, context)
        auxtopo = self.parameterAsBool(parameters, self.AUXTOPO, context)
        
        dem = DEM(dem_raster.source())
        fd = Flow(dem, auxtopo=auxtopo, filled=filled, verbose=verbose, verb_func=feedback.setProgressText)
        fac = fd.get_flow_accumulation()
        fd.save_gtiff(fd_raster)
        fac.save(fac_raster)
        results = {}
        results[self.OUTPUT_FD] = fd_raster
        results[self.OUTPUT_FAC] = fac_raster
        return results