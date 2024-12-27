import qupath.lib.gui.tools.MeasurementExporter
import qupath.lib.objects.PathDetectionObject
import qupath.lib.projects.ProjectImageEntry

def project = getProject()
def imagesToExport = project.getImageList()

def separator = "\t"
def columnsToInclude = new String[]{"Classification", "Centroid X µm", "Centroid Y µm", "Nucleus: Area"}
def exportType = PathDetectionObject.class
def outputDirectory = "D:/CD20//local_coord/" // output_path


new File(outputDirectory).mkdirs()


def exporter = new MeasurementExporter()
                  .separator(separator)            
                  .includeOnlyColumns(columnsToInclude) 
                  .exportType(exportType)          
                  .filter { obj -> obj.getPathClass() == getPathClass("Positive") } 


imagesToExport.each { imageEntry ->
    def imageName = imageEntry.getImageName() + '.txt'
    def outputPath = outputDirectory + imageName
    def outputFile = new File(outputPath)

    
    exporter.imageList([imageEntry])  
            .exportMeasurements(outputFile)  

    print "Exported measurements to ${outputPath}"
}

print "Done!"
