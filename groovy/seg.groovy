extension = '.png'  
def imageData = getCurrentImageData()
def outputDirectory = "D:/pdac" // output_directory
def name = GeneralTools.getNameWithoutExtension(imageData.getServer().getMetadata().getName())
def pathOutput = buildFilePath(outputDirectory, './tiles', name)
mkdirs(pathOutput)

double downsample = 4 
int outputSize = 400 // output_resolution
int overLap = 0

new TileExporter(imageData)
    .downsample(downsample)
    .imageExtension(extension) // often .tif, .jpg, '.png' or '.ome.tif'
    .tileSize(outputSize) // Define size of each tile, in pixels
    //.annotatedTilesOnly(true) // If you uncomment, only patches with annotations are exported
    //.annotatedCentroidTilesOnly(true)  // If true, only tiles with comments in the center will be exported (not tiles with only a few comments)
    .overlap(overLap)  
    .writeTiles(pathOutput)  

// Don't interrupt when script runs.
def dirOutput = new File(pathOutput)
for (def file in dirOutput.listFiles()) {
    if (!file.isFile() || file.isHidden() || !file.getName().endsWith(extension))
        continue
    def newName = file.getName().replaceAll("=","-").replaceAll("\\[","").replaceAll("\\]","").replaceAll(",","_")
    if (file.getName() == newName)
        continue
    def fileUpdated = new File(file.getParent(), newName)
    println("Renaming ${file.getName()} ---> ${fileUpdated.getName()}")
    file.renameTo(fileUpdated)
}

println('Done!')