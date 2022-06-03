using Images
include("./sff_rti.jl")
include("./sff-rti_utilities.jl")

### Begin test harness
base = "J:/Research/Generated Data/Blender/Statue du parc d'Austerlitz/SFFRTI/Restrained Z/"

innerFolderList = ["RTI_4_SFF_5", "RTI_4_SFF_20", "RTI_4_SFF_100", "RTI_20_SFF_5", "RTI_20_SFF_20"]
methodList = ["fvg", "max", "mean"]
kernelList = ["sml", "sobel"]

outputFolder = "J:/Image Out/Restrained Z/SSIM Tests/"

# Check file path endings to make sure they end in path separators
base = FixPathEnding(base)
outputFolder = FixPathEnding(outputFolder)

# Store path of ground truth depth map for comparison
gtPath = "J:/Research/Generated Data/Blender/Statue du parc d'Austerlitz/Ground Truth Restrained Z/Depth/Image0001.png"

ksize_list = [3, 5]

sff_rti(base, innerFolderList, methodList, kernelList; ksizeList=ksize_list, write_maps=false, write_csv=true, compute_psnr=true, compute_ssim=true, outputFolder=outputFolder, gtPath=gtPath)