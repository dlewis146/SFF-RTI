using Images
include("./sff_rti.jl")
include("./sff-rti_utilities.jl")

base = "F:/Research/Generated Data/Blender/Statue du parc d'Austerlitz/SFFRTI/Restrained Z/"

innerFolderList = ["RTI_4_SFF_5", "RTI_4_SFF_20", "RTI_4_SFF_100", "RTI_20_SFF_5", "RTI_20_SFF_20"]
methodList = ["fvg", "mean" ]
kernelList = ["sobel", "sml"]

outputFolder = "F:/Image Out/Restrained Z/"
write_maps_flag = true
compute_snr_flag = true

# Check file path endings to make sure they end in path separators
function FixPathEnding(f)  if last(f) !== '/' && last(f) !== '\\'; return f*"/" else return f end end
base = FixPathEnding(base)
outputFolder = FixPathEnding(outputFolder)

# Empty dictionaries "to gather RMSE and inverse RMSE 
rmseList = Dict()
irmseList = Dict()
snrDict = Dict()
psnrDict = Dict()

# Read in ground truth depth map for comparison
GT = Gray2Float64(Gray.(load("F:/Research/Generated Data/Blender/Statue du parc d'Austerlitz/Ground Truth Restrained Z/Depth/Image0001.png")))
# GT = Gray2Float64(Gray.(load("F:/Research/Generated Data/Blender/Statue du parc d'Austerlitz/RTI/Black Background - FullScale/Depth/Image0001.png")))

# outputStructList = Dict()
outputStructList = []

for f in innerFolderList
    for method in methodList
        for kernel in kernelList

            # Run SFF-RTI method
            Z, R, snrMean, psnrMean = sff_rti(base*f, method, kernel, outputFolder, write_maps_flag, compute_snr_flag)

            if outputFolder !== nothing
                push!(outputStructList, FileSet(Z,R,f,method,kernel))
            end

            # # Normalize computed depth map so that it's placed from 0-1
            Z_normalized = imageDisp01(Z)

            # Compute rmse and inverse RMSE then store all statistical measures in appropriate dictionaries
            rmseList[f,method,kernel]  = rmse(GT, Z_normalized)
            irmseList[f,method,kernel] = 1-rmseList[f,method,kernel]

            # snrDict[f,method,kernel] = 10*log10(mean(Z_normalized)/rmse(GT, Z_normalized))
            # psnrDict[f,method,kernel] = 10*log10(1/mse(GT, Z_normalized))
            snrDict[f,method,kernel] = snrMean
            psnrDict[f,method,kernel] = psnrMean

        end
    end
end

# ZMax, RMax = FindFileSetMax(outputStructList)
println("NOTE: Using SFF ZMax and RMax normalization coefficients")
ZMax = 8.5
RMax = 295.7640411457105

if outputFolder !== nothing
    WriteMaps(outputStructList, outputFolder, nothing, nothing)
end

WriteCSV(outputFolder*"/Ground truth comparison results.csv", innerFolderList, methodList, kernelList, rmseList, irmseList, ZMax, RMax, snrDict, psnrDict)
