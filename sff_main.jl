using Glob, Images, CSV, DataFrames

include("./sff.jl")
include("../Utilities/ImageUtilities.jl")

base = "C:/Users/dlewi/Downloads/BlenderData/Statue/SFF/25/"

# Get all interesting files
folderPath = base*"PNG/"
fileList = glob("*.png", folderPath)
csvPath = base*"/Image.csv"

# Compute focus maps and store in array
imageList = []
for file in fileList
    img = Gray.(load(file))

    focusMapX, focusMapY = FilterImage(img, "sobel")

    focusMap = abs.( focusMapX + focusMapY )

    push!(imageList, focusMap)
end

# Read z_cam out of CSV and then extract from array of CSV.Row objects into a vector of Float64's
focusList = CSV.File(csvPath; select=["z_cam"])
focusList = [row.z_cam for row in focusList]

Z = sff(imageList, focusList)

Gray.(imageNormalize(Z))