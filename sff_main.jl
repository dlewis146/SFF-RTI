using LinearAlgebra, Glob, Images, ImageView, Gtk, CSV, DataFrames, Printf

include("./sff.jl")
include("./ImageUtilities.jl")

# base = "D:/Research/Generated Data/Blender/Statue du parc d'Austerlitz/SFF/5"
# base = "D:/Research/Generated Data/Blender/Statue du parc d'Austerlitz/SFF/Black Background - 20_FullScale"
# base = "D:/Research/Generated Data/Blender/Statue du parc d'Austerlitz/SFF/20 with off center light"
base = "D:/Research/Generated Data/Blender/Statue du parc d'Austerlitz/SFF/Black Background - 50_FullScale/"

# Get all interesting files
folderPath = base*"/Renders/"
fileList = glob("*.png", folderPath)
csvPath = base*"/Image.csv"

# Compute focus maps and store in array
imageList = []
for file in fileList
    img = Gray.(load(file))

    # focusMapX, focusMapY = FilterImage(img, "sobel")

    # focusMap = focusMapX.^2 + focusMapY.^2
    # focusMap = sqrt.(focusMapX.^2 + focusMapY.^2)

    focusMap = SumModifiedLaplacian(img)

    # Average filter focus map
    focusMap = FilterImageAverage(focusMap, (3,3))

    push!(imageList, focusMap)
end

# Read z_cam out of CSV and then extract from array of CSV.Row objects into a vector of Float64's
focusList = CSV.File(csvPath; select=["z_cam"])

focusList = reverse([row.z_cam for row in focusList])

# Z = sff(imageList, focusList)
Z, R = sff(imageList, focusList, 2, true)

# Carve depth map with R
# R<20 is unreliable
ZC = Z
ZC[R.<20] .= minimum(focusList)

RC = ones(size(R))
RC[R.<20] .= 0
RC = Gray.(RC)

# Compute normals from depth map
normals = Depth2Normal(complement.(Z))

### STATISTICS

normalsX = shiftNormalsRange(normals[:,:,1])
normalsY = shiftNormalsRange(normals[:,:,2])
normalsZ = shiftNormalsRange(normals[:,:,3])

normalsColor = colorview(RGB, normalsX, normalsY, normalsZ)

save("D:\\Meeting Pictures 1-14\\R_50_SFF_SML.png", imageNormalize(R))
save("D:\\Meeting Pictures 1-14\\RC_50_SFF_SML.png", imageNormalize(RC))
save("D:\\Meeting Pictures 1-14\\Z_50_SFF_SML.png", imageNormalize(Z))
save("D:\\Meeting Pictures 1-14\\normalsColor_50_SFF_SML.png", normalsColor)
