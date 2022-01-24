using LinearAlgebra, Glob, Images, ImageView, CSV, DataFrames, Printf, PlotlyJS

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

# ROI = Nothing
# ROI = [(850, 150), (1175, 150), (850, 650), (1175, 650)]
# ROI = [(926, 199), (1096, 199), (926, 535), (1096, 535)]

# Compute focus maps and store in array
imageList = []
for file in fileList
    # img = Gray2UInt16(load(file))
    # img = Gray2Float64(load(file))
    img = Gray.(load(file))

    # if ROI != Nothing
        # img = img[ROI[1][2]:ROI[3][2], ROI[1][1]:ROI[2][1]]
    # end

    # focusMapX, focusMapY = FilterImage(img, "sobel")

    # focusMap = focusMapX.^2 + focusMapY.^2
    # focusMap = sqrt.(focusMapX.^2 + focusMapY.^2)

    focusMap = SumModifiedLaplacian(img)

    # Average filter focus map
    focusMap = FilterImageAverage(focusMap)

    push!(imageList, focusMap)
end

# Read z_cam out of CSV and then extract from array of CSV.Row objects into a vector of Float64's
focusList = CSV.File(csvPath; select=["z_cam"])

focusList = reverse([row.z_cam for row in focusList])

# Z = sff(imageList, focusList)
Z, R = sff(imageList, focusList, 2, true)

Z = FilterNaNs(Z, maximum(focusList))

# Carve depth map with R
# R<20 is unreliable
ZC = Z
ZC[R.<20] .= minimum(focusList)

RC = ones(size(R))
RC[R.<20] .= 0
# RC[R.<20] .= mapwindow(median, RC, (3,3))     

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
