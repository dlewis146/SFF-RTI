using Images, Glob, CSV, DataFrames, CoordinateTransformations, LinearAlgebra, ProgressMeter, Printf

include("./sff.jl")
include("./IlluminationInvariance.jl")
include("./ImageUtilities.jl")

base = "D:/Research/Generated Data/Blender/Statue du parc d'Austerlitz/SFFRTI/Black Background - RTI_20_SFF_20_FullScale/"
folderPath = base * "/Renders/"

csvPath = base * "/Image.csv"

outputFolder = "D:/Image Out/"

# Get all the XYZ camera positions in CSV.row structs
rowList = CSV.File(csvPath; select=["image", "x_lamp", "y_lamp", "z_lamp", "z_cam"])

# Get all Z positions and remove duplicates
zPosList = reverse(unique([row.z_cam for row in CSV.File(csvPath; select=["z_cam"])]))

fvgList = []

## PREP FOR WRITING OUT FVG IMAGES
# If directory doesn't exist, create it
# outputFolder = "./FVG_OUT"
# ispath(outputFolder) || mkpath(outputFolder)

prog = Progress(length(zPosList), "Computing full vector gradient images...")

# Compute FVG for each Z position
for idx in eachindex(zPosList)
    angleList = []
    fileList = []

    for row in rowList

        # Conditional to make sure we're only dealing with one Z position at a time
        if row.z_cam != zPosList[idx]
            continue
        end

        # Construct and store filename based off of CSV entries
        push!(fileList, folderPath*row.image*".png")

        # Compute spherical coordinates from Cartesian
        sph = SphericalFromCartesian()([row.x_lamp,row.y_lamp,row.z_lamp])

        # Place all coordinates (including spherical coordinate system angle) into LightAngle object list
        push!(angleList, LightAngle(row.x_lamp, row.y_lamp, row.z_lamp, sph.θ, sph.ϕ))

    end

    # Compute FVG
    imgFVG = ComputeFullVectorGradient(fileList, angleList, "sobel")

    # Average FVG image
    focusMap = FilterImageAverage(imgFVG)

    # Store FVG in list and iterate progress bar
    push!(fvgList, imgFVG)
    ProgressMeter.next!(prog; showvalues= [(:"Current Distance", zPosList[idx])])
end

@printf("\nComputing depth map...")
Z, R = sff(fvgList, zPosList, 2, true)

# Compute normals from depth map
@printf("\nComputing surface normals from depth map...")
normals = Depth2Normal(complement.(Z))

# Carve depth map with R
ZC = copy(Z)
ZC[R.<20] .= minimum(zPosList)

RC = copy(R)
RC[R.<20] .= 0
RC = Gray.(RC)

### VISUALIZATION & IO
normalsX = shiftNormalsRange(normals[:,:,1])
normalsY = shiftNormalsRange(normals[:,:,2])
normalsZ = shiftNormalsRange(normals[:,:,3])

normalsColor = colorview(RGB, normalsX, normalsY, normalsZ)

save(outputFolder*"/Z_20_20_FVG.png", imageNormalize(Z))
save(outputFolder*"/R_20_20_FVG.png", imageNormalize(R))
save(outputFolder*"/RC_20_20_FVG.png", imageNormalize(RC))
save(outputFolder*"/normalsColor_20_20_FVG.png", normalsColor)
