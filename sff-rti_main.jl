using Images, Glob, CSV, DataFrames, CoordinateTransformations, LinearAlgebra, ProgressMeter, Printf

include("./sff.jl")
include("./IlluminationInvariance.jl")
include("./ImageUtilities.jl")

# base = "D:/Research/Generated Data/Blender/Statue du parc d'Austerlitz/SFFRTI/Black Background - RTI_20_SFF_20/"
base = "D:/Research/Generated Data/Blender/Statue du parc d'Austerlitz/SFFRTI/Black Background - RTI_4_SFF_5/"
# base = "D:/Research/Generated Data/Blender/Statue du parc d'Austerlitz/SFFRTI/Black Background - RTI_20_SFF_20_FullScale/"
folderPath = base * "/Renders/"

# fileList = glob("*.png", folderPath)

csvPath = base * "/Image.csv"

# Get all the XYZ camera positions in CSV.row structs
rowList = CSV.File(csvPath; select=["image", "x_lamp", "y_lamp", "z_lamp", "z_cam"])

# Get all Z positions and remove duplicates
zPosList = unique([row.z_cam for row in CSV.File(csvPath; select=["z_cam"])])

fvgList = []

p = Progress(length(zPosList), "Computing full vector gradient images...")
# Compute FVG for each Z position
for idx in eachindex(zPosList)
    angleList = []
    fileList = []

    for row in rowList

        # Conditional to make sure we're only dealing with one Z position at a time
        if row.z_cam != zPosList[idx]
            continue
        end

        push!(fileList, folderPath*row.image*".png")

        # Compute spherical coordinates from Cartesian
        sph = SphericalFromCartesian()([row.x_lamp,row.y_lamp,row.z_lamp])

        # Place all coordinates (including spherical coordinate system angle) into LightAngle object list
        push!(angleList, LightAngle(row.x_lamp, row.y_lamp, row.z_lamp, sph.θ, sph.ϕ))

    end

    imgFVG = ComputeFullVectorGradient(fileList, angleList)

    push!(fvgList, imgFVG)
    ProgressMeter.next!(p; showvalues= [(:"Current Distance", zPosList[idx])])
end

@printf("\nComputing depth map...")
Z, R = sff(fvgList, reverse(zPosList)) #TEMP REVERSE
# Z, R = sff(fvgList, zPosList)

# @printf("\nBlurring depth map before computation of normals...")
# Z_blurred = imfilter(Z, Kernel.gaussian((2)))

# Compute normals from depth map
@printf("\nComputing surface normals from depth map...")
# normals = Depth2Normal(Z)
# normals = Depth2Normal(complement.(Z)*100)
# normals = Depth2Normal(Z*100)
normals = Depth2Normal(complement.(Z))

# Z = FilterNaNs(Z, maximum(zPosList))

# Carve depth map with R
ZC = Z
ZC[R.<20] .= minimum(zPosList)

RC[R.<20] .= 0
RC = Gray.(RC)

### STATISTICS & VISUALIZATION

# Print value ranges
@printf("\nRange of normalsX (Converted): [%f, %f]\n", minimum(normals[:,:,1]), maximum(normals[:,:,1]))
@printf("Range of normalsY (Converted): [%f, %f]\n", minimum(normals[:,:,2]), maximum(normals[:,:,2]))
@printf("Range of normalsZ (Converted): [%f, %f]\n", minimum(normals[:,:,3]), maximum(normals[:,:,3]))

# normals = normals * 255

normalsSFFRTIX = normals[:,:,1]
normalsSFFRTIY = normals[:,:,2]
normalsSFFRTIZ = normals[:,:,3]

normalsX = shiftNormalsRange(normals[:,:,1])
normalsY = shiftNormalsRange(normals[:,:,2])
normalsZ = shiftNormalsRange(normals[:,:,3])

normalsColor = colorview(RGB, normalsX, normalsY, normalsZ)
# normalsColorShifted = colorview(RGB, shiftNormalsRange(normalsX), shiftNormalsRange(normalsY), shiftNormalsRange(normalsZ))

@printf("\nRange of normalsX (Converted, normalized): [%f, %f]\n", minimum(normalsX), maximum(normalsX))
@printf("Range of normalsY (Converted, normalized): [%f, %f]\n", minimum(normalsY), maximum(normalsY))
@printf("Range of normalsZ (Converted, normalized): [%f, %f]\n", minimum(normalsZ), maximum(normalsZ))

@printf("\nVariable `normalsColor` ready to be written out in current format.")

# save("statueFVG.png", Gray.(imgFVG/maximum(imgFVG)))
