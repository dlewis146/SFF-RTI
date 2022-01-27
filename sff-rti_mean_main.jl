using Images, Glob, CSV, DataFrames, CoordinateTransformations, LinearAlgebra, ProgressMeter, Printf

include("./sff.jl")
include("./IlluminationInvariance.jl")
include("./ImageUtilities.jl")

function ComputeMeanImage(imageList)

    meanImageBuild = zero(imageList[1])

    for image in imageList
        meanImageBuild = meanImageBuild + image
    end

    meanImageBuild = meanImageBuild / length(imageList)

    return meanImageBuild

end

function ComputeSTDImage(imageList)

    stdImageStack = zeros(size(imageList[1], 1), size(imageList[1], 2), length(imageList))

    for idx in range(1, stop=length(imageList))
        stdImageStack[:,:,idx] = imageList[idx]
    end

    stdImageOut = zero(imageList[1])

    for y in range(1, stop=size(stdImageStack, 1))
        for x in range(1, stop=size(stdImageStack, 2))

            stdImageOut[y,x] = std(stdImageStack[y,x,:])
        end
    end

    return stdImageOut

end


# base = "D:/Research/Generated Data/Blender/Statue du parc d'Austerlitz/SFFRTI/Black Background - RTI_20_SFF_20/"
# base = "D:/Research/Generated Data/Blender/Statue du parc d'Austerlitz/SFFRTI/Black Background - RTI_4_SFF_5/"
base = "D:/Research/Generated Data/Blender/Statue du parc d'Austerlitz/SFFRTI/Black Background - RTI_20_SFF_20_FullScale/"
folderPath = base * "/Renders/"

# fileList = glob("*.png", folderPath)

csvPath = base * "/Image.csv"

# Get all the XYZ camera positions in CSV.row structs
rowList = CSV.File(csvPath; select=["image", "x_lamp", "y_lamp", "z_lamp", "z_cam"])

# Get all Z positions and remove duplicates
zPosList = reverse(unique([row.z_cam for row in CSV.File(csvPath; select=["z_cam"])]))

fvgList = []

## PREP FOR WRITING OUT FVG IMAGES
# If directory doesn't exist, create it
# outputFolder = "./STD_SML_OUT"
# outputFolder = "./MEAN_SML_OUT"

# ispath(outputFolder) || mkpath(outputFolder)

p = Progress(length(zPosList), "Computing output images...")    
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

    imageList = ReadImageList(fileList)

    gradientList = []
    for img in imageList

        imgMeanX, imgMeanY = FilterImage(img, "sobel")
        focusMap = sqrt.(imgMeanX.^2 + imgMeanY.^2)
        # focusMap = SumModifiedLaplacian(img)

        focusMap = FilterImageAverage(focusMap)
        push!(gradientList, focusMap)

    end

    # imgMean = ComputeMeanImage(gradientList)
    imgMean = ComputeSTDImage(gradientList)

    # save(outputFolder*string("/Z_", round(zPosList[idx], digits=3), ".png"), imageNormalize(imgMean))

    push!(fvgList, imgMean)

    ProgressMeter.next!(p; showvalues= [(:"Current Distance", zPosList[idx])])
end

@printf("\nComputing depth map...")
Z, R = sff(fvgList, zPosList, 2, true)

# Compute normals from depth map
@printf("\nComputing surface normals from depth map...")
normals = Depth2Normal(Z)
# normals = Depth2Normal(complement.(Z))

# Z = FilterNaNs(Z, maximum(zPosList))

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

@printf("\nVariable `normalsColor` ready to be written out in current format.")

save("D:\\Image Out\\R_20_20_STD_TENG.png", imageNormalize(R))
save("D:\\Image Out\\RC_20_20_STD_TENG.png", imageNormalize(RC))
save("D:\\Image Out\\Z_20_20_STD_TENG.png", imageNormalize(Z))
save("D:\\Image Out\\normalsColor_20_20_STD_TENG.png", normalsColor)
