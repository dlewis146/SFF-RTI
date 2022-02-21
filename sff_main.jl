using LinearAlgebra, Glob, Images, ImageView, CSV, DataFrames, Printf
include("./sff.jl")
include("./sff-rti_utilities.jl")
include("./ImageUtilities.jl")

global ACCEPTED_KERNELS = ["sml", "sobel"]

function sff_main(baseFolder, kernel="sml", outputFolder=nothing, write_maps=false)
       
    # Make sure that given kernel is compatible with methods
    if !(kernel in ACCEPTED_KERNELS)
        error(@printf("\nGiven kernel { %s } isn't compatible with current version of `sff_main`\n", kernel))
    end

    # Make sure there aren't multiple CSV files in the base path
    if length(glob("*.csv", baseFolder)) > 1
        error("\nMore than 1 CSV file found in given folder\n")
    end

    csvPath = glob("*.csv", baseFolder)[1]

    # Get all interesting files
    folderPath = baseFolder*"/Renders/"
    # fileList = glob("*.png", folderPath)

    # NOTE: Using this instead of just globbing all files so that they're read in and stored in the proper order. This doesn't happen by default due to a lack of zero-padding of frame numbers in file names
    rowList = CSV.File(csvPath; select=["image", "z_cam"])
    fileList = []
    for row in rowList
        push!(fileList, folderPath*row.image*".png")
    end

    # Compute focus maps and store in array
    imageList = []
    for idx in eachindex(fileList)

        img = Gray.(load(fileList[idx]))
        focusMap = FilterImageCombined(img, kernel)

        # Average filter focus map
        focusMap = FilterImageAverage(focusMap, (3,3))

        if write_maps == true
            outputFolderConcatenated = string(outputFolder, "/FOCUS MAPS/SFF_", basename(baseFolder), " ", uppercase(kernel), "/")

            # Make sure all nested folders exist and create them if not
            ispath(outputFolder*"/FOCUS MAPS/") || mkpath(outputFolder*"/FOCUS MAPS/")
            ispath(outputFolderConcatenated) || mkpath(outputFolderConcatenated)

            save(outputFolderConcatenated*string(idx)*".png", imageNormalize(focusMap))
        end

        push!(imageList, focusMap)
    end

    # Read z_cam out of CSV and then extract from array of CSV.Row objects into a vector of Float64's
    focusList = CSV.File(csvPath; select=["z_cam"])

    focusList = sort([row.z_cam for row in focusList])

    Z, R = sff(imageList, focusList, 2, true)

    # Normalize Z to between 0 and 1
    # Z_norm = (Z.-minimum(Z))/(maximum(Z)-minimum(Z))
    # Z_norm = imageCenterValues(Z_norm)

    # Carve depth map with R
    # R<20 is unreliable
    # ZC = Z
    # ZC[R.<20] .= minimum(focusList)

    # RC = ones(size(R))
    # RC[R.<20] .= 0
    # RC = Gray.(RC)

    # Compute normals from depth map
    # normals = Depth2Normal(complement.(Z_norm))

    ### STATISTICS

    # normalsX = shiftNormalsRange(normals[:,:,1])
    # normalsY = shiftNormalsRange(normals[:,:,2])
    # normalsZ = shiftNormalsRange(normals[:,:,3])
    # normalsColor = colorview(RGB, normalsX, normalsY, normalsZ)


    # if outputFolder !== nothing
    #     folderName = basename(baseFolder)

    #     outputFolderConcatenated = outputFolder*"/SFF "*uppercase(kernel)*"/"
    #     ispath(outputFolderConcatenated) || mkpath(outputFolderConcatenated)

    #     save(outputFolderConcatenated*"/Z_"*folderName*"_SFF_"*kernel*".png", complement.(imageDisp01(Z)))
    #     save(outputFolderConcatenated*"/R_"*folderName*"_SFF_"*kernel*".png", imageNormalize(R))
    #     save(outputFolderConcatenated*"/RC_"*folderName*"_SFF_"*kernel*".png", imageNormalize(RC))
    #     save(outputFolderConcatenated*"/normals_"*folderName*"_SFF_"*kernel*".png", normalsColor)
    # end

    return Z, R
    # return Z, R, RC, normalsColor

end


base = "F:/Research/Generated Data/Blender/Statue du parc d'Austerlitz/SFF/Restrained Z/"

innerFolderList = ["5", "20", "100"]
kernelList = ["sobel", "sml"]

outputFolder = "F:/Image Out/Restrained Z/"

# Empty dictionaries to gather RMSE and inverse RMSE 
rmseList = Dict()
irmseList = Dict()

# Read in ground truth depth map for comparison
# NOTE: Using RTI collection due to original GT collection having a different camera setup
GT = Gray2Float64(Gray.(load("F:/Research/Generated Data/Blender/Statue du parc d'Austerlitz/RTI/Black Background - FullScale/Depth/Image0001.png")))

outputStructList = []

for f in innerFolderList
    for kernel in kernelList

        println("Running ", kernel, " for ", f)

        # Run SFF method
        Z,R = sff_main(base*f, kernel, nothing, false)

        if outputFolder !== nothing
            push!(outputStructList, FileSet(Z,R,f,"SFF",kernel))
        end

        # Normalize computed depth map so that it's placed from 0-1
        Z_normalized = complement.(imageDisp01(Z))

        rmseList[f,kernel] = rmse(GT, Z_normalized)
        irmseList[f,kernel] = 1-rmseList[f,kernel]
    end
end

ZMax, RMax = FindFileSetMax(outputStructList)
# println("NOTE: Using SFF-RTI ZMax and RMax normalization coefficients")
# ZMax = 8.5
# RMax = 292.7648

if outputFolder !== nothing
    WriteMaps(outputStructList, outputFolder)
end

WriteCSV(outputFolder*"/Ground truth comparison results SFF.csv", innerFolderList, ["sff"], kernelList, rmseList, irmseList, ZMax, RMax)
