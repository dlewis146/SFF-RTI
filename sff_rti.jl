# using Images
using Glob, CSV, DataFrames, CoordinateTransformations, LinearAlgebra, ProgressMeter, Printf

include("./sff.jl")
include("./IlluminationInvariance.jl")
include("./ImageUtilities.jl")

global ACCEPTED_METHODS = ["fvg", "mean", "std"]
global ACCEPTED_KERNELS = ["sml", "sobel"]

function sff_rti(baseFolder, method="fvg", kernel="sobel", outputFolder=nothing, write_maps=false)

    # Make sure that `method` and `kernel` are lowercase for parsing compatibility
    method = lowercase(method)
    kernel = lowercase(kernel)

    # Make sure that given method and kernel are compatible with methods
    if !(method in ACCEPTED_METHODS)
        error(@printf("\nGiven method { %s } isn't compatible with current version of `sff_rti`\n", method))
    end
    if !(kernel in ACCEPTED_KERNELS)
        error(@printf("\nGiven kernel { %s } isn't compatible with current version of `sff_rti`\n", kernel))
    end

    # Make sure there aren't multiple CSV files in the base path
    if length(glob("*.csv", baseFolder)) > 1
        error("\nMore than 1 CSV file found in given folder\n")
    elseif length(glob("*.csv", baseFolder)) == 0
        error("\nCSV file not found in given folder\n")
    end

    csvPath = glob("*.csv", baseFolder)[1]
    folderPath = baseFolder * "/Renders/"

    # Get all the XYZ camera positions in CSV.row structs
    rowList = CSV.File(csvPath; select=["image", "x_lamp", "y_lamp", "z_lamp", "z_cam"])

    # Get all Z positions and remove duplicates
    # NOTE: Sort to make sure they're stored in ascending order
    zPosList = sort(unique([row.z_cam for row in CSV.File(csvPath; select=["z_cam"])]), rev=true)
    # zPosList = sort(unique([row.z_cam for row in CSV.File(csvPath; select=["z_cam"])]))
    println("NOTE: Reversing zPosList")

    fvgList = []
    global offsetAmount = 0

    # Not a NECESSARY addition, but it just makes it a bit clearer visually which method is currently being run
    if method == "fvg"
        prog = Progress(length(zPosList), "Computing full vector gradient images...")
    elseif method == "std"
        prog = Progress(length(zPosList), "Computing standard deviation images...")
    else
        prog = Progress(length(zPosList), "Computing "*method*" images...")
    end

    # Compute photometric-invariant image representation for each Z position
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

        # Compute desired focus map 
        focusMap = nothing
        if method == "fvg"
            # Compute FVG
            focusMap = ComputeFullVectorGradient(fileList, angleList, kernel)
        else
            gradientList = []
            for file in fileList
                img = Gray.(load(file))
                push!(gradientList, FilterImageCombined(img, kernel))
            end
            
            if method == "mean"
                focusMap = ComputeMeanImage(gradientList)
            elseif method == "std"
                focusMap = ComputeSTDImage(gradientList)
            end
        end

        # Average focus map
        focusMap = FilterImageAverage(focusMap)

        # Centering image values attempts to push dynamic range so that all values are positive and will no longer give errors when computing logs during Gaussian interpolation
        # focusMap = abs.(focusMap)
        # focusMap = imageCenterValues(focusMap)
        # NOTE: Let's just try to offset the image values by a constant amount so that they're all remaining on the same rough scale
        # focusMap = focusMap .+ 500
        # focusMap = focusMap .+ 255

        if minimum(focusMap) < offsetAmount
            global offsetAmount = minimum(focusMap)
        end

        # println(minimum(focusMap))

        # if write_maps == true

        #     outputFolderConcatenated = string(outputFolder, "/FOCUS MAPS/", basename(baseFolder), " ", uppercase(method), " ", uppercase(kernel), "/")

        #     # Make sure all nested folders exist and create them if not
        #     ispath(outputFolder*"/FOCUS MAPS/") || mkpath(outputFolder*"/FOCUS MAPS/")
        #     ispath(outputFolderConcatenated) || mkpath(outputFolderConcatenated)

        #     save(outputFolderConcatenated*string(idx)*".png", imageNormalize(focusMap))
        # end

        # Store focus map in list and iterate progress bar
        push!(fvgList, focusMap)
        ProgressMeter.next!(prog; showvalues= [(:"Current Distance", zPosList[idx]), (:"Offset Amount", offsetAmount)])
    end

    # Add offset amount to all focus maps
    # for focusMap in fvgList
    for idx in eachindex(fvgList)
        fvgList[idx] = fvgList[idx] .+ abs(offsetAmount)

        if write_maps == true

            outputFolderConcatenated = string(outputFolder, "/FOCUS MAPS/", basename(baseFolder), " ", uppercase(method), " ", uppercase(kernel), "/")

            # Make sure all nested folders exist and create them if not
            ispath(outputFolder*"/FOCUS MAPS/") || mkpath(outputFolder*"/FOCUS MAPS/")
            ispath(outputFolderConcatenated) || mkpath(outputFolderConcatenated)

            save(outputFolderConcatenated*string(idx)*".png", imageNormalize(fvgList[idx]))
        end

    end
    # [focusMap .+ offsetAmount for focusMap in fvgList]

    println("Computing depth map...")
    Z, R = sff(fvgList, zPosList, 2, true)

    # TEMP NORMALIZE 0-1
    # Z = (Z.-minimum(Z))/(maximum(Z)-minimum(Z))

    # Compute normals from depth map
    # println("Computing surface normals from depth map...")
    # normals = Depth2Normal(Z)
    # normals = Depth2Normal(complement.(Z))

    # Carve depth map with R
    # ZC = copy(Z)
    # ZC[R.<20] .= minimum(zPosList)

    # RC = copy(R)
    # RC[R.<20] .= 0
    # RC = Gray.(RC)

    ### VISUALIZATION & IO
    # normalsX = (normals[:,:,1])
    # normalsY = (normals[:,:,2])
    # normalsZ = (normals[:,:,3])

    # normalsX = shiftNormalsRange(normals[:,:,1])
    # normalsY = shiftNormalsRange(normals[:,:,2])
    # normalsZ = shiftNormalsRange(normals[:,:,3])

    # normalsColor = colorview(RGB, normalsX, normalsY, normalsZ)

    # If given output folder, write out images
    # if outputFolder !== nothing
    #     folderName = basename(baseFolder)

    #     outputFolderConcatenated = outputFolder*"/"*uppercase(method)*" "*uppercase(kernel)*"/"
    #     ispath(outputFolderConcatenated) || mkpath(outputFolderConcatenated)

    #     save(outputFolderConcatenated*"/Z_"*folderName*"_"*method*"_"*kernel*".png", imageNormalize(Z))
    #     save(outputFolderConcatenated*"/R_"*folderName*"_"*method*"_"*kernel*".png", imageNormalize(R))
    #     save(outputFolderConcatenated*"/RC_"*folderName*"_"*method*"_"*kernel*".png", imageNormalize(RC))
    #     save(outputFolderConcatenated*"/normals_"*folderName*"_"*method*"_"*kernel*".png", normalsColor)
    # end

    return Z, R
end
