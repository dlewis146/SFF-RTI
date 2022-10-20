using Glob, Statistics, ProgressMeter

include("./IlluminationInvariance.jl")
include("./MagnificationAdjust.jl")
# include("./ReadAcqFile.jl")

function ShapeFromFocus(folder_path, theta, phi, method="illumination", configPath=nothing, r=32)

    # Attempt to cast numbers that are actually integers back into integers so
    # that file searches don't get screwed up by the difference between 0 and 0.0 for example
    # if isinteger(theta)
        # theta = Int(theta)
    # end
    if isinteger(theta) theta=Int(theta) end
    if isinteger(phi) phi = Int(phi) end

    # Take folder_path, find all images (*.tif) that match given theta and phi
    searchPattern = string("*Theta_", theta, "*Phi_", phi, "*.tif")

    file_list = glob(searchPattern, folder_path)

    numCols = size(load(file_list[1]))[1]
    numRows = size(load(file_list[1]))[2]

    if configPath !== nothing
        println("Reading .acq file...")
        # Parse .acq file
        # acqTuple = ReadAcqFile(configPath)

        # Determine camera steps
        # ZStep = (acqTuple[3] - acqTuple[2]) / acqTuple[1]
        # ZStep = (Z2 - Z1) / ZNb
        # f = range(acqTuple[2], stop=acqTuple[3], step=ZStep)
    end

    # Compute reference distance and target size for resizing
    # refDistance = 0.0
    # targetRows = 0
    # targetCols = 0

    # if Z1 < Z2
    #     refDistance = r + Z1
    #
    #     magRatio = (r+Z1)/(r+Z2)
    #     # magRatio = (r+Z2)/(r+Z1)
    #     targetRows = floor(Int, numRows * magRatio)
    #     targetCols = floor(Int, numCols * magRatio)
    # end

    # Get image dimensions
    # imgSize = size(load(file_list[1]))
    # numCols = imgSize[1]
    # numRows = imgSize[2]

    # Create empty arrays for comparison
    # focus_compare = zeros(targetCols, targetRows)
    focus_compare = zeros(numCols, numRows)

    focus_idx_map = zero(focus_compare) # Similar functionality to 'zeros_like'

    # Create progress bar
    p = Progress(length(file_list), "Computing Shape from Focus map...")
    # p = Progress(length(file_list))
    # update!(p, 1)

    # magRatio = 0.0

    for idx in range(1, stop=length(file_list))
        # Read through images, compute focus maps, place in dictionary
        # @printf("\n Image: %d", idx)

        # magRatio = 0.0

        currentDistance = ""

        # Find current magnification ratio for cropping/resizing process
        stringArray = split(file_list[idx], "_")
        for j in range(1, stop=length(stringArray))
            if stringArray[j] == "Z"
                currentDistance = parse(Float64, stringArray[j+1])
                break
            end
        end

        if currentDistance == ""
            throw(ErrorException("Input images do not contain Z in their filename."))
        end

        # currentDistance = f[parse(Int, zPos)+1] + r
        # magRatio = refDistance / currentDistance

        focusMap = []
        if method == "sobel"
            img = load(file_list[idx])
            imgdx, imgdy = ComputeSobelImages(img)
            focusMap = imgdx + imgdy

        elseif method == "illumination"
            searchPatternLight = string("*Z_", currentDistance, "_*.tif")
            angle_file_list = glob(searchPatternLight, folder_path)
            focusMap = AssembleMatrices(angle_file_list)
        end

        # interestPts = findall(imgResized .> focus_compare)
        interestPts = findall(focusMap .> focus_compare)

        # focus_compare[interestPts] .= imgResized[interestPts]
        focus_compare[interestPts] .= focusMap[interestPts]
        focus_idx_map[interestPts] .= idx

        # Update progress bar
        ProgressMeter.next!(p; showvalues= [(:"Current Distance", round(currentDistance, digits=6))])
        # update!(p, idx)
    end

    return focus_idx_map
end
