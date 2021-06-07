using DelimitedFiles, Images, DataStructures

function MagnificationAdjust(file_list, configPath, r)

    # Read first image to get image shape
    numCols = size(load(file_list[1]))[1]
    numRows = size(load(file_list[1]))[2]

    # Read .acq file for acquisition parameters and find motor positions / distance
    ZNb = 0
    Z1 = 0.0
    Z2 = 0.0

    fileArray = readdlm(configPath)
    for rowIdx in range(1, stop=size(fileArray)[1])
        if fileArray[rowIdx, 1] == "Z_Nb"
            ZNb = fileArray[rowIdx, 2]
        end
        if fileArray[rowIdx, 1] == "Z_1"
            Z1 = fileArray[rowIdx, 2]
        end
        if fileArray[rowIdx, 1] == "Z_2"
            Z2 = fileArray[rowIdx, 2]
        end
    end


    ZStep = (Z2 - Z1)/(ZNb-1)
    f = range(Z1, stop=Z2, length=ZNb)
    # f = range(Z1, stop=Z2, step=ZStep)

    if length(f) != length(file_list)
        throw(ArgumentError("Number of expected images doesn't match number of input images."))
    end

    # println("r: ", r)
    # println("Object distance: ", r)
    # println("Z1: ", Z1)
    # println("Z2: ", Z2)
    # println("ZNb: ", ZNb)
    # println()

    refDistance = 0.0
    targetRows = 0
    targetCols = 0
    if Z1 < Z2
        refDistance = r + Z1
        magRatio = ((r+Z1)/(r+Z2))
        targetRows = floor(Int, numRows * magRatio)
        targetCols = floor(Int, numCols * magRatio)
    end


    p = Progress(length(file_list), "Adjusting images for magnification...")

    file_dict = SortedDict()
    for x in range(1, stop=length(f))

        currentDistance = f[x] + r

        # NOTE: Reading in image as grayscale
        img = Gray.(load(file_list[x]))

        magRatio = refDistance / currentDistance

        scaleRows = (numRows * magRatio)
        scaleCols = (numCols * magRatio)

        rowsDiff = numRows - scaleRows
        colsDiff = numCols - scaleCols

        imgCropped = img[1+floor(Int, colsDiff/2):numCols-floor(Int, colsDiff/2),
                         1+floor(Int, rowsDiff/2):numRows-floor(Int, rowsDiff/2)]

        file_dict[currentDistance] = imresize(imgCropped, (targetCols, targetRows))

        ProgressMeter.next!(p; showvalues=[(:"Magnification Ratio", magRatio)])

    end
    return file_dict
end
