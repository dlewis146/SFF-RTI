using Glob, ProgressMeter, ImageContrastAdjustment

# Include local files for functions
include("./IlluminationInvariance.jl")
include("./MagnificationAdjust.jl")
include("../SFF-RTI/sff.jl")
include("../Utilities/FilenameParser.jl")

base = "D:/Research Module/Generated Data/Blender/"
# base = "/media/david/USS Defiant/Research Module/Generated Data/Blender/"

folder_path = base * "Parade Shield/Shield_New/"
# folder_path = base * "Statue du parc d'Austerlitz/TIFF"
# acq_file = base * "Shield_25_New/configuration.acq"

method = "average"

if method == "illumination"
    # file_list = glob("LDR*Z_25*.png", folder_path)
    # file_list = glob("LDR*Z_00025*.png", folder_path)

    theta = 2.883
    phi = 12.326

    searchPattern = string("*Theta_", theta, "*Phi_", phi, "*.tif")
    # searchPattern = string("*Theta_", theta, "*Phi_", phi, "*.png")

    file_list = glob(searchPattern, folder_path)

    p = Progress(
        length(file_list),
        "Computing illumination invariant images...",
    )

    for idx in eachindex(file_list)
        zPos = ""

        stringArray = split(file_list[idx], "_")
        for j in range(1, stop = length(stringArray))
            if stringArray[j] == "Z"
                zPos = stringArray[j+1]
                break
            end
        end

        searchPatternLight = string("*Z_", zPos, "_*.tif")
        angle_file_list = glob(searchPatternLight, folder_path)

        focusMap = AssembleMatrices(angle_file_list)

        # filename = string("statue/lightGradient_Z_", zPos, ".csv")

        # CSV.write(filename, DataFrame(focusMap, :auto), header = false)

        ProgressMeter.next!(p; showvalues = [(:"Z Position", zPos)])

        folder_base = "paradeShieldNew/"

        imgeq = adjust_histogram(Equalization(), Gray.(focusMap / maximum(focusMap)), 256)

        # NOTE: Testing rpad for placing correct number of zeros at end of filenames
        save(folder_base * "Equalized/lightGradient_Z_" * rpad(zPos, 9, '0') * ".png" , imgeq)

        filename = string(folder_base, "lightGradient_Z_", rpad(zPos, 9, '0'), ".png")
        # save(filename, Gray.(focusMap / maximum(focusMap)))
        save(filename, focusMap ./ maximum(focusMap))

    end

elseif method == "magnification"

    r = 0 # [mm]
    # r = 32 # [mm]

    # file_list = glob("LDR*Z_00025*.png", folder_path)
    # file_list = glob("LDR*Z_25*.png", folder_path)
    file_list = glob("*.tif", folder_path)

    angleList = parseForAngles(file_list)

    for pair in angleList
        local searchPattern =
            "*Theta_" * string(pair[1]) * "*Phi_" * string(pair[2]) * "*"

        fileListAngular = glob(searchPattern, folder_path)

        imageDictNew = MagnificationAdjust(fileListAngular, acq_file, r)

        if !isdir(folder_path*"/Adjusted/")
            mkpath(folder_path*"/Adjusted/")
        end

        # for x in range(1, stop=length(imageListNew))
        idx = 1
        for (key, image) in imageDictNew

            filename = last(split(fileListAngular[idx], '/'))

            save(folder_path*"/Adjusted/"*filename, Gray.(image))

            idx += 1

        end

    end

elseif method == "shapefromfocusAll"

    # Get list of angle pairs (Θ, Φ) present in data set
    angleList = ParseForAngles(folder_path)

    # Get list of Z positions present in data set
    zPosList = ParseForZPositions(folder_path)

    p = Progress(length(angleList), "Computing SFF for all light positions...")
    for pair in angleList

        local imgList = []

        # Find images at current light position
        local searchPattern = string("*Theta_", pair[1], "_*Phi_", pair[2], ".tif")
        # searchPattern = string("*Theta_", pair[0], "_*Phi_", pair[1], "_*.tif")
        inputFileList = glob(searchPattern, folder_path)

        # Read in images and store in list
        # NOTE: Reading in grayscale
        for file in inputFileList
            img = Gray.(load(file))
            push!(imgList, img)
        end

        # Call SFF algorithm with focus stack of images
        local z = sff(imgList, zPosList)

        # Normalize image for range (0-1) since JuliaImages expects all image data
        # to fall within this range
        local zNorm = Gray.(z/maximum(z))

        # Create separate array with histogram equalization JUST FOR VISUALIZATION
        local zDISPLAY = adjust_histogram(zNorm, Equalization(256, 0, 1))
        # local zDISPLAY = adjust_histogram(Equalization(), Gray.(z/maximum(z)), 256)

        local filenameBase = string("sffOutput/depthMap_Theta_", pair[1], "_Phi_", pair[2])
        ext = ".png"

        # Save images
        save(filenameBase*ext, zNorm)
        save(filenameBase*"_DISPLAY"*ext, zDISPLAY)
        ProgressMeter.next!(p)
    end

elseif method == "shapefromfocus"

    # Get list of Z positions present in data set
    zPosList = ParseForZPositions(folder_path)

    imgList = []

    # For each Z position, compute illumination invariant image, convert to
    # `rawview`, then scale (0-255) and store in list. We want all the
    # images to be able to be compared to one another, hence the normalization.
    p = Progress(length(zPosList), "Computing illumination invariant images...")
    for zPos in zPosList

        # Find images at current Z position
        local searchPattern = string("*Z_", zPos, "_*.tif")
        inputFileList = glob(searchPattern, folder_path)

        # Compute illumination invariant image
        # focusMap = AssembleMatrices(inputFileList)
        focusMap = zero(Gray.(load(inputFileList[1])))
        for file in inputFileList

            img = Gray.(load(file))

            imgProcessedX, imgProcessedY = FilterImage(img, "sobel")

            imgProcessed = imgProcessedX + imgProcessedY

            focusMap = focusMap .+ imgProcessed
            # imgHold = imgHold .+ imgProcessed
        end


        # NOTE: Normalizing image and scaling (0-255)
        # focusMap = rawview(focusMap)
        # focusMap = focusMap ./ maximum(focusMap)
        # focusMap = focusMap * 255

        # Store image and iterate progress meter
        push!(imgList, focusMap)
        ProgressMeter.next!(p)
    end

    # Pass images into SFF function
    println("Calling Shape from Focus algorithm...")
    local z = sff(imgList, zPosList)

    # Normalize image for range (0-1) since JuliaImages expects all image data
    # to fall within this range
    local zNorm = Gray.(z/maximum(z))

    # Create separate array with histogram equalization JUST FOR VISUALIZATION
    local zDISPLAY = adjust_histogram(zNorm, Equalization(256, 0, 1))
    # local zDISPLAY = adjust_histogram(Equalization(), Gray.(z/maximum(z)), 256)

    local zInverted = Images.complement.(zNorm)

    filenameBase = "depthMapShield"

    # Save images
    save(filenameBase*".png", zNorm)
    save(filenameBase*"_DISPLAY.png", zDISPLAY)
    save(filenameBase*"_INVERTED.png", zInverted)

elseif method == "regionsegmentation"

    include("./RegionSegmentation.jl")

    type = "illumination"

    path_list = glob(
        "*.png",
        "/media/david/USS Defiant/julia/GramTest/images/lightGradients",
    )
    # path_list = glob("*_Theta_0_Phi_0_*.png", folder_path)

    imageSize = size(load(path_list[1]))

    focusSpaceArray = zeros(length(path_list), imageSize[1] * imageSize[2])

    # Iterate through all images, compute gradients, and place vectorized versions as array cols
    # for idx in range(1, stop=length(path_list))
    for idx in eachindex(path_list)
        img = load(path_list[idx])

        # Get X and Y derivatives and combine
        imgSobelX, imgSobelY = ComputeSobelImages(img)
        imgSobel = imgSobelX .+ imgSobelY

        # Flatten and place in array
        focusSpaceArray[idx, :] = vec(imgSobel)
    end

    # Pass focusSpaceArray to RegionSegmentation
    R = RegionSegmentation(focusSpaceArray)

elseif method == "average"

    fileList = glob("*.tif", base*"Parade Shield/Shield Single/TIFF")

    imgHold = zero(Gray.(load(fileList[1])))

    for file in fileList
        img = Gray.(load(file))

        imgProcessedX, imgProcessedY = FilterImage(img, "sobel")

        imgProcessed = imgProcessedX + imgProcessedY

        imgHold = imgHold .+ imgProcessed
    end

    imgHold = imgHold ./ length(fileList)

    imgHold = Gray.(imgHold ./ maximum(imgHold))

    save(base*"Parade Shield/Shield Single/imgAverage.png", abs.(imgHold))

end
