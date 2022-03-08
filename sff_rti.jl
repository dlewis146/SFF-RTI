# using Images
using Glob, CSV, DataFrames, CoordinateTransformations, LinearAlgebra, ProgressMeter, Printf

include("./sff.jl")
include("./IlluminationInvariance.jl")
include("./ImageUtilities.jl")

global ACCEPTED_METHODS = ["fvg", "mean", "std"]
global ACCEPTED_KERNELS = ["sml", "sobel"]

function FindSFFEquivalent(baseFolder)
    """
    This highly specific function takes in a given foldername for an SFF-RTI acquisition
    (named using standardized format), parses it for the number of Z (SFF) positions, 
    and searches for an equivalent SFF acquisition in the same folder that `SFFRTI` is
    located. Assumes that the SFF acquisitions are in folders named solely with their
    number of Z positions. 

    Example: SFF-RTI acquisition folder: `.../SFFRTI/RTI_4_SFF_20/`
             Equivalent SFF folder:      `.../SFF/20/`
    """

    _, numSFF = ParseFolderName(basename(baseFolder))

    # NOTE: Assumes path format of "{DRIVE}/Research/Generated Data/Blender/{OBJECT TITLE}/{METHOD}/
    folderName = ""

    if splitpath(baseFolder)[lastindex(splitpath(baseFolder))-1] != "SFFRTI"
        
        # Get base path up through the main object folder
        folderName = join(splitpath(baseFolder)[1:lastindex(splitpath(baseFolder))-3], "/")

        # Create expected SFF folder path
        folderName = folderName*"/SFF/"*splitpath(baseFolder)[lastindex(splitpath(baseFolder))-1]*"/"*numSFF*"/"
    else
        # Get base path up through the main object folder
        folderName = join(splitpath(baseFolder)[1:lastindex(splitpath(baseFolder))-2], "/")

        # Create expected SFF folder path
        folderName = folderName*"/SFF/"*numSFF*"/"
    end
   
    # Make sure there aren't multiple CSV files in the base path
    if isdir(folderName) == false
        error(@printf("\nCannot find equivalent SFF acquisition for `%s`\n", basename(baseFolder)))
    end

    return folderName
end

function sff_rti(baseFolder, method="fvg", kernel="sobel", outputFolder=nothing, write_maps=false, compute_snr=false)

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

    sffCSV = nothing
    sffFolderPath = ""
    snrDict = Dict()
    psnrDict = Dict()

    if compute_snr == true
        # Find equivalent SFF data collection
        sffFolderPath = FindSFFEquivalent(baseFolder)

        # Make sure there aren't multiple CSV files in the base path
        if length(glob("*.csv", sffFolderPath)) > 1
            error(@printf("\nMore than 1 CSV file found in equivalent SFF folder : `%s`\n", sffFolder))
        elseif length(glob("*.csv", sffFolderPath)) == 0
            error(@printf("\nCSV file not found in equivalent SFF folder : `%s`\n", sffFolder))
        end

        csvPathSFF = glob("*.csv", sffFolderPath)[1]
        sffCSV = CSV.File(csvPathSFF; select=["image", "z_cam"])
        # image = CSV.File(csvPathSFF; select=[""]
    end


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

        # Offsetting image values attempts to push dynamic range so that all values are positive and will no longer give errors when computing logs during Gaussian interpolation

        # Track maximum required offset value so that we can keep all images to be compared relative to one another. 
        global offsetAmount = minimum([minimum(focusMap), offsetAmount])

        if compute_snr == true

            for sffRow in sffCSV

                if sffRow.z_cam == zPosList[idx]

                    sffImg = Gray.(load(sffFolderPath*"/Renders/"*sffRow.image*".png"))

                    sffFocusMap = FilterImageCombined(sffImg, kernel)

                    # Determine what value to use for image normalization
                    snrNormalizationCoefficient = maximum([maximum(focusMap), maximum(sffFocusMap)])

                    # NOTE: Am I doing the right thing in normalizing these to be on then same scale???? Should I just be normalizing or perhaps not touching them at all?
                    snrDict[zPosList[idx]] = 10*log10(mean(focusMap./snrNormalizationCoefficient)/rmse(sffFocusMap./snrNormalizationCoefficient, focusMap./snrNormalizationCoefficient))

                    psnrDict[zPosList[idx]] = 10*log10(1/mse(sffFocusMap./snrNormalizationCoefficient, focusMap./snrNormalizationCoefficient))

                end
            end

        end



        # Store focus map in list and iterate progress bar
        push!(fvgList, focusMap)
        ProgressMeter.next!(prog; showvalues= [(:"Current Distance", zPosList[idx]), (:"Offset Amount", offsetAmount)])
    end 

    # Add offset amount to all focus maps
    # for focusMap in fvgList

    fvgList = [fvgList[idx] = fvgList[idx].+abs(offsetAmount) for idx in eachindex(fvgList)]
    focusMapNormCoeff = maximum(maximum.(fvgList))

    for idx in eachindex(fvgList)
        # fvgList[idx] = fvgList[idx] .+ abs(offsetAmount)

        if write_maps == true

            outputFolderConcatenated = string(outputFolder, "/FOCUS MAPS/", basename(baseFolder), " ", uppercase(method), " ", uppercase(kernel), "/")

            # Make sure all nested folders exist and create them if not
            ispath(outputFolder*"/FOCUS MAPS/") || mkpath(outputFolder*"/FOCUS MAPS/")
            ispath(outputFolderConcatenated) || mkpath(outputFolderConcatenated)

            save(outputFolderConcatenated*lpad(string(idx), length(string(length(fvgList))), "0")*".png", fvgList[idx]./focusMapNormCoeff)
        end

    end

    println("Computing depth map...")
    Z, R = sff(fvgList, zPosList, 2, true)

    if compute_snr == false
        snrDict = Dict(NaN=>NaN)
        psnrDict = Dict(NaN=>NaN)
    end

    return Z, R, mean(values(snrDict)), mean(values(psnrDict))
end