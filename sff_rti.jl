# using Images
using Glob, CSV, DataFrames, CoordinateTransformations, ProgressMeter, Printf

include("./sff.jl")
include("./IlluminationInvariance.jl")
include("../misc/ImageUtilities.jl")

global ACCEPTED_METHODS = ["fvg", "mean", "max"]
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
        folderName = folderName*"/SFF/"*splitpath(baseFolder)[lastindex(splitpath(baseFolder))-1]*"/"*string(numSFF)*"/"
    else
        # Get base path up through the main object folder
        folderName = join(splitpath(baseFolder)[1:lastindex(splitpath(baseFolder))-2], "/")

        # Create expected SFF folder path
        folderName = folderName*"/SFF/"*string(numSFF)*"/"
    end

    # Make sure directory exists
    if isdir(folderName) == false
        error(@printf("\nCannot find equivalent SFF acquisition for `%s`\n", basename(baseFolder)))
    end

    return folderName
end

function ComputeMultiLightGradients(baseFolder::String, method::String, kernel::String; ksize::Int=5, write_maps::Bool=false, outputFolder::String="", compute_psnr::Bool=false)

    # # Make sure there aren't multiple CSV files in the base path
    # if length(glob("*.csv", baseFolder)) > 1
    #     error("\nMore than 1 CSV file found in given folder\n")
    # elseif length(glob("*.csv", baseFolder)) == 0
    #     println(baseFolder)
    #     error("\nCSV file not found in given folder\n")
    # end

    println(ksize)

    # csvPath = glob("*.csv", baseFolder)[1]
    csvPath = baseFolder*"/Image.csv"
    folderPath = baseFolder * "/Renders/"

    # Get all the XYZ camera positions in CSV.row structs
    rowList = CSV.File(csvPath; select=["image", "x_lamp", "y_lamp", "z_lamp", "z_cam"])

    sffCSV = nothing
    sffFolderPath = ""
    psnrDict = Dict()

    # If computing PSNR, find equivalent SFF data collection and parse CSV for filenames
    if compute_psnr == true
        sffFolderPath = FindSFFEquivalent(baseFolder)

        csvPathSFF = sffFolderPath*"/Image.csv"
        sffCSV = CSV.File(csvPathSFF; select=["image", "z_cam"])
    end


    # Get all Z positions and remove duplicates
    # NOTE: Sort to make sure they're stored in ascending order
    zPosList = sort(unique([row.z_cam for row in CSV.File(csvPath; select=["z_cam"])]), rev=true)
    println("NOTE: Reversing zPosList")

    println("NOTE: Adding constant of 10 to x_lamp,y_lamp,z_lamp in attempt to keep all values about 0")

    gradientList = []
    global offsetAmount = 0
    global focusMapNormCoeff = 0

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

        # Use CSV to generate filelist and angle objects.
        for row in rowList

            # Conditional to make sure we're only dealing with one Z position at a time
            if row.z_cam != zPosList[idx]
                continue
            else
                # Construct and store filename based off of CSV entries
                push!(fileList, folderPath*row.image*".png")

                # Place all coordinates into LightAngle object list
                push!(angleList, LightAngle(row.x_lamp+10, row.y_lamp+10, row.z_lamp+10))
            end
        end

        # Compute desired focus map
        focusMap = nothing
        if method == "fvg"
            # Compute FVG
            focusMap = ComputeFullVectorGradient(fileList, angleList, kernel, ksize)
        else
            gradientListIntermediate = []
            for file in fileList
                img = Gray.(load(file))
                push!(gradientListIntermediate, FilterImageCombined(img, kernel, ksize))
            end

            if method == "mean"
                focusMap = ComputeMeanImage(gradientListIntermediate)
            elseif method == "max"
                focusMap = maximum(ImageList2Cube(gradientListIntermediate), dims=3)[:,:,1]
            end
        end

        # Average focus map
        focusMap = FilterImageAverage(focusMap)

        ##NOTE: Offsetting image values attempts to push dynamic range so that all values are positive and will no longer give errors when computing logs during Gaussian interpolation

        # Track maximum required offset value so that we can keep all images to be compared relative to one another.
        global offsetAmount = minimum([minimum(focusMap), offsetAmount])

        if compute_psnr == true
            for sffRow in sffCSV
                if sffRow.z_cam == zPosList[idx]

                    sffImg = Gray.(load(sffFolderPath*"/Renders/"*sffRow.image*".png"))

                    sffFocusMap = FilterImageCombined(sffImg, kernel, ksize)

                    # Determine what value to use for image normalization
                    snrNormalizationCoefficient = maximum([maximum(focusMap), maximum(sffFocusMap)])

                    # NOTE: Avoiding usage of built-in PSNR assessment function (assess_psnr) because it appears to be checking the peak values of the array to try and figure out the dynamic range for normalization and I'm not super convinced by that considering Julia is handling images a bit differently than other languages
                    psnrDict[zPosList[idx]] = 10*log10(1/mse(sffFocusMap./snrNormalizationCoefficient, focusMap./snrNormalizationCoefficient))

                end
            end

        end

        # Store focus map in list and iterate progress bar
        push!(gradientList, focusMap)
        ProgressMeter.next!(prog; showvalues= [(:"Current Distance", zPosList[idx]), (:"Offset Amount", offsetAmount), (:"Number of light angles", length(angleList)), (:"Number of Z positions", length(zPosList))])
    end

    # Add offset amount to all focus maps
    gradientList = [gradientList[idx] = gradientList[idx].+abs(offsetAmount) for idx in eachindex(gradientList)]
    focusMapNormCoeff = maximum(maximum.(gradientList))

    if write_maps == true
        for idx in eachindex(gradientList)
            outputFolderConcatenated = string(outputFolder, "/FOCUS MAPS/", basename(baseFolder), "/", uppercase(method), " ", uppercase(kernel), " KSIZE $ksize",  "/")

            # Make sure all nested folders exist and create them if not
            ispath(outputFolder*"/FOCUS MAPS/") || mkpath(outputFolder*"/FOCUS MAPS/")
            ispath(outputFolderConcatenated) || mkpath(outputFolderConcatenated)

            save(outputFolderConcatenated*lpad(string(idx), length(string(length(gradientList))), "0")*".png", gradientList[idx]./focusMapNormCoeff)
        end
    end

    return gradientList, zPosList, psnrDict

end

function sff_rti(baseFolder, innerFolderList::Array{String}, methodList::Array{String}, kernelList::Array{String}; ksizeList=[3], write_maps=false, write_csv=false, compute_psnr=false, compute_ssim=false, outputFolder="", gtPath="", ZMax=NaN)
    """
    Function that takes in a folder which contains a subfolder of simulated images and a CSV corresponding to the RTI and SFF parameters for those images. Multi-light integration methods are applied for each focus position (full vector gradient, mean gradient reponse, or maximum gradient response) and then are all passed into a Shape from Focus function to generate a depth map.

    Computation of MS-SSIM and PSNR are also possible for statistical comparison. The MS-SSIM requires the path to a ground truth depth map and the PSNR requires there to be an equivalent SFF acquisition (see `FindSFFEquivalent` called in the function `ComputeMultiLightGradients`).
    """

    # Read in ground truth for SSIM comparison if needed
    GT = nothing

    if compute_ssim == true && isfile(gtPath)
        GT = Gray2Float64(load(gtPath))
    elseif compute_ssim == true && !isfile(gtPath)
        throw(ArgumentError("Ground truth image filepath doesn't point to an existing file!\nInput: $gtPath"))
    end


    # Create empty dictionaries to gather MS-SIM, MS-SSIM (Just structure) and PSNR
    structureDict = Dict()
    msssimDict = Dict()
    psnrDict = Dict()

    # outputStructList = Dict()
    outputStructList = []
    for f in innerFolderList
        for method in methodList
            for kernel in kernelList
                for ksize in ksizeList

                    println()

                    # Make sure that `method` and `kernel` are lowercase for parsing compatibility
                    method = lowercase(method)
                    kernel = lowercase(kernel)

                    # Call function to compute multi-light gradients using desired method
                    gradientList, zPosList, psnrDictCurrent = ComputeMultiLightGradients(baseFolder*f, method, kernel; ksize=(ksize,ksize), write_maps=write_maps, compute_psnr=compute_psnr)

                    # Compute SFF
                    Z, R = sff(gradientList, zPosList; sampleStep=2, median=true)

                    psnrMean = 0.0
                    if compute_psnr
                        psnrMean = mean(values(psnrDictCurrent))
                    end

                    numRTI, numSFF = ParseFolderName(f)

                    if outputFolder !== nothing
                        push!(outputStructList, FileSet(copy(Z),copy(R),numRTI,numSFF,method,kernel))
                    end

                    # Normalize computed depth map so that it's placed from 0-1
                    Z_normalized = imageDisp01(Z)

                    if compute_ssim
                        # Compute structural comparison measure and MS-SSIM then store all statistical measures in appropriate dictionaries
                        """
                        Here, the first parameter is the kernel used to weight the neighbourhood of each pixel while calculating the SSIM locally, and defaults to KernelFactors.gaussian(1.5, 11). The second parameter is the set of weights (α, β, γ) given to the lunimance (L), contrast (C) and structure (S) terms while calculating the SSIM, and defaults to (1.0, 1.0, 1.0). Recall that SSIM is defined as Lᵅ × Cᵝ × Sᵞ.
                        Source: https://juliaimages.org/stable/examples/image_quality_and_benchmarks/structural_similarity_index/
                        """

                        iqi = MSSSIM(KernelFactors.gaussian(1.5,11), (0.0,0.0,1.0))
                        structureDict[f,method,kernel,ksize] = assess(iqi, GT, Z_normalized)

                        # ssim[f,method,kernel]  = assess_ssim(GT, Z_normalized)
                        msssimDict[f,method,kernel,ksize]  = assess_msssim(GT, Z_normalized)
                    elseif !compute_ssim
                        structureDict[f,method,kernel,ksize] = NaN
                        msssimDict[f,method,kernel,ksize] = NaN
                    end

                    psnrDict[f,method,kernel,ksize] = psnrMean

                end
            end
        end
    end

    # if isnan(ZMax)
    #     ZMax = FindFileSetMax(outputStructList)
    # else
    #     println("NOTE: Using given ZMax normalization coefficient")
    # end

    ZMax, RMax = FindFileSetMax(outputStructList)

    if write_maps == true
        for ksize in ksizeList
            WriteMaps(outputStructList, outputFolder*string(ksize), ZMax, RMax)
        end
    end

    if write_csv == true
        # csvPath = @printf("%s/Ground truth comparison results (%i,%i).csv", outputFolder, ksize[1], ksize[2])
        csvPath = outputFolder * "/Ground truth comparison results.csv"

        WriteCSV(csvPath, innerFolderList, methodList, kernelList, ksizeList, structureDict, msssimDict, psnrDict, ZMax, RMax)
    end

end


"""
Meant to be used for handling a single acquisition with a single set of input parameters. Maybe unnecessary, but with all the file handling and computation/writing of stats, this could be useful.

Copied on Nov. 26, 2022
"""
function sff_rti(folderPath, methodList, kernelList, ksizeList; write_maps=false, write_csv=false, compute_psnr=false, compute_ssim=false, outputFolder="", gtPath="", ZMax=NaN)
    """
    Function that takes in a folder which contains a subfolder of simulated images and a CSV corresponding to the RTI and SFF parameters for those images. Multi-light integration methods are applied for each focus position (full vector gradient, mean gradient reponse, or maximum gradient response) and then are all passed into a Shape from Focus function to generate a depth map.

    Computation of MS-SSIM and PSNR are also possible for statistical comparison. The MS-SSIM requires the path to a ground truth depth map and the PSNR requires there to be an equivalent SFF acquisition (see `FindSFFEquivalent` called in the function `ComputeMultiLightGradients`).
    """

    # Read in ground truth for SSIM comparison if needed
    GT = nothing

    if compute_ssim == true && isfile(gtPath)
        GT = Gray2Float64(load(gtPath))
    elseif compute_ssim == true && !(isfile(gtPath))
        throw(ArgumentError("Ground truth image doesn't exist at given path:  $gtPath"))
    end

    # Create empty dictionaries to gather MS-SIM, MS-SSIM (Just structure) and PSNR
    structureDict = Dict()
    msssimDict = Dict()
    psnrDict = Dict()

    # outputStructList = Dict()
    outputStructList = []
    ZMaxDict = Dict()
    RMaxDict = Dict()

    f = basename(folderPath)
    
    for method in methodList
        for kernel in kernelList
            for ksize in ksizeList
                println()

                # Make sure that `method` and `kernel` are lowercase for parsing compatibility
                method = lowercase(method)
                kernel = lowercase(kernel)

                println("Running SFF-RTI with variables:\nfolder: $f\nmethod: $method\nkernel: $kernel\nksize: $ksize\n")
                # Call function to compute multi-light gradients using desired method
                gradientList, zPosList, psnrDictCurrent = ComputeMultiLightGradients(folderPath, method, kernel; ksize=ksize, write_maps=write_maps, outputFolder=outputFolder, compute_psnr=compute_psnr)

                # Compute SFF
                Z, R = sff(gradientList, zPosList; sampleStep=2, median=true)

                psnrMean = 0.0
                if compute_psnr
                    psnrMean = mean(values(psnrDictCurrent))
                end

                numRTI, numSFF = ParseFolderName(f)

                if outputFolder !== nothing
                    WriteMapSingle(Z, R, numRTI, numSFF, method, kernel, ksize, outputFolder)
                    ZMaxDict[f,method,kernel,ksize] = maximum(Z)
                    RMaxDict[f,method,kernel,ksize] = maximum(R)

                    # push!(outputStructList, FileSet(Z,R,numRTI,numSFF,method,kernel))
                end

                # Normalize computed depth map so that it's placed from 0-1
                Z_normalized = imageDisp01(Z)

                if compute_ssim
                    # Compute structural comparison measure and MS-SSIM then store all statistical measures in appropriate dictionaries
                    """
                    Here, the first parameter is the kernel used to weight the neighbourhood of each pixel while calculating the SSIM locally, and defaults to KernelFactors.gaussian(1.5, 11). The second parameter is the set of weights (α, β, γ) given to the lunimance (L), contrast (C) and structure (S) terms while calculating the SSIM, and defaults to (1.0, 1.0, 1.0). Recall that SSIM is defined as Lᵅ × Cᵝ × Sᵞ.
                    Source: https://juliaimages.org/stable/examples/image_quality_and_benchmarks/structural_similarity_index/
                    """

                    iqi = MSSSIM(KernelFactors.gaussian(1.5,11), (0.0,0.0,1.0))
                    structureDict[f,method,kernel,ksize] = assess(iqi, GT, Z_normalized)

                    # ssim[f,method,kernel]  = assess_ssim(GT, Z_normalized)
                    msssimDict[f,method,kernel,ksize]  = assess_msssim(GT, Z_normalized)
                elseif !compute_ssim
                    structureDict[f,method,kernel,ksize] = NaN
                    msssimDict[f,method,kernel,ksize] = NaN
                end

                psnrDict[f,method,kernel,ksize] = psnrMean
            end
        end
    end

    # if isnan(ZMax)
    #     ZMax = FindFileSetMax(outputStructList)
    # else
    #     println("NOTE: Using given ZMax normalization coefficient")
    # end

    # ZMax, RMax = FindFileSetMax(outputStructList)

    # if write_maps == true
    #     for ksize in ksizeList
    #         WriteMaps(outputStructList, outputFolder*string(ksize), ZMax, RMax)
    #     end
    # end

    if write_csv == true
        # csvPath = @printf("%s/Ground truth comparison results (%i,%i).csv", outputFolder, ksize[1], ksize[2])
        csvPath = outputFolder * "/Ground truth comparison results.csv"

        WriteCSVSingles(csvPath, [f], methodList, kernelList, ksizeList, structureDict, msssimDict, psnrDict, ZMaxDict, RMaxDict)
        # WriteCSV(csvPath, [f], methodList, kernelList, ksizeList, structureDict, msssimDict, psnrDict, ZMax, RMax)
    end

end
