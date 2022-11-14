using LinearAlgebra, Glob, Images, ImageView, CSV, DataFrames, Printf
include("./sff.jl")
include("./sff-rti_utilities.jl")
include("./ImageUtilities.jl")

global ACCEPTED_KERNELS = ["sml", "sobel"]

function sff_main(baseFolder, kernel="sml"; ksize=(3,3), outputFolder=nothing, write_maps=false)

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
        focusMap = FilterImageAverage(focusMap, ksize)

        if write_maps == true
            outputFolderConcatenated = string(outputFolder, "/FOCUS MAPS/SFF_", basename(baseFolder), " ", uppercase(kernel), "/")

            # Make sure all nested folders exist and create them if not
            ispath(outputFolder*"/FOCUS MAPS/") || mkpath(outputFolder*"/FOCUS MAPS/")
            ispath(outputFolderConcatenated) || mkpath(outputFolderConcatenated)

            save(outputFolderConcatenated*lpad(string(idx),length(string(length(fileList))), "0")*".png", imageNormalize(focusMap))
        end

        push!(imageList, focusMap)
    end

    # Read z_cam out of CSV and then extract from array of CSV.Row objects into a vector of Float64's
    focusList = CSV.File(csvPath; select=["z_cam"])

    # NOTE: Reversing focusList for consistency with SFF-RTI and GT where point closest to camera is lowest pixel value and further is highest value
    focusList = sort([row.z_cam for row in focusList], rev=true)

    Z = sff(imageList, focusList; sampleStep=2, median=true)

    return Z
end

function sff_handler()
    base = "J:/Research/Generated Data/Blender/Statue du parc d'Austerlitz/SFF/"
    # base = "J:/Research/Generated Data/Blender/Statue du parc d'Austerlitz/SFF/Restrained Z/"

    innerFolderList = ["5", "100"]
    # kernelList = ["sml"]
    kernelList = ["sobel", "sml"]

    outputFolder = "J:/Image Out/Statue (5,5)/"

    # Create empty dictionaries to gather MS-SIM, MS-SSIM (Just structure) and PSNR
    msssimDict = Dict()
    structureDict = Dict()
    psnrDict = Dict()

    compute_ssim = true

    # Read in ground truth depth map for comparison
    # NOTE: Using RTI collection due to original GT collection having a different camera setup
    GT = Gray2Float64(Gray.(load("J:/Research/Generated Data/Blender/Statue du parc d'Austerlitz/RTI/Black Background - FullScale/Depth/Image0001.png")))

    outputStructList = []

    for f in innerFolderList
        for kernel in kernelList

            println("Running ", kernel, " for ", f)

            # Run SFF method
            Z = sff_main(base*f, kernel; ksize=(5,5), outputFolder=outputFolder, write_maps=true)

            if outputFolder !== nothing
                push!(outputStructList, FileSet(Z,0,parse(Int,f),"SFF",kernel))
            end

            # Normalize computed depth map so that it's placed from 0-1
            # Z_normalized = imageDisp01(Z)
            Z_normalized = complement.(imageDisp01(Z))

            if compute_ssim
                # Compute SSIM and MS-SSIM then store all statistical measures in appropriate dictionaries
                # TEMP: Trying to compute MS-SSIM and MS-SSIM with just structural comparison taken into account
                """
                Here, the first parameter is the kernel used to weight the neighbourhood of each pixel while calculating the SSIM locally, and defaults to KernelFactors.gaussian(1.5, 11). The second parameter is the set of weights (α, β, γ) given to the lunimance (L), contrast (C) and structure (S) terms while calculating the SSIM, and defaults to (1.0, 1.0, 1.0). Recall that SSIM is defined as Lᵅ × Cᵝ × Sᵞ.
                Source: https://juliaimages.org/stable/examples/image_quality_and_benchmarks/structural_similarity_index/
                """

                iqi = MSSSIM(KernelFactors.gaussian(1.5,11), (0.0,0.0,1.0))
                structureDict[f,method,kernel] = assess(iqi, GT, Z_normalized)

                # ssim[f,method,kernel]  = assess_ssim(GT, Z_normalized)
                msssimDict[f,method,kernel]  = assess_msssim(GT, Z_normalized)
            elseif !compute_ssim
                structureDict[f,method,kernel] = NaN
                msssimDict[f,method,kernel] = NaN
            end
            psnrDict[f,"sff",kernel] = NaN
        end
    end

    ZMax = FindFileSetMax(outputStructList)

    if outputFolder !== nothing
        WriteMaps(outputStructList, outputFolder, nothing)
    end

    WriteCSV(outputFolder*"/Ground truth comparison results SFF.csv", innerFolderList, ["sff"], kernelList, structureDict, msssimDict, ZMax, psnrDict)
end