using LinearAlgebra, Glob, Images, ImageView, CSV, DataFrames, Printf
include("./sff.jl")
include("./sff-rti_utilities.jl")
include("../misc/ImageUtilities.jl")

global ACCEPTED_KERNELS = ["sml", "sobel"]

function sff_main(baseFolder, kernel="sml"; ksize=(5,5), outputFolder=nothing, write_maps=false)
       
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

if abspath(PROGRAM_FILE) == @__FILE__

    base = "J:/Research/Generated Data/Blender/Statue du parc d'Austerlitz/SFF/Restrained Z/"

    innerFolderList = ["100"]
    kernelList = ["sml"]
    ksizeList = [5]

    outputFolder = "J:/Image Out/R TEST/"
    # outputFolder = "J:/Image Out/Restrained Z/Window Size Tests/(5,5)/"

    # Create empty dictionaries to gather MS-SIM, MS-SSIM (Just structure) and PSNR 
    msssimDict = Dict()
    structureDict = Dict()
    psnrDict = Dict()

    compute_ssim = true

    # Read in ground truth depth map for comparison
    # NOTE: Using RTI collection due to original GT collection having a different camera setup
    GT = Gray2Float64(Gray.(load("J:/Research/Generated Data/Blender/Statue du parc d'Austerlitz/Ground Truth Restrained Z/Depth/Image0001.png")))

    outputStructList = []

    for f in innerFolderList
        for ksize in ksizeList
            for kernel in kernelList

                println("Running ", kernel, " for ", f)

                # Run SFF method
                Z, R = sff_main(base*f, kernel; ksize=(ksize, ksize), outputFolder=outputFolder, write_maps=false)

                if outputFolder !== nothing
                    push!(outputStructList, FileSet(Z,R,0,parse(Int,f),"SFF",kernel))
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
                    structureDict[f,"sff",kernel] = assess(iqi, GT, Z_normalized)

                    # ssimDict[f,"sff",kernel, ksize]  = assess_ssim(GT, Z_normalized)
                    msssimDict[f,"sff",kernel, ksize]  = assess_msssim(GT, Z_normalized)
                elseif !compute_ssim
                    structureDict[f,"sff",kernel, ksize] = NaN
                    msssimDict[f,"sff",kernel, ksize] = NaN
                end

                # Not computing PSNR for SFF collections, but we create dummy entries for compatibility with `WriteCSV`
                psnrDict[f,"sff",kernel, ksize] = NaN
            end
        end
    end

    ZMax, RMax = FindFileSetMax(outputStructList)

    if outputFolder !== nothing
        for ksize in ksizeList
            WriteMaps(outputStructList, outputFolder*string(ksize), nothing, nothing)
        end
    end

    WriteCSV(outputFolder*"/Ground truth comparison results SFF (5,5).csv", innerFolderList, ["sff"], kernelList, [5], structureDict, msssimDict, psnrDict, ZMax, RMax)

end
