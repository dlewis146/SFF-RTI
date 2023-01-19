using LinearAlgebra, Glob, Images, CSV, DataFrames, Printf
include("./sff.jl")
include("./sff-rti_utilities.jl")
include("../misc/ImageUtilities.jl")

global ACCEPTED_KERNELS = ["sml", "sobel"]

function sff_main(baseFolder, kernel::String="sml"; ksize::Int=5, outputFolder::String=nothing, write_maps::Bool=false)

    # Make sure that given kernel is compatible with methods
    if !(kernel in ACCEPTED_KERNELS)
        error(@printf("\nGiven kernel { %s } isn't compatible with current version of `sff_main`\n", kernel))
    end

    # Make sure there aren't multiple CSV files in the base path
    # if length(glob("*.csv", baseFolder)) > 1
    #     error("\nMore than 1 CSV file found in given folder\n")
    # end

    csvPath = baseFolder*"/Image.csv"
    # csvPath = glob("*.csv", baseFolder)[1]

    # Get all interesting files
    folderPath = baseFolder*"/Renders/"

    # NOTE: Using this instead of just globbing all files so that they're read in and stored in the proper order. This doesn't happen by default due to a lack of zero-padding of frame numbers in file names
    rowList = CSV.File(csvPath; select=["image", "z_cam"])
    fileList = []
    zPosList = []

    for row in rowList
        push!(fileList, folderPath*row.image*".png")
        push!(zPosList, row.z_cam)
    end

    # NOTE: Reversing focusList for consistency with SFF-RTI and GT where point closest to camera is lowest pixel value and further is highest value
    # focusList = sort([row.z_cam for row in focusList], rev=true)


    p = sortperm(zPosList)
    fileList = reverse(fileList[p])
    zPosList = reverse(zPosList[p])


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


    Z, R = sff(imageList, zPosList; sampleStep=2, median=true)

    return Z, R
end

function sff_handler(folderPath, kernelList, ksizeList; write_maps=false, outputFolder="", compute_ssim=true, gtPath="", write_csv=true)


    # Create empty dictionaries to gather MS-SIM, MS-SSIM (Just structure) and PSNR
    msssimDict = Dict()
    structureDict = Dict()
    rmseDict = Dict()
    psnrDict = Dict()

    # Read in ground truth for SSIM comparison if needed
    GT = nothing

    if compute_ssim == true && isfile(gtPath)
        GT = Gray2Float64(load(gtPath))
    elseif compute_ssim == true && !(isfile(gtPath))
        throw(ArgumentError("Ground truth image doesn't exist at given path:  $gtPath"))
    end

   outputStructList = []
   ZMaxDict = Dict()
   RMaxDict = Dict()

    f = basename(folderPath)

    for ksize in ksizeList
        for kernel in kernelList

            println("Running SFF with variables:\nfolder: $f\nkernel: $kernel\nksize: $ksize\n")

            # Run SFF method
            Z, R = sff_main(folderPath, kernel; ksize=ksize, outputFolder=outputFolder, write_maps=write_maps)

            # Compute border mask using ground truth image and use that to mask the computed Z before doing anything else with it 
            borderMask = DetectBorders(GT)
            borderMask = abs.(borderMask .- 1.0)

            ZBorderVal = Z[1,1]
            Z = PaintMask(Z, borderMask; fillValue=ZBorderVal)

            save(outputFolder*"/mask.png", borderMask)

            if outputFolder !== nothing
                WriteMapSingle(Z, R, basename(f), kernel, ksize, outputFolder)
                ZMaxDict[f,"sff",kernel,ksize] = maximum(Z)
                RMaxDict[f,"sff",kernel,ksize] = maximum(R)
                # push!(outputStructList, FileSet(Z,R,0,parse(Int,f),"SFF",kernel))
            end

            # Normalize computed depth map so that it's placed from 0-1
            Z_normalized = imageDisp01(Z)
            Z_normalized = FillBorders(Z_normalized, 1.0)

            if compute_ssim
                # Compute SSIM and MS-SSIM then store all statistical measures in appropriate dictionaries
                # TEMP: Trying to compute MS-SSIM and MS-SSIM with just structural comparison taken into account
                """
                Here, the first parameter is the kernel used to weight the neighbourhood of each pixel while calculating the SSIM locally, and defaults to KernelFactors.gaussian(1.5, 11). The second parameter is the set of weights (α, β, γ) given to the lunimance (L), contrast (C) and structure (S) terms while calculating the SSIM, and defaults to (1.0, 1.0, 1.0). Recall that SSIM is defined as Lᵅ × Cᵝ × Sᵞ.
                Source: https://juliaimages.org/stable/examples/image_quality_and_benchmarks/structural_similarity_index/
                """

                iqi = MSSSIM(KernelFactors.gaussian(1.5,11), (0.0,0.0,1.0))
                structureDict[f,"sff",kernel, ksize] = assess(iqi, GT, Z_normalized)

                # ssimDict[f,"sff",kernel, ksize]  = assess_ssim(GT, Z_normalized)
                msssimDict[f,"sff",kernel, ksize]  = assess_msssim(GT, Z_normalized)
            elseif !compute_ssim
                structureDict[f,"sff",kernel, ksize] = 0
                msssimDict[f,"sff",kernel, ksize] = 0
            end

            rmseDict[f,"sff",kernel,ksize] = RMSE(GT, Z_normalized)      

            # Not computing PSNR for SFF collections, but we create dummy entries for compatibility with `WriteCSV`
            psnrDict[f,"sff",kernel, ksize] = 0
        end
    end

    # ZMax, RMax = FindFileSetMax(outputStructList)

    # if outputFolder !== nothing
    #     for ksize in ksizeList
    #         WriteMaps(outputStructList, outputFolder*string(ksize), nothing, nothing)
    #     end
    # end

    if write_csv == true
        # csvPath = @printf("%s/Ground truth comparison results (%i,%i).csv", outputFolder, ksize[1], ksize[2])
        csvPath = outputFolder * "/Ground truth comparison results.csv"

        WriteCSVSingles(csvPath, [f], ["sff"], kernelList, ksizeList, structureDict, msssimDict, psnrDict, rmseDict, ZMaxDict, RMaxDict)
        # WriteCSV(csvPath, [f], methodList, kernelList, ksizeList, structureDict, msssimDict, psnrDict, ZMax, RMax)
    end

end