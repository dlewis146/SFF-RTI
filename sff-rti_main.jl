using Images
include("./sff_rti.jl")
include("./sff-rti_utilities.jl")

function SFFRTIHandler(baseFolder, innerFolderList, methodList, kernelList; ksize=(3,3), write_maps=false, write_csv=false, compute_psnr=false, compute_ssim=false, outputFolder="", gtPath="", ZMax=NaN)

    # Read in ground truth for SSIM comparison if needed
    GT = nothing

    if compute_ssim == true && isfile(gtPath)
        GT = Gray2Float64(load(gtPath))
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

                # Run SFF-RTI method
                println()
                Z, psnrMean = sff_rti(baseFolder*f, method, kernel; ksize=ksize, outputFolder=outputFolder, write_maps=write_maps, compute_psnr=compute_psnr)

                numRTI, numSFF = ParseFolderName(f)

                if outputFolder !== nothing
                    push!(outputStructList, FileSet(Z,numRTI,numSFF,method,kernel))
                end

                # Normalize computed depth map so that it's placed from 0-1
                Z_normalized = imageDisp01(Z)

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

                psnrDict[f,method,kernel] = psnrMean

            end
        end
    end

    if isnan(ZMax) 
        ZMax = FindFileSetMax(outputStructList)
    else
        println("NOTE: Using given ZMax normalization coefficient")
    end

    if write_maps == true
        WriteMaps(outputStructList, outputFolder, ZMax)
    end

    if write_csv == true
        # csvPath = @printf("%s/Ground truth comparison results (%i,%i).csv", outputFolder, ksize[1], ksize[2])
        csvPath = outputFolder * "/Ground truth comparison results (" * string(ksize[1]) * "," * string(ksize[2]) * ").csv"

        WriteCSV(csvPath, innerFolderList, methodList, kernelList, structureDict, msssimDict, ZMax, psnrDict)
    end

end


function FixPathEnding(f)  if last(f) !== '/' && last(f) !== '\\'; return f*"/" else return f end end

### Begin test harness
function sffrti_main()
    """
    Used for dispatching SFF-RTI main function with different window sizes as well as handling input parameters in a cleaner fashion
    """

    base = "J:/Research/Generated Data/Blender/Statue du parc d'Austerlitz/SFFRTI/Restrained Z/"

    innerFolderList = ["RTI_4_SFF_5", "RTI_4_SFF_20", "RTI_4_SFF_100", "RTI_20_SFF_5", "RTI_20_SFF_20"]
    methodList = ["fvg", "mean", "max"]
    kernelList = ["sml"]

    outputFolder = "J:/Image Out/Restrained Z/SSIM Tests/"

    # Check file path endings to make sure they end in path separators
    base = FixPathEnding(base)
    outputFolder = FixPathEnding(outputFolder)

    # Store path of ground truth depth map for comparison
    gtPath = "J:/Research/Generated Data/Blender/Statue du parc d'Austerlitz/Ground Truth Restrained Z/Depth/Image0001.png"

    ksize_list = [3, 5, 7, 9, 11]

    for ksize in ksize_list
        SFFRTIHandler(base, innerFolderList, methodList, kernelList; ksize=(ksize,ksize), write_maps=true, write_csv=true, compute_psnr=true, compute_ssim=true, outputFolder=outputFolder*string(ksize)*"/", gtPath=gtPath)
    end
end