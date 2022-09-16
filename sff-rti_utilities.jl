struct FileSet
    Z::Array{Float64}
    R::Array{Float64}
    numRTI::Int
    numSFF::Int
    method::String
    kernel::String
end

function FixPathEnding(f)  if last(f) !== '/' && last(f) !== '\\'; return f*"/" else return f end end

function FindFileSetMax(structList)
    """
    Takes in list of FileSet structs and returns global maximums for Z and R. 
    Used for determining global normalization values.
    """
    zMaxGlobal = 0.0
    rMaxGlobal = 0.0

    for s in structList
        zMaxLocal = maximum(s.Z)
        rMaxLocal = maximum(s.R)

        if zMaxLocal > zMaxGlobal
            zMaxGlobal = zMaxLocal
        end

        if rMaxLocal > rMaxGlobal
            rMaxGlobal = rMaxLocal
        end
    end

    return zMaxGlobal, rMaxGlobal
end

function WriteMaps(structList, outputFolder, ZMax=nothing, RMax=nothing)
    """
    Takes in list of FileSet structs and an output folder path. Finds global normalization values for 
    depth maps (Z) in structList and then writes them out (along with computed surface normals) in appropriate folders inside of outputBase.
    """
        
    # Get new ZMax and/or RMax if not given
    if     ZMax === nothing && RMax !== nothing
        ZMax, _ = FindFileSetMax(structList)
    elseif ZMax !== nothing && RMax === nothing
        _, RMax = FindFileSetMax(structList)
    elseif ZMax === nothing && RMax === nothing
        ZMax, RMax = FindFileSetMax(structList)
    end

    for s in structList
        folderName = "RTI_"*string(s.numRTI)*"_SFF_"*string(s.numSFF)

        normals = Depth2Normal(s.Z)

        normalsX = shiftNormalsRange(normals[:,:,1])
        normalsY = shiftNormalsRange(normals[:,:,2])
        normalsZ = shiftNormalsRange(normals[:,:,3])
        normalsColor = colorview(RGB, normalsX, normalsY, normalsZ)

        outputFolder = FixPathEnding(outputFolder)

        outputFolderConcatenated = outputFolder*"/"*uppercase(s.method)*" "*uppercase(s.kernel)*"/"
        ispath(outputFolderConcatenated) || mkpath(outputFolderConcatenated)

        save(outputFolderConcatenated*"/Z_"*folderName*"_"*s.method*"_"*s.kernel*".tif", s.Z/ZMax)
        # save(outputFolderConcatenated*"/R_"*folderName*"_"*s.method*"_"*s.kernel*".tif", s.R)
        save(outputFolderConcatenated*"/R_"*folderName*"_"*s.method*"_"*s.kernel*".tif", s.R/RMax)
        save(outputFolderConcatenated*"/Znorm_"*folderName*"_"*s.method*"_"*s.kernel*".tif", imageDisp01(s.Z))
        save(outputFolderConcatenated*"/normals_"*folderName*"_"*s.method*"_"*s.kernel*".png", normalsColor)
    end
end

function WriteCSV(outputPath::String, folderList::Array{String}, methodList::Array{String}, kernelList::Array{String}, ksizeList::Array{Int}, ssimDict::Dict, msssimDict::Dict, psnrDict::Dict, ZMax::Float64, RMax::Float64)

    # Write out comparison results
    println("Writing results to text file...")
    open(outputPath, "w") do io

        # Write header
        write(io, "numRTI,numSFF,method,kernel,ms-ssim (just structure),ms-ssim,average psnr,ksize,Z normalization coefficient\n")

        for f in folderList
            for method in methodList
                for ksize in ksizeList

                    # Parse inputs for proper number of RTI and SFF based on method used
                    numRTI = 0
                    numSFF = 0

                    if lowercase(method) == "sff"
                        numRTI = 0
                        numSFF = parse(Int64, f)
                    elseif lowercase(method) != "sff"
                        numRTI, numSFF = ParseFolderName(f)
                    end

                    for kernel in kernelList
                        # Create and write line to CSV
                        lineOut = string(numRTI, ",", numSFF, ",", method, ",", kernel, ",", ssimDict[f,method,kernel,ksize], ",", msssimDict[f,method,kernel,ksize], ",", psnrDict[f,method,kernel,ksize], ",", ksize, ",", ZMax, ",", RMax, "\n")
                        write(io, lineOut)
                    end
            
                end
            end
        end
    end
end


function ParseFolderName(folderName)
    """
    Parses given string (folderName) to determine number of RTI and SFF positions based on following format:
    `RTI_{NUM RTI}_SFF_{NUM SFF}`
    """

    s = split(folderName, "_")

    numRTI = 0
    numSFF = 0

    for idx in eachindex(s)
        if uppercase(s[idx]) == "RTI"
            numRTI = parse(Int, s[idx+1])
        elseif uppercase(s[idx]) == "SFF"
            numSFF = parse(Int, s[idx+1])
        end
    end
    return numRTI, numSFF
end

function CompositeFromDepthMap(fileList, depthMap)
    """
    Takes in a list of image filepaths (so that RGB may be read) as well as a depth map whose pixel values correspond to the indices
    of the ordered list of images. This is to say that each pixel of the given depth map should have a value corresponding to 
    which image in the given list is considered to be focused for that point. 
    """

    imageList = ReadImageList(fileList, grayscale=false)

    outputImage = zeros(size(imageList[1],1), size(imageList[1],2), 3)

    for r in range(1, stop=size(depthMap, 1))
        for c in range(1, stop=size(depthMap, 2))

            z = Integer(round(depthMap[r,c]))

            outputImage[r,c,:] = imageList[z][r,c,:]

        end
    end

    return outputImage
end