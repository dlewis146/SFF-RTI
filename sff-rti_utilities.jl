struct FileSet
    Z::Array{Float64}
    R::Array{Float64}
    folderName::String
    method::String
    kernel::String
end

function FindFileSetMax(structList)
    """
    Takes in list of FileSet structs and returns global maximums for Z and R in that order. 
    Used for determining global normalization value.
    """

    zMaxGlobal = 0
    rMaxGlobal = 0

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
    Z and R maps in structList and then writes them out (along with computed surface normals) in appropriate folders inside of outputBase.
    """
        
    # Get new ZMax and/or RMax if not given
    if ZMax === nothing || RMax === nothing
        ZMaxNew, RMaxNew = FindFileSetMax(outputStructList)

        if ZMax === nothing; ZMax = ZMaxNew end
        if RMax === nothing; RMax = RMaxNew end
    end

    for s in structList
        # outputFolder = outputBase*s.folderName

        normals = zeros(size(s.Z, 1), size(s.Z, 2), 3)
        if s.method == "SFF"
            normals = Depth2Normal(complement.(s.Z))
        else
            normals = Depth2Normal(s.Z)
        end
        # normals = Depth2Normal(s.Z)
        normalsX = shiftNormalsRange(normals[:,:,1])
        normalsY = shiftNormalsRange(normals[:,:,2])
        normalsZ = shiftNormalsRange(normals[:,:,3])
        normalsColor = colorview(RGB, normalsX, normalsY, normalsZ)

        outputFolderConcatenated = outputFolder*"/"*uppercase(s.method)*" "*uppercase(s.kernel)*"/"
        ispath(outputFolderConcatenated) || mkpath(outputFolderConcatenated)

        save(outputFolderConcatenated*"/Z_"*s.folderName*"_"*s.method*"_"*s.kernel*".png", s.Z/ZMax)
        save(outputFolderConcatenated*"/Znorm_"*s.folderName*"_"*s.method*"_"*s.kernel*".png", imageDisp01(s.Z))
        save(outputFolderConcatenated*"/R_"*s.folderName*"_"*s.method*"_"*s.kernel*".png", s.R/RMax)
        save(outputFolderConcatenated*"/normals_"*s.folderName*"_"*s.method*"_"*s.kernel*".png", normalsColor)
    end
end

function WriteCSV(outputPath, folderList, methodList, kernelList, rmseList, irmseList, ZMax, RMax)

    # Write out comparison results
    println("Writing results to text file...")
    open(outputPath, "w") do io

        # Write header
        write(io, "numRTI,numSFF,method,kernel,rmse,inverse rmse,Z normalization coefficient,R normalization coefficient\n")

        for f in folderList
            for method in methodList

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
                    lineOut = string(numRTI, ",", numSFF, ",", method, ",", kernel, ",", rmseList[f,method,kernel], ",", irmseList[f,method,kernel], ",", ZMax, ",", RMax, "\n")
                    write(io, lineOut)
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
            numRTI = s[idx+1]
        elseif uppercase(s[idx]) == "SFF"
            numSFF = s[idx+1]
        end
    end
    return numRTI, numSFF
end

