using Distances, Printf, ProgressMeter

function FindMinimumVertices(filepath)

    f = open(filepath, "r")
    fileSize = countlines(f)
    close(f)

    line = 0  

    minimumDistanceVector = (9999.0, 9999.0, 9999.0)
    minimumDistance = 9999.0
 
    f = open(filepath, "r")

    verticesDoneFlag = true

    p = Progress(fileSize, "Reading file...")
    # read till end of file
    while ! eof(f) || verticesDoneFlag==false
    
        # read a new / next line for every iteration          
        s = readline(f)         
        line += 1
    
        string_parts = split(s)
    
        if length(string_parts) == 0; continue end
        if string_parts[1] == "f"; verticesDoneFlag = true end

        if string_parts[1] == "v"
            x = parse(Float64, string_parts[2])
            y = parse(Float64, string_parts[3])
            z = parse(Float64 ,string_parts[4])
            currentVec = Vector([x,y,z])
            currentDistance = euclidean(currentVec, Vector([0,0,0]))

            if currentDistance < minimumDistance
                minimumDistance = currentDistance
                minimumDistanceVector = currentVec
            end
        end

        ProgressMeter.next!(p)
    end
    
    close(f)

    @printf("\nxMin: %f\nyMin: %f\nzMin: %f", minimumDistanceVector[1], minimumDistanceVector[2], minimumDistanceVector[3])
    @printf("\nMinimum distance: %f", minimumDistance)

    return minimumDistanceVector

end

function CenterMesh(filepath)

    minVector = FindMinimumVertices(filepath)

    # Get filepath without extension 
    filepathOut = filepath[1:lastindex(filepath)-4] * " (Centered).obj"

    f = open(filepath, "r")
    fileSize = countlines(f)
    close(f)

    f = open(filepath)
    # fOut = open(filepathOut)

    linesOut = []

    verticesDoneFlag = true

    p = Progress(fileSize, "Centering mesh...")
    # read till end of file
    while ! eof(f) || verticesDoneFlag==false
    
        # read a new / next line for every iteration          
        s = readline(f)         
    
        string_parts = split(s)
    
        if length(string_parts) == 0; continue end
        if string_parts[1] == "f"; verticesDoneFlag = true end

        if string_parts[1] == "v"
            x = parse(Float64, string_parts[2])
            y = parse(Float64, string_parts[3])
            z = parse(Float64 ,string_parts[4])

            x2 = x - minVector[1]
            y2 = y - minVector[2]
            z2 = z - minVector[3]

            string_parts[2] = string(x2)
            string_parts[3] = string(y2)
            string_parts[4] = string(z2)

            sNew = "" * string_parts[1]
            for x in range(2, stop=length(string_parts))
                part = string_parts[x]
                sNew = sNew * " " * part
            end

            push!(linesOut, sNew)

        else
            push!(linesOut, s)
        end

        ProgressMeter.next!(p)
    end
    
    close(f)
    # close(fOut)

    file = open(filepathOut, "w")
    for line in linesOut
        write(file, line)
        write(file, "\n")
    end

    close(file)


end

filepath = "J:/Research/3D Models/Wilanow Data/Gemms/Gemma Wil.3083 2 II 2 167 - Ball Pivot.ply"

CenterMesh(filepath)
