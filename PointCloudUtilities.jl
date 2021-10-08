struct PointCloudData
    """
    Structure designed to hold all information from a parsed point cloud in separate variables for easy access.
    """
    rgbImage::Array{Float64}
    depthMap::Array{Float64}
    normalMap::Array{Float64}
end


function Depth2PointCloud(img, rgbMap=None, filepath="./point_cloud.ply ")
    """
    Takes a depth image and writes a simple XYZ point cloud to given filepath
    """

    # Find number of NaN in depth map
    numNaN = count(isnan.(img))

    lines = String[]

    push!(lines, "ply")
    push!(lines, "format ascii 1.0")
    push!(lines, string("element vertex ", length(img)-numNaN))
    push!(lines, "property float x")
    push!(lines, "property float y")
    push!(lines, "property float z")
    push!(lines, "property uint8 red")
    push!(lines, "property uint8 green")
    push!(lines, "property uint8 blue")
    push!(lines, "end_header")


    # Iterate across all points in image and push to pcd file
    for y in range(1, stop=size(img, 1))
        for x in range(1, stop=size(img, 2))

            # If NaN, skip point.
            if isnan(img[y,x])
                continue
            end

            z = img[y,x]
            # z = img[y,x] * 100
            red = Int(round(rgbMap[y,x,1]*255))
            green = Int(round(rgbMap[y,x,2]*255))
            blue = Int(round(rgbMap[y,x,3]*255))


            push!(lines, string(x, " ", y, " ", z, " ", red, " ", green, " ", blue))
        end
    end

    file = open(filepath, "w")
    for line in lines
        write(file, line)
        write(file, "\n")
    end

    close(file)
end


function Normals2PointCloud(img, rgpMap=None, filepath="./point_cloud.ply")

    """
    Takes a normal map and writes a simple XYZ point cloud to given filepath
    """

    lines = String[]

    push!(lines, "ply")
    push!(lines, "format ascii 1.0")
    push!(lines, string("element vertex ", length(img)-numNaN))
    push!(lines, "property float nx")
    push!(lines, "property float ny")
    push!(lines, "property float nz")
    push!(lines, "property uint8 red")
    push!(lines, "property uint8 green")
    push!(lines, "property uint8 blue")
    push!(lines, "end_header")


    # Iterate across all points in image and push to pcd file
    for y in range(1, stop=size(img, 1))
        for x in range(1, stop=size(img, 2))

            nx = img[y,x,1]
            ny = img[y,x,2]
            nz = img[y,x,3]
            red = Int(round(rgbMap[y,x,1]*255))
            green = Int(round(rgbMap[y,x,2]*255))
            blue = Int(round(rgbMap[y,x,3]*255))


            push!(lines, string(x, " ", y, " ", z, " ", red, " ", green, " ", blue))
        end
    end

    file = open(filepath, "w")
    for line in lines
        write(file, line)
        write(file, "\n")
    end

    close(file)
end


end


function ReadPointCloud(filepath)
    """
    Takes in filepath to ASCII XYZ point cloud and reads it into a structure designed to hold all modalities separately.
    """

    X = Float64[]
    Y = Float64[]
    Z = Float64[]

    Nx = Float64[]
    Ny = Float64[]
    Nz = Float64[]

    R = Float64[]
    G = Float64[]
    B = Float64[]

    ## NOTE: Assumes format for each line as follows:
    ## X Y Z normX normY normZ R G B

    # Open file for reading
    file = open(filepath)

    # Read through file until end of file is reached, pushing appropriate values into respective arrays    
    while ! eof(file)
        line = readline(file)

        try
            # Check to see if current line starts with a number. If not, we can assume it doesn't have image information to push into an array.
            if isa(parse(Float64, split(line, " ")[1]), Number) 
                charList = split(line, " ")

                push!(X,  parse(Float64, charList[1]))
                push!(Y,  parse(Float64, charList[2]))
                push!(Z,  parse(Float64, charList[3]))
                push!(Nx, parse(Float64, charList[4]))
                push!(Ny, parse(Float64, charList[5]))
                push!(Nz, parse(Float64, charList[6]))
                push!(R,  parse(Float64, charList[7]))
                push!(G,  parse(Float64, charList[8]))
                push!(B,  parse(Float64, charList[9]))
            end
        catch UndefVarError # Catch error if can't parse as Float64 AKA not a number
            continue
        end
    end

    close(file)

    # TODO: Change so that width and height only span object dimensions. Currently goes from image origin to far end of object.

    width = Int(maximum(X))
    height = Int(maximum(Y))

    # Create shaped arrays to write image information into
    rgbImage = zeros((height, width, 3))
    normalMap = zeros((height, width, 3))
    depthMap = zeros((height, width))

    # Create RGB image, normal map, and depth map
    for pt in eachindex(X)

        # NOTE: Cast to Int for indexing
        ptX = Int(X[pt])
        ptY = Int(Y[pt])

        # RGB image construction
        rgbVector = [ R[pt], G[pt], B[pt] ]
        rgbImage[ptY, ptX, :] = rgbVector

        # Normal map construction
        normalVector = [ Nx[pt], Ny[pt], Nz[pt] ]
        normalMap[ptY, ptX, :] = normalVector

        # Depth map construction
        depthMap[ptY, ptX] = Z[pt]
    end

    # Store images in data structure and return
    return PointCloudData(rgbImage, depthMap, normalMap)
end