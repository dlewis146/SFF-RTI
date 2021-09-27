using Images, Glob, CSV, DataFrames, GR

struct LP_Struct
    """
    Structure for containing information held in an LP file
    """
    X::Array{Float64}
    Y::Array{Float64}
    Z::Array{Float64}
end


function CreateFlatImageStack(fileList)
    """
    Takes in an array of file paths (as strings) and reads in all images to store as flattened columns. Used for pre-processing image data for computation of image normals
    """

    imageLength = size(load(fileList[1]), 1) * size(load(fileList[1]), 2)

    imageStack = zeros(length(fileList), imageLength)

    for idx in eachindex(fileList)
        img = Gray.(load(fileList[idx]))
        imageStack[idx, :] = vec(img)
    end
    return imageStack
end


function CreateLP(csvPath)
    """
    Reads in CSV from given path and parses for LP information.
    """

    lampData = CSV.File(csvPath; select=["x_lamp", "y_lamp", "z_lamp"])

    X = []
    Y = []
    Z = []

    for row in lampData
        push!(X, row.x_lamp)
        push!(Y, row.y_lamp)
        push!(Z, row.z_lamp)
    end

    LP = LP_Struct(X, Y, Z)
    return LP
end


function normalizeLostDynamic(normalMaps)
    """
    Translated from Matlab code written by Marvin Nurit.

    Takes in 3 x (m*n) array of normals computed by ComputeNormals and normalizes them to expected scales for visualization.
    """

    norm = sqrt.(sum(normalMaps.^2, dims=1))
    norm = GR.jlgr.repmat(norm, size(normalMaps, 1), 1)

    normalized = normalMaps ./ norm
    return normalized
end


function ComputeNormals(imageStack, LP)
    """
    Translated from Matlab code written by Marvin Nurit.

    Takes in 3 x (m*n) array of flattened images and LP struct to compute surface normals.
    """

    positions = [LP.X LP.Y LP.Z]
    positionsInverse = (transpose(positions) * positions) \ transpose(positions)
    
    normalMaps = positionsInverse * imageStack
    
    normalMaps = normalizeLostDynamic(normalMaps)

    return normalMaps
end


function NormalsPipeline(folderPath, csvPath)

    fileList = glob("*.png", folderPath)
    imageStack = CreateFlatImageStack(fileList)    

    LP = CreateLP(csvPath)

    normalMaps = ComputeNormals(imageStack, LP)

    imgSize = size(load(fileList[1]))

    normalX = reshape(normalMaps[1,:], imgSize)
    normalY = reshape(normalMaps[2,:], imgSize)
    normalZ = reshape(normalMaps[3,:], imgSize)

    return colorizeNormals(imageDisp01(normalX), imageDisp01(normalY), imageDisp01(normalZ))

    # return normalX, normalY, normalZ
end

function imageDisp01(img) img.+abs(minimum(img)) end

function colorizeNormals(normalX, normalY, normalZ) colorview(RGB, normalX, normalY, normalZ) end