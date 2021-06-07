using CSV

function ReadAcqFile(configPath::String)
    """
    Takes in path to .acq file, parses for focus variation sampling information,
    and returns values in tuple in order
        1. Number of positions sampled (Z_Nb)
        2. Bottom-most sampling boundary (Z_1)
        3. Uppper-most sampling boundary (Z_2)

    """

    # Read .acq file and set parameters for magnification compensation
    ZNb = 0
    Z1 = 0.0
    Z2 = 0.0

    fileArray = readdlm(configPath)
    for rowIdx in range(1, stop=size(fileArray)[1])
        if fileArray[rowIdx, 1] == "Z_Nb"
            ZNb = fileArray[rowIdx, 2]
        end
        if fileArray[rowIdx, 1] == "Z_1"
            Z1 = fileArray[rowIdx, 2]
        end
        if fileArray[rowIdx, 1] == "Z_2"
            Z2 = fileArray[rowIdx, 2]
        end
    end



    return (Znb, Z1, Z2)
end
