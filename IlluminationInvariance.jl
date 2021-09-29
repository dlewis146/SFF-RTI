using Images, LinearAlgebra, Printf

struct LightAngle
    """
    Structure for containing X,Y,Z coordinates of lamps as well as converted θ, ϕ angles.
    """
    x::Float64
    y::Float64
    z::Float64
    theta::Float64
    phi::Float64

end

function FilterImage(img, kernel="laplacian")

    imgX = zero(img)
    imgY = zero(img)

    if kernel == "laplacian"

        imgX = imfilter(img, Kernel.Laplacian((true, false)), "replicate")
        imgY = imfilter(img, Kernel.Laplacian((false, true)), "replicate")

    elseif kernel == "sobel"

        diffX, diffY = Kernel.sobel()

        imgX = imfilter(img, diffX)
        imgY = imfilter(img, diffY)

    end

    return imgX, imgY

end

function ComputeFullVectorGradient(file_list, angle_list, kernel="sobel")

    # Read first image to get image size
    numPixels = length(load(file_list[1]))
    outputSize = size(load(file_list[1]))

    numImages = length(file_list)

    # Initialize empty arrays
    G = zeros(numImages, numImages)
    dPdX = zeros(numPixels, numImages)
    # dPdX = zeros(numImages, numPixels)
    dPdY = zeros(numPixels, numImages)
    # dPdY = zeros(numImages, numPixels)

    for i in range(1, stop=numImages)

        str1 = file_list[i]
        splitVec1 = split(str1, ('_', '.'))

        theta1 = 0.0
        phi1 = 0.0

        # Read angles from respective input object
        angleObj = angle_list[i]
        theta1 = angleObj.theta
        phi1 = angleObj.phi

        # # Parse for theta and phi angles
        # for idx in range(1, stop=length(splitVec1))
        #     if splitVec1[idx] == "Theta"
        #         theta1 = deg2rad(parse(Float64, splitVec1[idx+1]))
        #     end
        #     if splitVec1[idx] == "Phi"
        #         phi1 = deg2rad(parse(Float64, splitVec1[idx+1]))
        #     end
        # end

        for j in range(1, stop=numImages)

            # If inverse coordinates aren't zero, it's therefore been set and,
            # as the inverse, can be copied to the current coordinates
            if G[j, i] != 0
                G[i, j] = G[j, i]
                continue
            end

            str2 = file_list[j]
            splitVec2 = split(str2, ('_', '.'))

            theta2 = 0.0
            phi2 = 0.0

            # Read angles from respective input object
            angleObj = angle_list[j]
            theta2 = angleObj.theta
            phi2 = angleObj.phi

            # Parse for theta and phi angles
            # for idx in range(1, stop=length(splitVec2))
            #     if splitVec2[idx] == "Theta"
            #         theta2 = deg2rad(parse(Float64, splitVec2[idx+1]))
            #     end
            #     if splitVec2[idx] == "Phi"
            #         phi2 = deg2rad(parse(Float64, splitVec2[idx+1]))
            #     end
            # end

            # Compute scalar product of image positions and place into G
            G[i, j] = dot([theta1, phi1], [theta2, phi2])

            # Check if we're on our first run through the images
            if i == 1

                # Compute focus maps
                imgX, imgY = FilterImage(Gray.(load(file_list[j])), kernel) # Make sure to read in image as grayscale
                # imgX, imgY = ComputeSobelImages(Gray.(load(file_list[j]))) # Make sure to read in image as grayscale

                Gray.(imgX)
                Gray.(imgY)
                # readline() // wait for 'enter' in terminal

                # Flatten images with vec() and then copy to appropriate
                # array rows
                dPdX[:,j] = vec(imgX)
                dPdY[:,j] = vec(imgY)
            end
        end
    end

    # G = (transpose(G) .* G)

    results = ComputeNorm(G, dPdX, dPdY, outputSize)

    return results
end

function ComputeNorm(G, dPdX, dPdY, outputSize)

    # Create empty vector to hold respective results for each file/light position
    results = zeros(outputSize)

    # Iterate through numImages (stored as either dimension of G)
    for idx in eachindex(results)
        dPdXRow = vec(dPdX[idx, :])
        dPdYRow = vec(dPdY[idx, :])

        results[idx] = dot(dPdXRow, G, dPdXRow) + dot(dPdYRow, G, dPdYRow)
    end
    return results
end
