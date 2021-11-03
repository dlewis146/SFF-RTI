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


function ComputeFullVectorGradient(file_list, angle_list, kernel="sobel")
    """
    Takes in a list of image filepaths as well as a list of LightAngle 
    objects (correlated with images by stored order). Returns a single
    image representing the full vector gradient of the given image space.
    """

    # Read first image to get image size
    numPixels = length(load(file_list[1]))
    outputSize = size(load(file_list[1]))

    numImages = length(file_list)

    # Initialize empty arrays
    G = zeros(numImages, numImages)
    dPdX = zeros(numPixels, numImages)
    dPdY = zeros(numPixels, numImages)

    for i in range(1, stop=numImages)

        theta1 = 0.0
        phi1 = 0.0

        # Read angles from respective input object (Assume degrees)
        angleObj = angle_list[i]
        theta1 = angleObj.theta
        phi1 = angleObj.phi

        for j in range(1, stop=numImages)

            # If inverse coordinates aren't zero, it's therefore been set and,
            # as the inverse, can be copied to the current coordinates
            if G[j, i] != 0
                G[i, j] = G[j, i]
                continue
            end

            theta2 = 0.0
            phi2 = 0.0

            # Read angles from respective input object (Assume degrees)
            angleObj = angle_list[j]
            theta2 = angleObj.theta
            phi2 = angleObj.phi

            # Compute scalar product of image positions and place into G
            G[i, j] = dot([theta1, phi1], [theta2, phi2])

            # Check if we're on our first run through the images
            if i == 1

                # Compute focus maps

                ## Make sure to normalize images with respect to brightness due to the changing of focus distance causing micro changes in aperture f-stop -- As per Subbarao & Choi
                imgX, imgY = FilterImage(brightnessNormalize(Gray.(load(file_list[j]))), kernel) # Make sure to read in image as grayscale
                # imgX, imgY = FilterImage(Gray.(load(file_list[j])), kernel) # Make sure to read in image as grayscale

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
    """
    For use by `ComputeFullVectorGradient`.
    """

    # Create empty vector to hold respective results for each file/light position
    results = zeros(outputSize)

    # Iterate through numImages (stored as either dimension of G)
    for idx in eachindex(results)
        dPdXRow = vec(dPdX[idx, :])
        dPdYRow = vec(dPdY[idx, :])

        # Compute dot product
        results[idx] = dot(dPdXRow, G, dPdXRow) + dot(dPdYRow, G, dPdYRow)
    end
    
    return results
end
