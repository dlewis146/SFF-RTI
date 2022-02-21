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

    numImages = length(file_list)

    # Initialize empty arrays
    G = zeros(numImages, numImages)
    dPdX = zeros(numPixels, numImages)
    dPdY = zeros(numPixels, numImages)

    for i in range(1, stop=numImages)

        # Read Cartesian coordinates from respective input object (Assume degrees)
        angleObj = angle_list[i]
        x1 = round(angleObj.x; digits=5)
        y1 = round(angleObj.y; digits=5)
        z1 = round(angleObj.z; digits=5)

        for j in range(1, stop=numImages)

            # If inverse coordinates aren't zero, it's therefore been set and,
            # as the inverse, can be copied to the current coordinates
            if G[j, i] != 0
                G[i, j] = G[j, i]
                continue
            end

            # Read Cartesian coordinates from respective input object (Assume degrees)
            angleObj = angle_list[j]
            x2 = round(angleObj.x; digits=5)
            y2 = round(angleObj.y; digits=5)
            z2 = round(angleObj.z; digits=5)

            # Compute cosine similarity between angles
            G[i,j] = 1 - acos( round((x1*x2 + y1*y2 + z1*z2) / (sqrt(x1^2 + y1^2 + z1^2) * sqrt(x2^2 + y2^2 + z2^2)); digits=5) )

            # Check if we're on our first run through the images
            if i == 1

                # Compute focus maps for image in `file_list` at `j`
                img = Gray.(load(file_list[j]))

                # Normalize images with respect to brightness due to the changing of focus distance causing micro changes in aperture f-stop -- As per Subbarao & Choi
                # img = brightnessNormalize(img)

                imgX, imgY = FilterImageSeparate(img, kernel)

                # Flatten images with vec() and then copy to appropriate
                # array rows
                dPdX[:,j] = vec(imgX)
                dPdY[:,j] = vec(imgY)
            end
        end
    end

    ### Compute FVG

    # Create empty vector to hold respective results for each file/light position
    results = zeros(size(load(file_list[1])))

    # Iterate through numImages (stored as either dimension of G)
    for idx in eachindex(results)
        dPdXRow = vec(dPdX[idx, :])
        dPdYRow = vec(dPdY[idx, :])

        # Compute and sum dot products
        results[idx] = dot(dPdXRow, G, dPdXRow) + dot(dPdYRow, G, dPdYRow)
    end

    return results
end