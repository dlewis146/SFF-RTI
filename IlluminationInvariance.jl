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

        # Read Cartesian coordinates from respective input object (assume degrees) and place in vector, rounding to 5 digits
        vec1 = vec([round(angle_list[i].x; digits=5),
                    round(angle_list[i].y; digits=5),
                    round(angle_list[i].z; digits=5)
                    ])

        for j in range(1, stop=numImages)

            # If inverse coordinates aren't zero, it's therefore been set and,
            # as the inverse, can be copied to the current coordinates
            if G[j, i] != 0
                G[i, j] = G[j, i]
                continue
            end

            # Read Cartesian coordinates from respective input object (assume degrees) and place in vector, rounding to five digits
            vec2 = vec([round(angle_list[j].x; digits=5),
                        round(angle_list[j].y; digits=5),
                        round(angle_list[j].z; digits=5)
                        ])

            # Compute cosine similarity between angles
            G[i,j] = round(dot(vec1,vec2) / (norm(vec1) * norm(vec2)); digits=5)

            # Check if we're on our first run through the images
            if i == 1

                # Compute focus maps for image in `file_list` at `j`
                img = Gray.(load(file_list[j]))

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