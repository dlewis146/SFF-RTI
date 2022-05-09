function sff(imageList, focusList, sampleStep=2, median=true)
    """
    Shape from Focus algorithm originally programmed in Matlab by Said Pertuz

    https://www.mathworks.com/matlabcentral/fileexchange/55103-shape-from-focus
    """

    M = size(imageList[1], 1)
    N = size(imageList[1], 2)
    P = length(imageList)

    imageStack = zeros(M,N,P)

    # Create focus stack of incoming images
    for idx in eachindex(imageList)

        # Get image from given list
        img = imageList[idx]

        # Place into stack
        imageStack[:,:,idx] = img
    end

    # Estimate depth map
    YMax, zi, InterpolatedFocusMeasures = GaussianInterpolation3Pt(focusList, imageStack, sampleStep)

    z = zi

    # NOTE:Casting YMax to Int for use as indices
    Ic = round.(Int, YMax)

    # Create and populate array to hold maximum focus value for each pixel
    fmax = zeros(size(Ic))
    fmax .= focusList[Ic]

    # z[isnan.(z)] .= maximum(focusList)
    z[isnan.(z)] = fmax[isnan.(z)]

    # Shift values beyond given focus limits to nearest boundary
    z[z.>maximum(focusList)] .= maximum(focusList)
    z[z.<minimum(focusList)] .= minimum(focusList)

    ## Median filter
    if median == true
        println("Running Z through median filter with size (3,3)...")
        z = mapwindow(median!, z, (3,3))
    end

    return z

end


function GaussianInterpolation3Pt(x, imageStack, Gstep=2)

    M,N,P = size(imageStack)

    # Find maximum value along Z dimension for each pixel
    coordsMax = argmax(imageStack, dims=3)
    
    idxMax = zeros(Int8, M, N)
    for idx in eachindex(coordsMax) idxMax[idx] = Int(coordsMax[idx][3]) end

    # Push Z values away from edges of focus stack
    idxMax[idxMax.<=Gstep] .= Gstep+1
    idxMax[idxMax.>=P-Gstep] .= P-Gstep

    # Make sure depth values are within boundaries of 3D interpolation window,
    # AKA make sure that if we Gstep from a max focus point to do interpolation,
    # make sure we're not going to try and access a point outside of our scope
    idxMaxPadded = padarray(idxMax, Pad(:replicate,Gstep+1, Gstep+1))

    interpolatedD = zeros(Float64, M, N)
    # interpolatedF = zeros(Float64, M, N)
    
    # For each pixel, interpolate Gaussian over three images in determined peak region
    for i in range(1, stop=M)
        for j in range(1, stop=N)

            z = idxMaxPadded[i,j]

            yLow  = imageStack[i,j,z-Gstep]
            yMid  = imageStack[i,j,z]
            yHigh = imageStack[i,j,z+Gstep]

            xLow  = x[z-Gstep]
            xMid  = x[z]
            xHigh = x[z+Gstep]

            # Compute Gaussian distribution parameters
            dBar = (( (log(yMid) - log(yHigh)) * (xMid^2 - xLow^2) ) / (2*(xMid-xLow) * ((log(yMid) - log(yLow)) + (log(yMid) - log(yHigh))) )) - (((log(yMid) - log(yLow)) * (xMid^2 - xHigh^2)) / (2*(xMid-xLow) * ((log(yMid) - log(yLow)) + (log(yMid) - log(yHigh)))))

            # sigmaF = -(((xMid^2 - xLow^2) + (xMid^2 - xHigh^2) ) / (2 * ( (log(yMid) - log(yLow)) + (log(yMid) - log(yHigh)))))

            # Interpolate focus value and place focus and depth values in respective maps
            # interpolatedF[i,j] = ( yMid / ( exp( (-1/2) * (((xMid-dBar)/sigmaF)^2) ) ) )
            interpolatedD[i,j] = dBar
        end
    end

    return idxMax, interpolatedD, interpolatedF
end


function Depth2Normal(img)
    """
    Takes in (m x n) image (depth map) and uses
    it to compute the surface normal at each point.

    Returns 3D array of normals where the channels are:
        1st channel - X-component of normal
        2nd channel - Y-component of normal
        3rd channel - Z-component of normal
    """

    # img = Gray2Float64(img)

    normalsOut = zeros(size(img, 1)+2, size(img, 2)+2, 3)

    for j in range(2, stop=size(img, 1)-1)
        for i in range(2, stop=size(img, 2)-1)

            dzdx = (Float64(img[j, i+1]) - Float64(img[j, i-1])) / 2.0
            dzdy = (Float64(img[j+1, i]) - Float64(img[j-1, i])) / 2.0
            # d = [-dzdx, -dzdy, Float64(1.0/255)]
            d = [-dzdx, -dzdy, 1.0/255]

            n = normalize(d)

            normalsOut[j,i,:] .= n

        end
    end

    return normalsOut[2:size(normalsOut, 1)-1, 2:size(normalsOut, 2)-1, :]

end
