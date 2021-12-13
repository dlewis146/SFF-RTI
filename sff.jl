function sffSimple(imageList, focusList)
    """
    Takes in list of focus maps and list of focal distances and computes 
    depth map by finding maximum image value for each pixel.

        imageList::Array{Array{Float64}} - List of computed focus maps for
                                         input images, ex) image gradients
        
        focusList::Array{Float64}      - List of focus distances for each 
                                         image. Correlated to imageList by 
                                         stored order.
    """

    # Get image dimensions
    M = size(imageList[1], 1)
    N = size(imageList[1], 2)
    P = length(imageList)

    # Create focus stack of incoming images

    imageStack = zeros(size(imageList[1])[1], size(imageList[1])[2], length(imageList))

    for idx in eachindex(imageList)

        # Get image from given list
        img = imageList[idx]

        # Place into stack
        imageStack[:,:,idx] = img
    end


    # Create empty depth map (background being NaN for parsing later)
    depthMap = fill(NaN, (M,N))

    for x in range(1, stop=size(imageStack, 1))
        for y in range(1, stop=size(imageStack, 2))

            v = imageStack[x,y,:]

            # If all points are the same, assume background and keep as NaN. Check this by comparing min and max of imageStack vector
            if (maximum(v) == minimum(v))
                continue
            end

            # Find image index where estimated focus is greatest
            idx = argmax(v)

            # Save focus value as pixel value in depth map
            depthMap[x,y] = focusList[idx]

        end
    end

    return depthMap
end


function sff(imageList, focusList)
    """
    Shape from Focus algorithm originally programmed in Matlab by Said Pertuz   

    https://www.mathworks.com/matlabcentral/fileexchange/55103-shape-from-focus
    """


    M = size(imageList[1], 1)
    N = size(imageList[1], 2)
    P = length(imageList)

    imageStack = zeros(size(imageList[1])[1], size(imageList[1])[2], length(imageList))

    # Create focus stack of incoming images
    for idx in eachindex(imageList)

        # Get image from given list
        img = imageList[idx]

        # Place into stack
        imageStack[:,:,idx] = img
    end

    # Estimate depth map
   YMax, zi, s, A = gauss3P(focusList, imageStack)

    z = zi;

    # Shift values beyond given focus limits to nearest boundary
    z[z.>maximum(focusList)] .= maximum(focusList)
    z[z.<minimum(focusList)] .= minimum(focusList)
    for idx in eachindex(z)
        if z[idx] > maximum(focusList)
            z[idx] = maximum(focusList)
        elseif z[idx] < minimum(focusList)
            z[idx] = minimum(focusList)
        end
    end

    # NOTE:Casting YMax to Int
    Ic = round.(Int, YMax)

    # Create and populate array to hold maximum focus value for each pixel
    # ORIG
    fmax = ones(size(Ic))
    fmax .= focusList[Ic]
    ## fmax .= focusList[YMax]

    z[isnan.(z)] = fmax[isnan.(z)]

    ## Median filter
    # @printf("\nRunning Z through median filter with size (3,3)...")
    # z = mapwindow(median, z, (3,3))

    ## Reliability measure

    # R = RMeasure(imageStack, focusList, zi, fmax)

    errorArray = zeros(M, N)
    
    for k in range(1, stop=size(focusList,1))
        errorArray = errorArray + abs.(imageStack[:,:,k] - A.*exp.(-(focusList[k] .- zi).^2 ./ (2*s.^2)))
    end

    errorArray = errorArray./(fmax.*size(focusList, 1))

    # @printf("\nRunning errorArray through averaging filter with size (3,3)...")
    # averagingKernel = fill(1/9, (3,3))
    # errorArray = imfilter(errorArray, averagingKernel)

    R = 20*log10.(P*fmax./errorArray)
    mask = isnan.(zi)
    
    R[mask] .= 0
    R[R.<0] .= 0
    R[isnan.(R)] .= 0

    return z, R

end


function GaussianInterpolation3Pt(x, imageStack, step=2)

    M,N,P = size(imageStack)

    # Find maximum value along Z dimension for each pixel
    idxMax = argmax(imageStack, dims=3)
    ZMax = imageStack[idxMax]

    # Make sure depth values are within boundaries of 3D interpolation window, 
    # AKA make sure that if we step from a max focus point to do interpolation, 
    # make sure we're not going to try and access a point outside of our scope
    ZMaxPadded = zero(ZMax)
    ZMaxPadded[ZMax.<=step] = step+1
    ZMaxPadded[ZMax.>=P-step] = p-step

    # NOTE:Casting ZMaxPadded to Int for use as indexing into imageStack
    ZMaxPadded = round.(Int, ZMaxPadded)

    interpolatedOut = zero(ZMax)
    errorOut = zero(ZMax)

    # For each pixel, interpolate Guassian over three images in determined peak region
    for i in range(1, stop=M)
        for j in range(1, stop=N) 

            z = ZMaxPadded[i,j]

            # Gather focus values and determined depth values
            yLow  = imageStack[i,j,z-step]
            yMid  = imageStack[i,j,z]
            yHigh = imageStack[i,j,z+step]

            xLow  = x[z-step]
            xMid  = x[z]
            xHigh = x[z+step]

            # Compute Gaussian distribution parameters
            dBar = ( (ln(yMid) - ln(yHigh)) * (xMid^2 - xLow^2) ) / ( 2*(xMid-xLow) * ( (ln(yMid) - ln(yLow)) + (ln(yMid) - ln(yHigh)) )  )
            dBar = dBar - ( ( (ln(yMid) - ln(yLow)) * (xMid^2 - xHigh^2) ) / ( 2*(xMid-xLow) * ( (ln(yMid) - ln(yLow)) + (ln(yMid) - ln(yHigh)) )  ) )

            sigmaF = - ( (xMid^2 - xLow^2) + (xMid^2 xHigh^2) ) / ( 2* ( (ln(yMid) - ln(yLow)) + (ln(yMid) - ln(yHigh)) ) ) 

            # Interpolate depth value and place in output depth map
            interpolatedOut[i,j] = ( yMid / ( e( (-1/2) * ((xMid-dBar)/sigmaF)^2 ) ) )

        end
    end

    return interpolatedOut
end


function gauss3P(x, Y)
    """
    For use by `sff`

    Interpolates focus measures to create a smoother depth map using a Guassian distribution sampled at three points in the peak region of the focus measure
    """

    STEP = 2

    M,N,P = size(Y)

    # Find maximum value INDEX for each pixel along 3rd dimension
    # YMax = maximum(Y,dims = 3)
    YMaxHold = argmax(Y, dims=3)

    YMax = zero(Y[:,:,1])

    for idx in eachindex(YMaxHold)
    # for px in YMaxHold
        YMax[idx] = YMaxHold[idx][3]
    end

    # Get vector of flattened YMax
    # Ic = YMax[:]

    #NOTE: Rounding YMax to Int because 2.0000044609612284 !< 2
    YMax = round.(Int, YMax)

    Ic = zero(YMax)
    # For all values of Ic that are <= STEP, set to STEP+1
    # Ic[Ic .<= STEP] .= (STEP+1)
    # For all values of Ic that are >= P-STEP, set to P-STEP
    # Ic[Ic .>= (P-STEP)] .= (P-STEP)
    for idx in eachindex(YMax)
        I_val = YMax[idx]

        if I_val <= STEP
            Ic[idx] = STEP+1
        elseif I_val >= P-STEP
            Ic[idx] = P-STEP
        else
            Ic[idx] = I_val
        end
    end


    # Get linear indices for Y
    YLinear = LinearIndices(Y)

    # Get cartesian indices for a single input image
    YCartesian = CartesianIndices(YLinear[:,:,1])

    # NOTE:Casting Ic to Int for indexing
    Ic = round.(Int, Ic)

    Index1 = zero(Ic)
    Index2 = zero(Ic)
    Index3 = zero(Ic)

    for idx in eachindex(Ic)
        cartesianPlanarCoords = YCartesian[idx]
        c = cartesianPlanarCoords[1]
        r = cartesianPlanarCoords[2]

        Ic_val = Ic[idx]

        Index1[idx] = YLinear[c, r, Ic_val-STEP]
        Index2[idx] = YLinear[c, r, Ic_val]
        Index3[idx] = YLinear[c, r, Ic_val+STEP]
    end

    x1 = zeros(M,N)
    x2 = zeros(M,N)
    x3 = zeros(M,N)

    for idx in eachindex(Ic)
        cartesianPlanarCoords = YCartesian[idx]
        c = cartesianPlanarCoords[1]
        r = cartesianPlanarCoords[2]

        Ic_val = Ic[idx]

        x1[idx] = x[Ic_val - STEP]
        x2[idx] = x[Ic_val]
        x3[idx] = x[Ic_val + STEP]

    end

    y1 = zeros(M,N)
    y2 = zeros(M,N)
    y3 = zeros(M,N)

    for idx in eachindex(y1)
        y1[idx] = log(Y[Index1[idx]])
        y2[idx] = log(Y[Index2[idx]])
        y3[idx] = log(Y[Index2[idx]])
    end

    c = ( (y1-y2).*(x2-x3)-(y2-y3).*(x1-x2) ) ./ ( (x1.^2-x2.^2).*(x2-x3)-(x2.^2-x3.^2).*(x1-x2) )
    b = ( (y2-y3)-c.*(x2-x3).*(x2+x3) )./(x2-x3)
    a = y1 - b.*x1 - c.*x1.^2

    # s = zero(c)
    s = zeros(ComplexF64, size(c))
    for idx in eachindex(c)
        # println(idx, " @ sqrt(-1/(2*", c[idx], ")")
        if isnan(c[idx])
            s[idx] = sqrt( -1 / (2*c[idx]))
        else
            s[idx] = sqrt(Complex(-1 / (2*c[idx])))
        end
    end
    # s[idx] = sqrt(-1 / (2*c[idx]))
    # println()
    # println(Y[459098])
    # println()
    # s = sqrt.(-1 ./ (2*c))
    # s = sqrt(-1./(2*c))
    u = b.*float32.(real.(s)).^2
    # u = b.*s.^2
    A = exp.(a + (u.^2)./(2*float32.(real.(s).^2)))
    # A = exp(a + u.^2./(2*s.^2))

    return YMax, u, s, A

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

    img = Gray2Float64(img)

    normalsOut = zeros(size(img, 1)+2, size(img, 2)+2, 3)

    for j in range(2, stop=size(img, 1)-1)
        for i in range(2, stop=size(img, 2)-1)

            dzdx = (Float64(img[j, i+1]) - Float64(img[j, i-1])) / 2.0
            dzdy = (Float64(img[j+1, i]) - Float64(img[j-1, i])) / 2.0
            d = [-dzdx, -dzdy, 1.0]

            # t = [ i  , j-1, Float64(img[j-1, i  ]) ]
            # f = [ i-1, j  , Float64(img[j  , i-1]) ]
            # c = [ i  , j  , Float64(img[j  , i  ]) ]
            # d = cross( (f-c), (t-c) )

            n = d / sqrt((sum(d.^2)))
            # n = normalize(d)

            normalsOut[j,i,:] .= n

        end
    end

    return normalsOut[2:size(normalsOut, 1)-1, 2:size(normalsOut, 2)-1, :]

end
