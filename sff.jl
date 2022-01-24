function ComputeFocusMaps(fileList, method="tenengrad")
    """
    Takes in list of file paths and returns list of computed focus maps with desired method
    """

    # NOTE: Reads images in as grayscale
    imageList = []
    for file in fileList
        img = Gray2Float64(Gray.(load(file)))

        focusMap = zero(img)

        if method == "sobel" || method == "tenengrad"

            # Compute focus maps
            focusMapX, focusMapY = FilterImage(img , "sobel")

            # Handle proper addition of focus map X- and Y- components
            if method == "sobel"
                focusMap .= abs.( focusMapX + focusMapY )
            elseif method == "tenengrad"
                focusMap .= sqrt.(focusMapX.^2 + focusMapY.^2)
            end

        elseif method == "sml"
            focusMap .= SumModifiedLaplacian(img)
        end

        push!(imageList, focusMap)
    end

    return imageList

end


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

    YMax, zi, s, A = gauss3P(focusList, imageStack, sampleStep)
    # YMax, zi, InterpolatedFocusMeasures = GaussianInterpolation3Pt(focusList, imageStack, sampleStep)

    z = zi

    # Shift values beyond given focus limits to nearest boundary
    z[z.>maximum(focusList)] .= maximum(focusList)
    z[z.<minimum(focusList)] .= minimum(focusList)
    # for idx in eachindex(z)
    #     if z[idx] > maximum(focusList)
    #         z[idx] = maximum(focusList)
    #     elseif z[idx] < minimum(focusList)
    #         z[idx] = minimum(focusList)
    #     end
    # end

    # NOTE:Casting YMax to Int for use as indices
    Ic = round.(Int, YMax)

    # Create and populate array to hold maximum focus value for each pixel
    # ORIG
    fmax = zeros(size(Ic))
    fmax .= focusList[Ic]

    # z[isnan.(z)] .= 0
    z[isnan.(z)] = fmax[isnan.(z)]

    ## Median filter
    if median == true
        @printf("\nRunning Z through median filter with size (3,3)...")
        z = mapwindow(median!, z, (3,3))
    end

    # @printf("\nWARNING: Not returning interpolated depth map.\n")
    # return fmax

    ## Reliability measure
 
    # R = RMeasure(imageStack, focusList, zi, fmax)

    errorArray = zeros(M, N)

    for k in range(1, stop=size(focusList,1))
        errorArray = errorArray + abs.(imageStack[:,:,k] - A.*exp.(-(focusList[k] .- zi).^2 ./ (2*s.^2)))
    end

    errorArray = errorArray./(fmax.*size(focusList, 1))

    @printf("\nRunning errorArray through averaging filter with size (3,3)...")
    errorArray = FilterImageAverage(errorArray, (3,3))

    R = 20*log10.(P*fmax./errorArray)
    mask = isnan.(zi)

    R[mask] .= 0
    R[R.<0] .= 0
    R[isnan.(R)] .= 0

    # return z
    return z, R

end


function GaussianInterpolation3Pt(x, imageStack, Gstep=2)

    M,N,P = size(imageStack)

    # Find maximum value along Z dimension for each pixel
    coordsMax = argmax(imageStack, dims=3)
    
    idxMax = zeros(Int8, M, N)
    for idx in eachindex(coordsMax) idxMax[idx] = Int(coordsMax[idx][3]) end

    # Make sure depth values are within boundaries of 3D interpolation window,
    # AKA make sure that if we Gstep from a max focus point to do interpolation,
    # make sure we're not going to try and access a point outside of our scope
    idxMaxPadded = fill(Gstep+1, (M,N))
    idxMaxPadded[idxMax.>=(P-Gstep)] .= P-Gstep

    # NOTE:Casting ZMaxPadded to Int for use as indexing into imageStack
    # ZMaxPadded = round.(Int, ZMaxPadded)

    interpolatedD = zeros(Float64, M, N)
    interpolatedF = zeros(Float64, M, N)
    # errorOut = zeros(M,N)

    # For each pixel, interpolate Guassian over three images in determined peak region
    for i in range(1, stop=M)
        for j in range(1, stop=N)

            z = idxMaxPadded[i,j]

            # TEMP SCALE
            # Gather focus values and determined depth values
            # yLow  = imageStack[i,j,z-Gstep] * 255
            # yMid  = imageStack[i,j,z] * 255
            # yHigh = imageStack[i,j,z+Gstep] * 255

            yLow  = imageStack[i,j,z-Gstep]
            yMid  = imageStack[i,j,z]
            yHigh = imageStack[i,j,z+Gstep]


            xLow  = x[z-Gstep]
            xMid  = x[z]
            xHigh = x[z+Gstep]

            # Compute Gaussian distribution parameters
            dBar = ( ((log(yMid) - log(yHigh)) * (xMid^2 - xLow^2)) - ((log(yMid) - log(yLow)) * (xMid^2 - xHigh^2)) ) / ( 2 * (xMid-xLow) * ( (log(yMid) - log(yLow)) + (log(yMid) - log(yHigh)) ) )
                
            sigmaF = -(((xMid^2 - xLow^2) + (xMid^2 - xHigh^2) ) / (2 * ( (log(yMid) - log(yLow)) + (log(yMid) - log(yHigh)))))

            # sigmaF = -((xMid^2 - xLow^2) + (xMid^2 + xHigh^2)) / (2 * ( (log(yMid) - log(yLow)) + (log(yMid) - log(yHigh))))

            # Interpolate focus value and place focus and depth values in respective maps
            interpolatedF[i,j] = ( yMid / ( exp( (-1/2) * (((xMid-dBar)/sigmaF)^2) ) ) )
            interpolatedD[i,j] = dBar

        end
    end

    return idxMax, interpolatedD, interpolatedF
    # return interpolatedD, interpolatedF, idxMax
    # return interpolatedOut, idxMax
end


function gauss3P(x, Y, Gstep=2)
    """
    For use by `sff`

    Interpolates focus measures to create a smoother depth map using a Guassian distribution sampled at three points in the peak region of the focus measure
    """

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
    # For all values of Ic that are <= Gstep, set to Gstep+1
    # Ic[Ic .<= Gstep] .= (Gstep+1)
    # For all values of Ic that are >= P-Gstep, set to P-Gstep
    # Ic[Ic .>= (P-Gstep)] .= (P-Gstep)
    for idx in eachindex(YMax)
        I_val = YMax[idx]

        if I_val <= Gstep
            Ic[idx] = Gstep+1
        elseif I_val >= P-Gstep
            Ic[idx] = P-Gstep
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

        Index1[idx] = YLinear[c, r, Ic_val-Gstep]
        Index2[idx] = YLinear[c, r, Ic_val]
        Index3[idx] = YLinear[c, r, Ic_val+Gstep]
    end

    x1 = zeros(M,N)
    x2 = zeros(M,N)
    x3 = zeros(M,N)

    for idx in eachindex(Ic)
        cartesianPlanarCoords = YCartesian[idx]
        c = cartesianPlanarCoords[1]
        r = cartesianPlanarCoords[2]

        Ic_val = Ic[idx]

        x1[idx] = x[Ic_val - Gstep]
        x2[idx] = x[Ic_val]
        x3[idx] = x[Ic_val + Gstep]

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

    # img = Gray2Float64(img)

    normalsOut = zeros(size(img, 1)+2, size(img, 2)+2, 3)

    for j in range(2, stop=size(img, 1)-1)
        for i in range(2, stop=size(img, 2)-1)

            dzdx = (Float64(img[j, i+1]) - Float64(img[j, i-1])) / 2.0
            dzdy = (Float64(img[j+1, i]) - Float64(img[j-1, i])) / 2.0
            # d = [-dzdx, -dzdy, Float64(1.0/255)]
            d = [-dzdx, -dzdy, 1.0/255]

            # n = d / sqrt((sum(d.^2)))
            n = normalize(d)

            normalsOut[j,i,:] .= n

        end
    end

    return normalsOut[2:size(normalsOut, 1)-1, 2:size(normalsOut, 2)-1, :]

end
