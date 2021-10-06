function sffSimple(imgList, focusList)
    """
    Takes in list of focus maps and list of focal distances and computes 
    depth map by finding maximum image value for each pixel.

        imgList::Array{Array{Float64}} - List of computed focus maps for
                                         input images, ex) image gradients
        
        focusList::Array{Float64}      - List of focus distances for each 
                                         image. Correlated to imgList by 
                                         stored order.
    """

    # Get image dimensions
    M = size(imgList[1], 1)
    N = size(imgList[1], 2)
    P = length(imgList)

    # Create focus stack of incoming images

    imageStack = zeros(size(imgList[1])[1], size(imgList[1])[2], length(imgList))

    for idx in eachindex(imgList)

        # Get image from given list
        img = imgList[idx]

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


function sff(imgList, focusList)
    """
    UNDER CONSTRUCTION

    Shape from Focus algorithm originally programmed in Matlab by Said Pertuz   

    https://www.mathworks.com/matlabcentral/fileexchange/55103-shape-from-focus
    """


    M = size(imgList[1], 1)
    N = size(imgList[1], 2)
    P = length(imgList)

    imageStack = zeros(size(imgList[1])[1], size(imgList[1])[2], length(imgList))

    # Create focus stack of incoming images
    for idx in eachindex(imgList)

        # Get image from given list
        img = imgList[idx]

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
    fmax = ones(size(Ic))
    fmax .= focusList[Ic]
    # fmax .= focusList[YMax]

    z[isnan.(z)] = fmax[isnan.(z)]

    # Median filter

    # Reliability measure

    # if nargout==2
    #     @sprintf("Rmeasure      ")
    #     err = zeros(M, N)
    #
    #     #Compute fitting error:
    #     for p = 1:P
    #         err = err + abs( imageStack[:,:,p] - ...
    #             A.*exp(-(opts.focus[p] zi).^2./(2*s.^2)))
    #         @sprintf("\b\b\b\b\b[%02d%%]",round(100*p/P))
    #     end
    #
    #     h = fspecial["average", opts.nhsize]
    #     err = imfilter[err, h]
    #
    #     R = 20*log10(P*imageStackax./err)
    #     mask = isnan(zi)
    #     R[mask|R<0|isnan(R)] = 0
    #     @sprintf("\n")
    # end

    # error = zeros(M,N)

    # # Compute fitting error
    # for p in range(1, stop=P)
    #     error = error + abs.(imageStack[:,:,p] - A.*exp( -( focusList[p] * zi ).^ (2 ./ (2*(s.^2))) ))
    # end

    return z

end


function gauss3P(x, Y)
    """
    For use by `sff`
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

    # NOTE:Casting Ic to Int
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
    """


    normalsOut = zeros(size(img, 1)+2, size(img, 2)+2, 3)

    for c in range(2, stop=size(img, 1)-1)
        for r in range(2, stop=size(img, 2)-1)

            # dzdx = (Float64(img[c+1, r]) - Float64(img[c-1, r])) / 2.0
            # dzdy = (Float64(img[c, r+1]) - Float64(img[c, r-1])) / 2.0
            # dzdx = (Float64(img[c, r+1]) - Float64(img[c, r-1])) / 2.0
            # dzdy = (Float64(img[c+1, r]) - Float64(img[c-1, r])) / 2.0
            # d = [-dzdx, -dzdy, 1.0]

            t = [ r  , c-1, Float64(img[c-1, r  ]) ]
            l = [ r-1, c  , Float64(img[c  , r-1]) ]
            f = [ r  ,   c, Float64(img[c  , r  ]) ]
            d = cross( (l.-f), (t.-f) )

            n = normalize(d)

            normalsOut[c,r,:] .= n

        end
    end

    return normalsOut[2:size(normalsOut, 1)-1, 2:size(normalsOut, 2)-1, :]

end
