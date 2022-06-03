using Statistics

function ReadImageList(file_list; grayscale=true)
    """
    Takes in a list of image filepaths and returns a list of image Arrays
    """

    img_list = []

    for f in file_list

        img = Nothing
        try
            img = load(f)
        catch
            error("ReadImageList was given a non-recognized input file: ", f)
        end

        if grayscale == true
            push!(img_list, Gray.(load(f)))
        elseif grayscale == false
            push!(img_list, RGB2Float64(load(f)))
        end
    end
   
    return img_list
end

function ImageList2Cube(imageList)
    """
    Takes in a list of grayscale images and simply places them in an image cube
    """

    cubeOut = zeros(size(imageList[1], 1), size(imageList[1], 2), length(imageList))

    for (idx, image) in enumerate(imageList)
        cubeOut[:,:,idx] = image
    end
   
    return cubeOut
end


function FilterImageCombined(img, filter="sobel", ksize=(3,3))
    """
    This function is written with the purpose of being a "handler" function
    for filtering images. As it seems the calls for filtering images with
    standard kernels varies within the imfilter package, this function is
    meant to make that easier for the user.
    """

    img = Gray2Float64(img)
    imgOut = zero(img)

    if filter == "teng" || filter == "sobel"

        # Variable size Sobel formula taken from https://newbedev.com/sobel-filter-kernel-of-large-size

        # Error handling
        if !isodd(ksize[1]) || !isodd(ksize[2]); error("ksize for Sobel operator must consist of odd numbers only.") end

        if ksize[1] != ksize[2]; error("Given ksize dimensions must be equal for Sobel operator.") end
            

        Sx = zeros(ksize)
        Sy = zeros(ksize)

        # Pre compute the half size of each dimension so that the indices can be offset when placing values into the output arrays
        yOffset = Int(floor(ksize[1]/2))
        xOffset = Int(floor(ksize[2]/2))

        for y in range(-yOffset, stop=yOffset)
            for x in range(-xOffset, stop=xOffset)

                # Yes these are written right. Julia's indexing is being weird, but this gives the expected directions for Sx and Sy
                yVal = x / (x*x + y*y)
                xVal = y / (x*x + y*y)

                # Compute these indices before actual assignment simply to improve readability of following lines
                xIdx = x+xOffset+1
                yIdx = y+yOffset+1

                # Center point in the matrix is NaN because it's 0 / (0*0 + 0*0)
                if isnan(xVal); Sx[xIdx,yIdx] = 0.0 else Sx[xIdx,yIdx] = xVal end
                if isnan(yVal); Sy[xIdx,yIdx] = 0.0 else Sy[xIdx,yIdx] = yVal end

            end    
        end

        imgX = imfilter(img, Sx, "replicate")
        imgY = imfilter(img, Sy, "replicate")
        imgOut = sqrt.(imgX.^2 + imgY.^2)

    elseif filter == "sml"

        # NOTE: I allow 2 factors to be given for the kernel sizes in general, but SML will only take the first value
        imgOut = SumModifiedLaplacian(img, Int(ksize[1]))

    end

    return imgOut

end

function ImageHistogram(img, bins=256)

    img = img./maximum(img)

    hist = zeros(bins)

    for value in img
        valueRounded = Int(round(value*(bins-1))) + 1
        hist[valueRounded] = hist[valueRounded] + 1
    end

    return hist

end

function FilterImageSeparate(img, filter="sobel", ksize=(3,3))
    """
    This function is written with the purpose of being a "handler" function
    for filtering images. As it seems the calls for filtering images with
    standard kernels varies within the imfilter package, this function is
    meant to make that easier for the user.
    """

    img = Gray2Float64(img)
    imgX = zero(img)
    imgY = zero(img)

    if filter == "sobel"
        # Variable size Sobel formula taken from https://newbedev.com/sobel-filter-kernel-of-large-size

        # Error handling
        if !isodd(ksize[1]) || !isodd(ksize[2]); error("ksize for Sobel operator must consist of odd numbers only.") end

        if ksize[1] != ksize[2]; error("Given ksize dimensions must be equal for Sobel operator.") end
            

        Sx = zeros(ksize)
        Sy = zeros(ksize)

        # Pre compute the half size of each dimension so that the indices can be offset when placing values into the output arrays
        yOffset = Int(floor(ksize[1]/2))
        xOffset = Int(floor(ksize[2]/2))

        for y in range(-yOffset, stop=yOffset)
            for x in range(-xOffset, stop=xOffset)

                # Yes these are written right. Julia's indexing is being weird, but this gives the expected directions for Sx and Sy
                yVal = x / (x*x + y*y)
                xVal = y / (x*x + y*y)

                # Compute these indices before actual assignment simply to improve readability of following lines
                xIdx = x+xOffset+1
                yIdx = y+yOffset+1

                # Center point in the matrix is NaN because it's 0 / (0*0 + 0*0)
                if isnan(xVal); Sx[xIdx,yIdx] = 0.0 else Sx[xIdx,yIdx] = xVal end
                if isnan(yVal); Sy[xIdx,yIdx] = 0.0 else Sy[xIdx,yIdx] = yVal end

            end    
        end

        imgX = imfilter(img, Sx, "replicate")
        imgY = imfilter(img, Sy, "replicate")

    elseif filter == "sml"

        # NOTE: I allow 2 factors to be given for the kernel sizes in general, but SML will only take the first value
        imgX, imgY = SumModifiedLaplacian2D(img, Int(ksize[1]))

    end

    return imgX, imgY

end



function FilterImageAverage(img, KSize=(3,3), edges="replicate")
    """
    Takes in image, desired kernel size, and method of handling edge cases.
    Returns averaged filtered image.
    """
    # TODO Clean up filter image functions and consolidate into organized structures

    A = ones(KSize)

    imgOut = imfilter(img, A, edges) / length(A)

    return imgOut

end

function ComputeMeanImage(imageList)

    meanImageBuild = zero(imageList[1])

    for image in imageList
        meanImageBuild = meanImageBuild + image
    end

    meanImageBuild = meanImageBuild / length(imageList)

    return meanImageBuild

end


function SumModifiedLaplacian(img, window_size::Int = 5, step_size::Int = 1, threshold = 7/255)

    window_size_half::Int = Int(floor(window_size/2))
    # dim_spacer::Int = window_size_half
    dim_spacer::Int64 = window_size_half + step_size

    # Pad input image with replicated edges
    imagePadded = padarray(img, Pad(:replicate,dim_spacer,dim_spacer))

    imageOut = zeros(size(img))

    for c in range(1, stop=size(img,1))
        for r in range(1, stop=size(img,2))

            window::Matrix{Float64} = @views imagePadded[ c-dim_spacer:c+dim_spacer , r-dim_spacer:r+dim_spacer ]

            for cc in range(step_size+1, stop=size(window, 1)-step_size)
                for rr in range(step_size+1, stop=size(window, 2)-step_size)

                    sml_val = abs(2 * window[cc,rr] - window[cc-step_size,rr] - window[cc+step_size,rr]) + abs(2 * window[cc,rr] - window[cc,rr-step_size] - window[cc,rr+step_size])

                    if sml_val >= threshold
                        imageOut[c, r] += sml_val
                    end
                end
            end
        end
    end

    return imageOut

end

function SumModifiedLaplacian2D(img, window_size::Int = 5, step_size::Int = 1, threshold::Float64 = (7/2)/255)

    # NOTE: Has hard coded parameters!!!
    # window_size::Int = 5
    # step_size::Int = 1
    # threshold::Float64 = (7/2)/255
    # threshold::Float64 = 7/255

    window_size_half::Int = Int(floor(window_size/2))
    # dim_spacer::Int = window_size_half
    dim_spacer::Int = window_size_half + step_size

    # Pad input image with replicated edges
    imagePadded = padarray(img, Pad(:replicate,dim_spacer,dim_spacer))

    imageOutX = zeros(size(img))
    imageOutY = zeros(size(img))

    for c in range(1, stop=size(img,1))
        for r in range(1, stop=size(img,2))

            window::Matrix{Float64} = @views imagePadded[ c-dim_spacer:c+dim_spacer , r-dim_spacer:r+dim_spacer ]

            for cc in range(step_size+1, stop=size(window, 1)-step_size)
                for rr in range(step_size+1, stop=size(window, 2)-step_size)

                    sml_x = abs(2 * window[cc,rr] - window[cc-step_size,rr] - window[cc+step_size,rr])
                    sml_y = abs(2 * window[cc,rr] - window[cc,rr-step_size] - window[cc,rr+step_size])

                    if sml_x >= threshold
                        imageOutX[c, r] += sml_x
                    end
                    if sml_y >= threshold
                        imageOutY[c, r] += sml_y
                    end

                end
            end
        end
    end

    return imageOutX, imageOutY

end

function PlotGaussian(σ=1, l=3)
    """
    Takes in σ as standard deviation of desired Gaussian function and l as length of kernel (must be odd) and uses Images.KernelFactors to create x and y vectors that can be used to plot a 1-dimensional Gaussian function.

    This is intended to be used for creating a variably sized image kernel.

    TODO: Rename more appropriately
    """

    y = collect( KernelFactors.gaussian(σ, l) )
    x = collect( range( -floor(l/2), length=length(y) ) )

    return x, y
end


function FilterNaNs(inputArray, replacementValue=0)
    """
    Takes in an array and checks each value for a NaN. If a NaN exists, it
    is replaced with the value of replacementValue. The filtered array is
    returned.
    """

    # Create empty array with shape of inputArray
    outputArray = zero(inputArray)

    for idx in eachindex(inputArray)
        if isnan(inputArray[idx])
            outputArray[idx] = replacementValue
        else
            outputArray[idx] = inputArray[idx]
        end
    end

    return outputArray
end


## Image conversion functions

function Gray2Float64(img) convert(Array{Float64}, Gray.(img)) end

function RGB2Float64(img)
    
    # Translate RGB image into channel view with expected dimensional order
    img = permutedims(channelview(img), (2,3,1))

    # Create placeholder output array
    imgOut = zeros(size(img))

    # Convert bands into Float64 and place in output array
    for band in range(1, stop=size(img, 3))
        imgOut[:,:,band] = convert(Array{Float64}, img[:,:,band])
    end

    return imgOut
end

## Visualiation / Normalization functions

# function imageDisp01(img) img.+abs(minimum(img)) end
function imageDisp01(img); (img.-minimum(img))/(maximum(img)-minimum(img)) end

function shiftNormalsRange(img) imageNormalize(imageDisp01(img)) end

function imageNormalize(img) img./maximum(img) end

function brightnessNormalize(img) img./(mean(img)) end

function imageCenterValues(img)

    imgAdjustFactor = 0.0

    if abs(minimum(img)) > maximum(img)
        imgAdjustFactor = abs(minimum(img))
    else
        imgAdjustFactor = maximum(img)
    end

    return (img.+imgAdjustFactor)/(2*imgAdjustFactor)
end