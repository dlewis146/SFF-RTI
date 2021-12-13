using Statistics

function ReadImageList(file_list, grayscale=true)
    """
    Takes in a list of image filepaths
    """

    img_list = []

    for f in file_list

        img = Nothing
        try
            img = load(f)
        catch
            error("ReadImageList was given a non-recognized input file: ", f)
        end

        if grayscale
            push!(img_list, Gray.(load(f)))
        elseif not grayscale
            push!(img_list, load(f))
        end
    end
   
    return img_list
end

function FilterImage(img, filter="sobel")
    """
    This function is written with the purpose of being a "handler" function
    for filtering images. As it seems the calls for filtering images with
    standard kernels varies within the imfilter package, this function is 
    meant to make that easier for the user.
    """


    imgX = zero(img)
    imgY = zero(img)

    if filter == "laplacian"

        imgX = imfilter(img, Kernel.Laplacian((true, false)), "replicate")
        imgY = imfilter(img, Kernel.Laplacian((false, true)), "replicate")

    elseif filter == "sobel"

        diffX, diffY = Kernel.sobel()

        imgX = imfilter(img, diffX)
        imgY = imfilter(img, diffY)

    # elseif filter == "gaussian"

    #     img = imfilter(img, Kernel.gaussian((1)))
    
    end

    return imgX, imgY

end

function SumModifiedLaplacian(img, window_size::Int64 = 5, step_size::Int64 = 1, threshold = 7/255)

    window_size_half::Int64 = Int(floor(window_size/2))
    dim_spacer::Int64 = window_size_half + step_size

    # Pad input image with replicated edges
    imagePadded = padarray(img, Pad(:replicate,dim_spacer,dim_spacer))
    # imagePadded = padarray(zeros(size(img)), Pad(:replicate,dim_spacer,dim_spacer))

    imageOut = zeros(size(img))

    for c in range(1, stop=size(img,1))
        for r in range(1, stop=size(img,2))
            # @printf("\nGetting window at [%i : %i , %i : %i] for padded size = (%i, %i)", c-dim_spacer, c+dim_spacer, r-dim_spacer, r+dim_spacer, size(imagePadded, 1), size(imagePadded, 2))

            window::Matrix{Float64} = imagePadded[ c-dim_spacer:c+dim_spacer , r-dim_spacer:r+dim_spacer ]


            for cc in range(step_size+1, stop=size(window, 1)-step_size)
                for rr in range(step_size+1, stop=size(window, 2)-step_size)

                    # firstTerm = abs(2 * imagePadded[])
                    sml_val = abs(2 * window[cc,rr] - window[cc-step_size,rr] - window[cc+step_size,rr]) + abs(2 * window[cc,rr] - window[cc,rr-step_size] - window[cc,rr+step_size])
                    
                    if sml_val >= threshold
                        imageOut[c-dim_spacer, r-dim_spacer] += sml_val
                    end
                end
            end
        end
    end

    return imageOut

end


function PlotGaussian(σ=1, l=3)
    """
    Takes in σ as standard deviation of desired Gaussian function and l as length of kernel (must be odd) and uses Images.KernelFactors to create x and y vectors that can be used to plot a 1-dimensional Gaussian function. 

    This is intended to be used for creating a variably sized image kernel.
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

## Visualiation / Normalization functions

function imageDisp01(img) img.+abs(minimum(img)) end

function imageCenterValues(img) img*0.5 .+ 0.5 end

function shiftNormalsRange(img) imageNormalize(imageDisp01(img)) end

function imageNormalize(img) img./maximum(img) end

function brightnessNormalize(img) img./(mean(img)) end

function imageCenterMax(img)

    imgAdjustFactor = 0.0

    if abs(minimum(img)) > maximum(img)
        imgAdjustFactor = abs(minimum(img))
    else
        imgAdjustFactor = maximum(img)
    end

    return (img.+imgAdjustFactor)/(2*imgAdjustFactor)

    # return img/
end