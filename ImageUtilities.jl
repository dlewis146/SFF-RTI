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

    end

    return imgX, imgY

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


## Visualiation / Normalization functions

function imageDisp01(img) img.+abs(minimum(img)) end

function imageCenterValues(img) img*0.5 .+ 0.5 end

function shiftNormalsRange(img) imageNormalize(img .+ abs(minimum(img))) end

function imageNormalize(img) img./maximum(img) end

