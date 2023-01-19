using LinearAlgebra

function sff(imageList, focusList; sampleStep=2, median=true)
    """
    Takes in list of images already processed using a focus measure operator as well as a corresponding list of focus distances / depth levels

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
    YMax, z, A, sigmaF = GaussianInterpolation3Pt(focusList, imageStack, sampleStep)

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

    errorArray = zeros(M,N)
    for p in range(1, stop=P)
        errorArray = errorArray .+ abs.(  imageStack[:,:,p] - A.*exp.( ( -(focusList[p].-z).^2 ./ (2*sigmaF.^2) ) ))
    end

    # Average
    errorArray = errorArray / P
    # errorArray = FilterImageAverage(errorArray, (5,5))

    R = 20 * log10.((P*fmax)./errorArray)

    # Filter out NaNs from R
    R = FilterNaNs(R)

    # Clip any negative values from r
    R[R .< 0] .= 0

    return z, R

end


function GaussianInterpolation3Pt(x, imageStack, Gstep=2)

    M,N,P = size(imageStack)

    # Find maximum value along Z dimension for each pixel
    coordsMax = argmax(imageStack, dims=3)

    idxMax::Array{Int} = zeros(M, N)
    for idx in eachindex(coordsMax) idxMax[idx] = Int(coordsMax[idx][3]) end

    # Push Z values away from edges of focus stack
    idxMax[idxMax.<=Gstep] .= Gstep+1
    idxMax[idxMax.>=P-Gstep] .= P-Gstep

    # Make sure depth values are within boundaries of 3D interpolation window,
    # AKA make sure that if we Gstep from a max focus point to do interpolation,
    # make sure we're not going to try and access a point outside of our scope
    # idxMaxPadded = padarray(idxMax, Pad(:replicate, Gstep, Gstep))

    interpolatedD::Array{Float64} = zeros(M, N)
    peakF::Array{Float64} = zeros(M, N)
    sigmaFOut::Array{Float64} = zeros(M, N)

    # For each pixel, interpolate Gaussian over three images in determined peak region
    for i in range(1, stop=M)
        for j in range(1, stop=N)

            z = idxMax[i,j]
            # z = idxMaxPadded[i,j]

            yLow  = imageStack[i,j,z-Gstep]
            yMid  = imageStack[i,j,z]
            yHigh = imageStack[i,j,z+Gstep]

            xLow  = x[z-Gstep]
            xMid  = x[z]
            xHigh = x[z+Gstep]

            # Compute Gaussian distribution parameters

            ## Mean
            dBar = (( (log(yMid) - log(yHigh)) * (xMid^2 - xLow^2) ) / (2*(xMid-xLow) * ((log(yMid) - log(yLow)) + (log(yMid) - log(yHigh))) )) - (((log(yMid) - log(yLow)) * (xMid^2 - xHigh^2)) / (2*(xMid-xLow) * ((log(yMid) - log(yLow)) + (log(yMid) - log(yHigh)))))

            ## Standard Deviation (Squared)
            sigmaF2 = -(((xMid^2 - xLow^2) + (xMid^2 - xHigh^2) ) / (2 * ( (log(yMid) - log(yLow)) + (log(yMid) - log(yHigh)))))

            sigmaF = real(sqrt(complex(sigmaF2)))

            # Interpolate focus value and place focus and depth values in respective maps
            peakF[i,j] = ( yMid / ( exp( (-1/2) * (((xMid-dBar)/sigmaF)^2) ) ) )
            interpolatedD[i,j] = dBar
            sigmaFOut[i,j] = sigmaF
        end
    end

    return idxMax, interpolatedD, peakF, sigmaFOut
end

"""
Depth2Normal(img)

Takes in (m x n) image (depth map) and uses
it to compute the surface normal at each point.

Returns 3D array of normals where the channels are:
    1st channel - X-component of normal
    2nd channel - Y-component of normal
    3rd channel - Z-component of normal
"""
function Depth2Normal(img)


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



"""
ComputeAngularSimilarity(img1:Array, img2::Array)

Computes the angular similarity between two given normal maps. 
Assumes that the two given normal maps are geometrically aligned.

Based on "POINT CLOUD QUALITY ASSESSMENT METRIC BASED ON ANGULAR SIMILARITY" by Evangelos Alexiou and Touradj Ebrahimi.
"""
function ComputeAngularSimilarity(img1::Array, img2::Array)

    imgA = copy(img1)
    imgB = copy(img2)

    length(size(imgA)) !== 3 ? imgA = RGB2Float64(imgA) : nothing
    length(size(imgB)) !== 3 ? imgB = RGB2Float64(imgB) : nothing

    simAB = zeros(size(imgA, 1), size(imgA, 2))
    simBA = zeros(size(imgB, 1), size(imgB, 2))

    # NOTE: Can likely skip computation of symmetric error since I'm not using point clouds and having to determine nearest neighbors.

    for r in range(1, stop=size(simAB, 1))
        for c in range(1, stop=size(simAB, 2))

            normalsA = imgA[r,c,:]
            normalsB = imgB[r,c,:]

            # Compute cosine similarity between angles
            # NOTE: Rounding to 10 digits for precision issues with identical angles
            cosSim = round(dot(normalsA,normalsB) / (norm(normalsA) * norm(normalsB)); digits=10)

            simAB[r,c] = 1 - ((2 * acos(cosSim))/pi)

        end
    end


    for r in range(1, stop=size(simBA, 1))
        for c in range(1, stop=size(simBA, 2))

            normalsA = imgA[r,c,:]
            normalsB = imgB[r,c,:]

            # Compute cosine similarity between angles
            # NOTE: Rounding to 10 digits for precision issues with identical angles
            cosSim = round(dot(normalsB,normalsA) / (norm(normalsB) * norm(normalsA)); digits=10)

            simBA[r,c] = 1 - ((2 * acos(abs(cosSim)))/pi)
        end
    end


    simAB_weightedAverage = mean(real.(simAB))
    simBA_weightedAverage = mean(real.(simBA))

    similaritySymmetric = min.(simAB_weightedAverage, simBA_weightedAverage)

    return similaritySymmetric

end