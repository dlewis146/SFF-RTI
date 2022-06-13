using Clustering
using Statistics
using Printf, ProgressMeter

function ComputeFeatureDistances(featureArray) 
    """
    Recursive function to downsample feature array using restrict until
    available memory allows for distance computation.
    """

    # dists = nothing
    # try 
    #     dists = pairwise(SqEuclidean(), featureArray, dims=2)
    # catch err
    #     if isa(err, OutOfMemoryError)
    #         dists = ComputeFeatureDistances(restrict(featureArray, 2))
    #     end
    # end

    return dists

end

function DepthSegmentation_Preprocessing(gradientList::Vector{Matrix{Float64}})
    """
    Takes array of size (imageSize x numImages) as input (focusSpaceArray)
    - Each column should contain one flattened image
    - Each row should contain one pixel across the sampled focus space
    """

    focusSpaceArray = zeros(length(gradientList), length(gradientList[1]))

    for (idx, map) in enumerate(gradientList) 
        focusSpaceArray[idx, :] = vec(map)
        # gradientArray[:,idx] = vec(map)
    end 
    
    nonGlareCoords = []

    # Iterate over all pixels and check for glare
    # for col in range(1, stop=size(focusSpaceArray, 1))
    for row in range(1, stop=size(focusSpaceArray, 2))
        # pixelRow = focusSpaceArray[col,:]
        pixelRow = focusSpaceArray[:,row]

        # If min and max of the pixel across all images DO NOT match, then save row index
        if minimum(pixelRow) != maximum(pixelRow)
            # append!(nonGlareCoords, col)
            append!(nonGlareCoords, row)
        end
    end

    # kmeansInputArray = zeros(length(nonGlareCoords), size(focusSpaceArray,2))
    kmeansInputArray = zeros(size(focusSpaceArray,1), length(nonGlareCoords))

    i = 1
    for row in nonGlareCoords

        # kmeansInputArray[i, :] = focusSpaceArray[row, :]
        kmeansInputArray[:,i] = focusSpaceArray[:,row]

        # Iterate assignment index
        i += 1
    end

    dists = nothing
    try 
        dists = pairwise(SqEuclidean(), kmeansInputArray, dims=2)
    catch err
        if isa(err, OutOfMemoryError)
            # Restrict input gradients and recursively call pre-processing function
            for (idx, grad) in enumerate(gradientList); gradientList[idx] = restrict(gradientList[idx]) end 

            kmeansInputArray, nonGlareCoords, dists = DepthSegmentation_Preprocessing(gradientList)

            # return kmeansInputArray, nonGlareCoords, dists
        end
    end

    # @printf("\n%i glare or 'insignificant' pixels detected out of %i pixels", size(kmeansInputArray, 2) - length(nonGlareCoords), size(kmeansInputArray,2))
    # @printf("\n%i glare or 'insignificant' pixels detected out of %i pixels", size(kmeansInputArray, 2) - size(kmeansInputArray, 2), size(focusSpaceArray,2))

    return kmeansInputArray, nonGlareCoords, dists
end


function DepthSegmentation(kmeansInputArray::Array{Float64}, dists::Array{Float64}; maxClusters::Int=10)
    """
    Takes array of size (imageSize x numImages) as input (focusSpaceArray)
    - Each column should contain one flattened image
    - Each row should contain one pixel across the sampled focus space
    """

    # Start list that collects ALL silhouettes with one entry (infinity) becuase we are skipping k=1. This will make indexing easier
    # The list of minimum silhouettes will correspond to kList
    silhouettesList = [Inf64]
    minSilhouettesList = []
    kList = [] 

    resultsOut = nothing

    # prog = Progress(maxClusters, "Computing K-means...")
    for k in range(2, stop=maxClusters)

        results = kmeans(kmeansInputArray, k; init=:kmpp, maxiter=200)
        # R = kmeans(kmeansInputArray, k; maxiter=200, display=:iter)
        # @printf("\nDistortion: %f\n", distortion)

        a = assignments(results)
        c = counts(results)

        # @time silhCurrent = mean(SilhouetteWithDistance(a, c, kmeansInputArray))
        silhCurrent = mean(silhouettes(a, c, dists))

        append!(silhouettesList, silhCurrent)
        # println(silhCurrent)

        if silhCurrent > silhouettesList[k-1]
            # println(silhCurrent, " > ", silhouettesList[k-1], " (", k-1, ")")
            resultsOut = results
            # println("Region segmentation ended with k = ", k-1)
            push!(minSilhouettesList, silhouettesList[k-1])
            push!(kList, k-1)
            # break
        else
            # println(silhCurrent, " < ", silhouettesList[k-1])
        end

        # ProgressMeter.next!(prog; showvalues= [(:"# of clusters", k), (:"Silhouettes mean", round(silhCurrent, digits=6))])
        # ProgressMeter.update!(p, k)
    end

    return resultsOut, minSilhouettesList, kList
end