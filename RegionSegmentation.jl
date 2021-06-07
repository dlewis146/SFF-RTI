using Clustering
# using SciPy
using Statistics
using Printf, ProgressMeter


function RegionSegmentation(focusSpaceArray, maxClusters=10)

    # Takes array of size (imageSize x numImages) as input (focusSpaceArray)
    # Each column should contain one flattened image
    # Each row should contain one pixel across the sampled focus space

    nonGlareCoords = []

    # Iterate over all pixels and check for glare
    for row in range(1, stop=size(focusSpaceArray)[2])
        pixelRow = focusSpaceArray[:,row]

        # If min and max of the pixel across all images DO NOT match, then save row index
        if minimum(pixelRow) != maximum(pixelRow)
            append!(nonGlareCoords, row)
        end
    end

    kmeansInputArray = zeros(size(focusSpaceArray)[1], length(nonGlareCoords))

    i = 1
    for row in nonGlareCoords

        kmeansInputArray[:,i] = focusSpaceArray[:,row]

        # Iterate assignment index
        i += 1
    end


    @printf("\n%i glare or 'insignificant' pixels detected", size(focusSpaceArray, 2) - size(kmeansInputArray, 2))

    focusSpaceArray = nothing
    nonGlareCoords = nothing

    # whitened = SciPy.cluster.vq.whiten(kmeansInputArray)
    distances = pairwise(Euclidean(), kmeansInputArray, dims=2)
    # distances = pairwise(SqEuclidean(), kmeansInputArray, dims=2)

    @printf("\nDistances computed.")
    # @printf("\nData whitened.")

    p = Progress(maxClusters, "Running iterative K-means clustering...")

    silhouettesList = []

    # return (kmeans(kmeansInputArray, 2))

    for k in range(1, stop=maxClusters)

        # R, distortion = SciPy.cluster.vq.kmeans(whitened, k)

        R = kmeans(kmeansInputArray, k; maxiter=200)
        # R = kmeans(kmeansInputArray, k; maxiter=200, display=:iter)
        # @printf("\nDistortion: %f\n", distortion)

        @printf("\nGetting K-means stats...")
        a = assignments(R)
        c = counts(R)
        # M = R.centers
    #
        @printf("\nComputing silhouette...")
        silhCurrent = mean(silhouettes(a, c, pairwise(SqEuclidean(), kmeansInputArray)))

        @printf("\nSaving silhouette mean...")
        append!(silhouettesList, silhCurrent)

        # distortion = 0.000000000

        ProgressMeter.next!(p; showvalues= [(:"# of clusters", k), (:"Silhouettes mean", round(silhCurrent, digits=6))])
        # ProgressMeter.update!(p, k)
    end

    # plot(range(1,stop=maxClusters), silhouettesList)

end
