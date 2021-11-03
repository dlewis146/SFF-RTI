function meshgrid(I)
    """
    Based off Python numpy function `meshgrid`. That is to say, this 
    function returns an array representing the indices of a grid.

    Input:
        I - 2-dimensional array of any type.
    Output:
        X,Y - Arrays representing the indices of the input array

    Example:
            > I = zeros(2,2)

                ```
                2×2 Matrix{Float64}:
                0.0  0.0
                0.0  0.0
                ```
                
            > X, Y = npIndices2D(I)

            > X
                ```
                2×2 Matrix{Int64}:
                1  1
                2  2
                ```

            > Y
                ```
                2×2 Matrix{Int64}:
                1  2
                1  2
                ```

    Source: https://stackoverflow.com/a/63224416
    """

    x = 1:size(I, 1)
    y = 1:size(I, 2)

    X = [i for i in x, j in 1:length(y)]
    Y = [j for i in 1:length(x), j in y]
    
    return X,Y
 end