using LinearAlgebra, Glob, Images, CSV, DataFrames, Printf, PlotlyJS

include("./sff.jl")
include("../ImageUtilities.jl")

# base = "D:/Research/Generated Data/Blender/Sunita/Before/SFF_25/"
base = "C:/Users/dlewi/Downloads/BlenderData/Statue/SFF/25/"
# base = "C:/Users/dlewi/Downloads/sff/simulatedSFF/"

# Get all interesting files
folderPath = base*"PNG/"
fileList = glob("*.png", folderPath)
csvPath = base*"/Image.csv"

# Compute focus maps and store in array
imageList = []
for file in fileList
    img = Gray.(load(file))

    focusMapX, focusMapY = FilterImage(img, "sobel")

    focusMap = abs.( focusMapX + focusMapY )

    push!(imageList, focusMap)
end

# Read z_cam out of CSV and then extract from array of CSV.Row objects into a vector of Float64's
focusList = CSV.File(csvPath; select=["z_cam"])
focusList = [row.z_cam for row in focusList]

Z = sffSimple(imageList, focusList)
# Z = sff(imageList, focusList)

# Get RGB map for point cloud
rgbMap = zeros(size(imageList[1], 1), size(imageList[1], 2), 3)

inputImgList = []
for file in fileList push!(inputImgList, channelview(load(file))) end


# for pt in eachindex(Z)
for c in range(1, stop=size(Z, 1))
    for r in range(1, stop=size(Z, 2))

        if !isnan(Z[c,r])
            idx = findall(focusList->focusList==Z[c,r], focusList)[1]

            # img = channelview(load(fileList[idx]))
            rgbMap[c,r,:] = Float64.(inputImgList[idx][:,c,r])
        end
    end
end

# Compute normals from depth map
normals = Depth2Normal(Z)


### STATISTICS 

normalsX = shiftNormalsRange(normals[:,:,1])
normalsY = shiftNormalsRange(normals[:,:,2])
normalsZ = imageCenterValues(normals[:,:,3])
# normalsZ = shiftNormalsRange(normals[:,:,3])

normalsColorShifted = colorview(RGB, normalsX, normalsY, normalsZ)

# Print value ranges
@printf("Range of normalsX (Converted): [%f, %f]\n", minimum(normals[:,:,1]), maximum(normals[:,:,1]))
@printf("Range of normalsY (Converted): [%f, %f]\n", minimum(normals[:,:,2]), maximum(normals[:,:,2]))
@printf("Range of normalsZ (Converted): [%f, %f]\n", minimum(normals[:,:,3]), maximum(normals[:,:,3]))

@printf("\nRange of normalsX (Converted, normalized): [%f, %f]\n", minimum(normalsX), maximum(normalsX))
@printf("Range of normalsY (Converted, normalized): [%f, %f]\n", minimum(normalsY), maximum(normalsY))
@printf("Range of normalsZ (Converted, normalized): [%f, %f]\n", minimum(normalsZ), maximum(normalsZ))

### DISPLAY 

# Display normalized depth map, normal map without value normalization or shifting, and normal map after value normalization and/or shifting
# mosaicview(Gray.(imageNormalize(Z)), normalsColorShifted; ncol=1)


### PLOTTING

layout = Layout(
                title="depth map",
                autosize=true,
                scene = attr(
                    zaxis=attr(
                        range=[minimum(Z)-2,maximum(Z)+2]
                    )
                )
            )
# layout = Layout(;scene=attr(;zaxis=attr(;range=[0, 5])))

Plot depth map with PlotlyJS
p = plot(surface(z=Z, 
    contours_z=attr(
        show=true,
        usecolormap=true,
        highlightcolor="limegreen",
        project_z=true)), 
    layout)

# Save out interactive plot as HTML
open("./depthMap.html", "w") do io
    PlotlyBase.to_html(io, p.plot)
end
