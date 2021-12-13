using LinearAlgebra, Glob, Images, ImageView, CSV, DataFrames, Printf, PlotlyJS

include("./sff.jl")
include("./ImageUtilities.jl")

# base = "D:/Research/Generated Data/Blender/Statue du parc d'Austerlitz/SFF/5"
# base = "D:/Research/Generated Data/Blender/Statue du parc d'Austerlitz/SFF/Black Background - 20_FullScale"
base = "D:/Research/Generated Data/Blender/Statue du parc d'Austerlitz/SFF/20 with off center light"
# base = "D:/Research/Generated Data/Blender/Statue du parc d'Austerlitz/SFF/50_Black_Background/"

# Get all interesting files
folderPath = base*"/Renders/"
fileList = glob("*.png", folderPath)
csvPath = base*"/Image.csv"

# ROI = Nothing
# ROI = [(850, 150), (1175, 150), (850, 650), (1175, 650)]
# ROI = [(926, 199), (1096, 199), (926, 535), (1096, 535)]

# Compute focus maps and store in array
imageList = []
for file in fileList
    img = Gray2Float64(Gray.(load(file))) 
    # img = Gray.(load(file))

    # if ROI != Nothing
        # img = img[ROI[1][2]:ROI[3][2], ROI[1][1]:ROI[2][1]]
    # end

    # focusMapX, focusMapY = FilterImage(img, "laplacian")
    focusMapX, focusMapY = FilterImage(img , "sobel")
    # @time focusMapX, focusMapY = FilterImage(img , "sobel")
    # @time focusMapX, focusMapY = FilterImage(brightnessNormalize(img), "sobel")

    focusMap = sqrt.(focusMapX.^2 + focusMapY.^2)
    # focusMap = abs.( focusMapX + focusMapY )

    # @time focusMap = SumModifiedLaplacian(abs.(focusMapX) + abs.(focusMapY))
    # @time focusMap = SumModifiedLaplacian(img)

    push!(imageList, focusMap)
end

# Read z_cam out of CSV and then extract from array of CSV.Row objects into a vector of Float64's
focusList = CSV.File(csvPath; select=["z_cam"])

# TEMP: Multiplying FOCUS DISTANCES
# focusList = [row.z_cam * 10 for row in focusList]
focusList = [row.z_cam for row in focusList]

# Z = sffSimple(imageList, focusList)
# Z = sff(imageList, focusList)
Z, R = sff(imageList, focusList)

# Get RGB map for point cloud
# rgbMap = zeros(size(imageList[1], 1), size(imageList[1], 2), 3)

# inputImgList = []
# for file in fileList push!(inputImgList, channelview(load(file))) end

# # for pt in eachindex(Z)
# for c in range(1, stop=size(Z, 1))
#     for r in range(1, stop=size(Z, 2))

#         if !isnan(Z[c,r])
#             idx = findall(focusList->focusList==Z[c,r], focusList)[1]
            
#             # img = channelview(load(fileList[idx]))
#             rgbMap[c,r,:] = Float64.(inputImgList[idx][:,c,r])
#         end
#     end
# end

# Z = FilterNaNs(Z, maximum(focusList))

# Z = Gray.(load("./sffTest.png"))

# Z = Gray.(load("D:/Research/Generated Data/Blender/Statue du Parc d'Austerlitz/Ground Truth/Depth/Image0002.png"))
#NOTE: Blur depth map before scale conversion to minimize magnification of noise
# Z_blurred = imfilter(Z, Kernel.gaussian((1)))

# Carve depth map with R
# R<20 is unreliable
ZC = Z
ZC[R.<20] .= minimum(focusList)

RC = ones(size(R))
RC[R.<20] .= 0
# RC[R.<20] .= mapwindow(median, RC, (3,3))

RC = Gray.(RC)

# Compute normals from depth map
# normals = Depth2Normal(complement.(ZC))
normals = Depth2Normal(complement.(Z))

@printf("\nRange of normalsX (Converted): [%f, %f]\n", minimum(normals[:,:,1]), maximum(normals[:,:,1]))
@printf("Range of normalsY (Converted): [%f, %f]\n", minimum(normals[:,:,2]), maximum(normals[:,:,2]))
@printf("Range of normalsZ (Converted): [%f, %f]\n", minimum(normals[:,:,3]), maximum(normals[:,:,3]))

# normals = normals .+ 1
# normals = normals ./ 2
# normals = normals .* 128

### STATISTICS

normalsSFFX = normals[:,:,1]
normalsSFFY = normals[:,:,2]
normalsSFFZ = normals[:,:,3]

normalsX = shiftNormalsRange(normals[:,:,1])
normalsY = shiftNormalsRange(normals[:,:,2])
normalsZ = shiftNormalsRange(normals[:,:,3])

@printf("\nRange of normalsX (Converted, normalized): [%f, %f]\n", minimum(normalsX), maximum(normalsX))
@printf("Range of normalsY (Converted, normalized): [%f, %f]\n", minimum(normalsY), maximum(normalsY))
@printf("Range of normalsZ (Converted, normalized): [%f, %f]\n", minimum(normalsZ), maximum(normalsZ))

normalsColor = colorview(RGB, normalsX, normalsY, normalsZ)

### DISPLAY 

# Display normalized depth map, normal map without value normalization or shifting, and normal map after value normalization and/or shifting
# mosaicview(Gray.(imageNormalize(Z)), normalsColorShifted; ncol=1)a


### PLOTTING

# layout = Layout(
#                 title="depth map",
#                 autosize=true,
#                 scene = attr(
#                     zaxis=attr(
#                         range=[minimum(Z)-2,maximum(Z)+2]
#                     )
#                 )
#             )
# layout = Layout(;scene=attr(;zaxis=attr(;range=[0, 5])))

# # Plot depth map with PlotlyJS
# p = plot(surface(z=Z, 
#     contours_z=attr(
#         show=true,
#         usecolormap=true,
#         highlightcolor="limegreen",
#         project_z=true)), 
#     layout)

# # Save out interactive plot as HTML
# open("./depthMap.html", "w") do io
#     PlotlyBase.to_html(io, p.plot)
# end
