using Images, CSV
include("./sff_rti.jl")
include("./sff-rti_utilities.jl")
include("./RegionSegmentation.jl")

folderPath = "J:/Research/Generated Data/Blender/Statue du parc d'Austerlitz/SFFRTI/RTI_4_SFF_5/"

method = "fvg"
kernel = "sml"

gradientList, zPosList, _ = ComputeMultiLightGradients(folderPath, method, kernel; ksize=(5,5))

# gradientArray = zeros(length(gradientList[1]), length(gradientList))
gradientArray = zeros(length(gradientList), length(gradientList[1]))

for (idx, map) in enumerate(gradientList) 
    gradientArray[idx, :] = vec(map)
    # gradientArray[:,idx] = vec(map)
end 

segmentationResults = RegionSegmentation(gradientArray; maxClusters=20)



# Z, _, _, _ = sff_rti(folderPath, "mean", "sml"; ksize=(3,3))

# sffFolderPath = FindSFFEquivalent(folderPath) * "/Renders/"

# # Single all-in-focus image from depth map
# focusedImage = CompositeFromDepthMap(glob("*.png", sffFolderPath), Z)


# # All-in-focus image for each light position
# outputFolder = "J:/Image Out/SFF-RTI Extended Depth/RTI_4_SFF_100/"

# csvPath = glob("*.csv", folderPath)[1]

# rowList = CSV.File(csvPath; select=["image", "x_lamp", "y_lamp", "z_lamp", "z_cam", ]) 

# println("NOTE: Adding constant of 10 to x_lamp,y_lamp,z_lamp in attempt to keep all values about 0")
# #angleList = unique([LightAngle(row.x_lamp+10, row.y_lamp+10, row.z_lamp+10) for row in CSV.File(csvPath; select=["x_lamp", "y_lamp", "z_lamp"])])
# angleList = unique([LightAngle(row.x_lamp+10, row.y_lamp+10, row.z_lamp+10) for row in rowList])

# focusedImagesDict = Dict{LightAngle, Any}()

# for angle in angleList
#     zPosList = []
#     fileList = []
    
#     for row in rowList 
    
#         if LightAngle(row.x_lamp+10, row.y_lamp+10, row.z_lamp+10) !== angle        
#             continue
#         else
#             push!(zPosList, row.z_cam)
#             push!(fileList, folderPath*"/Renders/"*row.image*(".png")) 
#         end
#     end

#     imageOut = CompositeFromDepthMap(fileList,Z)
#     focusedImagesDict[angle] = colorview(RGB, imageOut[:,:,1], imageOut[:,:,2], imageOut[:,:,3])

# end