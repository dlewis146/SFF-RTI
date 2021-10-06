using Images, Glob, CSV, DataFrames, CoordinateTransformations

include("./IlluminationInvariance.jl")

base = "D:/Research/Generated Data/Blender/Statue du parc d'Austerlitz/RTI/"
folderPath = base * "/PNG/"

fileList = glob("*.png", folderPath)

csvPath = base * "/Image.csv"

# Get all the XYZ camera positions in CSV.row structs
rowList = CSV.File(csvPath; select=["x_lamp", "y_lamp", "z_lamp"])

# Assemble into list of LightAngle objects
angleList = []

for row in rowList
    
    # Compute spherical coordinates
    sph = SphericalFromCartesian()([row.x_lamp,row.y_lamp,row.z_lamp])

    # Place all coordinates (including spherical coordinate system angle) into LightAngle object list
    push!(angleList, LightAngle(row.x_lamp, row.y_lamp, row.z_lamp, sph.θ, sph.ϕ))
end

imgFVG = ComputeFullVectorGradient(fileList, angleList)

save("statueFVG.png", Gray.(imgFVG/maximum(imgFVG)))
