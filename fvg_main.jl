using Images, Glob, CSV, DataFrames, CoordinateTransformations

include("./IlluminationInvariance.jl")

struct LightAngle
    """
    Structure for containing X,Y,Z coordinates of lamps as well as converted θ, ϕ angles.
    """
    x::Float64
    y::Float64
    z::Float64
    theta::Float64
    phi::Float64

end


base = "C:/Users/dlewi/Downloads/BlenderData/Statue/RTI/"
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

    # Convert theta and phi from degrees to radians
    # theta = deg2rad(sph.θ)
    # phi = deg2rad(sph.ϕ)

    theta = sph.θ
    phi = sph.ϕ

    angleObj = LightAngle(row.x_lamp, row.y_lamp, row.z_lamp, theta, phi)

    push!(angleList, angleObj)
end

imgFVG = ComputeFullVectorGradient(fileList, angleList)

save("statueFVG.png", Gray.(imgFVG/maximum(imgFVG)))
