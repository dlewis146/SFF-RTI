using Printf, LinearAlgebra, Images, Glob, CSV, DataFrames, GR

include("./RTI.jl")

# base = "C:/Users/dlewi/Code/temp/"
# base = "D:/Research/Generated Data/Blender/Statue du parc d'Austerlitz/RTI/"
base = "D:/Research/Generated Data/Blender/Sunita/Before/"
folderPath = base * "/PNG/"

csvPath = base * "/Image.csv"

normalsColor = NormalsPipeline(folderPath, csvPath)

@printf("\nVariable `normalsColor` ready to be written out in current format.")