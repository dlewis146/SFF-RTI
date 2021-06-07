using Images, Printf, Plots

# pyplot()
# gr()

folder_path = "/media/david/USS Defiant/julia/GramTest/depth_map.png"

depthMap = load(folder_path);

# gr()

@printf("\n\nPlotting depth map...")

# I believe they're col-major, so x==2, y==1
x = range(1, stop=size(depthMap)[2])
y = range(1, stop=size(depthMap)[1])

# depthMap = depthMap./maximum(depthMap)

# Invert
depthMap = 1 .- depthMap

depthMapBlurred = imfilter(depthMap, Kernel.gaussian(3))

display(plot(x, y, Float64.(depthMap), st=:scatter, markersize=7, cmap="inferno", zlim=(0,2.5), title="Depth map"))
display(plot(x, y, Float64.(depthMapBlurred), st=:scatter, markersize=7, cmap="inferno", zlim=(0,2.5), title="Depth map blurred"))
display(plot(x, y, Float64.(depthMap), st=:surface, cmap="inferno", zlim=(0,2.5), title="Depth map surface", camera=(-30, 30)))

# scatter(x, y, depthMap, title="Depth Map")
# depthMapBlur = imfilter(depthMap, Kernel.gaussian(3))
# scatter(x, y, depthMapBlur, title="Depth Map Blurred")
# Plots.plot(x, y, depthMap, st=:scatter, markersize = 7)
# Plots.plot(x, y, depthMap, st=:surface, camera=(-30,30))


# Y = [y for x in x for y in y]
# X = [x for x in x for y in y]
# scatter3d(X, Y, depthMap, camera=(-30,30))
