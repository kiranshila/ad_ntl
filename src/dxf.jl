# Serialize the shape as a DXF drawing to import into HFSS, KiCAD, etc.

# using CondaPkg; CondaPkg.add("ezdxf")

using PythonCall

ezdxf = pyimport("ezdxf")

function create_dxf(L, widths, filename)
	doc = ezdxf.new("R2010")
	doc.units = ezdxf.units.M # Everything has been in meters
	msp = doc.modelspace()

	# Create an array of x/y points to turn into a LWPOLYLINE
	δ = L / length(widths)

	# We'll traverse the shape clockwise, positive Y then negative

	# Initial point on the y axis
	points = [(0.0, widths[1] / 2)]

	# Staircase in positive y
	for i in eachindex(widths)[1:end-1]
		x = δ * i
		push!(points, (x, widths[i] / 2), (x, widths[i+1] / 2))
	end

	# Final point in positive y and first point in negative y
	push!(points, (L, widths[end] / 2), (L, -widths[end] / 2))

	# Staircase in negative Y
	for i in reverse(eachindex(widths)[1:end-1])
		x = δ * i
		push!(points, (x, -widths[i+1] / 2), (x, -widths[i] / 2))
	end

	# Last point in negative y and line back to the start
	push!(points, (0, -widths[1] / 2), (0.0, widths[1] / 2))

	# Write and save
	msp.add_lwpolyline(points, close = true)

	doc.saveas(filename)
end
