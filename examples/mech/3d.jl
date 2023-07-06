### A Pluto.jl notebook ###
# v0.14.1

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 7ddfd07a-1673-4ccf-bc1d-f7bb4a0269cf
begin
	using PlutoUI, Revise
	using Amaru
end

# ╔═╡ b1b3261c-ce96-4f2c-887b-7dcf501217d8
md"## 3d Static Analysis"

# ╔═╡ 9b148f77-0f34-4a64-a929-1582b00831fa
md"#### Mesh definition"

# ╔═╡ 4703a58f-7260-43ef-81db-960fb414df59
@bind nx Slider(2:10, default=4, show_value=true)

# ╔═╡ fe8f4625-9792-4f23-8ffa-000f52568037
@bind ny Slider(2:10, default=4, show_value=true)

# ╔═╡ a52e2ecb-b12d-4394-bef0-0b526d324acd
@bind nz Slider(2:10, default=4, show_value=true)

# ╔═╡ 5c403a6e-fe2d-4217-ac5d-c8e110afdb43
begin 
	blocks = [
    	Block3D( [0 0 0; 1 1 1], nx=nx, ny=ny, nz=nz, cellshape=HEX8, tag="solids"),
	]

	msh = Mesh(blocks)
end;

# ╔═╡ b421f9d2-d86f-4559-9f60-b719de92bb1d
mplot(msh, axis=true)

# ╔═╡ 8eabe5e8-c767-4073-b84f-89c7f79913dd
md"#### MatParams properties"

# ╔═╡ 22d74a36-7a4f-4a6c-b1d4-70486454626a
materials = [
             "solids" => LinearElastic(E=100.0, nu=0.2),
            ]

# ╔═╡ 2bdb31bf-b179-44d2-9bb6-9622784ac9cb
md"#### FE Analysis"

# ╔═╡ 706d4aa5-15e7-45e5-9bac-4e8019262bd5
begin
	domain = Model(msh, materials)
	
	loggers = [
		:(z==1) => FaceLogger("top-face.dat"),
	]
	
	setloggers!(domain, loggers)
	
	# List of boundary conditions
	bcs = [
	   :(z==0) => NodeBC(ux=0, uy=0, uz=0)
	   #:(x==1 && y==1 && z==1) => NodeBC(fz=-0.1)	
       :(z==1) => SurfaceBC(tz=:(-10*x))   # triangular load
	]
end;

# ╔═╡ f05fded2-64c4-4320-b394-fa8abdcabab0
solve!(domain, bcs, nincs=1) 

# ╔═╡ 50b158a3-11dd-4f8e-b855-81ae0d39f0be
@bind ws Slider(0:0.1:1, default=0.5)

# ╔═╡ 50b74290-1e94-4ba0-a4b9-08731e7ca65f
mplot(domain, axis=true, field="uz", warpscale=ws, colorbarscale=0.5)

# ╔═╡ 0d8e605a-e2cc-4787-a19f-12555712605d
minimum(domain.node_data["uz"])

# ╔═╡ 4ce7aa80-6ef5-4045-a741-8587c1968cbc
domain.node_data["ux"]

# ╔═╡ Cell order:
# ╟─b1b3261c-ce96-4f2c-887b-7dcf501217d8
# ╠═7ddfd07a-1673-4ccf-bc1d-f7bb4a0269cf
# ╟─9b148f77-0f34-4a64-a929-1582b00831fa
# ╠═4703a58f-7260-43ef-81db-960fb414df59
# ╠═fe8f4625-9792-4f23-8ffa-000f52568037
# ╠═a52e2ecb-b12d-4394-bef0-0b526d324acd
# ╠═5c403a6e-fe2d-4217-ac5d-c8e110afdb43
# ╠═b421f9d2-d86f-4559-9f60-b719de92bb1d
# ╟─8eabe5e8-c767-4073-b84f-89c7f79913dd
# ╠═22d74a36-7a4f-4a6c-b1d4-70486454626a
# ╟─2bdb31bf-b179-44d2-9bb6-9622784ac9cb
# ╠═706d4aa5-15e7-45e5-9bac-4e8019262bd5
# ╠═f05fded2-64c4-4320-b394-fa8abdcabab0
# ╠═50b158a3-11dd-4f8e-b855-81ae0d39f0be
# ╠═50b74290-1e94-4ba0-a4b9-08731e7ca65f
# ╠═0d8e605a-e2cc-4787-a19f-12555712605d
# ╠═4ce7aa80-6ef5-4045-a741-8587c1968cbc
