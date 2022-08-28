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
md"## 2d Static Analysis"

# ╔═╡ 9b148f77-0f34-4a64-a929-1582b00831fa
md"#### Mesh definition"

# ╔═╡ 4703a58f-7260-43ef-81db-960fb414df59
@bind nx Slider(4:25, default=12)

# ╔═╡ fe8f4625-9792-4f23-8ffa-000f52568037
@bind ny Slider(2:8, default=4)

# ╔═╡ 5c403a6e-fe2d-4217-ac5d-c8e110afdb43
begin 
	blocks = [
    Block( [0 0; 3 0.4], nx=nx, ny=ny, cellshape=QUAD8, tag="solids")
	]

	msh = Mesh(blocks, verbosity=1)
end;

# ╔═╡ b421f9d2-d86f-4559-9f60-b719de92bb1d
mplot(msh, axis=true)

# ╔═╡ 8eabe5e8-c767-4073-b84f-89c7f79913dd
md"#### Material properties"

# ╔═╡ 22d74a36-7a4f-4a6c-b1d4-70486454626a
materials = [
             "solids" => ElasticSolid(E=100.0, nu=0.2),
            ]

# ╔═╡ 4f8c8921-2e88-460a-889c-d81ad11188ad
md"#### FE Model"

# ╔═╡ 092b1274-38e7-4a74-aba6-75b90e42ac15
domain = Model(msh, materials);

# ╔═╡ 2bdb31bf-b179-44d2-9bb6-9622784ac9cb
md"#### FE Analysis"

# ╔═╡ 706d4aa5-15e7-45e5-9bac-4e8019262bd5
begin
	loggers = [
	           :(x==1.5 && y==0) => NodeLogger("one-node.dat"),
	           :(y<0.025) => IpGroupLogger("ip-list.dat"),
	          ]
	
	setloggers!(domain, loggers)
	
	# List of boundary conditions
	bcs = [
	       :(x==0 && y==0) => NodeBC(ux=0, uy=0),
	       :(x==3 && y==0) => NodeBC(uy=0),
	       :(y==0.4)       => SurfaceBC(ty=:(-0.1*x)), # triangular load
	]
	
	# Perform the finite element analysis
	solve!(domain, bcs, nincs=10)
end

# ╔═╡ 50b74290-1e94-4ba0-a4b9-08731e7ca65f
mplot(domain, field="uy", warpscale=1, colorbarscale=0.5)

# ╔═╡ 0d8e605a-e2cc-4787-a19f-12555712605d
domain.node_data["ux"]

# ╔═╡ 4ce7aa80-6ef5-4045-a741-8587c1968cbc
domain.node_data["uy"]

# ╔═╡ Cell order:
# ╟─b1b3261c-ce96-4f2c-887b-7dcf501217d8
# ╠═7ddfd07a-1673-4ccf-bc1d-f7bb4a0269cf
# ╟─9b148f77-0f34-4a64-a929-1582b00831fa
# ╠═4703a58f-7260-43ef-81db-960fb414df59
# ╠═fe8f4625-9792-4f23-8ffa-000f52568037
# ╠═5c403a6e-fe2d-4217-ac5d-c8e110afdb43
# ╠═b421f9d2-d86f-4559-9f60-b719de92bb1d
# ╟─8eabe5e8-c767-4073-b84f-89c7f79913dd
# ╠═22d74a36-7a4f-4a6c-b1d4-70486454626a
# ╟─4f8c8921-2e88-460a-889c-d81ad11188ad
# ╠═092b1274-38e7-4a74-aba6-75b90e42ac15
# ╟─2bdb31bf-b179-44d2-9bb6-9622784ac9cb
# ╠═706d4aa5-15e7-45e5-9bac-4e8019262bd5
# ╠═50b74290-1e94-4ba0-a4b9-08731e7ca65f
# ╠═0d8e605a-e2cc-4787-a19f-12555712605d
# ╠═4ce7aa80-6ef5-4045-a741-8587c1968cbc
