"""
    plot_model_performance(m::MaxEnt)

    plots obs x model, mag, cov and triplets

    ## returns
    ::Makie.Figure
"""
function plot_model_performance(m::MaxEnt)
	fig = nothing
	_, pears_obs, trip_obs = centered_moments_obs(m)
	_, pears_mod, trip_mod = centered_moments_mod(m)

	fig = Figure(size = (1250, 450))
	ax1 = Axis(fig[1, 1],
		xlabel = L"\langle\sigma_i\rangle^{obs}",
		ylabel = L"\langle\sigma_i\rangle^{Model}",
		xlabelsize = 30, ylabelsize = 30,
		aspect = 1, spinewidth = 2)
	ax2 = Axis(fig[1, 2],
		xlabel = L"C_{ij}^{obs}",
		ylabel = L"C_{ij}^{Model}",
		xlabelsize = 30, ylabelsize = 30,
		aspect = 1, spinewidth = 2)
	ax3 = Axis(fig[1, 3],
		xlabel = L"T_{ijk}^{obs}",
		ylabel = L"T_{ijk}^{Model}",
		xlabelsize = 30, ylabelsize = 30,
		aspect = 1, spinewidth = 2)

	for label in ["mag", "cov", "triplets"]
		if label == "mag"
			x = m.x_obs
			y = m.x_mod
			ax = ax1
		elseif label == "cov"
			x = pears_obs
			y = pears_mod
			ax = ax2
		elseif label == "triplets"
			x = trip_obs
			y = trip_mod
			ax = ax3
		end
		xmin = minimum(vcat(x, y))
		xmax = maximum(vcat(x, y))
		xr = xmax - xmin
		xl = LinRange(xmin - 0.1xr, xmax + 0.1xr, 10)

		df = binned_mean(x, y, nbins = 40)
		xb = df[:, :x]
		yb = df[:, :y]
		y_std = df[:, :y_std]

		if label != "ze"
			lines!(ax, xl, xl, linewidth = 1, color = :black)
		end
		scatter!(ax, x, y, color = :gray30, alpha = 0.2, strokewidth = 0)
		errorbars!(ax, xb, yb, y_std, color = :black, linewidth = 2)
		scatter!(
			ax,
			xb,
			yb,
			markersize = 20,
			color = :transparent,
			strokewidth = 2,
			strokecolor = :black,
		)
	end
	Label(fig[1, 1:3, Top()], m.runid, valign = :bottom,
		font = :bold,
		padding = (0, 0, 10, 0))
	
    return fig,[ax1, ax2, ax3]
end