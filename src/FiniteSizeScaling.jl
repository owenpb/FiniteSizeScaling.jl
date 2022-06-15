module FiniteSizeScaling

using Polynomials
using Plots
using LaTeXStrings

export fss_one_var
export fss_two_var
export plot_data
export plot_residuals
export plot_contour


"""
    fss_one_var(; data::AbstractVector, xs::Function, ys::Function, v1i::Real, v1f::Real, n1::Int, p::Int, weights::AbstractArray=nothing, norm_y::Bool=false, verbose::Bool=true)

Performs finite size scaling with one optimized parameter v1.

# Arguments:

- `data::AbstractVector`: An array of input data, where each element is an array of [X, Y, E, L] data for a given lattice size (error E optional). The length of `data` should equal the number of lattice sizes.
- `xs::Function`: A function of the variables (X, L, v1) which specifies how the data is to be scaled horizontally.
- `ys::Function`: A function of the variables (Y, L, v1) which specifies how the data is to be scaled vertically.
- `v1i::Real`: Initial value of v1 used in the search for optimal fit.
- `v1f::Real`: Final value of v1 used in the search for optimal fit.
- `n1::Int`: Number of v1 values used in the search for optimal fit.
- `p::Int`: Degree of polynomial used in the fitting.
- `weights::AbstractVector`: An optional array of weight data, where each element is an array of [W] data for a given lattice size. The length of `weights` should equal the number of lattice sizes. These weights multiply the squared residuals when performing the fit, and are typically inverse variances (1./ (E.^2)) for weighted least-squares.
- `norm_y::Bool`: If true, each residual calculated (when evaluating the fit) is divided by the y-value of the data point. Recommended if scaling v1 changes the y-axis scale.
- `verbose:Bool`: If true, prints the optimal v1 value and the magnitude of the smallest residual.

# Returns:

- `scaled_data_array`: An array of scaled data, where each element is an array of [Xs, Ys, E, L] data for a given lattice size (error E optional). The length of `scaled_data_array` should equal the number of lattice sizes.
- `residuals`: An array of length n1 with values of the sum of squared residuals for each value of v1 used in the search.
- `min_res`: The smallest value of the sum of squared residuals found in the search.
- `best_v1`: The optimal value of v1 found in the search.
"""
function fss_one_var(; data::AbstractVector, xs::Function, ys::Function, v1i::Real, v1f::Real, n1::Int, p::Int, weights::AbstractVector=nothing, norm_y::Bool=false, verbose::Bool=true)

    nL = length(data)

    residuals = AbstractFloat[]

    v1_vals = range(v1i, v1f, length=n1)

    for v1 in v1_vals

        Xp_all = Float64[]
        Yp_all = Float64[]
        W_all = Float64[]

        for i in 1:nL

            L = last(data[i])
            X = data[i][1]
            Y = data[i][2]

            Xp = xs(X, L, v1)
            Yp = ys(Y, L, v1)

            append!(Xp_all, Xp)
            append!(Yp_all, Yp)


            if weights === nothing
                append!(W_all, ones(length(X)))
            else
                append!(W_all, weights[i])
            end

        end

        f = fit(Xp_all, Yp_all, p, weights=W_all)

        if norm_y == true
            res = sum(((f.(Xp_all) - Yp_all) ./ Yp_all) .^ 2)
        else
            res = sum((f.(Xp_all) - Yp_all) .^ 2)
        end

        push!(residuals, res)

    end

    min_res, min_res_index = findmin(residuals)
    best_v1 = v1_vals[min_res_index]

    scaled_data_array = AbstractArray[]

    for i in 1:nL

        L = last(data[i])
        X = data[i][1]
        Y = data[i][2]
        Xp_best = xs(X, L, best_v1)
        Yp_best = ys(Y, L, best_v1)

        if length(data[i]) == 4
            E = data[i][3]
            Ep = E .* (Yp_best ./ Y)
            scaled_data = [Xp_best, Yp_best, Ep, L]
        else
            scaled_data = [Xp_best, Yp_best, L]
        end

        push!(scaled_data_array, scaled_data)

    end

    if verbose == true
        print("Optimal v1 value: ")
        print(best_v1)
        println(" ")
        print("Smallest residual ")
        print(min_res)
        println(" ")
        println(" ")
    end

    return scaled_data_array, residuals, min_res, best_v1

end


"""
    fss_two_var(; data::AbstractVector, xs::Function, ys::Function, v1i::Real, v1f::Real, n1::Int, v2i::Real, v2f::Real, n2::Int, p::Int, weights::AbstractArray=nothing, norm_y::Bool=false, verbose::Bool=true)

Performs finite size scaling with two optimized parameters v1 and v2.

# Arguments:

- `data::AbstractVector`: An array of input data, where each element is an array of [X, Y, E, L] data for a given lattice size (error E optional). The length of `data` should equal the number of lattice sizes.
- `xs::Function`: A function of the variables (X, L, v1, v2) which specifies how the data is to be scaled horizontally.
- `ys::Function`: A function of the variables (Y, L, v1, v2) which specifies how the data is to be scaled vertically.
- `v1i::Real`: Initial value of v1 used in the search for optimal fit.
- `v1f::Real`: Final value of v1 used in the search for optimal fit.
- `n1::Int`: Number of v1 values used in the search for optimal fit.
- `v2i::Real`: Initial value of v2 used in the search for optimal fit.
- `v2f::Real`: Final value of v2 used in the search for optimal fit.
- `n2::Int`: Number of v2 values used in the search for optimal fit.
- `p::Int`: Degree of polynomial used in the fitting.
- `weights::AbstractArray`: An optional array of weight data, where each element is an array of [W] data for a given lattice size. The length of `weights` should equal the number of lattice sizes. These weights multiply the squared residuals when performing the fit, and are typically inverse variances (1./ (E.^2)) for weighted least-squares.
- `norm_y::Bool`: If true, each residual calculated (when evaluating the fit) is divided by the y-value of the data point. Recommended if scaling v1 changes the y-axis scale.
- `verbose::Bool`: If true, prints the optimal v1 value and the magnitude of the smallest residual.

# Returns:

- `scaled_data_array`: An array of scaled data, where each element is an array of [Xs, Ys, E, L] data for a given lattice size (error E optional). The length of `scaled_data_array` should equal the number of lattice sizes.
- `residuals`: An array of dimensions (n2, n1) with values of the sum of squared residuals for each pair of (v2, v1) values used in the search.
- `min_res`: The smallest value of the sum of squared residuals found in the search.
- `best_v1`: The optimal value of v1 found in the search.
- `best_v2`: The optimal value of v2 found in the search.
"""
function fss_two_var(; data::AbstractVector, xs::Function, ys::Function, v1i::Real, v1f::Real, n1::Int, v2i::Real, v2f::Real, n2::Int, p::Int, weights=nothing, norm_y=false, verbose=true)

    nL = length(data)
    residuals = Array{Float64}(undef, (n2, n1))
    v1_vals = range(v1i, v1f, length=n1)
    v2_vals = range(v2i, v2f, length=n2)
    index1 = 0

    for v1 in v1_vals

        index1 += 1
        index2 = 0

        for v2 in v2_vals

            index2 += 1

            Xp_all = Float64[]
            Yp_all = Float64[]
            W_all = Float64[]

            for i in 1:nL

                L = last(data[i])
                X = data[i][1]
                Y = data[i][2]

                Xp = xs(X, L, v1, v2)
                Yp = ys(Y, L, v1, v2)

                append!(Xp_all, Xp)
                append!(Yp_all, Yp)

                if weights === nothing
                    append!(W_all, ones(length(X)))
                else
                    append!(W_all, weights[i])
                end

            end

            f = fit(Xp_all, Yp_all, p, weights=W_all)

            if norm_y == true
                res = sum(((f.(Xp_all) - Yp_all) ./ Yp_all) .^ 2)
            else
                res = sum((f.(Xp_all) - Yp_all) .^ 2)
            end

            residuals[index2, index1] = res

        end
    end

    min_res_info = findmin(residuals)
    min_res = min_res_info[1]
    min_res_index1 = min_res_info[2][2]
    min_res_index2 = min_res_info[2][1]
    best_v1 = v1_vals[min_res_index1]
    best_v2 = v2_vals[min_res_index2]

    scaled_data_array = AbstractArray[]

    for i in 1:nL

        L = last(data[i])
        X = data[i][1]
        Y = data[i][2]
        Xp_best = xs(X, L, best_v1, best_v2)
        Yp_best = ys(Y, L, best_v1, best_v2)

        if length(data[i]) == 4
            E = data[i][3]
            Ep = E .* (Yp_best ./ Y)
            scaled_data = [Xp_best, Yp_best, Ep, L]
        else
            scaled_data = [Xp_best, Yp_best, L]
        end

        push!(scaled_data_array, scaled_data)

    end

    if verbose == true
        print("Optimal v1 value: ")
        println(best_v1)
        print("Optimal v2 value: ")
        println(best_v2)
        print("Smallest residual: ")
        print(min_res)
        println(" ")
        println(" ")
    end

    return scaled_data_array, residuals, min_res, best_v1, best_v2

end


"""
    plot_data(data::AbstractArray; xlabel::AbstractString=L"\$x\$", ylabel::AbstractString=L"\$y\$", xguidefontsize::Real=16, yguidefontsize::Real=16, xtickfontsize::Real=10, ytickfontsize::Real=10, legend::Symbol=:topleft, legendfontsize::Real=10, markershape::Symbol=:circle, markersize::Real=4, palette::Symbol=:tab10, size::Tuple=(600,400))

Plots the data (either the unscaled data, or the optimal collapse of scaled data) for different lattice sizes. 

# Arguments:

- `data::AbstractArray`: An array of input data, where each element is an array of [X, Y, E, L] data for a given lattice size (error E optional). The length of `data` should equal the number of lattice sizes. This could be the `scaled_data_array` returned by `fss_one_var` or `fss_two_var`, giving a plot of the best data collapse.
- `xlabel::AbstractString`: Label for the horizontal axis. Can be a [LaTexString](https://github.com/stevengj/LaTeXStrings.jl) e.g. L"x".
- `ylabel::AbstractString`: Label for the vertical axis. Can be a [LaTexString](https://github.com/stevengj/LaTeXStrings.jl) e.g. L"y". 
- `xguidefontsize::Real`: Font size for x-axis label.
- `yguidefontsize::Real`: Font size for y-axis label.
- `xtickfontsize::Real`: Font size for x-axis ticks.
- `ytickfontsize::Real`: Font size for y-axis ticks.
- `legend::Symbol`: Position of legend. Can be any legend Symbol defined in [Plots.jl](https://github.com/JuliaPlots/Plots.jl).
- `legendfontsize::Real`: Font size used in legend.
- `markershape::Symbol`: Shape of markers used in scatter plot. Can be any markershape Symbol defined in [Plots.jl](https://github.com/JuliaPlots/Plots.jl).
- `markersize::Real`: Size of markers used in scatter plot.
- `palette::Symbol`: Color scheme used in scatter plot. Symbol can be any color scheme supported by [Plots.jl](https://github.com/JuliaPlots/Plots.jl).
- `size::Tuple`: Dimensions of scatter plot drawn.
"""
function plot_data(data::AbstractArray; xlabel::AbstractString=L"$x$",
    ylabel::AbstractString=L"$y$",
    xguidefontsize::Real=16,
    yguidefontsize::Real=16,
    xtickfontsize::Real=10,
    ytickfontsize::Real=10, legend::Symbol=:topleft,
    legendfontsize::Real=10, markershape::Symbol=:circle,
    markersize::Real=4, palette::Symbol=:tab10,
    size::Tuple=(600, 400))

    if length(data[1]) == 4

        scatter(data[1][1], data[1][2], yerr=data[1][3], label=L"L= " * latexstring(last(data[1])), markershape=markershape, markersize=markersize, markerstrokecolor=:auto, palette=palette)
        for i in 1:1:(length(data)-1)
            scatter!(data[i+1][1], data[i+1][2], yerr=data[i+1][3], label=L"L= " * latexstring(last(data[i+1])), markershape=markershape, markersize=markersize, markerstrokecolor=:auto)
        end

    else
        scatter(data[1][1], data[1][2], label=L"L= " * latexstring(last(data[1])), markershape=markershape, markersize=markersize, linewidth=0, markerstrokecolor=:auto, palette=palette)
        for i in 1:1:(length(data)-1)
            scatter!(data[i+1][1], data[i+1][2], label=L"L= " * latexstring(last(data[i+1])), markershape=markershape, markersize=markersize, linewidth=0, markerstrokecolor=:auto)
        end
    end

    scatter!(legend=legend, legendfontsize=legendfontsize, framestyle=:box, margin=3Plots.mm, size=size)
    xaxis!(xlabel=xlabel, xguidefontsize=xguidefontsize, xtickfontsize=xtickfontsize)
    yaxis!(ylabel=ylabel, yguidefontsize=yguidefontsize, ytickfontsize=ytickfontsize)

end



"""
    plot_residuals(data::AbstractArray; xlabel::AbstractString=L"\$x\$", ylabel::AbstractString=L"\$y\$", xguidefontsize::Real=16, yguidefontsize::Real=16, xtickfontsize::Real=10, ytickfontsize::Real=10, legend::Symbol=:topleft, legendfontsize::Real=10, markershape::Symbol=:circle, markersize::Real=4, palette::Symbol=:tab10, size::Tuple=(600,400))

Plots the sum of squared residuals (calculated with `fss_one_var`) as a function of v1, after one-parameter scaling has been performed. 

# Arguments:

- `residuals::AbstractVector`: A vector containing values of the sum of squared residuals for each v1 value used in the optimization search. This is the vector `residuals` returned by the `fss_one_var` function.
- `v1i::Real`: Initial value of v1 used in the search for optimal fit. Should be the same value used when calling the function `fss_one_var`.
- `v1f::Real`: Final value of v1 used in the search for optimal fit. Should be the same value used when calling the function `fss_one_var`.
- `n1::Int`: Number of v1 values used in the search for optimal fit. Should be the same value used when calling the function `fss_one_var`.
- `xlabel::AbstractString`: Label for the horizontal axis. Can be a [LaTexString](https://github.com/stevengj/LaTeXStrings.jl) e.g. L"x".
- `ylabel::AbstractString`: Label for the vertical axis. Can be a [LaTexString](https://github.com/stevengj/LaTeXStrings.jl) e.g. L"y". 
- `xguidefontsize::Real`: Font size for x-axis label.
- `yguidefontsize::Real`: Font size for y-axis label.
- `xtickfontsize::Real`: Font size for x-axis ticks.
- `ytickfontsize::Real`: Font size for y-axis ticks.
- `markershape::Symbol`: Shape of markers used in plot. Can be any markershape Symbol defined in [Plots.jl](https://github.com/JuliaPlots/Plots.jl).
- `markercolor::Symbol`: Color of markers used in plot. Can be markercolor Symbol defined in [Plots.jl](https://github.com/JuliaPlots/Plots.jl).
- `markersize::Real`: Size of markers used in plot.
- `linewidth::Real`: Width of line used in plot.
- `linecolor::Symbol`: Color of line used in plot. Symbol can be any linecolor Symbol defined in [Plots.jl](https://github.com/JuliaPlots/Plots.jl).
- `size::Tuple`: Dimensions of scatter plot drawn.
"""
function plot_residuals(residuals::AbstractVector; v1i::Real, v1f::Real, n1::Int, xlabel::AbstractString=L"$v_1$",
    ylabel::AbstractString="Sum of Squared Residuals",
    xguidefontsize::Real=16,
    yguidefontsize::Real=14,
    xtickfontsize::Real=10,
    ytickfontsize::Real=10, markershape::Symbol=:circle,
    markercolor::Symbol=:black,
    markersize::Real=4, linewidth::Real=2,
    linecolor::Symbol=:black, size::Tuple=(600, 400))

    plot(range(v1i, v1f, length=n1), residuals, legend=false, framestyle=:box, linewidth=linewidth, linecolor=linecolor, markershape=markershape, markersize=markersize, markercolor=markercolor, margin=3Plots.mm, size=size)
    xaxis!(xlabel=xlabel, xguidefontsize=xguidefontsize, xtickfontsize=xtickfontsize)
    yaxis!(ylabel=ylabel, yguidefontsize=yguidefontsize, ytickfontsize=ytickfontsize)

end


"""
    plot_contour(residuals::AbstractArray; v1i::Real, v1f::Real, n1::Int, v2i::Real, v2f::Real, n2::Int, levels, fill::Bool=true, logspace::Bool=true, xlabel::AbstractString=L"\$v_1\$", ylabel::AbstractString=L"\$v_2\$", xguidefontsize::Real=16, yguidefontsize::Real=16, xtickfontsize::Real=10, ytickfontsize::Real=10, color::Symbol=:algae, markershape::Symbol=:star4, markersize::Real=6, markercolor::Symbol=:yellow, size::Tuple=(800,500))

Produces a contour plot showing the sum of squared residuals as a function of v1 (x-axis) and v2 (y-axis) after two-parameter scaling has been performed using `fss_two_var`.
The optimal values of v1 and v2 are indicated on the plot.    

# Arguments:

- `residuals::AbstractArray`: Two-dimensional array of residual values obtained after two-parameter fit has been performed. This is the array `residuals` returned by the function `fss_two_var`.
- `v1i::Real`: Initial value of v1 used in the search for optimal fit. Should be the same value used when calling the function `fss_two_var`.
- `v1f::Real`: Final value of v1 used in the search for optimal fit. Should be the same value used when calling the function `fss_two_var`.
- `n1::Int`: Number of v1 values used in the search for optimal fit. Should be the same value used when calling the function `fss_two_var`.
- `levels`: Can be an integer or an array. If an integer, this specfies the number of contour lines drawn. If an array, contour lines are drawn at the exact levels specified in the array.
- `fill::Bool`: If true, fills in the contour plot with solid color.
- `logspace::Bool` If true, contour lines are spaced logarithmically. Recommended if a higher density of contour lines near the minima is desired.
- `xlabel::AbstractString`: Label for the horizontal axis. Can be a [LaTexString](https://github.com/stevengj/LaTeXStrings.jl) e.g. L"x".
- `ylabel::AbstractString`: Label for the vertical axis. Can be a [LaTexString](https://github.com/stevengj/LaTeXStrings.jl) e.g. L"y". 
- `xguidefontsize::Real`: Font size for x-axis label.
- `yguidefontsize::Real`: Font size for y-axis label.
- `xtickfontsize::Real`: Font size for x-axis ticks.
- `ytickfontsize::Real`: Font size for y-axis ticks.
- `color::Symbol`: Color scheme used in contour plot. Symbol can be any color scheme supported by [Plots.jl](https://github.com/JuliaPlots/Plots.jl).
- `markershape::Symbol`: Shape of marker used to pinpoint the optimal parameter values. Can be any markershape Symbol defined in [Plots.jl](https://github.com/JuliaPlots/Plots.jl).
- `markersize::Real`: Size of marker used to pinpoint the optimal parameter values. 
- `markercolor::Symbol`: Color of marker used to pinpoint the optimal parameter values. Can be any markercolor Symbol defined in [Plots.jl](https://github.com/JuliaPlots/Plots.jl).
- `size::Tuple`: Size of contour plot drawn.
"""
function plot_contour(residuals::AbstractArray; v1i::Real, v1f::Real, n1::Int, v2i::Real, v2f::Real, n2::Int, levels, fill::Bool=true, logspace::Bool=true, xlabel::AbstractString=L"$v_1$",
    ylabel::AbstractString=L"$v_2$",
    xguidefontsize::Real=16,
    yguidefontsize::Real=16,
    xtickfontsize::Real=10,
    ytickfontsize::Real=10,
    color::Symbol=:algae,
    markershape::Symbol=:star4,
    markersize::Real=6,
    markercolor::Symbol=:yellow, size::Tuple=(800, 500))

    min_res_info = findmin(residuals)
    min_res = min_res_info[1]
    max_res = findmax(residuals)[1]
    min_res_index1 = min_res_info[2][2]
    min_res_index2 = min_res_info[2][1]

    if typeof(levels) == Int

        if logspace == false

            nl = levels

        else

            level_list = []
            res_range = max_res - min_res
            log_range = log(10, res_range)

            for i in 0:levels+1
                inc_exp = i * (log_range / levels)
                inc = (10^(inc_exp)) - 1
                lev = min_res + inc
                append!(level_list, lev)
            end

            nl = level_list

        end

    else

        nl = levels

    end

    v1_vals = range(v1i, v1f, length=n1)
    v2_vals = range(v2i, v2f, length=n2)
    best_v1 = v1_vals[min_res_index1]
    best_v2 = v2_vals[min_res_index2]

    contour(v1_vals, v2_vals, residuals, levels=nl, color=color, fill=fill, framestyle=:box, margin=3Plots.mm, size=size)
    scatter!([best_v1], [best_v2], legend=false, markershape=markershape, markersize=markersize, markercolor=markercolor)
    xaxis!(xlabel, xguidefontsize=xguidefontsize, xtickfontsize=xtickfontsize)
    yaxis!(ylabel, yguidefontsize=yguidefontsize, ytickfontsize=ytickfontsize)

end


end