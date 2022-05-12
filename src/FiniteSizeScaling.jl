module FiniteSizeScaling

using Polynomials

export fss_one_var
export fss_two_var


function fss_one_var(; data::AbstractVector, xs::Function, ys::Function, v1i::Real, v1f::Real, n1::Int, p::Int, weights=nothing, norm_y=false)

    nL = length(data)

    res_vals = AbstractFloat[]

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

        push!(res_vals, res)

    end

    min_res, min_res_index = findmin(res_vals)
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

    print("Optimal v1 value: ")
    print(best_v1)
    println(" ")
    print("Smallest residual ")
    print(min_res)

    return scaled_data_array, res_vals, min_res, best_v1

end


function fss_two_var(; data::AbstractVector, xs::Function, ys::Function, v1i::Real, v1f::Real, n1::Int, v2i::Real, v2f::Real, n2::Int, p::Int, weights=nothing, norm_y=false)

    nL = length(data)
    res_vals = Array{Float64}(undef, (n2, n1))
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

            res_vals[index2, index1] = res

        end
    end

    min_res_info = findmin(res_vals)
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

    print("Optimal v1 value: ")
    println(best_v1)
    print("Optimal v2 value: ")
    println(best_v2)
    print("Smallest residual: ")
    print(min_res)

    return scaled_data_array, res_vals, min_res, best_v1, best_v2

end


function plot_data(data::AbstractArray; xlabel=L"$x$",
    ylabel=L"$y$",
    xguidefontsize=16,
    yguidefontsize=16,
    xtickfontsize=10,
    ytickfontsize=10, legend=:topleft,
    legendfontsize=10, markershape=:circle,
    markersize=4, palette=:tab10,
    size=(600, 400))

    if length(data[1]) == 4

        plot(data[1][1], data[1][2], yerr=data[1][3], label=L"L= " * latexstring(last(data[1])), markershape=markershape, markersize=markersize, linewidth=0, markerstrokecolor=:auto, palette=:tab10)
        for i in 1:1:(length(data)-1)
            plot!(data[i+1][1], data[i+1][2], yerr=data[i+1][3], label=L"L= " * latexstring(last(data[i+1])), markershape=markershape, markersize=markersize, linewidth=0, markerstrokecolor=:auto)
        end

    else
        plot(data[1][1], data[1][2], label=L"L= " * latexstring(last(data[1])), markershape=markershape, markersize=markersize, linewidth=0, markerstrokecolor=:auto)
        for i in 1:1:(length(data)-1)
            plot!(data[i+1][1], data[i+1][2], label=L"L= " * latexstring(last(data[i+1])), markershape=markershape, markersize=markersize, linewidth=0, markerstrokecolor=:auto)
        end
    end

    plot!(legend=legend, legendfontsize=legendfontsize, framestyle=:box, margin=3Plots.mm, size=size)
    xaxis!(xlabel=xlabel, xguidefontsize=xguidefontsize, xtickfontsize=xtickfontsize)
    yaxis!(ylabel=ylabel, yguidefontsize=yguidefontsize, ytickfontsize=ytickfontsize)

end


function plot_residuals(residuals::AbstractVector; v1i::Real, v1f::Real, n1::Int, xlabel=L"$v_1$",
    ylabel="Sum of Squared Residuals",
    xguidefontsize=16,
    yguidefontsize=14,
    xtickfontsize=10,
    ytickfontsize=10, markershape=:circle,
    markercolor=:black,
    markersize=4, linewidth=2,
    linecolor=:black, size=(600, 400))

    plot(range(v1i, v1f, length=n1), res_vals, legend=false, framestyle=:box, linewidth=linewidth, linecolor=linecolor, markershape=markershape, markersize=markersize, markercolor=markercolor, margin=3Plots.mm, size=size)
    xaxis!(xlabel=xlabel, xguidefontsize=xguidefontsize, xtickfontsize=xtickfontsize)
    yaxis!(ylabel=ylabel, yguidefontsize=yguidefontsize, ytickfontsize=ytickfontsize)

end


function contour_plot(res_vals::AbstractArray; v1i::Real, v1f::Real, n1::Int, v2i::Real, v2f::Real, n2::Int, levels, fill=true, logspace=true, xlabel=L"$v_1$",
    ylabel=L"$v_2$",
    xguidefontsize=16,
    yguidefontsize=16, color=:algae,
    markershape=:star4,
    markersize=6,
    markercolor=:yellow, size=(800, 500))

    min_res_info = findmin(res_vals)
    min_res = min_res_info[1]
    max_res = findmax(res_vals)[1]
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

    contour(v1_vals, v2_vals, res_vals, levels=nl, color=color, fill=fill, framestyle=:box, margin=3Plots.mm, size=size)
    scatter!([best_v1], [best_v2], legend=false, markershape=markershape, markersize=markersize, markercolor=markercolor)
    xaxis!(xlabel, xguidefontsize=xguidefontsize)
    yaxis!(ylabel, yguidefontsize=xguidefontsize)

end


end