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

end