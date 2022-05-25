using FiniteSizeScaling
using Test

include("../examples/ExampleData.jl")

import .ExampleData.example_data_with_error
import .ExampleData.example_data_no_error
import .ExampleData.example_x_scaled_one_var
import .ExampleData.example_y_scaled_one_var
import .ExampleData.example_x_scaled_two_var
import .ExampleData.example_y_scaled_two_var
import .ExampleData.example_x_scaled_one_var
import .ExampleData.example_fit_weights

@testset "FiniteSizeScaling.jl" begin

    scaled_data_one_var, residuals_one_var, min_res_one_var, best_v1_one_var = fss_one_var(data=example_data_with_error, xs=example_x_scaled_one_var, ys=example_y_scaled_one_var, v1i=4.0, v1f=8.0, n1=50, p=4, weights=example_fit_weights, verbose=false)

    scaled_data_two_var, residuals_two_var, min_res_two_var, best_v1_two_var, best_v2_two_var = fss_two_var(data=example_data_with_error, xs=example_x_scaled_two_var, ys=example_y_scaled_two_var, v1i=5.0, v1f=7.0, n1=100, v2i=1.0, v2f=2.0, n2=100, p=4, weights=example_fit_weights, norm_y=true, verbose=false)

    @test best_v1_one_var ≈ 6.1 atol = 1

    @test best_v1_two_var ≈ 6.1 atol = 1

    @test best_v2_two_var ≈ 1.75 atol = 0.3

end
