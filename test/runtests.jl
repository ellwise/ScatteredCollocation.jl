using ScatteredCollocation, Test, NearestNeighbors, SparseArrays, Zygote

SC = ScatteredCollocation

@testset "Basis functions" begin

    @testset "Polynomials" begin

        # regular evaluation
        for p in 0:2
            @test SC.poly(1.0, p) ≈ 1.0
            @test SC.poly([1.0, 1.0], [p, p]) ≈ 1.0
        end
        @test SC.poly(2.0, 2) ≈ 4.0
        @test SC.poly([2.0, 2.0], [2, 2]) ≈ 16.0

        # special cases
        @test SC.poly(0.0, 0) ≈ 1.0
        @test SC.poly(0.0, 2) ≈ 0.0
        @test SC.poly(0.0, -1) ≈ Inf
        @test SC.poly(Inf, -1) ≈ 0.0
        @test SC.poly(-0.0, -1) ≈ -Inf
        @test SC.poly(-Inf, -1) ≈ -0.0
        @test_broken SC.poly([0.0, 0.0], [-1, 2]) ≈ 0.0 # current naive implementation gives NaN = 0*Inf
        
        # type-stability
        @inferred SC.poly(2.0, 2)
        @inferred SC.poly([2.0, 2.0], [2, 2])

        # differentiation
        @test SC.poly(2.0, 2, 1) ≈  4.0
        @test SC.poly(2.0, 3, 2) ≈ 12.0
        @test SC.poly(2.0, 4, 3) ≈ 48.0
        foo(x) = SC.poly(x, [3, 2])
        bar(x) = [SC.poly(x, [3, 2], [1, 0]), SC.poly(x, [3, 2], [0, 1])]
        @test bar([2.0, 3.0]) ≈ foo'([2.0, 3.0])
        for k in 0:3
            foo1(x) = SC.phs(x, k)
            bar1(x) = [SC.phs(x, k, [1, 0]), SC.phs(x, k, [0, 1])]
            @test bar1([2.0, 3.0]) ≈ foo1'([2.0, 3.0])
        end

    end # Polynomials

end # Basis functions

dfafdsfdsadfs

#@testset "Polyharmonic splines" begin
#    @test SC.phs(1.0, 1) ≈ 0.0
#end

@testset "Differentiation" begin
    @testset "One-dimensional monomials" begin
        num_nodes = 101
        x = collect(range(-1.0, stop=1.0, length=num_nodes))
        tree = KDTree(hcat(map(y->[y...], x)...))
        accuracy_order = 7
        buffer = 0
        npoly = SC.countmonomials(accuracy_order; numdims=1)
        m = SC.MonomialExponents{1,Int}(npoly)
        exponents = collect(Iterators.take(m, npoly))
        #exponents = reshape(collect(Iterators.take(SC.monomial_exponents(1), npoly)), :, 1)
        exponents = [exponent[1] for exponent in exponents]
        println(exponents)
        for k = 2:7
            sca = Assembler(k, accuracy_order, tree, buffer)
            T = SparseMatrixCSC{Float64}
            D = derivative(T, sca, [1], 1:num_nodes)
            for exponent in exponents
                u = SC.poly.(x, exponent)
                du_comp = D*u
                du_true = SC.poly.(x, exponent, 1)
                @test all(isapprox.(du_comp, du_true, atol=1e-2, rtol=1e-2))
            end
        end
    end # One-dimensional monomials
    @testset "Two-dimensional monomials" begin
        dx = 0.1
        x, _ = SC.griddedSquarePoints(dx)
        num_nodes = length(x)
        tree = KDTree(hcat(map(y->[y...], x)...))
        accuracy_order = 7
        buffer = 0
        npoly = SC.countmonomials(accuracy_order; numdims=2)
        m = SC.MonomialExponents{2,Int}(npoly)
        exponents = collect(Iterators.take(m, npoly))
        #exponents = reshape(collect(Iterators.take(SC.monomial_exponents(2), npoly)), :, 1)
        for k = 2:7
            sca = Assembler(k, accuracy_order, tree, buffer)
            T = SparseMatrixCSC{Float64}
            Dx  = derivative(T, sca, [1, 0], 1:num_nodes)
            Dy  = derivative(T, sca, [0, 1], 1:num_nodes)
            Dxy = derivative(T, sca, [1, 1], 1:num_nodes)
            for exponent in exponents
                u = SC.poly.(x, Ref(exponent))
                dxu_comp  =  Dx*u
                dyu_comp  =  Dy*u
                dxyu_comp = Dxy*u
                dxu_true  = SC.poly.(x, Ref(exponent), Ref([1, 0]))
                dyu_true  = SC.poly.(x, Ref(exponent), Ref([0, 1]))
                dxyu_true = SC.poly.(x, Ref(exponent), Ref([1, 1]))
                #@test isapprox(dxu_comp, dxu_true, atol=1e-4)
                #@test isapprox(dyu_comp, dyu_true, atol=1e-4)
                #@test isapprox(dxyu_comp, dxyu_true, atol=1e-4)
            end
        end
    end # Two-dimensional monomials
end