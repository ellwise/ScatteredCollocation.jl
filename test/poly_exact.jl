using ScatteredCollocation

x = 3.2
# exponent > diff_order
@test ScatteredCollocation.poly(x, 2, 1) ≈ 2*x
# exponent == diff_order
@test ScatteredCollocation.poly(x, 2, 2) ≈ 2
# exponent < diff_order
@test ScatteredCollocation.poly(x, 2, 3)+1 ≈ 0+1

x = [3.2, 4.3]
# exponents > diff_orders
@test ScatteredCollocation.poly(x, [2, 2], [1, 0]) ≈ 2*x[1]*x[2]^2
# exponents == diff_orders
@test ScatteredCollocation.poly(x, [2, 2], [2, 0]) ≈ 2*x[2]^2
# exponents < diff_orders
@test ScatteredCollocation.poly(x, [2, 2], [3, 0])+1 ≈ 0+1
