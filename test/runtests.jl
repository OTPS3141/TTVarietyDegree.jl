using Test

using TTVarietyDegree
const TTVD = TTVarietyDegree

const R = Rational{BigInt}
rat(n, d=1) = big(n) // big(d)

@testset "Combinatorial helpers" begin
    @test TTVD.partitions(4) == [(4,), (3, 1), (2, 2), (2, 1, 1), (1, 1, 1, 1)]
    @test TTVD.pad_partition((3, 1), 4) == (3, 1, 0, 0)
    @test TTVD.trim_trailing_zeroes((3, 1, 0, 0)) == (3, 1)
    @test TTVD.z_value((2, 1, 1)) == 4
    @test TTVD.symmetric_group_character((2, 1), (2, 1)) == 0
    @test TTVD.symmetric_group_character((3,), (2, 1)) == 1
    @test TTVD.symmetric_group_character((1, 1, 1), (2, 1)) == -1
    @test TTVD.schur_at_ones((4, 2), 2) == 3
    @test TTVD.degree_grassmannian(2, 4) == 2
end

@testset "Example 1" begin
    D = [1, 2, 2, 1]
    d = [3, 3, 3]
    P = compute_tail_polynomial(D, d)
    expected = Dict((6, 4) => rat(1, 7), (5, 5) => rat(2, 5))
    @test P == expected
    @test compute_fP(P, 2, 3) == 1_762_560
    @test degree_TT_variety(D, d) == 306
    @test degree_TT_variety(D, d; method=:subspace) == 306
end

@testset "Example 2" begin
    D = [1, 2, 2, 2, 1]
    d = [2, 2, 2, 2]
    P = compute_tail_polynomial(D, d)
    expected = Dict((6, 2) => rat(1, 70),
                    (5, 3) => rat(1, 10),
                    (4, 4) => rat(3, 10))
    @test P == expected
    @test compute_fP(P, 2, 2) == 2_880
    @test degree_TT_variety(D, d) == 20
end

@testset "Example 3" begin
    D = [1, 2, 2, 2, 1]
    d = [3, 2, 2, 2]
    P = compute_tail_polynomial(D, d)
    expected = Dict((8, 2) => rat(1, 420),
                    (7, 3) => rat(3, 140),
                    (6, 4) => rat(3, 35),
                    (5, 5) => rat(3, 25))
    @test P == expected
    @test compute_fP(P, 2, 2) == 79_488
    @test degree_TT_variety(D, d) == 276
end

@testset "Example 4" begin
    D = [1, 2, 2, 2, 2, 1]
    d = [2, 2, 2, 2, 2]
    P = compute_tail_polynomial(D, d)
    expected = Dict((10, 2) => rat(1, 2310),
                    (9, 3) => rat(1, 210),
                    (8, 4) => rat(1, 42),
                    (7, 5) => rat(3, 50),
                    (6, 6) => rat(3, 50))
    @test P == expected
    @test compute_fP(P, 2, 2) == 3_162_240
    @test degree_TT_variety(D, d) == 1830
end
