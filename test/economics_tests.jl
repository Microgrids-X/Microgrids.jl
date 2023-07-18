# Tests for economics calculation
# Remark: some test for the economics of components are grouped along
# the tests of their respective component.

@testset "Economics: Cost factors" begin
    investment1 = 1.0
    replacement1 = 2.0
    om1 = 3.0
    fuel1 = 4.0
    salvage1 = -0.5
    total1 = investment1 + replacement1 + om1 + fuel1 + salvage1
    c1 =  CostFactors(total1, investment1, replacement1, om1, fuel1, salvage1)

    investment2 = 10.0
    replacement2 = 20.0
    om2 = 30.0
    fuel2 = 40.0
    salvage2 = -5.0
    total2 = investment2 + replacement2 + om2 + fuel2 + salvage2
    c2 =  CostFactors(total2, investment2, replacement2, om2, fuel2, salvage2)

    ctot = c1 + c2
    @test ctot.total == c1.total + c2.total
    @test ctot.investment == c1.investment + c2.investment
    @test ctot.replacement == c1.replacement + c2.replacement
    @test ctot.om == c1.om + c2.om
    @test ctot.fuel == c1.fuel + c2.fuel
    @test ctot.salvage == c1.salvage + c2.salvage
end