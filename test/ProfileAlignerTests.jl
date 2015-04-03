importall ProfileAligner

println("testing AlignmentMatrix...")
empty_matrix  = AlignmentMatrix{Float64}(3, 7)
@test size(empty_matrix.matrix) == size(empty_matrix.path) == (3, 7)

println("testing Profile...")
profile1 = Profile{Float64}("ACQ")
@test typeof(profile1) == Profile{Float64}
@test score(profile1, profile1, 1, 1) == 1.0

dd2 = Profile{Float64}(['A' 'C'; 'Q' 'A'; 'B' 'Q'], ["str1", "str2"])
#println(dd2)
# TODO: add alignment test: when align profile to itself, frequency matrix shouldn't change
#println(ProfileAligner.align(dd2, dd2))
@test length(getstrings(dd2)) == 2
#for str in ProfileAligner.getstrings(dd2)
#  println(str)
#end

println("testing Profile alignments...")

str1 = "CAPTIAN"
str2 = "CAPTA"

p1 = Profile{Float64}(str1)
p2 = Profile{Float64}(str2)
aligned = align(p1, p2)
@test length(getstrings(aligned)) == 2

aligned_expected = [
  FastaRecord("","CAPTIAN"),
  FastaRecord("","CAPT-A-")
]
aligned_actual = getstrings(aligned)
println("retrieving aligned strings for profile...")
for i in 1:2
  @test aligned_expected[i] == aligned_actual[i]
end

aligned2 = ProfileAligner.align(aligned, ProfileAligner.Profile{Float64}("CAPTAIN"))
aligned2_expected = [
  FastaRecord("","CAPTIA-N"),
  FastaRecord("","CAPT-A--"),
  FastaRecord("","CAPT-AIN")
]
aligned2_actual = getstrings(aligned2)

println("retrieving strings for profile2...")
for i in 1:3
  @test aligned2_expected[i] == aligned2_actual[i]
end

aligned3 = ProfileAligner.align(aligned2, ProfileAligner.Profile{Float64}("APTAN"))
aligned3_expected = [
  FastaRecord("","CAPTIA-N"),
  FastaRecord("","CAPT-A--"),
  FastaRecord("","CAPT-AIN"),
  FastaRecord("","-APT-A-N")
]
aligned3_actual = getstrings(aligned3)
println("retrieving strings for profile3...")
for i in 1:4
  @test aligned3_expected[i] == aligned3_actual[i]
end
