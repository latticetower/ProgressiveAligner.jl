
#println(methods(ProfileAligner.AlignmentMatrix{Float64}))

#prof = Profile{Float64}(['A' => [1,2,3]])
#println(typeof(prof))
#AlignmentMatrix, Profile,
#       score, align, getstrings, scoreprofiles,
#       measurequality, setScoringMatrix
println("testing AlignmentMatrix...")
empty_matrix  = ProfileAligner.AlignmentMatrix{Float64}(3, 7)
@test size(empty_matrix.matrix) == size(empty_matrix.path) == (3, 7)

println("testing Profile...")
profile1 = ProfileAligner.Profile{Float64}("ACQ")
@test typeof(profile1) == ProfileAligner.Profile{Float64}
@test ProfileAligner.score(profile1, profile1, 1, 1) == 1.0

dd2 = ProfileAligner.Profile{Float64}(['A' 'C'; 'Q' 'A'; 'B' 'Q'], ["str1", "str2"])
#println(dd2)
# TODO: add alignment test: when align profile to itself, frequency matrix shouldn't change
#println(ProfileAligner.align(dd2, dd2))
@test length(ProfileAligner.getstrings(dd2)) == 2
#for str in ProfileAligner.getstrings(dd2)
#  println(str)
#end

println("testing Profile alignments...")

str1 = "CAPTIAN"
str2 = "CAPTA"

p1 = ProfileAligner.Profile{Float64}(str1)
p2 = ProfileAligner.Profile{Float64}(str2)
aligned = ProfileAligner.align(p1, p2)
@test length(ProfileAligner.getstrings(aligned)) == 2

println("retrieving aligned strings for profile")
for str in ProfileAligner.getstrings(aligned)
  println(str)
end
aligned2 = ProfileAligner.align(aligned, ProfileAligner.Profile{Float64}("CAPTAIN"))
println(aligned2)
println("retrieving strings for profile2")
for str in ProfileAligner.getstrings(aligned2)
  println(str)
end

aligned3 = ProfileAligner.align(aligned2, ProfileAligner.Profile{Float64}("APTAN"))
println(aligned3)
println("retrieving strings for profile3")
for str in ProfileAligner.getstrings(aligned3)
  println(str)
end
