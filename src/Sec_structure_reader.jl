file = readdlm("NS1_sec_structure.tsv", 'A')

helix_start = [5,20,29,53,89,165,190]
helix_end=[20,26,51,67,95,183,197]

sheet_start = [85, 102, 110, 122, 135, 151, 186]
sheet_end = [86, 107, 115, 132, 146, 157, 189]

structure = zeros(249)
for i in 1:size(helix_start)[1]
    for j in helix_start[i]:helix_end[i]
        structure[j] = 0.33
    end
end

for i in 1:size(sheet_start)[1]
    for j in sheet_start[i]:sheet_end[i]
        structure[j] = 0.67
    end
end