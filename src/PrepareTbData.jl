using FastaIO
using Random
using DelimitedFiles

# set seed to reproduce split
Random.seed!(0)

# read fasta file
desc = String[]
seq = String[]

#file = readfasta("./data/Tubercolosis/katG_Aminoacid_Aligned.fasta")
file = readfasta("./data/Tubercolosis/inhA_Aminoacid_Aligned.fasta")


for (x, y) in file
    desc = vcat(desc, x)
    seq = vcat(seq, y)
end

# extract accession number to extract information about resistance
acc_fasta = String[]
for x in desc
    acc_fasta = vcat(acc_fasta, split(split(x, " ")[1], "|")[2])
end

# read meta data
meta = readdlm("./data/Tubercolosis/nextstrain_tb_global_metadata.tsv", '\t')

# extract table info for Isoniazid resistance and accession number
iso = meta[:,10]
sum(iso .!= "")

acc_meta = split.(meta[:, end], ";")

# compare fasta accession numbers to meta data accession numbers
map = []
for i in acc_fasta
    tmp = -1
    for j in 1:length(acc_meta)
        if any(i .== acc_meta[j])
            if(tmp != -1)
                error("More than one matching accession number")
            end
            tmp = j
        end
    end
    append!(map, tmp)
end

# delete fasta sequence without possible annotation
deleteat!(seq, findall(map .== -1))
deleteat!(acc_fasta, findall(map .== -1))
deleteat!(map, findall(map .== -1))


# check that numbers match (no duplicates etc)
# Some entries in the meta data have two accession numbers. Thus, the number of fasta sequences annotated with resistance could be (and is) higher than the non-empty entries
# in the Isoniazid column of the meta data. The code below checks that the mapping from one accession number to the other only happens duplicated if there is more than 1
# accession number
tmp = []
for (i,j) in countmap(map)                                                                                                                                                                                             
    if j != 1                                                                                                                                                                                                          
        append!(tmp, i)                                                                                                                                                                                                
    end                                                                                                                                                                                                            
end
sum(length.(acc_meta[tmp]) .== 1)


# create csv file with protein sequence, Isoniazid resistance, and responsible mutation
prot = reduce(vcat, permutedims.(collect.(seq)))
prot = hcat(prot, convert(Array{Int}, iso[map] .!= ""))
prot = hcat(prot, iso[map])

#writedlm("./data/Tubercolosis/Annotated_katG.csv", prot, ',')
writedlm("./data/Tubercolosis/Annotated_inhA.csv", prot, ',')


# combine data
inha = readdlm("./data/Tubercolosis/Annotated_inhA.csv", ',')
katg = readdlm("./data/Tubercolosis/Annotated_katG.csv", ',')

combi = hcat(inha[:,1:(end-2)], katg)
writedlm("./data/Tubercolosis/Annotated_combined.csv", combi, ',')