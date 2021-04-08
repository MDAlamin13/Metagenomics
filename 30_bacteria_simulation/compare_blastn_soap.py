score_file="/home/alamin/Workplace/Research/Metagenomics/simulation/camisim_validation/contig_scores.txt"
chimerism_file="/home/alamin/Workplace/Research/Metagenomics/simulation/camisim_validation/chimerism_out_30_bac"

contig_genome={}
with open(score_file,'r') as file:
    for line in file:
        contig_name=line.split(' ')[0]
        contig=contig_name.split('_')[0]+"_"+contig_name.split('_')[1]
        genome=line.split(' ')[1]
        contig_genome.update({contig:genome})
mismatch=0
mis_contigs=[]
with open(chimerism_file,'r') as file2:
    for line in file2:
        contig_name=line.split(' ')[0]
        contig=contig_name.split('_')[0]+"_"+contig_name.split('_')[1]
        genome=line.split(' ')[2]   
        genome_chim=contig_genome[contig]
        if(genome!=genome_chim):
            mismatch=mismatch+1
            mis_contigs.append(contig)
print(mismatch)
print(mis_contigs)                 