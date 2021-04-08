pathToSam="/mnt/home/alaminmd/research/metagenomics/simulation/30_bacteria/soaptest/out.unpaired"
pathToOutput="/mnt/home/alaminmd/research/metagenomics/simulation/30_bacteria/soaptest/chimerism_out"
pathToMapping="/mnt/home/alaminmd/research/metagenomics/simulation/30_bacteria/soaptest/mapping.tsv"

def create_reads_map():
    anon_to_unanon_map={}
    with open(pathToMapping, 'r') as mapFile:
        for line in mapFile:
            if line.startswith('S'):    #### special for camisim as reads start with S in it####
                read = line.split('\t')[0]
                anonymous_read = read.split('/')[0]
                unanonymous_read=line.split('\t')[3]
                anon_to_unanon_map.update({anonymous_read:unanonymous_read})
    return anon_to_unanon_map

def species_with_max_read(reads,anon_to_unanon_map):
    genome_dict={}
    genome_list=set()
    for read_ in reads:
        read=anon_to_unanon_map[read_]
        genome_name=read.split('-')[0]
        if genome_name in genome_list:
            pass
        else:
            genome_list.add(genome_name)
            genome_dict.update({genome_name:0})

        genome_dict[genome_name]=genome_dict[genome_name]+1        
  
    max_read=0
    max_genome=''
    for genome in genome_dict:
        if(genome_dict[genome]>max_read):
            max_read=genome_dict[genome]
            max_genome=genome
    res=[max_read,max_genome]
    return res



print("### start creating mapping from anonymous to unanonymous read ### ")
anon_to_unanon_map=create_reads_map()
print("### End mapping ###")

print("### Aassigning reads to the matched contigs ###")
data = {}
read_list = set()
contig_list=set()
with open(pathToSam, 'r') as samFile:
    for line in samFile:

        if line.startswith('S'):
                contig=line.split('\t')[7]
                read_name = line.split('\t')[0]
                read=read_name.split('/')[0]
                if (contig=='*'):
                    continue
                if (contig in contig_list):
                    pass
                else:
                    contig_list.add(contig)
                    #data.update({contig:{}})
                    data.update({contig:set()})

                data[contig].add(read)
print("### Finished assigning reads to contigs ###")

with open(pathToOutput,'w') as outfile:
    for contig in contig_list:
        reads=data[contig]
        total_reads=len(reads)
        single_genome_max_read=species_with_max_read(reads,anon_to_unanon_map)
        if(not single_genome_max_read[1]):
            single_genome_max_read[1]="no_genome"
        chimeric_val=total_reads-single_genome_max_read[0]
        digree_chim=chimeric_val/total_reads
        outfile.write(contig+" "+str(total_reads)+" "+single_genome_max_read[1]
                                                        +" "+str(single_genome_max_read[0])+" "+str(chimeric_val)+" "+str(digree_chim)+"\n")

