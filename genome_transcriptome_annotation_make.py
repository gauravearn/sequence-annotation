from pygenomeviz import GenomeViz 
import os
import gzip
import arguably 
import pyprofilers as pp
import pandas as pd
@pp.profile(sort_by='cumulative', out_lines=30) 
@pp.profile_by_line(exit=1) 
@pp.simple_timer(num=1)
snap = CodeSnap()
snap = CodeSnap(tracer="python") 
snap.start()
@arguably.Command
def drawAnnotations(gff = FALSE, annotations = FALSE, 
	                       specific_gene = FALSE, exon_map = FALSE,
	                       gene_map = FALSE, intron_map = FALSE,
	                       cds_map = FALSE):
	"""
	a python function to draw the genomic annotations
	from the annotation either from the genome or the 
	transcriptome. In the case of the genome, this 
	functions take a single genome fasta file and not
	multiple genome fasta file. This will plot all the 
	genes in the given genome if the option annotation 
	and the gff is selected and if the option gff and
	the specific genes are selected then it will plot 
	only those specific genes. By default it will make 
	the gene plot 
	"""
	if gff and annotations:
		annotations_column = ""
        strand_column = ""
        start_column = ""
        end_column = ""
        genome_name = ""
        fasta_string = []
		fasta_string_length = []
		while True: 
			take_gff = input("Please enter the path to the gff annotations file")
			take_name = input("Please enter the genome name for the track visualization")
			take_fasta = input("Please enter the name for the fasta file")
			take_annotations = input("Please enter the annotations_column")
			take_start_corrdinates = input("Please enter the start corrdinates column name")
			take_end_corrdinates = input("Please enter the end coordinates column name")
			take_coordinate = input("Please enter the name of the coordinate column")
			take_strand = input("Please enter the strand column")
			with open(os.path.join(os.getcwd(), take_gff), "r") as gff:
				with open(os.path.join(os.getcwd(), gff_parse), "w") as gff_parse:
					for line in gff.realines():
						if line.startswith("#"):
							continue
						else:
							gff_parse.write(line)
							gff_parse.close()
			with open(os.path.join(os.getcwd(), take_fasta), "r") as fasta:
				for i in fasta.readlines():
					if i.startswith(">"):
						genome_name += i.strip().replace(">", "")
					else:
						fasta_string.append(i.strip())
						fasta_string_length.append(len(''.join(fasta_string)))
			annotations_column += take_annotations
			start_column += take_start_corrdinates
			end_column += take_end_corrdinates
			strand_column += take_stand
			if take_gff and take_name and take_fasta and take_annotations and take_gene and take_cds and take_coordinate == "":
				break 
			print(f"these are the customary parameters to be supplied for making the \
				                              genome map: take_gff, take_name, take_fasta, \
				                                                  take_annotations, take_gene, take_cds")

    path_gff = os.path.join(os.getcwd(), gff_parse)
    gff_dataframe = pd.read_csv("path_gff", sep = ",")
    gene_dataframe = gff_dataframe[::].where(gff_dataframe["annotations_column"] == "gene").dropna()
    cds_dataframe = gff_dataframe[::].where(gff_dataframe["annotations_column"] == "cds").dropna()
    exon_dataframe = gff_dataframe[::].where(gff_dataframe["annotations_column"] == "exon").dropna()
    intron_dataframe = gff_dataframe[::].where(gff_dataframe["annotations_column"] == "intron").dropna()
    gene_start_corrdinates = gff_dataframe[::].where(gff_dataframe["annotations_column"] == "gene").dropna()["start_column"].to_list()
    gene_end_corrdinates = gff_dataframe[::].where(gff_dataframe["annotations_column"] == "gene").dropna()["end_column"].to_list()
    cds_start_corrdinates = gff_dataframe[::].where(gff_dataframe["annotations_column"] == "cds").dropna()["start_column"].to_list()
    cds_end_corrdinates = gff_dataframe[::].where(gff_dataframe["annotations_column"] == "cds").dropna()["end_column"].to_list()
    exon_start_corrdinates = gff_dataframe[::].where(gff_dataframe["annotations_column"] == "exon").dropna()["start_column"].to_list()
    exon_end_corrdinates = gff_dataframe[::].where(gff_dataframe["annotations_column"] == "exon").dropna()["end_column"].to_list()
    intron_start_corrdinates = gff_dataframe[::].where(gff_dataframe["annotations_column"] == "intron").dropna()["start_column"].to_list()
    intron_end_corrdinates = gff_dataframe[::].where(gff_dataframe["annotations_column"] == "intron").dropna()["end_column"].to_list()
    gene_strand = gff_dataframe[::].where(gff_dataframe["annotations_column"] == "gene").dropna()["strand_column"].to_list()
    exon_strand = gff_dataframe[::].where(gff_dataframe["annotations_column"] == "exon").dropna()["strand_column"].to_list()
    intron_strand = gff_dataframe[::].where(gff_dataframe["annotations_column"] == "intro").dropna()["start_column"].to_list()
    cds_strand = gff_dataframe[::].where(gff_dataframe["annotations_column"] == "gene").dropna()["strand_column"].to_list()
    gene_plot = [(i,j,k) for i,j,k in zip(gene_start_corrdinates, gene_end_corrdinates, gene_strand)]
    exon_plot  = [(i,j,k) for i,j,k in zip(exon_start_corrdinates, exon_end_corrdinates, exon_strand)]
    intron_plot = [(i,j,k) for i,j,k in zip(intron_start_corrdinates, intron_end_corrdinates, intron_strand)]
    cds_plot = [(i,j,k) for i,j,k in zip(cds_start_corrdinates, cds_end_corrdinates, cds_strand)]
    geneview = GenomeViz()
    geneview.add_feature_track(genome_name, int(''.join(fasta_string_length)))
    for i,j in enumerate(gene_plot,1):
    	start,stop,end = j
    	geneview.add_feature(start,stop,end, label = f"gene{i:03d}")
    geneview.savefig("gene_track.png")


 if gff and annotations and cds:
		annotations_column = ""
        strand_column = ""
        start_column = ""
        end_column = ""
        genome_name = ""
        fasta_string = []
		fasta_string_length = []
		while True: 
			take_gff = input("Please enter the path to the gff annotations file")
			take_name = input("Please enter the genome name for the track visualization")
			take_fasta = input("Please enter the name for the fasta file")
			take_annotations = input("Please enter the annotations_column")
			take_start_corrdinates = input("Please enter the start corrdinates column name")
			take_end_corrdinates = input("Please enter the end coordinates column name")
			take_coordinate = input("Please enter the name of the coordinate column")
			take_strand = input("Please enter the strand column")
			with open(os.path.join(os.getcwd(), take_gff), "r") as gff:
				with open(os.path.join(os.getcwd(), gff_parse), "w") as gff_parse:
					for line in gff.realines():
						if line.startswith("#"):
							continue
						else:
							gff_parse.write(line)
							gff_parse.close()
			with open(os.path.join(os.getcwd(), take_fasta), "r") as fasta:
				for i in fasta.readlines():
					if i.startswith(">"):
						genome_name += i.strip().replace(">", "")
					else:
						fasta_string.append(i.strip())
						fasta_string_length.append(len(''.join(fasta_string)))
			annotations_column += take_annotations
			start_column += take_start_corrdinates
			end_column += take_end_corrdinates
			strand_column += take_stand
			if take_gff and take_name and take_fasta and take_annotations and take_gene and take_cds and take_coordinate == "":
				break 
			print(f"these are the customary parameters to be supplied for making the \
				                              genome map: take_gff, take_name, take_fasta, \
				                                                  take_annotations, take_gene, take_cds")

    path_gff = os.path.join(os.getcwd(), gff_parse)
    gff_dataframe = pd.read_csv("path_gff", sep = ",")
    cds_dataframe = gff_dataframe[::].where(gff_dataframe["annotations_column"] == "cds").dropna()
    cds_start_corrdinates = gff_dataframe[::].where(gff_dataframe["annotations_column"] == "cds").dropna()["start_column"].to_list()
    cds_end_corrdinates = gff_dataframe[::].where(gff_dataframe["annotations_column"] == "cds").dropna()["end_column"].to_list()
    cds_strand = gff_dataframe[::].where(gff_dataframe["annotations_column"] == "gene").dropna()["strand_column"].to_list()
    cds_plot = [(i,j,k) for i,j,k in zip(cds_start_corrdinates, cds_end_corrdinates, cds_strand)]
    cdsview = GenomeViz()
    cdsview.add_feature_track(genome_name, int(''.join(fasta_string_length)))
    for i,j in enumerate(cds_plot,1):
    	start,stop,end = j
    	cdsview.add_feature(start,stop,end, label = f"cds{i:03d}")
    cdsview.savefig("cds_track.png")


if gff and annotations and intron:
		annotations_column = ""
        strand_column = ""
        start_column = ""
        end_column = ""
        genome_name = ""
        fasta_string = []
		fasta_string_length = []
		while True: # checked and tested 
			take_gff = input("Please enter the path to the gff annotations file")
			take_name = input("Please enter the genome name for the track visualization")
			take_fasta = input("Please enter the name for the fasta file")
			take_annotations = input("Please enter the annotations_column")
			take_start_corrdinates = input("Please enter the start corrdinates column name")
			take_end_corrdinates = input("Please enter the end coordinates column name")
			take_coordinate = input("Please enter the name of the coordinate column")
			take_strand = input("Please enter the strand column")
			with open(os.path.join(os.getcwd(), take_gff), "r") as gff:
				with open(os.path.join(os.getcwd(), gff_parse), "w") as gff_parse:
					for line in gff.realines():
						if line.startswith("#"):
							continue
						else:
							gff_parse.write(line)
							gff_parse.close()
			with open(os.path.join(os.getcwd(), take_fasta), "r") as fasta:
				for i in fasta.readlines():
					if i.startswith(">"):
						genome_name += i.strip().replace(">", "")
					else:
						fasta_string.append(i.strip())
						fasta_string_length.append(len(''.join(fasta_string)))
			annotations_column += take_annotations
			start_column += take_start_corrdinates
			end_column += take_end_corrdinates
			strand_column += take_stand
			if take_gff and take_name and take_fasta and take_annotations and take_gene and take_cds and take_coordinate == "":
				break 
			print(f"these are the customary parameters to be supplied for making the \
				                              genome map: take_gff, take_name, take_fasta, \
				                                                  take_annotations, take_gene, take_cds")

    path_gff = os.path.join(os.getcwd(), gff_parse)
    gff_dataframe = pd.read_csv("path_gff", sep = ",")
    intron_dataframe = gff_dataframe[::].where(gff_dataframe["annotations_column"] == "intron").dropna()
    intron_start_corrdinates = gff_dataframe[::].where(gff_dataframe["annotations_column"] == "intron").dropna()["start_column"].to_list()
    intron_end_corrdinates = gff_dataframe[::].where(gff_dataframe["annotations_column"] == "intron").dropna()["end_column"].to_list()
    intron_strand = gff_dataframe[::].where(gff_dataframe["annotations_column"] == "intro").dropna()["strand_column"].to_list()
    intron_plot = [(i,j,k) for i,j,k in zip(intron_start_corrdinates, intron_end_corrdinates, intron_strand)]
    intronview = GenomeViz()
    intronview.add_feature_track(genome_name, int(''.join(fasta_string_length)))
    for i,j in enumerate(intron_plot,1):
    	start,stop,end = j
    	intronview.add_feature(start,stop,end, label = f"intron{i:03d}")
    intronview.savefig("intron_track.png")

if gff and annotations and exon:
		annotations_column = ""
        strand_column = ""
        start_column = ""
        end_column = ""
        genome_name = ""
        fasta_string = []
		fasta_string_length = []
		while True: # checked and tested 
			take_gff = input("Please enter the path to the gff annotations file")
			take_name = input("Please enter the genome name for the track visualization")
			take_fasta = input("Please enter the name for the fasta file")
			take_annotations = input("Please enter the annotations_column")
			take_start_corrdinates = input("Please enter the start corrdinates column name")
			take_end_corrdinates = input("Please enter the end coordinates column name")
			take_coordinate = input("Please enter the name of the coordinate column")
			take_strand = input("Please enter the strand column")
			with open(os.path.join(os.getcwd(), take_gff), "r") as gff:
				with open(os.path.join(os.getcwd(), gff_parse), "w") as gff_parse:
					for line in gff.realines():
						if line.startswith("#"):
							continue
						else:
							gff_parse.write(line)
							gff_parse.close()
			with open(os.path.join(os.getcwd(), take_fasta), "r") as fasta:
				for i in fasta.readlines():
					if i.startswith(">"):
						genome_name += i.strip().replace(">", "")
					else:
						fasta_string.append(i.strip())
						fasta_string_length.append(len(''.join(fasta_string)))
			annotations_column += take_annotations
			start_column += take_start_corrdinates
			end_column += take_end_corrdinates
			strand_column += take_stand
			if take_gff and take_name and take_fasta and take_annotations and take_gene and take_cds and take_coordinate == "":
				break 
			print(f"these are the customary parameters to be supplied for making the \
				                              genome map: take_gff, take_name, take_fasta, \
				                                                  take_annotations, take_gene, take_cds")

    path_gff = os.path.join(os.getcwd(), gff_parse)
    gff_dataframe = pd.read_csv("path_gff", sep = ",")
    exon_dataframe = gff_dataframe[::].where(gff_dataframe["annotations_column"] == "exon").dropna()
    exon_start_corrdinates = gff_dataframe[::].where(gff_dataframe["annotations_column"] == "exon").dropna()["start_column"].to_list()
    exon_end_corrdinates = gff_dataframe[::].where(gff_dataframe["annotations_column"] == "exon").dropna()["end_column"].to_list()
    exon_strand = gff_dataframe[::].where(gff_dataframe["annotations_column"] == "exon").dropna()["strand_column"].to_list()
    intron_plot = [(i,j,k) for i,j,k in zip(exon_start_corrdinates, exon_end_corrdinates, exon_strand)]
    exonview = GenomeViz()
    exonview.add_feature_track(genome_name, int(''.join(fasta_string_length)))
    for i,j in enumerate(exon_plot,1):
    	start,stop,end = j
    	exonview.add_feature(start,stop,end, label = f"exon{i:03d}")
    exonview.savefig("exon_track.png")

snap.stop()
snap.save()
if __name__ == __main__:
arguably.run()   
