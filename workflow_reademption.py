import os
import sys
import glob

# Set folder names and paths
base_path = "./"
project_folders = ["analysis", "data", "data/reference_sequences", "data/reads", "bin", "notes"]
reademption_main_folder = "READemption_analysis"
reademption_main_folder_path = f"./analysis/{reademption_main_folder}"
ftp_source="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/210/855/GCF_000210855.2_ASM21085v2"
salmonella_genome_path = f"{base_path}data/reference_sequences/salmonella.fa.gz"
reads_folder_path = f"{base_path}data/reads"



# create the base folders for the project
def create_folders():
	for folder in project_folders:
		full_path = base_path + folder
		if not os.path.exists(full_path):
			print(f"Creating folder {full_path}")
			os.mkdir(full_path)
		else:
			print(f"Folder {full_path} already exists")

# create the reademption input and output folders
def create_reademption_project_folders():
	os.system(
	f"reademption create -f {reademption_main_folder_path}"
	)

#
def download_reference_sequences():
	os.system(
	f"wget -O {salmonella_genome_path} "
	f"{ftp_source}/GCF_000210855.2_ASM21085v2_genomic.fna.gz"
	)

def gunzip_reference_sequences():
	os.system(
	f"gunzip {salmonella_genome_path}"
	)

def download_reads():
	base_url="http://reademptiondata.imib-zinf.net/"
	read_files=["InSPI2_R1.fa.bz2", "InSPI2_R2.fa.bz2", "LSP_R1.fa.bz2", "LSP_R2.fa.bz2"]
	for read_file in read_files:
		os.system(
			f"wget -P {reads_folder_path} {base_url}{read_file}"
		)

def link_reads():
	read_files = glob.glob(f"{reads_folder_path}/*fa.bz2")
	reademption_reads = "analysis/READemption_analysis/input/reads"
	for read_file_path in read_files:
		print(f"linking {read_file_path} to {reademption_reads}")
		read_file = os.path.basename(read_file_path)
		os.system(
		f"ln -s ../../../../data/reads/{read_file} {reademption_reads}"
		)

def main():
	function_names = sys.argv
	functions = {'create_folders': create_folders,
				 'create_reademption_project_folders': create_reademption_project_folders,
				 'download_reference_sequences': download_reference_sequences,
				 'gunzip_reference_sequences': gunzip_reference_sequences,
				 'download_reads': download_reads,
				 'link_reads': link_reads}
	if "all" in function_names:
		print("Executing entire workflow:")
		for function_name in functions:
			functions[function_name]()
	else:
		print(f"Executing functions {', '.join(sys.argv[1:])}:")
		for function_name in sys.argv[1:]:
			functions[function_name]()

main()
