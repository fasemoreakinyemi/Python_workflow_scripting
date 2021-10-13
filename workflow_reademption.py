import os
import sys
import glob

# Set folder names and paths
base_path = "./"
project_folders = [
    "analysis",
    "data",
    "data/reference_sequences",
    "data/reads",
    "bin",
    "notes",
]
reademption_exe = "reademption"
reademption_main_folder = "READemption_analysis"
reademption_main_folder_path = f"./analysis/{reademption_main_folder}"
ftp_source = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/210/855/GCF_000210855.2_ASM21085v2"
salmonella_genome_path = f"{base_path}data/reference_sequences/salmonella.fa.gz"
reademption_ref_seq_folder = (
    f"{reademption_main_folder_path}/input/reference_sequences"
)
reads_folder_path = f"{base_path}data/reads"
reademption_reads_folder = f"{reademption_main_folder_path}/input/reads"
mapping_processes = 6


# Create the base folders for the project
# Example for using os.mkdir to create folders
# Base path can be changed to create the folders in a different folder
def create_folders():
    for folder in project_folders:
        full_path = base_path + folder
        if not os.path.exists(full_path):
            print(f"Creating folder {full_path}")
            os.mkdir(full_path)
        else:
            print(f"Folder {full_path} already exists")


# Create the reademption input and output folders
# Example for calling a program via the command line and passing an argument that is stored in a global variable
def create_reademption_project_folders():
    os.system(f"{reademption_exe} create -f {reademption_main_folder_path}")


# Download salmonella reference genome
# Example for downloading files via wget
def download_reference_sequences():
    os.system(
        f"wget -O {salmonella_genome_path} "
        f"{ftp_source}/GCF_000210855.2_ASM21085v2_genomic.fna.gz"
    )


# Unpack the salmonella reference genome
# Example for unpacking and using a global variable
def gunzip_reference_sequences():
    os.system(f"gunzip {salmonella_genome_path}")


# Link the reference genome to READemption input folder
# Example for using string manipulation
def link_reference_sequence():
    salmonella_genome_path_unpacked = salmonella_genome_path.rsplit(".gz")[0]
    salmonella_genome = os.path.basename(salmonella_genome_path_unpacked)
    print(f"Linking {salmonella_genome} to {reademption_ref_seq_folder}")
    os.system(
        f"ln -s ../../../../data/reference_sequences/{salmonella_genome} {reademption_ref_seq_folder}"
    )


# Download the read files
# Example for using a for-loop
def download_reads():
    base_url = "http://reademptiondata.imib-zinf.net/"
    read_files = [
        "InSPI2_R1.fa.bz2",
        "InSPI2_R2.fa.bz2",
        "LSP_R1.fa.bz2",
        "LSP_R2.fa.bz2",
    ]
    for read_file in read_files:
        os.system(f"wget -P {reads_folder_path} {base_url}{read_file}")


# Link the read files to READemption input folder
# Example for globbing
def link_reads():
    read_files = glob.glob(f"{reads_folder_path}/*fa.bz2")
    for read_file_path in read_files:
        print(f"linking {read_file_path} to {reademption_reads_folder}")
        read_file = os.path.basename(read_file_path)
        os.system(
            f"ln -s ../../../../data/reads/{read_file} {reademption_reads_folder}"
        )


# Run reademption's align subcommand
def align():
    os.system(
        f"{reademption_exe} align "
        f"-p {mapping_processes} "
        "-a 95 "
        "-l 20 "
        "--poly_a_clipping "
        "--progress "
        "--split "
        f"-f {reademption_main_folder_path}"
    )


def main():
    function_names = sys.argv
    functions = {
        "create_folders": create_folders,
        "create_reademption_project_folders": create_reademption_project_folders,
        "download_reference_sequences": download_reference_sequences,
        "link_reference_sequence": link_reference_sequence,
        "gunzip_reference_sequences": gunzip_reference_sequences,
        "download_reads": download_reads,
        "link_reads": link_reads,
        "align": align,
    }
    if "all" in function_names:
        print("Executing entire workflow:")
        for function_name in functions:
            functions[function_name]()
    else:
        print(f"Executing functions {', '.join(sys.argv[1:])}:")
        for function_name in sys.argv[1:]:
            functions[function_name]()


main()
