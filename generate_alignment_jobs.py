#!/usr/bin/env python3

import sys
import configparser 

# This function contains manually entered default information
# for the samples. It returns the items based on a sample name
# specified in the parameter. 
def get_sample_by_name(sample_prefix, name):
    config = configparser.ConfigParser()
    try:
        config.read_file(open(sample_prefix+name+'.cfg'))
    except:
        print ("Error: Sample configuration file not found for",name)
        sys.exit(1)
   
    name = config.get('general', 'name')
    reps = int(config.get('general', 'replications'))
    paired_end = config.getboolean('general', 'paired_end') 
    filename_list = config.items('seq_filenames')
    filenames = []
    for keyname, filename in filename_list:
        filenames.append(filename)
    
    sample = Sample(name, reps, filenames, paired_end)

    return sample

class Sample:
    def __init__(self, name, reps, filenames, paired_end):
        self.reps = reps
        self.name = name
        self.filenames = filenames
        self.paired_end = paired_end
    
    def get_name(self):
        return self.name

    def set_name(self, name):
        self.name = name

    def get_reps(self):
        return self.reps

    def set_reps(self, reps):
        self.reps = reps

    def get_filenames(self):
        return self.filenames

    def set_filenames(self, filenames):
        self.filenames = filenames

    def get_paired_end(self):
        return self.paired_end

class GlobalOptions:
        def __init__(self, genome_location, gtf_location, output_directory, sample_configuration_directory, header, thread_number, STAR_location):
            self.genome_location = genome_location
            self.gtf_location = gtf_location
            self.output_directory = output_directory
            self.sample_configuration_directory = sample_configuration_directory
            self.header = header
            self.thread_number = thread_number
            self.STAR_location = STAR_location

        def get_genome_location(self):
            return self.genome_location
        
        def get_gtf_location(self):
            return self.gtf_location

        def get_output_directory(self):
            return self.output_directory

        def get_sample_configuration_directory(self):
            return self.sample_configuration_directory

        def get_header(self):
            return self.header

        def get_thread_number(self):
            return self.thread_number

        def get_STAR_location(self):
            return self.STAR_location

def get_global_options(configuration_prefix):
    config = configparser.ConfigParser()
    global_config_filepath = configuration_prefix+"global.cfg"
    try:
        config.read_file(open(global_config_filepath))
    except:
        print ("Error: global configuration file not found at",global_config_filepath)
        sys.exit(1)

    genome_location = config.get('global', 'genome_directory')
    gtf_location = config.get('global', 'gtf_location')
    output_directory = config.get('global', 'output_directory')
    sample_configurations_directory = config.get('global', 'sample_configurations_directory')
    
    header_list = config.items('header') 
    header=[]
    for keyname, header_line in header_list:
        header.append(header_line)

    thread_number = int(config.get('global', 'thread_number'))
    STAR_location = config.get('global', 'STAR_location')
    global_options = GlobalOptions(genome_location, gtf_location, output_directory, sample_configurations_directory, header, thread_number, STAR_location)
    return global_options

def create_job(sample, global_options):
    print("Creating job for", sample.get_name())

    f = open(sample.get_name()+'_align.sh', 'w')
    
    # header
    for header_line in global_options.get_header():
        f.write(header_line + '\n')

    # date stamp
    f.write('DATE_STAMP=$(date +%Y%m%d_%H%M)\n')
    
    # run stamp
    f.write('RUN_STAMP=${DATE_STAMP}_'+sample.get_name()+'\n')

    # thread number
    f.write("THREAD_NUM=" + str(global_options.get_thread_number()) + "\n")

    # genome location
    f.write("GENOME_DIR="+global_options.get_genome_location()+'\n')

    # GTF file
    f.write("GTF_FILE="+global_options.get_gtf_location()+'\n')
    # create output directory
    f.write('OUT_DIR='+global_options.get_output_directory()+'/"$RUN_STAMP"/'+'\n' )
    f.write('mkdir -p $OUT_DIR\n')
    f.write('\n')

    # for every rep, output align command
    display_rep = 1
    if sample.get_paired_end():
        for rep in range(0, sample.get_reps()*2,2):
            if rep == (sample.get_reps()*2-2):
                f.write(global_options.get_STAR_location() + ' --runThreadN $THREAD_NUM --genomeDir $GENOME_DIR --readFilesCommand zcat --readFilesIn ' + ' '.join(sample.get_filenames()[rep:rep+2]) + ' --outFileNamePrefix ' + '"$OUT_DIR""$RUN_STAMP"_R' + str(display_rep) + '_' + ' --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --genomeLoad NoSharedMemory --sjdbGTFfile $GTF_FILE' + '\n\n')
            else:
                f.write(global_options.get_STAR_location() + ' --runThreadN $THREAD_NUM --genomeDir $GENOME_DIR --readFilesCommand zcat --readFilesIn ' + ' '.join(sample.get_filenames()[rep:rep+2]) + ' --outFileNamePrefix ' + '"$OUT_DIR""$RUN_STAMP"_R' + str(display_rep) + '_' +  ' --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --genomeLoad NoSharedMemory --sjdbGTFfile $GTF_FILE' +  ' &&' + '\n\n')
            display_rep += 1
    else:
        for rep in range(sample.get_reps()):
            if rep == sample.get_reps()-1:
                f.write(global_options.get_STAR_location() + ' --runThreadN $THREAD_NUM --genomeDir $GENOME_DIR --readFilesCommand zcat --readFilesIn ' + sample.get_filenames()[rep] + ' --outFileNamePrefix ' + '"$OUT_DIR""$RUN_STAMP"_R' + str(display_rep) + '_' +  ' --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --genomeLoad NoSharedMemory --sjdbGTFfile $GTF_FILE' + '\n\n')
            else:
                f.write(global_options.get_STAR_location() + ' --runThreadN $THREAD_NUM --genomeDir $GENOME_DIR --readFilesCommand zcat --readFilesIn ' + sample.get_filenames()[rep] + ' --outFileNamePrefix ' + '"$OUT_DIR""$RUN_STAMP"_R' + str(display_rep) + '_' +  ' --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --genomeLoad NoSharedMemory --sjdbGTFfile $GTF_FILE' +  ' &&' + '\n\n')
            display_rep += 1

    return

def main():
    config_prefix="./configuration_files/"
    global_options = get_global_options(config_prefix)
    sample_config_prefix = global_options.get_sample_configuration_directory()

    sample1 = get_sample_by_name(sample_config_prefix, 'leaf_control')

    create_job(sample1, global_options)

    sample2 = get_sample_by_name(sample_config_prefix, 'leaf_short_drought')
    create_job(sample2, global_options)

if __name__ == '__main__':
    main()
