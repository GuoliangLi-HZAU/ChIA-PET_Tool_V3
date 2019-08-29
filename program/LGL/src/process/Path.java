package process;

import java.io.File;

public class Path {
	
	String MODE;
	String Fastq_file_1;
	String Fastq_file_2;
	String linker;
	String minimum_linker_alignment_score;
	String GENOME_INDEX;
	String GENOME_LENGTH;
	String CHROM_SIZE_INFO;
	String CYTOBAND_DATA;
	String SPECIES;
	
	String PROGRAM_DIRECTORY;
	
	String START_STEP = "1";
	String OUTPUT_DIRECTORY;
	String OUTPUT_PREFIX = "out";
	String minimum_tag_length = "18";
	String maximum_tag_length = "1000";
	String minSecondBestScoreDiff = "3";
	String output_data_with_ambiguous_linker_info = "1";
	String NTHREADS = "1";
	String MAPPING_CUTOFF = "30";
	String MERGE_DISTANCE = "2";
	String SELF_LIGATION_CUFOFF = "8000";
	String EXTENSION_LENGTH = "500";
	String MIN_COVERAGE_FOR_PEAK = "5";
	String PEAK_MODE = "2";
	String MIN_DISTANCE_BETWEEN_PEAK = "500";
	String GENOME_COVERAGE_RATIO = "0.8";
	String PVALUE_CUTOFF_PEAK = "0.00001";
	String INPUT_ANCHOR_FILE = "null";
	String PVALUE_CUTOFF_INTERACTION = "0.05";
	
	/**
	 * @author tksun
	 * @function set parameters
	 * @param args parameter array
	 */
	public void setParameter(String[] args) {
		PROGRAM_DIRECTORY = new Path().getClass().getProtectionDomain().getCodeSource().getLocation().getPath();
		if(PROGRAM_DIRECTORY.endsWith(".jar")) {
			PROGRAM_DIRECTORY = PROGRAM_DIRECTORY.substring(0, PROGRAM_DIRECTORY.lastIndexOf("/") + 1);
        }
		OUTPUT_DIRECTORY = PROGRAM_DIRECTORY;
		int necessary = 0;
		for (int i = 0; i < args.length; i += 2) {
			if (args[i].equals("--start_step")) {
				START_STEP = args[i + 1];
			} else if (args[i].equals("--output")) {
				OUTPUT_DIRECTORY = args[i + 1];
			} else if (args[i].equals("--prefix")) {
				OUTPUT_PREFIX = args[i + 1];
			} else if (args[i].equals("--mode")) {
				MODE = args[i + 1];
				if (!(args[i + 1].equals("0") || args[i + 1].equals("1"))) {
					System.out.println("Error: mode " + args[i + 1] + " is incorrect");
					System.exit(0);
				}
				necessary++;
			} else if (args[i].equals("--fastq1")) {
				Fastq_file_1 = args[i + 1];
				checkPath(args[i + 1]);
				necessary++;
			} else if (args[i].equals("--fastq2")) {
				Fastq_file_2 = args[i + 1];
				checkPath(args[i + 1]);
				necessary++;
			} else if (args[i].equals("--linker")) {
				linker = args[i + 1];
				checkPath(args[i + 1]);
				necessary++;
			} else if (args[i].equals("--minimum_linker_alignment_score")) {
				minimum_linker_alignment_score = args[i + 1];
				necessary++;
			} else if (args[i].equals("--minimum_tag_length")) {
				minimum_tag_length = args[i + 1];
			} else if (args[i].equals("--maximum_tag_length")) {
				maximum_tag_length = args[i + 1];
			} else if (args[i].equals("--minSecondBestScoreDiff")) {
				minSecondBestScoreDiff = args[i + 1];
			} else if (args[i].equals("--output_data_with_ambiguous_linker_info")) {
				output_data_with_ambiguous_linker_info = args[i + 1];
				if (!(args[i + 1].equals("0") || args[i + 1].equals("1"))) {
					System.out.println("Error: output_data_with_ambiguous_linker_info " + args[i + 1] + " is incorrect");
					System.exit(0);
				}
			} else if (args[i].equals("--thread")) {
				NTHREADS = args[i + 1];
			} else if (args[i].equals("--GENOME_INDEX")) {
				GENOME_INDEX = args[i + 1];
				String prefix = new File(args[i + 1]).getName();
				File parentDir = new File(args[i + 1]).getParentFile();
				File[] files = parentDir.listFiles();
				int n = 0;
				for (File f : files) {
					if (f.getName().equals(prefix + ".amb") || f.getName().equals(prefix + ".ann") || f.getName().equals(prefix + ".bwt") || 
							f.getName().equals(prefix + ".pac") || f.getName().equals(prefix + ".sa")) {
						n++;
					}
				}
				if (n != 5) {
					System.out.println("Error: please check genome index (*.amb, *.ann, *.bwt, *pac and *.sa, these files are necessary)");
					System.exit(0);
				}
				necessary++;
			} else if (args[i].equals("--MAPPING_CUTOFF")) {
				MAPPING_CUTOFF = args[i + 1];
			} else if (args[i].equals("--MERGE_DISTANCE")) {
				MERGE_DISTANCE = args[i + 1];
			} else if (args[i].equals("--SELF_LIGATION_CUFOFF")) {
				SELF_LIGATION_CUFOFF = args[i + 1];
			} else if (args[i].equals("--EXTENSION_LENGTH")) {
				EXTENSION_LENGTH = args[i + 1];
			} else if (args[i].equals("--MIN_COVERAGE_FOR_PEAK")) {
				MIN_COVERAGE_FOR_PEAK = args[i + 1];
			} else if (args[i].equals("--PEAK_MODE")) {
				PEAK_MODE = args[i + 1];
				if (!(args[i + 1].equals("1") || args[i + 1].equals("2"))) {
					System.out.println("Error: PEAK_MODE " + args[i + 1] + " is incorrect");
					System.exit(0);
				}
			} else if (args[i].equals("--MIN_DISTANCE_BETWEEN_PEAK")) {
				MIN_DISTANCE_BETWEEN_PEAK = args[i + 1];
			} else if (args[i].equals("--GENOME_LENGTH")) {
				GENOME_LENGTH = args[i + 1];
				necessary++;
			} else if (args[i].equals("--GENOME_COVERAGE_RATIO")) {
				GENOME_COVERAGE_RATIO = args[i + 1];
			} else if (args[i].equals("--PVALUE_CUTOFF_PEAK")) {
				PVALUE_CUTOFF_PEAK = args[i + 1];
			} else if (args[i].equals("--CHROM_SIZE_INFO")) {
				CHROM_SIZE_INFO = args[i + 1];
				checkPath(args[i + 1]);
				necessary++;
			} else if (args[i].equals("--INPUT_ANCHOR_FILE")) {
				INPUT_ANCHOR_FILE = args[i + 1];
			} else if (args[i].equals("--PVALUE_CUTOFF_INTERACTION")) {
				PVALUE_CUTOFF_INTERACTION = args[i + 1];
			} else if (args[i].equals("--CYTOBAND_DATA")) {
				CYTOBAND_DATA = args[i + 1];
				checkPath(args[i + 1]);
				necessary++;
			} else if (args[i].equals("--SPECIES")) {
				SPECIES = args[i + 1];
				if (!(args[i + 1].equals("1") || args[i + 1].equals("2") || args[i + 1].equals("3"))) {
					System.out.println("Error: species " + args[i + 1] + " is incorrect");
					System.exit(0);
				}
				necessary++;
			}
		}
		if (necessary < 10) {
			System.out.println("Error: necessary options are not enough");
			System.exit(0);
		}
	}
	
	/**
	 * @author sun
	 * @function check the path of file
	 * @param str the path of file
	 */
	public void checkPath(String str) {
		File f = new File(str);
		if (!f.exists()) {
			System.out.println("Error: " + str + " doesn't exist");
			System.exit(0);
		}
	}

	public String getPROGRAM_DIRECTORY() {
		return PROGRAM_DIRECTORY;
	}
	
	public String getOUTPUT_DIRECTORY() {
		return OUTPUT_DIRECTORY;
	}

	public String getOUTPUT_PREFIX() {
		return OUTPUT_PREFIX;
	}
	
	public String getFastq_file_1() {
		return Fastq_file_1;
	}

	public String getFastq_file_2() {
		return Fastq_file_2;
	}

	public String getLinker() {
		return linker;
	}
	
	public String getMinimum_linker_alignment_score() {
		return minimum_linker_alignment_score;
	}

	public String getMinimum_tag_length() {
		return minimum_tag_length;
	}

	public String getMaximum_tag_length() {
		return maximum_tag_length;
	}

	public String getMinSecondBestScoreDiff() {
		return minSecondBestScoreDiff;
	}

	public String getOutput_data_with_ambiguous_linker_info() {
		return output_data_with_ambiguous_linker_info;
	}

	public String getNTHREADS() {
		return NTHREADS;
	}

	public String getMAPPING_CUTOFF() {
		return MAPPING_CUTOFF;
	}

	public String getMERGE_DISTANCE() {
		return MERGE_DISTANCE;
	}
}