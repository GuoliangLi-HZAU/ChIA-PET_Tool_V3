package process;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.Calendar;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class Path {
	
	String MODE="1";
	String Fastq_file_1="";
	String Fastq_file_2="";
	String linker;
	String minimum_linker_alignment_score="14";
	String GENOME_INDEX;
	String GENOME_LENGTH;
	String CHROM_SIZE_INFO;
	String CYTOBAND_DATA;
	String SPECIES;
	
	String PROGRAM_DIRECTORY;
	
	String START_STEP = "1";
	String STOP_STEP = "100";
	String OUTPUT_DIRECTORY="./";
	String OUTPUT_PREFIX = "out";
	String minimum_tag_length = "18";
	String maximum_tag_length = "1000";
	String minSecondBestScoreDiff = "3";
	String output_data_with_ambiguous_linker_info = "1";
	String NTHREADS = "1";
	String MAPPING_CUTOFF = "20";
	String MERGE_DISTANCE = "2";
	String SELF_LIGATION_CUFOFF = "8000";
	String provide_slc = "N";
	String EXTENSION_LENGTH = "500";
	String MIN_COVERAGE_FOR_PEAK = "5";
	String PEAK_MODE = "2";
	String MIN_DISTANCE_BETWEEN_PEAK = "500";
	String GENOME_COVERAGE_RATIO = "0.8";
	String PVALUE_CUTOFF_PEAK = "0.00001";
	String INPUT_ANCHOR_FILE = "null";
	String macs2 = "N";
	String PVALUE_CUTOFF_INTERACTION = "0.05";
	String FQMODE = "paired-end";
	String XOR_cluster = "N";
	public String ALLMAP = "false";
	public String MAP2Linker = "false";
	String MAPMEM = "false";
	String printreadID = "N";
	public String search_all_linker = "N";
	public int printallreads = 0; //0 print all, 1 print known +-
	String restrictionsiteFile = "None";
	//public boolean hichipmode = false;
	String hichipM = "N";
	String keeptemp = "N";
	String ligation_site = "-";
	String[] ligation_sites;
	String genomefile = "";
	int Ngenome = 0;
	HashMap<String, Integer> chrMAP = new HashMap<>();
	HashMap<Integer, String> chrMAP_r = new HashMap<>();
	int minfragsize = 100;
	int maxfragsize = 1000000;
	int minInsertsize = 1;
	int maxInsertsize = 1000;
	String zipbedpe = "N";
	String zipsam = "N";
	String deletesam = "N";
	String fastp = "";
	String skipmap = "N";
	String MAP_ambiguous = "N";
	public String AutoLinker = "true";
	public String kmerlen = "9";
	public String removeResblock = "Y";
	//HindIII，则为AAGCTAGCTT；如果是MboI则序列为GATCGATC; dpnII GATCGATC
	/*
	 * "mboi": ["^GATC"],
     * "dpnii": ["^GATC"],
     * "bglii": ["A^GATCT"],
     * "hindiii": ["A^AGCTT"]
     * Hinf1 G^ANTC
     * NlaIII CATG^
     * AluI AG^CT
	 * */
	
	public String getLigationSite(String ligation_site) {
		String site = "";
		if(ligation_site.equalsIgnoreCase("HindIII")) {
			site = "A^AGCTT";
		}else if(ligation_site.equalsIgnoreCase("MboI")) {
			site = "^GATC";
		}else if(ligation_site.equalsIgnoreCase("BglII")) {
			site = "A^GATCT";
		}else if(ligation_site.equalsIgnoreCase("DpnII")) {
			site = "^GATC";
		}else if(ligation_site.equalsIgnoreCase("Sau3AI")) {
			site = "^GATC";
		}else if(ligation_site.equalsIgnoreCase("NlaIII")) {
			site = "CATG^";
		}else if(ligation_site.equalsIgnoreCase("Hinf1")) {
			site = "G^ANTC";
		}else if(ligation_site.equalsIgnoreCase("AluI")) {
			site = "AG^CT";
		}else if(ligation_site.contains("^")) {
			System.out.println("Error: the restriction position has to be specified using '^'\n"
					+ "Please, use '^' to specify the cutting position\n"
					+ "i.e A^GATCT for HindIII digestion.\n");
			System.exit(0);
		}else {
			site = ligation_site;
		}
		for(int j=0; j<site.length(); j++) {
			if(!(site.charAt(j) == 'C' || site.charAt(j) == 'c' || site.charAt(j) == 'T' 
					|| site.charAt(j) == 't' || site.charAt(j) == 'G' || site.charAt(j) == 'g'
					|| site.charAt(j) == 'A' || site.charAt(j) == 'a' || site.charAt(j) == '^' || 
					site.charAt(j) == 'N' || site.charAt(j) == 'n') ) {
				System.out.println("Error: \n"
						+ "Please print HindIII/MboI/BglII/DpnII or "
						+ "restriction site with '^' and contains 'ATCG' without other character!!! "
						+ site + " : " + site.charAt(j));
				System.exit(0);
			}
		}
		return site;
	}
	public String[] replaceN(String[] ligs) {
		//String[] ligs = ligation.split(",");
		int Nlig = 0;
		for(int i=0; i< ligs.length; i++) {
			if(ligs[i].contains("N") || ligs[i].contains("n")) {
				Nlig+=3;
			}
		}
		String[] ligs2 = new String[Nlig+ligs.length];
		if(Nlig>0) {
			int k=0;
			for(int i=0; i< ligs.length; i++) {
				if(ligs[i].contains("N")) {
					ligs2[k] = ligs[i].replace('N', 'A');
					ligs2[k+1] = ligs[i].replace('N', 'T');
					ligs2[k+2] = ligs[i].replace('N', 'C');
					ligs2[k+3] = ligs[i].replace('N', 'G');
					k=k+4;
				}else if(ligs[i].contains("n")) {
					ligs2[k] = ligs[i].replace('n', 'A');
					ligs2[k+1] = ligs[i].replace('n', 'T');
					ligs2[k+2] = ligs[i].replace('n', 'C');
					ligs2[k+3] = ligs[i].replace('n', 'G');
					k=k+4;
				}else {
					ligs2[k] = ligs[i];
					k++;
				}
			}
		}else {
			ligs2 = ligs;
		}
		for(String str:ligs2) {
			System.out.println("[Enzyme] " + str);
		}

        List mylist = Arrays.asList(ligs2);
        Set myset = new HashSet(mylist);
        String [] ligs3=(String [])myset.toArray(new String[0]);
		return ligs3;
	}
	
	/**
	 * @author tksun
	 * @function set parameters
	 * @param args parameter array
	 * @throws IOException 
	 */
	public void setParameter(String[] args) throws IOException {
		PROGRAM_DIRECTORY = new Path().getClass().getProtectionDomain().getCodeSource().getLocation().getPath();
		if(PROGRAM_DIRECTORY.endsWith(".jar")) {
			PROGRAM_DIRECTORY = PROGRAM_DIRECTORY.substring(0, PROGRAM_DIRECTORY.lastIndexOf("/") + 1);
        }
		OUTPUT_DIRECTORY = PROGRAM_DIRECTORY;
		int necessary = 0;
		for (int i = 0; i < args.length; i += 2) {
			if (args[i].equals("--start_step")) {
				START_STEP = args[i + 1];
			} else if (args[i].equals("--stop_step")) {
				STOP_STEP = args[i + 1];
			} else if (args[i].equals("--map_ambiguous")) {
				MAP_ambiguous = args[i + 1];
			} else if (args[i].equals("--skipmap")) {
				skipmap = args[i + 1];
			} else if (args[i].equals("--printreadID")) {
				printreadID = args[i + 1];
			} else if (args[i].equals("--zipbedpe")) {
				zipbedpe = args[i + 1];
			} else if (args[i].equals("--zipsam")) {
				zipsam = args[i + 1];
			} else if (args[i].equals("--deletesam")) {
				deletesam = args[i + 1];
			} else if (args[i].equals("--printallreads")) {
				printallreads =  Integer.parseInt(args[i + 1]);
			} else if (args[i].equals("--search_all_linker")) {
				search_all_linker =  args[i + 1];
			} else if (args[i].equals("--XOR_cluster")) {
				XOR_cluster =  args[i + 1];
			} else if (args[i].equals("--hichip")) {
				hichipM = args[i + 1];
				if(provide_slc.equalsIgnoreCase("N")) {
					SELF_LIGATION_CUFOFF = "1000";
				}
			} else if (args[i].equals("--keeptemp")) {
				keeptemp = args[i + 1];
			} else if (args[i].equals("--ResRomove")) {
				removeResblock = args[i + 1];
			}else if (args[i].equals("--restrictionsiteFile")) {
				restrictionsiteFile = args[i + 1];
			} else if (args[i].equals("--ligation_site")) {
				ligation_site = args[i + 1];
			    String[] ligations = ligation_site.split(",");
				String[] ligation_sites1 = ligation_site.split(",");
				for(int k=0; k<ligations.length; k++) {
					ligation_sites1[k] = getLigationSite(ligations[k]);
				}
				ligation_sites = replaceN(ligation_sites1);
				System.out.println("[Enzyme site number] " + ligation_sites.length);
				//System.out.println("XXXX " + ligation_sites.toString());
			} else if (args[i].equals("--genomefile")) {
				genomefile = args[i + 1];
			} else if (args[i].equals("--minfragsize")) {
				minfragsize = Integer.parseInt(args[i + 1]);
			} else if (args[i].equals("--maxfragsize")) {
				maxfragsize = Integer.parseInt(args[i + 1]);
			} else if (args[i].equals("--minInsertsize")) {
				minInsertsize = Integer.parseInt(args[i + 1]);
			} else if (args[i].equals("--maxInsertsize")) {
				maxInsertsize = Integer.parseInt(args[i + 1]);
			} else if (args[i].equals("--output")) {
				OUTPUT_DIRECTORY = args[i + 1];
			} else if(args[i].equals("--mapall")) {
				ALLMAP = args[i + 1];
			} else if(args[i].equals("--map2linker")) {
				MAP2Linker = args[i + 1];
			} else if(args[i].equals("--mapmem")) {
				MAPMEM = args[i + 1];
			} else if (args[i].equals("--prefix")) {
				OUTPUT_PREFIX = args[i + 1];
			} else if (args[i].equals("--mode")) {
				MODE = args[i + 1];
				if (!(args[i + 1].equals("0") || args[i + 1].equals("1"))) {
					System.out.println("Error: mode " + args[i + 1] + " is incorrect");
					System.exit(0);
				}
				necessary++;
			} else if (args[i].equals("--kmerlen")) {
				kmerlen = args[i + 1];
			}else if (args[i].equals("--fastp")) {
				fastp = args[i + 1];
				if(!fastp.equals("fastp")) {
					checkPath(args[i + 1]);
				}
			} else if (args[i].equals("--fastq1")) {
				Fastq_file_1 = args[i + 1];
				//checkPath(args[i + 1]);
				necessary++;
			} else if (args[i].equals("--fastq2")) {
				Fastq_file_2 = args[i + 1];
				//checkPath(args[i + 1]);
				necessary++;
			} else if (args[i].equals("--fqmode")) {
				FQMODE = args[i + 1];
				necessary++;
			}else if (args[i].equals("--autolinker")) {
				AutoLinker = args[i + 1]; //true false
				necessary++;
			} else if (args[i].equals("--linker")) {
				linker = args[i + 1];
				//checkPath(args[i + 1]);
				AutoLinker="false";
				necessary++;
			} else if (args[i].equals("--minimum_linker_alignment_score")) {
				minimum_linker_alignment_score = args[i + 1];
				necessary++;
			} else if (args[i].equals("--minimum_tag_length")) {
				minimum_tag_length = args[i + 1];
				if(Integer.parseInt(minimum_tag_length)==0) {
					minimum_tag_length = "1";
				}
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
				provide_slc = "Y";
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
				// for buffer write
				BufferedReader reader = new BufferedReader(new FileReader(CHROM_SIZE_INFO));
    			String line = reader.readLine();
    			while(line != null) {
    				chrMAP.put(line.split("[ \t]+")[0], Ngenome);
    				chrMAP_r.put(Ngenome, line.split("[ \t]+")[0]);
    				Ngenome++;
    				line = reader.readLine();
    			}
			} else if (args[i].equals("--INPUT_ANCHOR_FILE")) {
				INPUT_ANCHOR_FILE = args[i + 1];
			} else if (args[i].equals("--macs2")) {
				macs2 = args[i + 1];
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
			} else {
				System.out.println("Error: unexpected paramater: " + args[i]);
				System.exit(0);
			}
		}
		if (necessary < 5 && hichipM.equals("N")) {
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
