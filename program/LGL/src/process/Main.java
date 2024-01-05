package process;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.zip.GZIPInputStream;

import LGL.align.LocalAlignment;
import LGL.data.FastQ;
import multhread.BigFileProcess;

public class Main {
	public static boolean isGZipped(File f) {
		int magic = 0;
	
		try {
			RandomAccessFile raf = new RandomAccessFile(f, "r");
			magic = raf.read() & 0xff | ((raf.read() << 8) & 0xff00);
			raf.close();
		} catch (Throwable e) {
			e.printStackTrace(System.err);
		}
		
		return magic == GZIPInputStream.GZIP_MAGIC;
	}
	public static void writeFile(String f, String s, boolean append) {
		File file = new File(f);
		try {
	        if(!file.exists()){
			    file.createNewFile();
	        }
		    FileWriter fw = new FileWriter(file, append);
		    fw.write(s+"\n");
		    fw.close();
		} catch (IOException e) {
		    e.printStackTrace();  
		}
	}
	
	public static void shellrun(String command, String runcase) {
		Process process = null;
		
		try {
			process = Runtime.getRuntime().exec(command);
			BufferedReader reader = new BufferedReader(
			 			new InputStreamReader(
			 					process.getInputStream()));
			String data = "";
			while((data = reader.readLine()) != null) {
				System.out.println(data);
			}
			
			int exitValue = process.waitFor();
			
			if(exitValue != 0 && !runcase.equals("rm")) {
				System.out.println("error");
			}
		} catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	public static int findstartinres(String res) {
		int i = -1;
		for(i=0; i<res.length(); i++) {
			if(res.charAt(i)=='^') {
				return i;
			}
		}
		return i;
	}
	
	public static int[] findstartinres_arry(String[] resS) {
		int[] sites = new int[resS.length];
		for(int i=0; i<resS.length; i++) {
			sites[i] = findstartinres(resS[i]);
			if(sites[i] == -1) {
				System.out.println("Error: can not find '^' in ligation site!!!!\n");
				System.exit(0);
			}
		}
		return sites;
	}
	
	public static int[] findsmIndex(StringBuilder chrseqbuilder, String[] resS, int start, int[] resloci, int chrlen) {
		int newindex = chrlen+1, residx = -1;
		
		for(int i=0;i<resS.length; i++) {
			if(resloci[i]<=start) {
			    resloci[i] = chrseqbuilder.indexOf(resS[i], start);
			}
		}
		
		for(int i=0; i<resloci.length; i++) {
		    if(resloci[i]>=0 && resloci[i] < newindex) {
		    	newindex = resloci[i];
		    	residx = i;
		    }
		}
		int[] inarray = new int[2];
		inarray[0] = newindex;
		if(residx>=0) {
			resloci[residx] = chrseqbuilder.indexOf(resS[residx], start);
		}else {
			residx=0;
		}
		inarray[1] = residx;
		return inarray;
	}
	
	@SuppressWarnings("null")
	public static void split_genome_byressite(String genomefile, String[] resS, String outfile, String ressite) throws IOException {
	    //replace and find postion of '^' in res
		int[] posoff = findstartinres_arry(resS);
		int[] reslen = new int[resS.length];
		for(int i=0; i<resS.length; i++) {
		    resS[i] = resS[i].replace("^", "").toUpperCase();
		    reslen[i] = resS[i].length();
		}
		
		Calendar rightNow = Calendar.getInstance();
		
		//out res file point
		BufferedWriter resBufferedWriter = new BufferedWriter(new FileWriter(outfile));
		// for genoem file
		BufferedReader reader = new BufferedReader(new FileReader(genomefile));
		String line = reader.readLine();
		String chrom = "";
		//String chrseq = "";
		StringBuilder chrseqbuilder = new StringBuilder();
		String newline = "";
		//int index = -1;
		int chrlen = 0;
		int pstart = 0, pend = 0;
		int Nblock = 0;
		int[] index_residx = new int[2];
		index_residx[0] = -1;
		int[] resloci = new int[resS.length];
		while( line != null) {
			if(line.startsWith(">")) {
				if(!chrom.equals("") && !chrseqbuilder.equals("")) {
					//System.out.println("111 "+ index_residx[0] + " " + index_residx[1] + " " + resloci[0]);
					chrlen = chrseqbuilder.length();
					index_residx = findsmIndex(chrseqbuilder, resS, 0, resloci,chrlen); //[index, reslen]
					while(index_residx[0]>=0 && index_residx[0]<chrlen-1) {
						//print res site in outfile
						if(index_residx[0] + posoff[index_residx[1]]>chrlen) {
							pend = chrlen;
						}else {
							pend = index_residx[0]+posoff[index_residx[1]];
						}
						if(index_residx[0]>0) {
							newline = chrom + "\t" + pstart + "\t" + pend;
							pstart = pend;
							//System.out.println(newline);
							resBufferedWriter.write(newline);
							resBufferedWriter.newLine();
							Nblock ++;
						}
						
						//index = chrseqbuilder.indexOf(res, index+reslen);
						index_residx = findsmIndex(chrseqbuilder, resS, index_residx[0] + 1, resloci,chrlen); //reslen[index_residx[1]] cause many site for multipe eny
						//System.out.println("xxx "+ index_residx[0] + " " + index_residx[1] + " " + pend);
					}
					// last block in this chrom
					if(pend<chrlen) { //pend!=0 && pend ==0 means no res in chr
						//pend++; // not add 1
						newline = chrom + "\t" + pend + "\t" + chrlen;
						//System.out.println(newline);
						resBufferedWriter.write(newline);
						resBufferedWriter.newLine();
						Nblock ++;
					}
					//init
					pstart = 0;
					rightNow = Calendar.getInstance();
					System.out.println("[" + rightNow.getTime().toString() +"] Split " + chrom + " to " + Nblock + " blocks by " + ressite);
					Nblock = 0;
				}
				
			    //chrom = line.substring(1).split(" ")[0].split("\t")[0];
				chrom = line.substring(1).split("[ \t]+")[0];
			    //System.out.println("Spliting genome file, "
			    //		+ chrom);
			    //init chrseq
			    chrseqbuilder = chrseqbuilder.delete( 0, chrlen);
			}else {
				// whether need replace '\r\n' in java
				chrseqbuilder.append(line.toUpperCase());
			}
			line = reader.readLine();
		}
		// print last chrom
		chrlen = chrseqbuilder.length();
		index_residx = findsmIndex(chrseqbuilder, resS, 0, resloci, chrlen);
		Nblock=0;
		while(index_residx[0]>=0 && index_residx[0]<chrlen-1) {
			//print res site in outfile
			if(index_residx[0] + posoff[index_residx[1]]>chrlen) {
				pend = chrlen;
			}else {
				pend = index_residx[0]+posoff[index_residx[1]];
			}
			if(index_residx[0]>0) {
				newline = chrom + "\t" + pstart + "\t" + pend;
				pstart = pend;
				resBufferedWriter.write(newline);
				resBufferedWriter.newLine();
				Nblock ++;
			}
			//index = chrseqbuilder.indexOf(res, index+);
			index_residx = findsmIndex(chrseqbuilder, resS, index_residx[0] + 1, resloci, chrlen); //reslen[index_residx[1]]
		}
		if(pend<chrlen) { //pend!=0 && 
			//pend++;
			newline = chrom + "\t" + pend + "\t" + chrlen;
			resBufferedWriter.write(newline);
			resBufferedWriter.newLine();
			Nblock ++;
		}
		rightNow = Calendar.getInstance();
		System.out.println("[" + rightNow.getTime().toString() +"] Split " + chrom + " to " + Nblock + " blocks splited by " + ressite
				+ "\n[" + rightNow.getTime().toString() +"] Split genome done!");
		//init
		chrseqbuilder = chrseqbuilder.delete( 0, chrlen);
		resBufferedWriter.close();
	}
	
	public static void checkPath(String str) {
		File f = new File(str);
		if (!f.exists()) {
			System.out.println("Error: " + str + " doesn't exist");
			System.exit(0);
		}
	}
	
    public static String[] CleanFastq(String fastq1, String fastq2, String Nthread, String outdir, String prefix) throws IOException {
    	String outPrefix = outdir + "/" + prefix + "/" + prefix;
		String line = "";
		BufferedWriter cleanfqF = new BufferedWriter(new FileWriter(outPrefix + 
				".cleanfastp.sh", false));
	    
		String newfq1="", newfq2="";
    	String[] fastq1s = fastq1.split(",");
    	if(!fastq2.equals("")) {
    		String[] fastq2s = fastq2.split(",");
    		if(fastq1s.length!=fastq2s.length) {
        		System.out.println("Error: fastq1 and fastq2 with different number of files!!!");
        		System.exit(0);
        	}
    		for(int i=0;i<fastq1s.length;i++) {
    			line="fastp -M 10 -Q -w "+ Nthread + " -i " + fastq1s[i] + " -o " + outPrefix+".clean." +i+"_R1.fq.gz " + " -I " + fastq2s[i] + " -O "+outPrefix+".clean."+
    		       i+"_R2.fq.gz";
    			cleanfqF.write(line);
    			cleanfqF.newLine();
    			if(i==0) {
    				newfq1=outPrefix+".clean." +i+"_R1.fq.gz";
    				newfq2=outPrefix+".clean." +i+"_R2.fq.gz";
    			}else {
    				newfq1=newfq1+","+outPrefix+".clean." +i+"_R1.fq.gz";
    				newfq2=newfq2+","+outPrefix+".clean." +i+"_R2.fq.gz";
    			}
    		}
    	}else {
    		for(int i=0;i<fastq1s.length;i++) {
    			line="fastp -M 10 -Q -w "+ Nthread + " -i " + fastq1s[i] + " -o "+outPrefix+".clean." +i+"_R1.fq.gz ";
    			cleanfqF.write(line);
    			cleanfqF.newLine();
    			if(i==0) {
    				newfq1=outPrefix+".clean." +i+"_R1.fq.gz";
    			}else {
    				newfq1=newfq1+","+outPrefix+".clean." +i+"_R1.fq.gz";
    			}
    		}
    	}
    	
    	cleanfqF.close();
    	Shell shell = new Shell();
		shell.runShell(outPrefix + ".cleanfastp.sh");
		new File(outPrefix + ".cleanfastp.sh").delete();
		String[] retString = new String[] {newfq1, newfq2};
		return retString;
    	
    }
	
	public static void main(String []args) throws IOException, InterruptedException {
		if (args.length < 5) {
			System.out.println("Error: please set the necessary parameters");
			System.out.println("Usage: java -jar <path of ChIA_PET.jar> [options]");
			System.out.println("Necessary options:");
			System.out.println("    --fastq1\tpath of read1 fastq file");
			System.out.println("    --fastq2\tpath of read2 fastq file");
			System.out.println("    --autolinker\tdetect linker by our program, true [default] or false, then no need provide --linker and --mode paramater. default: true");
			System.out.println("      When the parameter --stop_step 0 is present, only the automatic detection linker program will be run.");
			System.out.println("      When the parameter --stop_step 1 is present, only the automatic detection linker AND linker filter program will be run.");
			System.out.println("    --mode\tmode of tool, 0: short read; 1: long read, need for ChIA-PET data");
			System.out.println("    --linker\tpath of linker file, need for ChIA-PET mode");
			System.out.println("    --fastp\tfastp path, strong suggest for ChIA-PET data.");
			System.out.println("    --skipheader\tskip header N reads for detect linker, default 1000000.");
			System.out.println("    --linkerreads\tN reads used for detect linker, default 100000.");
			
			System.out.println("    --hichip\tY(es) or N(o)[default] or O(nly print restriction site file without run other step), need for hichip data");
			System.out.println("    --ligation_site\tIt can be the name of restriction enzyme, such as HindIII, MboI, DpnII, Bglii, Sau3AI, Hinf1, NlaIII, AluI "
					+ "\n                or the site of enzyme digestion, A^AGCTT, ^GATC, ^GATC, A^GATCT, G^ANTC, CATG^, AG^CT or others."
					+ "\n                multipe restriction enzyme can be seperated by comma, such as G^ANTC,^GATC."
					+ "\n                restriction site with '^' and contains 'ATCG' without other character!!! "
					+ "\n                if the genomic enzyme digestion file --restrictionsiteFile is provided,"
					+ "\n                this parameter does not need to be provided. "
					+ "\n                only needed for hichip data");
			System.out.println("    --ResRomove\tY or N, whether remove PET in same restriction contig. default: Y");
			System.out.println("    --restrictionsiteFile\trestriction site file, can be genarated while has --ligation_site and without this paramater or "
					+ "\n                provide restriction enzyme information with --ligation_site, we will automatically generate the file."
					+ "\n                only needed for hichip data");
			System.out.println("    --genomefile\tgenome fasta file path, needed for with --ligation_site and without --restrictionsiteFile"
					+ "\n                only needed for hichip data");
			System.out.println("    --minfragsize\tMinimum restriction fragment length to consider, default 20"); //100
			System.out.println("    --maxfragsize\tMaximum restriction fragment length to consider, default 1000000");
			System.out.println("    --minInsertsize\tMinimum restriction fragment skip of mapped reads to consider, default 1");
			//System.out.println("    --maxInsertsize\tMaximum restriction fragment skip of mapped reads to consider, default 1000");
			
			System.out.println("    --fqmode\tsingle-end or paired-end (default), only required --fastq1 when single-end mode for ChIA-PET data");
			System.out.println("    --minimum_linker_alignment_score\tminimum alignment score");
			System.out.println("    --GENOME_INDEX\tthe path of BWA index");
			System.out.println("    --GENOME_LENGTH\tthe number of base pairs in the whole genome");
			System.out.println("    --CHROM_SIZE_INFO\tthe file that contains the length of each chromosome, example file is in ChIA-PET_Tool_V3/chrInfo,"
					         + "\n                     \tthis is necessary for > step 2 analysis."
					         + "\n                     \tNote. please make sure chromosome name in this file is same as name in genome file!!!");
			System.out.println("    --CYTOBAND_DATA\tthe ideogram data used to plot intra-chromosomal peaks and interactions, example file is in "
					+ "ChIA-PET_Tool_V3/chrInfo");
			System.out.println("    --SPECIES\t1: human; 2: mouse; 3: others");
			System.out.println("Other options:");
			System.out.println("    --start_step\tstart with which step, 1: linker filtering; 2: mapping to genome; 3: removing redundancy; 4: categorization of "
					+ "PETs; 5: peak calling; 6: interaction calling; 7: visualizing, default: 1");
			System.out.println("    --stop_step\tstop with which step, 1: linker filtering; 2: mapping to genome; 3: removing redundancy; 4: categorization of "
					+ "PETs; 5: peak calling; 6: interaction calling; 7: visualizing, default: 100, should be bigger than --start_step");
			System.out.println("    --output\tpath of output, default: ChIA-PET_Tool_V3/output");
			System.out.println("    --prefix\tprefix of output files, default: out");
			System.out.println("    --minimum_tag_length\t minimum tag length, default: 18");
			System.out.println("    --maximum_tag_length\t maximum tag length, default: 1000");
			System.out.println("    --minSecondBestScoreDiff\tthe score difference between the best-aligned and the second-best aligned linkers, default: 3");
			System.out.println("    --output_data_with_ambiguous_linker_info\twhether to print the linker-ambiguity PETs, 0: not print; 1: print, default: 1");
			//add new
			System.out.println("    --printreadID\t write read ID to bedpe file, default: N");
			System.out.println("    --printallreads\t print all reads no matter strand, default: 0[print all]; 1, only print valid strand reads.");
			System.out.println("    --search_all_linker\t search all linkers in reads or just search one time, default: N");
			
			System.out.println("    --thread\tthe number of threads used in linker filtering and mapping to genome, default: 1");
			System.out.println("    --MAPPING_CUTOFF\tcutoff of mapping quality score for filtering out low-quality or multiply-mapped reads, default: 20");
			System.out.println("    --MERGE_DISTANCE\tthe distance limit to merge the PETs with similar mapping locations, default: 2");
			System.out.println("    --SELF_LIGATION_CUFOFF\tthe distance threshold between self-ligation PETs and intra- chromosomal inter-ligation PETs, "
					+ "default: 8000 for ChIA, and 1000 for HiChIP");
			System.out.println("    --EXTENSION_LENGTH\tthe extension length from the location of each tag, default: 500, 1500 suggest for single-end mode");
			System.out.println("    --MIN_COVERAGE_FOR_PEAK\tthe minimum coverage to define peak regions, default: 5");
			System.out.println("    --PEAK_MODE\t1: peak region mode, which takes all the overlapping PET regions above the coverage threshold as peak "
					+ "regions; 2: peak summit mode, which takes the highest coverage of overlapping regions as peak regions, default: 2");
			System.out.println("    --MIN_DISTANCE_BETWEEN_PEAK\tthe minimum distance between two peaks, default: 500");
			System.out.println("    --GENOME_COVERAGE_RATIO\tthe estimated proportion of the genome covered by the reads, default: 0.8");
			System.out.println("    --PVALUE_CUTOFF_PEAK\tp-value to filter peaks that are not statistically significant, default: 0.00001");
			System.out.println("    --INPUT_ANCHOR_FILE\ta file which contains user-specified anchors for interaction calling, default: null");
			System.out.println("    --PVALUE_CUTOFF_INTERACTION\tp-value to filter false positive interactions, default: 0.5");
			System.out.println("    --zipbedpe\tgzip bedpe related file, after analysis done. default: N. Y for gzip, N for not.");
			System.out.println("    --zipsam\tConvert sam file to bam, after analysis done. default: N");
			System.out.println("    --deletesam\tDelete sam files. default: N");
			System.out.println("    --keeptemp\tKeep temp sam and bedpe file. default: N");
			System.out.println("    --map_ambiguous\tAlso mapping ambiguous reads without linker. default: N");
			System.out.println("    --skipmap\tSkip mapping read1 and read2, start from paired R1.sam and R2.sam, only valid in HiChIP mode now. default: N");
			System.out.println("    --macs2\tmacs2 path, using macs2 callpeak to detect anchor peak with alignment file. default: N");
			System.out.println("    --nomodel\tmacs2 parameter, Whether or not to build the shifting model in macs2. default: N");
			System.out.println("    --shortestP\textend and keep shorest peak length longer than N for loop calling, suggest 1500, user can set 0 to skip this step. default: 1500");
			System.out.println("    --shortestA\textend and keep shorest anchor length longer than N for loop calling, user can set 0 to skip this step. default: 0");
			System.out.println("    --XOR_cluster\tWhether keep loops if only one side of anchor is overlap with peak. default: N");
			System.out.println("    --addcluster\tKeep all regions with more than 2 count reads as potential anchor for calling loop. default: N. if peaks number of macs2 smaller than 10000, this paramater will work automaticly.");
			System.exit(0);
		}
		
		Path p = new Path();
		p.setParameter(args);
		
		if(Integer.valueOf(p.STOP_STEP)< Integer.valueOf(p.START_STEP)) {
			System.out.println("Stop step must bigger than start step!!!");
			System.exit(0);
		}
		
		File file = new File(p.OUTPUT_DIRECTORY+"/"+p.OUTPUT_PREFIX);
		if (!file.exists()) {
			file.mkdirs();
		}
		// rm 
		String rmcmd = "rm "+p.OUTPUT_DIRECTORY +"/"+ p.OUTPUT_PREFIX +".*.bedpe.txt";
		shellrun(rmcmd, "rm");
		
		Calendar rightNow = Calendar.getInstance();
		System.out.println("[" + rightNow.getTime().toString() +"] start ChIA-PET analysis");
		
		if(p.skipmap.equalsIgnoreCase("Y")) {
			if(Integer.valueOf(p.START_STEP) <= 1) {
				p.START_STEP = "2";
			}
			
			//Total PETs
			long nPETs_hichip = 0;
			String[] fastqs = p.Fastq_file_1.split(",");
		    for(int jk = 0; jk < fastqs.length; jk++) {		    	
		    	BufferedReader fastqFileIn;
		    	File fastq = new File(fastqs[jk]).getCanonicalFile();
		    	if(isGZipped(fastq)) {
					System.out.println("[" + rightNow.getTime().toString() +"] Counting total pets with gzip fastq file, " + fastqs[jk]);
					fastqFileIn  = new BufferedReader(
			                new InputStreamReader(new GZIPInputStream(new FileInputStream(fastq))));
				}else {
					System.out.println("[" + rightNow.getTime().toString() +"] Counting total pets with fastq file, " + fastqs[jk]);
					fastqFileIn = new BufferedReader(new FileReader(fastq));
				}
		    	
		        String readline = fastqFileIn.readLine();
		
		        while (readline != null) {
		        	nPETs_hichip++;
		            readline = fastqFileIn.readLine();
		        }
		        fastqFileIn.close();
		    }
		    //BufferedWriter localPrintWriter = new BufferedWriter(new FileWriter(p.OUTPUT_DIRECTORY + "/" + p.OUTPUT_PREFIX + ".basic_statistics.txt", false));
			BufferedWriter localPrintWriter = new BufferedWriter(new FileWriter(p.OUTPUT_DIRECTORY + "/" +p.OUTPUT_PREFIX+"/" + p.OUTPUT_PREFIX + ".basic_statistics.txt", false));
		    localPrintWriter.write("Total PETs\t" + nPETs_hichip/4);
		    localPrintWriter.newLine();
		    localPrintWriter.close();
		}
		
		long time = System.currentTimeMillis() / 1000;
		
		if(!p.fastp.isEmpty() && Integer.valueOf(p.START_STEP) <= 1) {
			System.out.println("[" + rightNow.getTime().toString() +"] Step0: fastq mode " + p.FQMODE);
			if(!p.fastp.equals("fastp")) {
				File file_tmp = new File(p.fastp);
				if(!file_tmp.exists()) {
					System.out.println("Uncorrect fastp path: " + p.fastp);
					System.exit(0);
				}
			}
		    time = System.currentTimeMillis() / 1000;
		    System.out.println("[" + rightNow.getTime().toString() +"] Step1: clean fastq with " + p.fastp);
		    String[] newfqs = CleanFastq(p.Fastq_file_1, p.Fastq_file_2, p.NTHREADS, p.OUTPUT_DIRECTORY, p.OUTPUT_PREFIX);
		    p.Fastq_file_1 = newfqs[0];
		    p.Fastq_file_2 = newfqs[1];
		    //System.out.println("TTT " + p.Fastq_file_1 + "  " + p.Fastq_file_2);
		}
		
		if(Integer.valueOf(p.START_STEP) <= 1 && p.hichipM.equals("N")) {
			if(p.AutoLinker.equals("true")) {
				String outPrefix = p.OUTPUT_DIRECTORY+"/"+p.OUTPUT_PREFIX+"/"+p.OUTPUT_PREFIX;
				if(p.FQMODE.equals("single-end")) {
					String mycmd = "java -cp "+p.PROGRAM_DIRECTORY+"/ChIA-PET.jar LGL.util.kmer "+
							" --fq1 " + p.Fastq_file_1 +" -k " + p.kmerlen + " -N " + p.linkerreads + " -s 0 -e 0 --score 12 " +
							"--prefix " + outPrefix + ".k" + p.kmerlen + " --minlinker 13 "+ " --skip " + p.skipheader + " -p 0.05";
					shellrun(mycmd, "2");
				}else {
					String mycmd = "java -cp "+p.PROGRAM_DIRECTORY+"/ChIA-PET.jar LGL.util.kmer "+
							" --fq1 " + p.Fastq_file_1 +" --fq2 " + p.Fastq_file_2 +" -k " + p.kmerlen + " -N " + p.linkerreads + " -s 0 -e 0 --score 12 " +
							"--prefix " + outPrefix + ".k" + p.kmerlen + " --minlinker 13 "+ " --skip " + p.skipheader + " -p 0.05";
					shellrun(mycmd, "2");
				}
				//prefix.k9.linker.txt
				/*
				 * Linker: 0
				 * 0_0
				 * */
				//p.linker
				BufferedReader linkerBufferedReader  = new BufferedReader(new FileReader( outPrefix + ".k" + p.kmerlen + ".linker.txt"));
			    String str;String newmode = "";float readlenaver =0;int Nlinker=0;
			    while ((str = linkerBufferedReader.readLine()) != null) {
			    	if (str.startsWith("Linker_mode:")) {
			    		newmode = str.split(" ")[1];
			    	}else if(str.startsWith("Seq_len:")) {
			    		readlenaver = Float.parseFloat(str.split(" ")[1]);
			    	}else {
			    		Nlinker++;
			    	}
			    }
			    //System.out.println("\n----xxxxaa11----- " + newmode + ";");
			    linkerBufferedReader.close();
			    if(!newmode.equalsIgnoreCase("")) {
				    if(p.MODE.equalsIgnoreCase("0")) {
		    			if(newmode.equalsIgnoreCase("1") && readlenaver>80 ) {
		    				//System.out.println("\n----111-----\n");
		    				p.MODE = "1";
		    			}
		    		}else {
		    			if(newmode.equalsIgnoreCase("0") && readlenaver<=80) {
		    				//System.out.println("\n----111.1-----\n");
		    				p.MODE = "0";
		    			}else if(newmode.equalsIgnoreCase("0") && readlenaver>80) {
		    				//System.out.println("\n----111.22-----\n");
		    				p.MAP2Linker = "true"; // 切linker时会将只有一端swlinker序列根据这个输出1-2，2-1还是1-1，2-2
		    				//p.ALLMAP = "true";
		    			}
		    		}
			    }
			    if(Nlinker==1) { // AA
			    	if(readlenaver>80) {
			    		//System.out.println("\n----111.32-----\n");
			    		p.MODE = "1";
			    		p.MAP2Linker = "true";
			    	}else {
			    		//System.out.println("\n----111.xxx-----\n");
			    		p.MODE = "0";
			    	}
			    }
			    if(readlenaver<80) {
			    	p.EXTENSION_LENGTH = "1500";
			    	p.MAPPING_CUTOFF = "20";
			    }
			    p.linker = outPrefix + ".k" + p.kmerlen + ".linker.txt";
			    time = System.currentTimeMillis() / 1000;
				System.out.println("[" + rightNow.getTime().toString() +"] Step0: Linker mode " + p.MODE);
			}
		}
		
		if(Integer.valueOf(p.STOP_STEP) <= 0) {
			System.exit(0);
		}
		
		if( Integer.valueOf(p.START_STEP) <= 1 && (p.hichipM.equals("Y") || p.hichipM.equals("O") ) ) {
			p.MERGE_DISTANCE = "-1";
			System.out.println("[" + rightNow.getTime().toString() +"] HiChIP mode " + p.hichipM);
			System.out.println("[" + rightNow.getTime().toString() +"] Step1: Checking and processing restriction file ...");
			if(p.removeResblock.equals("Y")) {
				if(p.restrictionsiteFile.equals("None")) {
					if(!p.ligation_site.equals("-")) {
						System.out.println("[" + rightNow.getTime().toString() +"] Without restriction file, generating by ligation site!!!");
						//generate res file by ligation site
						if(p.genomefile.equals("")) {
							p.genomefile = p.GENOME_INDEX;
						}
						File Gefile = new File(p.genomefile); 
						if (!Gefile.exists()) {
							System.out.println("Warning: Unvalid genome file path!!! "
									+ p.genomefile);
							System.exit(0);
						}
						p.restrictionsiteFile = p.OUTPUT_DIRECTORY + "/" + p.OUTPUT_PREFIX + "/" + p.OUTPUT_PREFIX +".res.bed";
						split_genome_byressite(p.genomefile, p.ligation_sites, p.restrictionsiteFile, p.ligation_site);
						
						//System.exit(0);
					}else {
						System.out.println("Here is HiChiP mode, please defind --restrictionsiteFile or --ligation_site!!!!");
						System.exit(0);
					}
				}else {
					//check file path is valid
					File resfile = new File(p.restrictionsiteFile); 
					if (!resfile.exists()) {
						System.out.println("Warning: Unvalid restriction file path!!! "
								+ p.restrictionsiteFile + ", please define valid path or "
										+ "\njust provide --ligation_site and remove --restrictionsiteFile paramater,"
										+ "\nthen pragram will generate restriction site file !!!\n");
						System.exit(0);
					}else {
						System.out.println("[" + rightNow.getTime().toString() +"] Valid restriction file path, "
								+ p.restrictionsiteFile );
						if(p.hichipM.equals("O")) {
						    System.out.println("[" + rightNow.getTime().toString() +"] We will not regenerate again. If you want to regenerate a new enzyme digestion site file, "
										+ "\n                please remove the --restrictionsitefile parameter and run the program again. "
										+ "\n                And then we will generate the corresponding enzyme digestion file in the specified output directory.");
						}
					}
				}
			}
			//Total PETs
			long nPETs_hichip = 0;
			String[] fastqs = p.Fastq_file_1.split(",");
		    for(int jk = 0; jk < fastqs.length; jk++) {
		    	BufferedReader fastqFileIn;
		    	File fastq = new File(fastqs[jk]);
		    	if(isGZipped(fastq)) {
					System.out.println("[" + rightNow.getTime().toString() +"] Counting total pets with gzip fastq file, " + fastqs[jk]);
					fastqFileIn  = new BufferedReader(
			                new InputStreamReader(new GZIPInputStream(new FileInputStream(fastq))));
				}else {
					System.out.println("[" + rightNow.getTime().toString() +"] Counting total pets with fastq file, " + fastqs[jk]);
					fastqFileIn = new BufferedReader(new FileReader(fastq));
				}
		    	
		        String readline = fastqFileIn.readLine();
		
		        while (readline != null) {
		        	nPETs_hichip++;
		            readline = fastqFileIn.readLine();
		        }
		        fastqFileIn.close();
		     }
		    //BufferedWriter localPrintWriter = new BufferedWriter(new FileWriter(p.OUTPUT_DIRECTORY + "/" + p.OUTPUT_PREFIX + ".basic_statistics.txt", false));
			BufferedWriter localPrintWriter = new BufferedWriter(new FileWriter(p.OUTPUT_DIRECTORY + "/" +p.OUTPUT_PREFIX+"/" + p.OUTPUT_PREFIX + ".basic_statistics.txt", false));
		    localPrintWriter.write("Total PETs\t" + nPETs_hichip/4);
		    localPrintWriter.newLine();
		    localPrintWriter.close();
		    
			p.MODE = "1";
			if(p.hichipM.equals("O")) {
				System.out.println("\n[" + rightNow.getTime().toString() +"] Note. Here is HiChiP mode, cause you specified that "
						+ "\n                only enzyme digestion genome processing should be performed, "
						+ "\n                so the program exited now. If you want to run all programs, please modify the --hichip parameter.");
				System.exit(0);
			}
		}
		else if(p.FQMODE.equals("single-end")) {
			//fastq_file linker_file minimum_linker_alignment_score
			//minimum_tag_length maximum_tag_length minSecondBestScoreDiff
			//output_data_with_ambiguous_linker_info
			if(Integer.valueOf(p.START_STEP) <= 1) {
				checkPath(p.linker);
				String mycmd = "java -cp "+p.PROGRAM_DIRECTORY+"/ChIA-PET.jar LGL.chiapet.LinkerFiltering_FastQ_singleEnd"+
					" " + p.Fastq_file_1 +" "+ p.linker +" "+ p.minimum_linker_alignment_score+" "+  
					p.minimum_tag_length +" "+ p.maximum_tag_length +" "+ p.minSecondBestScoreDiff+" "+ 
					p.output_data_with_ambiguous_linker_info +" "+ p.OUTPUT_DIRECTORY + "/" + 
					p.OUTPUT_PREFIX +" "+ p.OUTPUT_PREFIX;
				shellrun(mycmd, "2");
			}
		}else {
			if (Integer.valueOf(p.START_STEP) == 1) {
				checkPath(p.linker);
				LinkerFiltering lf = new LinkerFiltering(p);
				lf.reset(Integer.valueOf(p.START_STEP));
				lf.writeFile(p.OUTPUT_DIRECTORY+"/"+p.OUTPUT_PREFIX+"/"+p.OUTPUT_PREFIX+".time.txt", String.valueOf(time), false);
				System.out.println("[" + rightNow.getTime().toString() +"] Step1: Linker filtering ...");
				lf.run();
			}
		}
		if(p.FQMODE.equals("single-end")) {
			p.MODE = "0";
			//p.EXTENSION_LENGTH = "1500";
	    	//p.MAPPING_CUTOFF = "20";
		}
		
		
		if( Integer.valueOf(p.START_STEP) > 1 && Integer.valueOf(p.START_STEP) <=3 && (p.hichipM.equals("Y") || p.hichipM.equals("O") ) ) {
			p.MERGE_DISTANCE = "-1";
			System.out.println("[" + rightNow.getTime().toString() +"] HiChIP mode " + p.hichipM);
			if(p.removeResblock.equals("Y")) {
				if(p.restrictionsiteFile.equals("None")) {
					if(!p.ligation_site.equals("-")) {
						System.out.println("[" + rightNow.getTime().toString() +"] generating restriction file by ligation site!!!");
						//generate res file by ligation site
						if(p.genomefile.equals("")) {
							p.genomefile = p.GENOME_INDEX;
						}
						File Gefile = new File(p.genomefile); 
						if (!Gefile.exists()) {
							System.out.println("Warning: Unvalid genome file path!!! "
									+ p.genomefile);
							System.exit(0);
						}
						p.restrictionsiteFile = p.OUTPUT_DIRECTORY + "/" + p.OUTPUT_PREFIX + "/" + p.OUTPUT_PREFIX +".res.bed";
					    split_genome_byressite(p.genomefile, p.ligation_sites, p.restrictionsiteFile, p.ligation_site);
						
						//System.exit(0);
					}else {
						System.out.println("Here is HiChiP mode, please defind --restrictionsiteFile or --ligation_site!!!!");
						System.exit(0);
					}
				}
			}
		}
		
		//if(p.hichipM.equals("N")) {
		//    System.exit(0);
		//}
		
		if(Integer.valueOf(p.STOP_STEP) <= 1) {
			System.exit(0);
		}
		
		time = System.currentTimeMillis() / 1000;
		writeFile(p.OUTPUT_DIRECTORY+"/"+p.OUTPUT_PREFIX+"/"+p.OUTPUT_PREFIX+".time.txt", String.valueOf(time), true);
		if (Integer.valueOf(p.START_STEP) <= 2 && p.hichipM.equals("NNNNNNN") ) {
			//check linker information, reverse linker
			BufferedReader localBufferedReader  = new BufferedReader(new FileReader(p.linker));
		    ArrayList<String> localArrayList = new ArrayList<String>();
		    String str;
		    while ((str = localBufferedReader.readLine()) != null) {
		    	if(str.startsWith("Lin") || str.startsWith("Seq")) {
		    		continue;
		    	}
		    	if (str.length() > 1) {
		    		localArrayList.add(str.trim());
		    	}
		    }
		    localBufferedReader.close();
		    if(localArrayList.size()==1) {
		    	//if(p.MODE)
		    }
		    /*
		    if(localArrayList.size()>1) {
		    	String seq1 = localArrayList.get(0);
		    	String seq2 = LGL.util.SeqUtil.revComplement(localArrayList.get(1));
		    	int maxlenl = seq1.length();
		    	if(maxlenl != seq2.length()) {
		    		p.MAP2Linker = "true";
		    	}else {
			    	if(maxlenl < seq2.length()) {
			    		maxlenl = seq2.length();
			    	}
			    	LocalAlignment localAligner = new LocalAlignment(maxlenl, maxlenl);
			    	localAligner.align(seq1, seq2, 0);
			        if(localAligner.getMaxScore()<maxlenl-1) {
			        	p.MAP2Linker = "true";
			        }
		    	}
		    	rightNow = Calendar.getInstance();
				System.out.println("[" + rightNow.getTime().toString() +"] Step1: 2 linker mode " + p.MAP2Linker);
		    }
		    */
		}
		// for mapping
		if (Integer.valueOf(p.START_STEP) <= 2) {
			rightNow = Calendar.getInstance();
			System.out.println("[" + rightNow.getTime().toString() +"] Step2: Mapping to genome ...");
			Mapping mapping = new Mapping(p);
			mapping.Map();
		}
		
		if(Integer.valueOf(p.STOP_STEP) <= 2) {
			System.exit(0);
		}
		
		//step3 .bedpe
		time = System.currentTimeMillis() / 1000;
		writeFile(p.OUTPUT_DIRECTORY+"/"+p.OUTPUT_PREFIX+"/"+p.OUTPUT_PREFIX+".time.txt", String.valueOf(time), true);
		if (Integer.valueOf(p.START_STEP) <= 3) {
			rightNow = Calendar.getInstance();
			System.out.println("[" + rightNow.getTime().toString() +"] Step3: Removing redundancy ...");
			Purifying purifying = new Purifying(p);
			purifying.Purify();
			purifying.combiningData();
		}
		
		time = System.currentTimeMillis() / 1000;
		writeFile(p.OUTPUT_DIRECTORY+"/"+p.OUTPUT_PREFIX+"/"+p.OUTPUT_PREFIX+".time.txt", String.valueOf(time), true);
		if(Integer.valueOf(p.STOP_STEP) <= 3) {
			System.exit(0);
		}
		
		//step4 .bedpe.selected.unique.txt .bedpe.selected.pet.txt
		if (Integer.valueOf(p.START_STEP) <= 4) {
			rightNow = Calendar.getInstance();
			System.out.println("[" + rightNow.getTime().toString() +"] Step4: Categorization of PETs ...");
			DividePets dvidePets = new DividePets(p);
			dvidePets.dividePets();
		}
		time = System.currentTimeMillis() / 1000;
		writeFile(p.OUTPUT_DIRECTORY+"/"+p.OUTPUT_PREFIX+"/"+p.OUTPUT_PREFIX+".time.txt", String.valueOf(time), true);
		if(Integer.valueOf(p.STOP_STEP) <= 4) {
			System.exit(0);
		}
		
		//step5 .bedpe.selected.unique.txt ipet spet opet
		if (Integer.valueOf(p.START_STEP) <= 5) {
			//macs2
			if(!p.macs2.equalsIgnoreCase("N")) {
				if(!p.INPUT_ANCHOR_FILE.equals("null")) {
					System.out.println("[** Warning **]: you simutately provide --INPUT_ANCHOR_FILE " + p.INPUT_ANCHOR_FILE + " and --macs2 " + p.macs2 + ", we will"
							+ " run anchor mode with your provide --INPUT_ANCHOR_FILE file instead of macs2 call peak mode.");
				}else {
					rightNow = Calendar.getInstance();
					System.out.println("[" + rightNow.getTime().toString() +"] Step5: Calling peak with macs2 ...");
					String outPrefix = p.OUTPUT_DIRECTORY + "/" + p.OUTPUT_PREFIX + "/" + p.OUTPUT_PREFIX;
					BufferedWriter macs2F = new BufferedWriter(new FileWriter(outPrefix + 
							".callpeak.macs2.sh", false));
					macs2F.write("## macs2 callpeak");
					macs2F.newLine();
					//String runcmd = "macs2 callpeak -t " + outPrefix+ "*.pe.sam -f SAM --nomodel --extsize 147 --keep-dup all -q 0.01 -g "+
					//		p.GENOME_LENGTH + " -n " + p.OUTPUT_PREFIX;// + " --extsize 147"; //.merge.srt.sam
					String runcmd = "awk -v OFS=\"\\t\" '{if($9==\"+\"){newstrand=\"-\"}else{newstrand=\"+\"}; print $1,$2,$3,\".\\t.\",newstrand}' " + outPrefix + ".bedpe.selected.unique.txt > " + outPrefix + ".allValidPairs.bed";
					macs2F.write(runcmd);
					macs2F.newLine();
					runcmd = "awk -v OFS=\"\\t\" '{if($10==\"+\"){newstrand=\"-\"}else{newstrand=\"+\"}; print $4,$5,$6,\".\\t.\",newstrand}' " + outPrefix + ".bedpe.selected.unique.txt >> " + outPrefix + ".allValidPairs.bed";
					macs2F.write(runcmd);
					macs2F.newLine();
					runcmd = "macs2 callpeak -t " + outPrefix + ".allValidPairs.bed --keep-dup all -g " + 
					      p.GENOME_LENGTH + " -f BED -n " + outPrefix + " " + p.broadpeak;// -B  --verbose 1
					if(p.nomodel.equalsIgnoreCase("Y")) {
						runcmd = runcmd + " --nomodel --extsize 147 -q 0.01 "; 
					}
					macs2F.write(runcmd);
					macs2F.newLine();
	
					macs2F.close();
			    	Shell shell = new Shell();
					shell.runShell(outPrefix + ".callpeak.macs2.sh");
					//new File(outPrefix + ".callpeak.macs2.sh").delete();
					p.INPUT_ANCHOR_FILE = outPrefix+ "_peaks.narrowPeak";
					if(!p.broadpeak.equals("")) {
						p.INPUT_ANCHOR_FILE = outPrefix+ "_peaks.broadPeak";
					}
					
					File file_a = new File(p.INPUT_ANCHOR_FILE);
					if(file_a.length() == 0 && !p.nomodel.equalsIgnoreCase("Y")) {
						BufferedWriter macs2F_nomodel = new BufferedWriter(new FileWriter(outPrefix + 
								".callpeak.macs2-nomodel.sh", false));
						macs2F_nomodel.write("## macs2 callpeak with no model");
						macs2F_nomodel.newLine();
						runcmd = "macs2 callpeak -t " + outPrefix + ".allValidPairs.bed --keep-dup all -g " + 
							      p.GENOME_LENGTH + " -f BED -n " + outPrefix + " --nomodel --extsize 147 -q 0.01" + " " + p.broadpeak;
						macs2F_nomodel.write(runcmd);
						macs2F_nomodel.newLine();
						macs2F_nomodel.close();
				    	//Shell shell = new Shell();
						shell.runShell(outPrefix + ".callpeak.macs2-nomodel.sh");
					}

				}
			}

			if(!p.INPUT_ANCHOR_FILE.equals("null")) {
				String panchor = "N";
				if(!p.macs2.equals("N")) {
					BufferedReader reader = new BufferedReader(new FileReader(p.INPUT_ANCHOR_FILE));
					String data; int Npeak = 0;
					while((data = reader.readLine()) != null) {
						Npeak++;
					}
					if(Npeak < p.peakcutoff || p.addcluster.equalsIgnoreCase("Y")) {
						//find cluster as loop
						String temp_len = p.EXTENSION_LENGTH;
						p.EXTENSION_LENGTH = "500";
						findpeak findpeak = new findpeak(p);
						findpeak.run();
						p.EXTENSION_LENGTH = temp_len;
						panchor = "Y";
					}
				}
				
				String outPrefix = p.OUTPUT_DIRECTORY + "/" + p.OUTPUT_PREFIX + "/" + p.OUTPUT_PREFIX;
				BufferedWriter anchor2bedpe = new BufferedWriter(new FileWriter(outPrefix + 
						".anchor2bedpe.sh", false));
				anchor2bedpe.write("## anchor to bedpe");
				anchor2bedpe.newLine();
				String runcmd="";
				if( Integer.parseInt(p.EXTENSION_LENGTH) > 0) {
				    runcmd="bedtools slop -i " + p.INPUT_ANCHOR_FILE + " -g " + p.CHROM_SIZE_INFO + 
						" -b " + p.EXTENSION_LENGTH + " > "+ outPrefix + ".cpt.exd.peak";
				}else {
					runcmd="cp " + p.INPUT_ANCHOR_FILE + " " + outPrefix + ".cpt.exd.peak";
				}
				anchor2bedpe.write(runcmd);
				anchor2bedpe.newLine();
				
				if(p.shortestPeak>0) {
					runcmd="awk -v OFS=\"\\t\" -v minlen=" + p.shortestPeak + " '{if($3-$2>=minlen){print $0}else{middle=int(minlen/2); summit=int(($3+$2)/2); if(summit>middle){m1=summit-middle}else{m1=0}; m2=summit+middle; $2=m1; $3=m2; print $0}}' " + outPrefix+ ".cpt.exd.peak | sort -k1,1 -k2,2n -k3,3n > " + outPrefix + ".cpt.exd.srt.peak";	
					anchor2bedpe.write(runcmd);
					anchor2bedpe.newLine();
					
					runcmd="rm "+outPrefix+".cpt.exd.peak";
					anchor2bedpe.write(runcmd);
					anchor2bedpe.newLine();
				}else {
					runcmd="mv "+outPrefix+".cpt.exd.peak "+outPrefix+".cpt.exd.srt.peak";
					anchor2bedpe.write(runcmd);
					anchor2bedpe.newLine();
				}
				
				if(panchor.equals("Y")) {
					runcmd = "sort -k1,1V -k2,2n -k3,3n " + outPrefix+ ".panchor.bed | uniq > tmp.bed; mv tmp.bed " + outPrefix+ ".panchor.bed";
					anchor2bedpe.write(runcmd);
					anchor2bedpe.newLine();
					runcmd="awk -v OFS=\"\\t\" '{print $1,$2,$3}' " + outPrefix+".cpt.exd.srt.peak " + outPrefix+ ".panchor.bed | sort -k1,1V -k2,2n -k3,3n | uniq > tmp.bed; mv tmp.bed " + outPrefix+".cpt.exd.srt.peak";
					anchor2bedpe.write(runcmd);
					anchor2bedpe.newLine();
				}
				
				runcmd= "bedtools merge -d 256 -i "+outPrefix+".cpt.exd.srt.peak > " + outPrefix + ".cpt.merge.peak";
				anchor2bedpe.write(runcmd);
				anchor2bedpe.newLine();
				
				runcmd="rm "+outPrefix+".cpt.exd.srt.peak";
				anchor2bedpe.write(runcmd);
				anchor2bedpe.newLine();
				
				if(p.shortestAnchor>0) {
					runcmd="awk -v OFS=\"\\t\" -v minlen=" + p.shortestAnchor + " '{if($3-$2>=minlen){print $0}else{middle=int(minlen/2); summit=int(($3+$2)/2); if(summit>middle){m1=summit-middle}else{m1=0}; m2=summit+middle; $2=m1; $3=m2; print $0}}' "  + outPrefix  + ".cpt.merge.peak | bedtools merge -i - > "+outPrefix+".cpt_peaks.tmp.bed; mv "+outPrefix+".cpt_peaks.tmp.bed " + outPrefix  + ".cpt.merge.peak";	
					anchor2bedpe.write(runcmd);
					anchor2bedpe.newLine();
				}
				
			    runcmd="awk 'BEGIN{OFS=\"\\t\";i=1}{print $1,$2,$3,\"peak_\"i;i=i+1}' " +  outPrefix  + ".cpt.merge.peak > "+outPrefix+"_peaks.slopPeak";
				anchor2bedpe.write(runcmd);
				anchor2bedpe.newLine();
				//runcmd="mv cpt_peaks.tmp.bed " + outPrefix  + "_peaks.slopPeak";
				runcmd="rm "+outPrefix+".cpt.merge.peak";
				anchor2bedpe.write(runcmd);
				anchor2bedpe.newLine();
				
				
				runcmd="awk '$1!=$4 || $6-$2>=" + p.SELF_LIGATION_CUFOFF + "' " + outPrefix + ".bedpe.selected.unique.txt > " + outPrefix + ".bedpe.selected.unique.validpet.txt";
				anchor2bedpe.write(runcmd);
				anchor2bedpe.newLine();
				//interact
				runcmd="pairToBed -a " + outPrefix + ".bedpe.selected.unique.validpet.txt -b " + outPrefix + "_peaks.slopPeak -type both > " + outPrefix + ".rmdup.bedpe.tmp";
				anchor2bedpe.write(runcmd);
				anchor2bedpe.newLine();
												
				anchor2bedpe.close();
				Shell shell = new Shell();
				shell.runShell(outPrefix + ".anchor2bedpe.sh");
				
				peak2interaction.main(new String[]{outPrefix + ".rmdup.bedpe.tmp", outPrefix + ".cluster.filtered"});
				p.INPUT_ANCHOR_FILE = outPrefix + "_peaks.slopPeak";
				p.EXTENSION_LENGTH = "150"; //cause ipet is only keep start
				//p.EXTENSION_MODE = "1"; //0 only extend to downstream.
				
				//delete
				if(p.keeptemp.equals("N")) {
					new File(outPrefix + ".bedpe.selected.unique.validpet.txt").delete();
					new File(outPrefix + ".rmdup.bedpe.tmp").delete();
				}
			}
			
			//cluster mode
			if(p.INPUT_ANCHOR_FILE.equals("null")) {
				rightNow = Calendar.getInstance();
				System.out.println("[" + rightNow.getTime().toString() +"] Step5: Interaction Calling ...");
				InteractionCalling interactionCalling = new InteractionCalling(p);
				interactionCalling.run();
			}
			
			String outPrefix = p.OUTPUT_DIRECTORY + "/" + p.OUTPUT_PREFIX + "/" + p.OUTPUT_PREFIX;
			BufferedWriter anchor2bedpe = new BufferedWriter(new FileWriter(outPrefix + 
					".validanchor.sh", false));
			anchor2bedpe.write("## get unique anchor");
			anchor2bedpe.newLine();
			String runcmd="awk -v OFS=\"\\t\" '{print $1,$2,$3; print $7,$8,$9}' " + outPrefix + ".cluster.filtered | sort -k1,1 -k2,2n -k3,3n | uniq > " + outPrefix + ".cluster.filtered.validunique.anchor";
			anchor2bedpe.write(runcmd);
			anchor2bedpe.newLine();
											
			anchor2bedpe.close();
			Shell shell = new Shell();
			shell.runShell(outPrefix + ".validanchor.sh");
			
			System.out.println("[" + rightNow.getTime().toString() +"] Pvalue calculation ...");
			Pvalues pValues = new Pvalues(p);
			pValues.calculation(); //ipet one-read count in anchor. ipet to aln, count tag
			pValues.globalTag();
			pValues.calculate();
		}
		
		time = System.currentTimeMillis() / 1000;
		writeFile(p.OUTPUT_DIRECTORY+"/"+p.OUTPUT_PREFIX+"/"+p.OUTPUT_PREFIX+".time.txt", String.valueOf(time), true);
		if(Integer.valueOf(p.STOP_STEP) <= 5) {
			System.exit(0);
		}
		
		PeakCalling peakCalling = new PeakCalling(p);
		if (Integer.valueOf(p.START_STEP) <= 6) {
			rightNow = Calendar.getInstance();
			System.out.println("[" + rightNow.getTime().toString() +"] Step6: Peak Calling ...");
			peakCalling.run();
			peakCalling.localTag();
			peakCalling.globalSpet();
			peakCalling.calculatePvalue();
			peakCalling.chromosomalPet();
		}
		time = System.currentTimeMillis() / 1000;
		writeFile(p.OUTPUT_DIRECTORY+"/"+p.OUTPUT_PREFIX+"/"+p.OUTPUT_PREFIX+".time.txt", String.valueOf(time), true);
		if(Integer.valueOf(p.STOP_STEP) <= 6) {
			System.exit(0);
		}
		
		if (Integer.valueOf(p.START_STEP) <= 7) {
			rightNow = Calendar.getInstance();
			System.out.println("[" + rightNow.getTime().toString() +"] Step7: Visualizing ...");
			peakCalling.runningTime();
			peakCalling.summary();
			peakCalling.statisticsReport();
			peakCalling.genomicBrowser();
			//new File(p.PROGRAM_DIRECTORY+"/"+p.OUTPUT_PREFIX+".peakcalling2.sh").delete();
			new File(p.OUTPUT_DIRECTORY+"/"+p.OUTPUT_PREFIX+"/"+p.OUTPUT_PREFIX+".peakcalling2.sh").delete();
		}
		if(Integer.valueOf(p.STOP_STEP) <= 7) {
			System.exit(0);
		}
		
		if(p.zipbedpe.equals("Y")) {
			rightNow = Calendar.getInstance();
			System.out.println("[" + rightNow.getTime().toString() +"] gzip .bedpe* file");
			String outPrefix = p.OUTPUT_DIRECTORY+"/"+p.OUTPUT_PREFIX+"/"+p.OUTPUT_PREFIX;
			String line = "";
			BufferedWriter zipbedF = new BufferedWriter(new FileWriter(outPrefix + 
					".zipbedpe.sh", false));
			
		    line = "gzip -f "+outPrefix+".bedpe";
		    zipbedF.write(line);
		    zipbedF.newLine();
		    line = "gzip -f "+outPrefix+".bedpe.selected*";
		    zipbedF.write(line);
		    zipbedF.newLine();
		    line = "gzip -f "+outPrefix+".bedpe.filter.byres*";
		    zipbedF.write(line);
		    zipbedF.newLine();
		    zipbedF.close();
		    Shell shell = new Shell();
			shell.runShell(outPrefix + ".zipbedpe.sh");
			new File(outPrefix + ".zipbedpe.sh").delete();
		}
		if(p.zipsam.equals("Y")) {
			rightNow = Calendar.getInstance();
			System.out.println("[" + rightNow.getTime().toString() +"] convert sam file 2 bam file");
			String outPrefix = p.OUTPUT_DIRECTORY+"/"+p.OUTPUT_PREFIX+"/"+p.OUTPUT_PREFIX;
			String line = "";
			BufferedWriter sam2bam = new BufferedWriter(new FileWriter(outPrefix + 
					".sam2bam.sh", false));
			if(p.hichipM.equals("Y")) {
		    	//HiChIP mode, R1.sam R2.sam
				if(p.Fastq_file_1.contains(",")) {
					int Nfq = p.Fastq_file_1.split(",").length;
					for(int i=0;i<Nfq;i++) {
						line = "samtools view -bSh -o " + outPrefix + "." + i + ".R1.bam " + outPrefix + "." + i + ".R1.sam";
				    	sam2bam.write(line);
				    	sam2bam.newLine();
				    	line = "samtools view -bSh -o " + outPrefix + "." + i + ".R2.bam " + outPrefix + "." + i + ".R2.sam";
				    	sam2bam.write(line);
				    	sam2bam.newLine();
				    	line = "rm " + outPrefix + "." + i + ".R1.sam " + outPrefix + "." + i + ".R2.sam";
				    	sam2bam.write(line);
				    	sam2bam.newLine();
					}
				}else {
			    	line = "samtools view -bSh -o " + outPrefix + ".R1.bam " + outPrefix + ".R1.sam";
			    	sam2bam.write(line);
			    	sam2bam.newLine();
			    	line = "samtools view -bSh -o " + outPrefix + ".R2.bam " + outPrefix + ".R2.sam";
			    	sam2bam.write(line);
			    	sam2bam.newLine();
			    	line = "rm " + outPrefix + ".R1.sam " + outPrefix + ".R2.sam";
			    	sam2bam.write(line);
			    	sam2bam.newLine();
				}
			}else {
				String header = "for x in 1_2 2_1";
	    		if(p.ALLMAP.equalsIgnoreCase("true")) {
	    			header = "for x in 1_2 2_1 1_1 2_2";
	    		}
	    		if(p.MAP2Linker.equalsIgnoreCase("true")) {
	    			header = "for x in 1_1 2_2";
	    		}
	    		sam2bam.write(header);
	    		sam2bam.newLine();
				line = "do";
				sam2bam.write(line);
				sam2bam.newLine(); //.1_2.R8.sam
				  header = "for y in {1..8}";
				  sam2bam.write(header);
				  sam2bam.newLine();
				  line="do";
				  sam2bam.write(line);
				  sam2bam.newLine();
				    line = "    samtools view -bSh -o "+outPrefix +".${x}.R${y}.bam "+outPrefix+".${x}.R${y}.sam && rm "+outPrefix+".${x}.R${y}.sam";
				    sam2bam.write(line);
				    sam2bam.newLine();
				  line="done";
			      sam2bam.write(line);
				  sam2bam.newLine();
				line="done";
				sam2bam.write(line);
				sam2bam.newLine();
			}
			sam2bam.close();
			Shell shell = new Shell();
			shell.runShell(outPrefix + ".sam2bam.sh");
			new File(outPrefix + ".sam2bam.sh").delete();
		}
		
		// rm sam and bam
		if(p.deletesam.equals("Y")) {
			rightNow = Calendar.getInstance();
			System.out.println("[" + rightNow.getTime().toString() +"] delete sam file");
			String outPrefix = p.OUTPUT_DIRECTORY+"/"+p.OUTPUT_PREFIX+"/"+p.OUTPUT_PREFIX;
			String line = "";
			BufferedWriter sam2bam = new BufferedWriter(new FileWriter(outPrefix + 
					".rmsam.sh", false));
			if(p.hichipM.equals("Y")) {
		    	//HiChIP mode, R1.sam R2.sam
				if(p.Fastq_file_1.contains(",")) {
					int Nfq = p.Fastq_file_1.split(",").length;
					for(int i=0;i<Nfq;i++) {
				    	line = "rm " + outPrefix + "." + i + ".R1.sam " + outPrefix + "." + i + ".R2.sam";
				    	sam2bam.write(line);
				    	sam2bam.newLine();
					}
				}else {
			    	line = "rm " + outPrefix + ".R1.sam " + outPrefix + ".R2.sam";
			    	sam2bam.write(line);
			    	sam2bam.newLine();
				}
			}else {
				String header = "for x in 1_2 2_1";
	    		if(p.ALLMAP.equalsIgnoreCase("true")) {
	    			header = "for x in 1_2 2_1 1_1 2_2";
	    		}
	    		if(p.MAP2Linker.equalsIgnoreCase("true")) {
	    			header = "for x in 1_1 2_2";
	    		}
	    		sam2bam.write(header);
	    		sam2bam.newLine();
				line = "do";
				sam2bam.write(line);
				sam2bam.newLine(); //.1_2.R8.sam
				  header = "for y in {1..8}";
				  sam2bam.write(header);
				  sam2bam.newLine();
				  line="do";
				  sam2bam.write(line);
				  sam2bam.newLine();
				    line = "    rm "+outPrefix+".${x}.R${y}.sam";
				    sam2bam.write(line);
				    sam2bam.newLine();
				  line="done";
			      sam2bam.write(line);
				  sam2bam.newLine();
				line="done";
				sam2bam.write(line);
				sam2bam.newLine();
			}
			sam2bam.close();
			Shell shell = new Shell();
			shell.runShell(outPrefix + ".rmsam.sh");
			new File(outPrefix + ".rmsam.sh").delete();
		}
		rightNow = Calendar.getInstance();
		System.out.println("[" + rightNow.getTime().toString() +"] finish ChIA-PET analysis");
	}	
	
}


