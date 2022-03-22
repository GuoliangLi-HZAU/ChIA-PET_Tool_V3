package LGL.chiapet;

/**
 * @author update by sun
 * 
 */
import multhread.BigFileProcess;
import process.Path;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

public class LinkerFiltering_FastQ_PET_longread {
	
	int debug_level = 1;
	String fastQFile1 = null;
	String fastQFile2 = null;
	String linker = null;
	String inputInfoFile = null;
	String outputFolder = null;
	String outputPrefix = "out";
	String[] linkers = null;
	int nLinkers = 0;
	int[] barCodeStart = null;
    int[] barCodeLength = null;
    int maxLinkerLength;
    int minLinkerLength;
    int minimum_linker_alignment_score = 14;
    int minimum_tag_length = 18;
    int minSecondBestScoreDiff = 3;
    int output_data_with_ambiguous_linker_info = 1;
    int flip_head_tag = 0;
    int flip_tail_tag = 0;
    int[][] scoreDistribution;
    int[][] secondBestScoreDiffDistribution;
    int maxSecondBestScoreDiff = 0;
    int maximum_tag_length = 1000;
    int maxRealTagLength = 0;
    int[][] tagLengthDistribution;
    int[][] linkerCompositionDistribution;
    int nAmbiguousLinkerComposition = 0;
    BufferedWriter debug_output = null;
    String[] letter = { "A", "B", "C", "D", "E", "F" };
    
    private long nPETs = 0;
    private int threadNum = 0;
    static int maxThreadNum = 16;
    private Path p;
	public boolean search_all_linker = false;
	int printallreads = 0; //0 print all, 1 print known +-
	public String AutoLinker = "true";
	public String MAP2Linker = "false";
  
    public LinkerFiltering_FastQ_PET_longread(Path path) throws IOException {
    	p = path;
    	this.outputFolder = p.getOUTPUT_DIRECTORY()+"/"+p.getOUTPUT_PREFIX();
        this.outputPrefix = p.getOUTPUT_PREFIX();
        this.threadNum = Integer.valueOf(p.getNTHREADS());
        if(p.search_all_linker.equals("Y")) {
        	this.search_all_linker = true;
        }else {
        	this.search_all_linker = false;
        }
        
        this.printallreads = p.printallreads;
        this.AutoLinker = p.AutoLinker;
        this.MAP2Linker = p.MAP2Linker;
        
        if (this.threadNum > maxThreadNum) {
        	System.out.println("Error: parameter error. The maximun numbers of threads must <= " + maxThreadNum);
            System.exit(0);
        }
    
        readInputInformation();
    
        testInputInformation();
    
        readLinkers();
    
	    outputRunningInfo();
	    
	    this.tagLengthDistribution = new int[2][this.maximum_tag_length];
	    Arrays.fill(this.tagLengthDistribution[0], 0);
	    Arrays.fill(this.tagLengthDistribution[1], 0);
	    
	    this.scoreDistribution = new int[2][this.maxLinkerLength + 1];
	    Arrays.fill(this.scoreDistribution[0], 0);
	    Arrays.fill(this.scoreDistribution[1], 0);
	    this.secondBestScoreDiffDistribution = new int[2][this.maxLinkerLength * 2 + 1];
	    Arrays.fill(this.secondBestScoreDiffDistribution[0], 0);
	    Arrays.fill(this.secondBestScoreDiffDistribution[1], 0);
	    
	    filterSequenceByLinker();
	    printDistribution();
    }
  
    public void readInputInformation() throws IOException {
    	this.fastQFile1 = p.getFastq_file_1();
    	this.fastQFile2 = p.getFastq_file_2();
    	this.linker = p.getLinker();
    	this.minimum_linker_alignment_score = Integer.parseInt(p.getMinimum_linker_alignment_score());
    	this.minimum_tag_length = Integer.parseInt(p.getMinimum_tag_length());
    	this.maximum_tag_length = Integer.parseInt(p.getMaximum_tag_length());
    	this.minSecondBestScoreDiff = Integer.parseInt(p.getMinSecondBestScoreDiff());
    	this.output_data_with_ambiguous_linker_info = Integer.parseInt(p.getOutput_data_with_ambiguous_linker_info());
	}

	public void testInputInformation() throws IOException {
		boolean requiredInputMissed = false;
	    if (this.fastQFile1 == null) {
	    	System.out.println("No FASTQ file 1 input!");
	    	requiredInputMissed = true;
	    }
	    if (this.fastQFile2 == null) {
	    	System.out.println("No FASTQ file 2 input!");
	    	requiredInputMissed = true;
	    }
	    if (this.linker == null) {
            System.out.println("No linker input!");
            requiredInputMissed = true;
        }
	    if (requiredInputMissed == true) {
	    	System.exit(0);
	    }
	    if (this.maximum_tag_length < this.minimum_tag_length) {
	    	System.out.println("maximum_tag_length is smaller than the minimum_tag_length! Stop ...");
	    	System.exit(0);
	    }
	    if (this.outputFolder == null) {
	    	System.out.println("No output folder! Set it as \"out\"");
	      	this.outputFolder = "out";
	    }
	    File localFile = new File(this.outputFolder);
	    if (!localFile.exists()) {
	    	localFile.mkdir();
	    } else if (!localFile.isDirectory()) {
	    	System.out.println(this.outputFolder + " exists, but it is not a directory!!! please check...");
	    	System.exit(0);
	    }
	    this.outputFolder = localFile.getPath();
	}
	
	public void readLinkers() throws IOException {
		BufferedReader localBufferedReader  = new BufferedReader(new FileReader(this.linker));
	    ArrayList<String> localArrayList = new ArrayList<String>();
	    String str;
	    while ((str = localBufferedReader.readLine()) != null) {
	    	if(str.startsWith("Linker") || str.startsWith("Seq")) {
        		continue;
        	}
	    	if (str.length() > 1) {
	    		localArrayList.add(str.trim());
	    	}
	    }
	    localBufferedReader.close();
	    if(localArrayList.size()==1) {
	    	String seq2 = LGL.util.SeqUtil.revComplement(localArrayList.get(0));
	    	localArrayList.add(seq2);
	    }
	    
		this.linkers = new String[localArrayList.size()];
        this.maxLinkerLength = 0;
        this.minLinkerLength = Integer.MAX_VALUE;
	    for (int i = 0; i < this.linkers.length; i++) {
	    	this.linkers[i] = localArrayList.get(i);
	    	if (this.minLinkerLength > this.linkers[i].length()) {
	    		this.minLinkerLength = this.linkers[i].length();
	    	}
	    	if (this.maxLinkerLength < this.linkers[i].length()) {
	    		this.maxLinkerLength = this.linkers[i].length();
	    	}
	    }
	    this.nLinkers = this.linkers.length;
	    if (this.nLinkers <= 0) {
	    	System.out.println("No linker sequence information.\tStop!!!");
	    	System.exit(0);
	    }
	    if (this.nLinkers > 100) {
	    	System.out.println("Too many linkers. Please check...\tStop!!!");
	    	System.exit(0);
	    }
	    if (this.minLinkerLength < 5) {
	    	System.out.println("the shortest linker is less than 5bp.\tStop!!!");
	    	System.exit(0);
	    }
	    this.linkerCompositionDistribution = new int[this.nLinkers + 1][this.nLinkers + 1];
	    for (int i = 0; i <= this.nLinkers; i++) {
	    	Arrays.fill(this.linkerCompositionDistribution[i], 0);
	    }
	}

	public void outputRunningInfo() throws IOException {
		BufferedWriter localPrintWriter = new BufferedWriter(new FileWriter(this.outputFolder + "/" + this.outputPrefix + 
				".runningInformation.LinkerFiltering_FastQ_PET.txt", false));
	    localPrintWriter.write("Fastq file 1\t" + new File(this.fastQFile1).getAbsolutePath());
	    localPrintWriter.newLine();
	    localPrintWriter.write("Fastq file 2\t" + new File(this.fastQFile2).getAbsolutePath());
	    localPrintWriter.newLine();
	    for (int i = 0; i < this.linkers.length; i++) {
	    	localPrintWriter.write("Linker" + this.letter[i] + "\t" + this.linkers[i]);
	    	localPrintWriter.newLine();
	    }
	    localPrintWriter.write("Output folder\t" + new File(this.outputFolder).getAbsolutePath());
	    localPrintWriter.newLine();
	    localPrintWriter.write("Output prefix\t" + this.outputPrefix);
	    localPrintWriter.newLine();
	    localPrintWriter.write("Minimum linker alignment score\t" + this.minimum_linker_alignment_score);
	    localPrintWriter.newLine();
	    localPrintWriter.write("Minimum tag length\t" + this.minimum_tag_length);
	    localPrintWriter.newLine();
	    localPrintWriter.write("Maximum tag length\t" + this.maximum_tag_length);
	    localPrintWriter.newLine();
	    localPrintWriter.write("Minimum SecondBestScore difference\t" + this.minSecondBestScoreDiff);
	    localPrintWriter.newLine();
	    localPrintWriter.write("Output data with ambiguous linker info\t" + this.output_data_with_ambiguous_linker_info);
	    localPrintWriter.close();
	}
	
	public void filterSequenceByLinker() throws IOException {
		BufferedWriter[] arrayOfBufferedWriter = new BufferedWriter[this.nLinkers * this.nLinkers * 8 + 2];
	    for (int i = 0; i < this.nLinkers; i++) {
	    	for (int j = 0; j < this.nLinkers; j++) {
	    		arrayOfBufferedWriter[((i * this.nLinkers + j) * 8)] = new BufferedWriter(new FileWriter(this.outputFolder + "/" + this.outputPrefix + "." + 
	    		(i + 1) + "_" + (j + 1) + ".R1.fastq", false));
	    		arrayOfBufferedWriter[((i * this.nLinkers + j) * 8 + 1)] = new BufferedWriter(new FileWriter(this.outputFolder + "/" + this.outputPrefix + "." + 
	    		(i + 1) + "_" + (j + 1) + ".R2.fastq", false));
	    		arrayOfBufferedWriter[((i * this.nLinkers + j) * 8 + 2)] = new BufferedWriter(new FileWriter(this.outputFolder + "/" + this.outputPrefix + "." + 
	    		(i + 1) + "_" + (j + 1) + ".R3.fastq", false));
	    		arrayOfBufferedWriter[((i * this.nLinkers + j) * 8 + 3)] = new BufferedWriter(new FileWriter(this.outputFolder + "/" + this.outputPrefix + "." + 
		        (i + 1) + "_" + (j + 1) + ".R4.fastq", false));
	    		arrayOfBufferedWriter[((i * this.nLinkers + j) * 8 + 4)] = new BufferedWriter(new FileWriter(this.outputFolder + "/" + this.outputPrefix + "." + 
		        (i + 1) + "_" + (j + 1) + ".R5.fastq", false));
	    		arrayOfBufferedWriter[((i * this.nLinkers + j) * 8 + 5)] = new BufferedWriter(new FileWriter(this.outputFolder + "/" + this.outputPrefix + "." + 
		        (i + 1) + "_" + (j + 1) + ".R6.fastq", false));
	    		arrayOfBufferedWriter[((i * this.nLinkers + j) * 8 + 6)] = new BufferedWriter(new FileWriter(this.outputFolder + "/" + this.outputPrefix + "." +
		        (i + 1) + "_" + (j + 1) + ".R7.fastq", false));
	    		arrayOfBufferedWriter[((i * this.nLinkers + j) * 8 + 7)] = new BufferedWriter(new FileWriter(this.outputFolder + "/" + this.outputPrefix + "." + 
		        (i + 1) + "_" + (j + 1) + ".R8.fastq", false));
	    	}
	    }
	    if (this.output_data_with_ambiguous_linker_info == 1) {
	    	arrayOfBufferedWriter[(this.nLinkers * this.nLinkers * 8)] = new BufferedWriter(new FileWriter(this.outputFolder + "/" + this.outputPrefix + 
		    		  ".ambiguous.R1.fastq", false));
	    	arrayOfBufferedWriter[(this.nLinkers * this.nLinkers * 8 + 1)] = new BufferedWriter(new FileWriter(this.outputFolder + "/" + this.outputPrefix + 
		    		  ".ambiguous.R2.fastq", false));
	    }
	    if (this.debug_level >= 2) {
	    	this.debug_output = new BufferedWriter(new FileWriter(this.outputFolder + "/" + this.outputPrefix + "debug.LinkerFiltering_FastQ_PET.txt", false));
	    }
	    
	    //File fastq1 = new File(this.fastQFile1);
        //File fastq2 = new File(this.fastQFile2);
	    String[] fastq1s = this.fastQFile1.split(",");
	    String[] fastq2s = this.fastQFile2.split(",");
	    for(int jk = 0; jk < fastq1s.length; jk++) {
	    	File fastq1 = new File(fastq1s[jk]);
	        File fastq2 = new File(fastq2s[jk]);
	        BigFileProcess bfp = new BigFileProcess(fastq1, fastq2, this, arrayOfBufferedWriter, threadNum);
	        bfp.start();
	        bfp.join();
	    }
        
	    for (int k = 0; k < this.nLinkers; k++) {
	    	for (int m = 0; m < this.nLinkers; m++) {
	    		arrayOfBufferedWriter[((k * this.nLinkers + m) * 8)].close();
	    		arrayOfBufferedWriter[((k * this.nLinkers + m) * 8 + 1)].close();
	    		arrayOfBufferedWriter[((k * this.nLinkers + m) * 8 + 2)].close();
	    		arrayOfBufferedWriter[((k * this.nLinkers + m) * 8 + 3)].close();
	    		arrayOfBufferedWriter[((k * this.nLinkers + m) * 8 + 4)].close();
	    		arrayOfBufferedWriter[((k * this.nLinkers + m) * 8 + 5)].close();
	    		arrayOfBufferedWriter[((k * this.nLinkers + m) * 8 + 6)].close();
	    		arrayOfBufferedWriter[((k * this.nLinkers + m) * 8 + 7)].close();
	    	}
	    }
	    if (this.output_data_with_ambiguous_linker_info == 1) {
	    	arrayOfBufferedWriter[(this.nLinkers * this.nLinkers * 8)].close();
	    	arrayOfBufferedWriter[(this.nLinkers * this.nLinkers * 8 + 1)].close();
	    }
	    BufferedWriter localPrintWriter = new BufferedWriter(new FileWriter(this.outputFolder + "/" + this.outputPrefix + ".basic_statistics.txt", false));
	    localPrintWriter.write("Total PETs\t" + nPETs);
	    localPrintWriter.newLine();
	    localPrintWriter.close();
	    if (this.debug_level >= 2) {
	    	this.debug_output.close();
	    }
    }

	private void printDistribution() throws IOException {
		BufferedWriter localPrintWriter = new BufferedWriter(new FileWriter(this.outputFolder + "/" + this.outputPrefix + 
				".linker_alignment_score_distribution.txt", false));
	    for (int i = 0; i < this.scoreDistribution[0].length; i++) {
	    	localPrintWriter.write(i + "\t" + this.scoreDistribution[0][i] + "\t" + this.scoreDistribution[1][i]);
	    	localPrintWriter.newLine();
	    }
	    localPrintWriter.close();
	    localPrintWriter = new BufferedWriter(new FileWriter(this.outputFolder + "/" + this.outputPrefix + ".linker_alignment_score_difference_distribution.txt", 
	    		false));
	    for (int i = 0; i < this.secondBestScoreDiffDistribution[0].length && i <= this.maxSecondBestScoreDiff; i++) {
	    	localPrintWriter.write(i + "\t" + this.secondBestScoreDiffDistribution[0][i] + "\t" + this.secondBestScoreDiffDistribution[1][i]);
	    	localPrintWriter.newLine();
	    }
	    localPrintWriter.close();
	    localPrintWriter = new BufferedWriter(new FileWriter(this.outputFolder + "/" + this.outputPrefix + ".tag_length_distribution.txt", false));
	    for (int i = 0; i < this.tagLengthDistribution[0].length && i <= this.maxRealTagLength; i++) {
	    	localPrintWriter.write(i + "\t" + this.tagLengthDistribution[0][i] + "\t" + this.tagLengthDistribution[1][i]);
	    	localPrintWriter.newLine();
	    }
	    localPrintWriter.close();
	    localPrintWriter = new BufferedWriter(new FileWriter(this.outputFolder + "/" + this.outputPrefix + ".linker_composition_distribution.txt", false));
	    String str1 = null;
	    String str2 = null;
	    for (int j = 0; j <= this.nLinkers; j++) {
	    	for (int k = 0; k <= this.nLinkers; k++) {
		        if (j == 0) {
		          str1 = "A";
		        } else if (j == 1) {
		          str1 = "B";
		        } else if (j == 2) {
		          str1 = "X";
		        }
		        if (k == 0) {
		          str2 = "A";
		        } else if (k == 1) {
		          str2 = "B";
		        } else if (k == 2) {
		          str2 = "X";
		        }
		        if ((str2 != "X") || (str1 != "X")) {
		          localPrintWriter.write(str1 + "_" + str2 + "\t");
		        }
	    	}
	    }
	    localPrintWriter.write("Ambiguous\tTotal");
	    localPrintWriter.newLine();
	    int j = 0;
	    for (int k = 0; k <= this.nLinkers; k++) {
	    	for (int m = 0; m <= this.nLinkers; m++) {
	    		if (k == this.nLinkers && m == this.nLinkers) {
	    			break;
	    		}
		        localPrintWriter.write(this.linkerCompositionDistribution[k][m] + "\t");
		        j += this.linkerCompositionDistribution[k][m];
	    	}
	    }
	    j += this.nAmbiguousLinkerComposition;
	    localPrintWriter.write(this.nAmbiguousLinkerComposition + "\t" + j);
	    localPrintWriter.newLine();
	    for (int k = 0; k <= this.nLinkers; k++) {
	    	for (int m = 0; m <= this.nLinkers; m++) {
	    		if (k == this.nLinkers && m == this.nLinkers) {
	    			break;
	    		}
		        localPrintWriter.write(String.format("%.2f", new Object[] {
		        		Double.valueOf(100.0D * this.linkerCompositionDistribution[k][m] / j) 
		        		}
		        ) + "%" +"\t");
	    	}
	    }
	    localPrintWriter.write(String.format("%.2f", new Object[] {
	    		Double.valueOf(100.0D * this.nAmbiguousLinkerComposition / j)
	    		}
	    ) + "%" + "\t" + "100%");
	    localPrintWriter.close();
    }
	//
	public void setnPETs(long nPETs) {
    	this.nPETs = nPETs;
    }
    
    public long getnPETs() {
    	return this.nPETs;
    }
    
    public int getnLinkers() {
    	return nLinkers;
    }
    
    public int getoutput_data_with_ambiguous_linker_info() {
    	return output_data_with_ambiguous_linker_info;
    }
    
    public int getdebug_level() {
    	return debug_level;
    }

    public BufferedWriter getdebug_output() {
    	return debug_output;
    }
    // processOnePET
    public int getminimum_linker_alignment_score() {
    	return minimum_linker_alignment_score;
    }
    
    public int getminimum_tag_length() {
    	return minimum_tag_length;
    }
    
    public int getmaximum_tag_length() {
    	return maximum_tag_length;
    }
    
    public int getminSecondBestScoreDiff() {
    	return minSecondBestScoreDiff;
    }

    public void setlinkerCompositionDistribution(int i, int j, int n) {
    	this.linkerCompositionDistribution[i][j] = n;
    }
    
    public int getlinkerCompositionDistribution(int i, int j) {
    	return linkerCompositionDistribution[i][j];
    }
    
    public int getflip_head_tag() {
    	return flip_head_tag;
    }
    
    public int getflip_tail_tag() {
    	return flip_tail_tag;
    }
    
    public void setnAmbiguousLinkerComposition(int nAmbiguousLinkerComposition) {
    	this.nAmbiguousLinkerComposition = nAmbiguousLinkerComposition;
    }
    
    public int getnAmbiguousLinkerComposition() {
    	return nAmbiguousLinkerComposition;
    }   
    // processOneSequence
    public String[] getlinkers() {
    	return linkers;
    }
    
    public int getbarCodeStart(int i) {
    	return barCodeStart[i];
    }
    
    public int getbarCodeLength(int i) {
    	return barCodeLength[i];
    }

    public int getmaxLinkerLength() {
    	return maxLinkerLength;
    }
    
    public void setscoreDistribution(int i, int j, int n) {
    	this.scoreDistribution[i][j] = n;
    }
    
    public int getscoreDistribution(int i, int j) {
    	return scoreDistribution[i][j];
    }
    
    public void setsecondBestScoreDiffDistribution(int i, int j, int n) {
    	this.secondBestScoreDiffDistribution[i][j] = n;
    }
    
    public int getsecondBestScoreDiffDistribution(int i, int j) {
    	return secondBestScoreDiffDistribution[i][j];
    }
    
    public void setmaxSecondBestScoreDiff(int maxSecondBestScoreDiff) {
    	this.maxSecondBestScoreDiff = maxSecondBestScoreDiff;
    }
    
    public int getmaxSecondBestScoreDiff() {
    	return maxSecondBestScoreDiff;
    }
    
    public void settagLengthDistribution(int i, int j, int n) {
    	this.tagLengthDistribution[i][j] = n;
    }
    
    public int gettagLengthDistribution(int i, int j) {
    	return tagLengthDistribution[i][j];
    }
    
    public void setmaxRealTagLength(int maxRealTagLength) {
    	this.maxRealTagLength = maxRealTagLength;
    }
    
    public int getmaxRealTagLength() {
    	return maxRealTagLength;
    }

	public int printallreads() {
		// TODO Auto-generated method stub
		return printallreads;
	}
}
