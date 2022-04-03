package LGL.chiapet;

/**
 * @author update by sun
 */
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

import multhread.BigFileProcess;
import process.Path;

public class LinkerFiltering_FastQ_PET {
    
    static int MAX_Linker_Number = 100;
    static int MIN_Linker_Length = 5;
    int debug_level = 1;
    String fastQFile1 = null;
    String fastQFile2 = null;
    String linker = null;
    String outputFolder = null;
    String outputPrefix = "out";
    String[] linkers = null;
    int nLinkers = 0;
    int[] barCodeStart  = null;
    int[] barCodeLength = null; // default: 4
    int maxLinkerLength;
    int minLinkerLength;
    int minimum_linker_alignment_score = 14; // minimum_linker_alignment score: default 14
    int minimum_tag_length = 18; // default: 18
    int minSecondBestScoreDiff = 3;
    int output_data_with_ambiguous_linker_info = 1; // default: output the data with ambiguous linker info
    int flip_head_tag = 0; // default: Not flip head tag
    int flip_tail_tag = 0; // default: Not flip tail tag  
    int[][] scoreDistribution;
    int[][] secondBestScoreDiffDistribution;
    int maxSecondBestScoreDiff = 0;   
    int maximum_tag_length = 1000; // default: 1000bp
    int maxRealTagLength = 0; // used to control the output length for the tag length distribution
    int[][] tagLengthDistribution;
    int[][] linkerCompositionDistribution;
    int nAmbiguousLinkerComposition = 0;
    BufferedWriter debug_output = null;
    
    String[] letter = {"A", "B", "C", "D", "E", "F"};
    
    private long nPETs = 0;
    private int nOutput = 10;
    private int iOutput = 0;
    private int threadNum = 0;
    static int maxThreadNum = 16;
    private Path p;
	public boolean search_all_linker = false;
	public String AutoLinker = "true";
    
    public LinkerFiltering_FastQ_PET(Path path) throws IOException {
    	p = path;
    	this.outputFolder = p.getOUTPUT_DIRECTORY()+"/"+p.getOUTPUT_PREFIX();
        this.outputPrefix = p.getOUTPUT_PREFIX();
        this.threadNum = Integer.valueOf(p.getNTHREADS());
        if(p.search_all_linker.equals("Y")) {
        	this.search_all_linker = true;
        }else {
        	this.search_all_linker = false;
        }
        this.AutoLinker = p.AutoLinker;
        
        if (this.threadNum > maxThreadNum) {
        	System.out.println("Error: parameter error. The maximun numbers of threads must <= " + maxThreadNum);
            System.exit(0);
        }

        readInputInformation();
        testInputInformation();
        readLinkers();
        outputRunningInfo();
        
        // initiate the distribution arrays
        // tag length distribution
        tagLengthDistribution = new int[2][maximum_tag_length];
        Arrays.fill(tagLengthDistribution[0], 0);
        Arrays.fill(tagLengthDistribution[1], 0);

        scoreDistribution = new int[2][this.maxLinkerLength + 1];
        Arrays.fill(scoreDistribution[0], 0);
        Arrays.fill(scoreDistribution[1], 0);
        
        secondBestScoreDiffDistribution = new int[2][this.maxLinkerLength * 2 + 1];
        Arrays.fill(secondBestScoreDiffDistribution[0], 0);
        Arrays.fill(secondBestScoreDiffDistribution[1], 0);
        
        filterSequenceByLinker();
        printDistribution();
    }

    public void readInputInformation() throws IOException {// 读取输入文件信息
    	this.fastQFile1 = p.getFastq_file_1();
    	this.fastQFile2 = p.getFastq_file_2();
    	this.linker = p.getLinker();
    	this.minimum_linker_alignment_score = Integer.parseInt(p.getMinimum_linker_alignment_score());
    	this.minimum_tag_length = Integer.parseInt(p.getMinimum_tag_length());
    	this.maximum_tag_length = Integer.parseInt(p.getMaximum_tag_length());
    	this.minSecondBestScoreDiff = Integer.parseInt(p.getMinSecondBestScoreDiff());
    	this.output_data_with_ambiguous_linker_info = Integer.parseInt(p.getOutput_data_with_ambiguous_linker_info());
    }
    
    public void testInputInformation() throws IOException {// 测试输入
        boolean requiredInputMissed = false;
        // test required input parameters
        if (this.fastQFile1==null) {
            System.out.println("No FASTQ file 1 input!");
            requiredInputMissed = true;
        }
        if (this.fastQFile2==null) {
            System.out.println("No FASTQ file 2 input!");
            requiredInputMissed = true;
        }
        if (this.linker == null) {
            System.out.println("No linker input!");
            requiredInputMissed = true;
        }
        // if required inputs are missed, exit
        if (requiredInputMissed == true) {
            System.exit(0);
        }
        // compare the maximum_tag_length and minimum_tag_length
        if (maximum_tag_length < minimum_tag_length) {
            System.out.println("maximum_tag_length is smaller than the minimum_tag_length! Stop ...");
            System.exit(0);
        }
        // test optional input parameters
        if (this.outputFolder==null) {
            System.out.println("No output folder! Set it as \"out\"");
            this.outputFolder = "out";
        }
        // test whether outputFolder exists
        File fOutputFolder = new File(this.outputFolder);
        if (!fOutputFolder.exists()) { // no this folder
            fOutputFolder.mkdir(); // create this folder
        } else if (!fOutputFolder.isDirectory()) { // a file exists with the same name, but not a folder
            System.out.println(this.outputFolder + " exists, but it is not a directory!!! please check...");
            System.exit(0);
        }
        this.outputFolder = fOutputFolder.getPath();
    }
    
    public void outputRunningInfo() throws IOException {// 输出信息
    	BufferedWriter fileOut = new BufferedWriter(new FileWriter(this.outputFolder + "/" + this.outputPrefix +
    			".runningInformation.LinkerFiltering_FastQ_PET.txt", false));
    	fileOut.write("Fastq file 1\t" + new File(this.fastQFile1).getAbsolutePath());
    	fileOut.newLine();
        fileOut.write("Fastq file 2\t" + new File(this.fastQFile2).getAbsolutePath());
        fileOut.newLine();
        for (int i = 0; i < linkers.length; i++) {
            fileOut.write("Linker" + letter[i] + "\t" + this.linkers[i] + " " + this.barCodeStart[i] + " " + this.barCodeLength[i]);
            fileOut.newLine();
        }
        fileOut.write("Output folder\t" + new File(this.outputFolder).getAbsolutePath());
        fileOut.newLine();
        fileOut.write("Output prefix\t" + this.outputPrefix);
        fileOut.newLine();
        fileOut.write("Minimum linker alignment score\t" + this.minimum_linker_alignment_score);
        fileOut.newLine();
        fileOut.write("minimum tag length\t" + this.minimum_tag_length);
        fileOut.newLine();
        fileOut.write("maximum tag length\t" + this.maximum_tag_length);
        fileOut.newLine();
        fileOut.write("Minumun SecondBestScore difference\t" + this.minSecondBestScoreDiff);
        fileOut.newLine();
        fileOut.write("Output data with ambiguous linker info\t" + this.output_data_with_ambiguous_linker_info);
        /*fileOut.newLine();
        fileOut.write("flip_head_tag\t" + this.flip_head_tag);
        fileOut.newLine();
        fileOut.write("flip_tail_tag\t" + this.flip_tail_tag);
        fileOut.newLine();
        fileOut.write("debug_level\t" + this.debug_level);*/
        fileOut.close();
    }
    
    public void readLinkers() throws IOException {
    	BufferedReader fileIn = new BufferedReader(new FileReader(this.linker));
        ArrayList<String> tempLinkers = new ArrayList<String>();
        String line;
        while ((line = fileIn.readLine()) != null) {
        	if(line.startsWith("Linker") || line.startsWith("Seq")) {
        		continue;
        	}
            if (line.length() <= 1) // skip the short lines
            {
                continue;
            }
            tempLinkers.add(line.trim());
        }
        fileIn.close();
        if(tempLinkers.size()==1) {
	    	String seq2 = LGL.util.SeqUtil.revComplement(tempLinkers.get(0));
	    	tempLinkers.add(seq2);
	    }
        
        this.linkers = new String[tempLinkers.size()];
        this.maxLinkerLength = 0;
        this.minLinkerLength = Integer.MAX_VALUE;
        this.barCodeStart  = new int[tempLinkers.size()];
        this.barCodeLength = new int[tempLinkers.size()];
        for (int i = 0; i < linkers.length; i++) {
        	this.linkers[i] = tempLinkers.get(i);
            // shortest linker length
            if (this.minLinkerLength > this.linkers[i].length()) {
                this.minLinkerLength = this.linkers[i].length();
            }
            // longest linker length
            if (this.maxLinkerLength < this.linkers[i].length()) {
                this.maxLinkerLength = this.linkers[i].length();
            }
        }
        barCodeStart[0] = -1; 
        barCodeLength[0] = 0;
        for (int i = 0, j = 0; i < linkers[0].length() && j < linkers[1].length(); i++, j++) {
        	char base1 = linkers[0].charAt(i);
			char base2 = linkers[1].charAt(j);
			if (base1 != base2) {
				if (barCodeStart[0] == -1) {
					barCodeStart[0] = i + 1;
					barCodeStart[1] = barCodeStart[0];
				}
				barCodeLength[0]++;
				barCodeLength[1]++;
			}
        }
        nLinkers = linkers.length;
        if (nLinkers <= 0) {
            System.out.println("No linker sequence information.\tStop!!!");
            System.exit(0);
        }
        if (nLinkers > MAX_Linker_Number) {
            System.out.println("More than " + MAX_Linker_Number + " linkers. Too many linkers. Please check...\t" + "Stop!!!");
            System.exit(0);
        }
        if (this.minLinkerLength < MIN_Linker_Length) {
            System.out.println("the shortest linker is less than " + MIN_Linker_Length + "bp.\tStop!!!");
            System.exit(0);
        }
        // create linker composition distribution matrix, and initiate the values to 0
        linkerCompositionDistribution = new int[nLinkers][nLinkers];
        for (int i = 0; i < nLinkers; i++) {
            Arrays.fill(linkerCompositionDistribution[i], 0);
        }
    }
    
    public void filterSequenceByLinker() throws IOException {
        // for every linker pair, output a pair of FastQ files
    	BufferedWriter[] fileOut = new BufferedWriter[nLinkers * nLinkers * 2 + 2];
        for (int i = 0; i < nLinkers; i++) {
            for (int j = 0; j < nLinkers; j++) {
                fileOut[(i * nLinkers + j) * 2] = new BufferedWriter(new FileWriter(this.outputFolder + "/" + this.outputPrefix + "." + (i + 1) + "_" + (j + 1) +
                		".R1.fastq", false));
                fileOut[(i * nLinkers + j) * 2 + 1] = new BufferedWriter(new FileWriter(this.outputFolder + "/" + this.outputPrefix + "." + (i + 1) + "_" + 
                		(j + 1) + ".R2.fastq", false));
            }
        }
        if (this.output_data_with_ambiguous_linker_info == 1) {
            fileOut[nLinkers * nLinkers * 2] = new BufferedWriter(new FileWriter(this.outputFolder + "/" + this.outputPrefix + ".ambiguous.R1.fastq", false));
            fileOut[nLinkers * nLinkers * 2 + 1] = new BufferedWriter(new FileWriter(this.outputFolder + "/" + this.outputPrefix + ".ambiguous.R2.fastq", false));
        }
        if (debug_level >= 2) {
            debug_output = new BufferedWriter(new FileWriter(this.outputFolder + "/" + this.outputPrefix + "debug.LinkerFiltering_FastQ_PET.txt", false));
        }
  
        //File fastq1 = new File(this.fastQFile1);
        //File fastq2 = new File(this.fastQFile2);
        String[] fastq1s = this.fastQFile1.split(",");
	    String[] fastq2s = this.fastQFile2.split(",");
	    for(int jk = 0; jk < fastq1s.length; jk++) {
	    	File fastq1 = new File(fastq1s[jk]);
	        File fastq2 = new File(fastq2s[jk]);
	        BigFileProcess bfp = new BigFileProcess(fastq1, fastq2, this, fileOut, threadNum);
	        bfp.start();
	        bfp.join();
	    }
        

        for (int i = 0; i < nLinkers; i++) {
            for (int j = 0; j < nLinkers; j++) {
                fileOut[(i * nLinkers + j) * 2].close();
                fileOut[(i * nLinkers + j) * 2 + 1].close();
            }
        }
        if (this.output_data_with_ambiguous_linker_info == 1) {
            fileOut[nLinkers * nLinkers * 2].close();
            fileOut[nLinkers * nLinkers * 2 + 1].close();
        }
        new File(this.outputFolder + "/" + this.outputPrefix + ".basic_statistics.txt").delete();
        BufferedWriter basicOut = new BufferedWriter(new FileWriter(this.outputFolder + "/" + this.outputPrefix + ".basic_statistics.txt", false));
        basicOut.write("Total PETs\t" + nPETs);
        basicOut.newLine();
        basicOut.close();
        if (debug_level >= 2) {
            debug_output.close();
        }
    }
    
    private void printDistribution() throws IOException {
    	BufferedWriter fileOut = new BufferedWriter(new FileWriter(this.outputFolder + "/" + this.outputPrefix + ".linker_alignment_score_distribution.txt", 
    			false));
        //fileOut.println("Distribution of linker alignment scores");
        for (int i = 0; i < scoreDistribution[0].length; i++) {
            fileOut.write(i + "\t" + scoreDistribution[0][i] + "\t" + scoreDistribution[1][i]);
            fileOut.newLine();
        }
        fileOut.close();
        fileOut = new BufferedWriter(new FileWriter(this.outputFolder + "/" + this.outputPrefix + ".linker_alignment_score_difference_distribution.txt", false));
        //fileOut.println("\nDistribution of alignment score differences between best alignment and second best
        //alignment");
        for (int i = 0; i < secondBestScoreDiffDistribution[0].length && i <= maxSecondBestScoreDiff; i++) {
            fileOut.write(i + "\t" + secondBestScoreDiffDistribution[0][i] + "\t" + secondBestScoreDiffDistribution[1][i]);
            fileOut.newLine();
        }
        fileOut.close();
        fileOut = new BufferedWriter(new FileWriter(this.outputFolder + "/" + this.outputPrefix + ".tag_length_distribution.txt", false));
        //fileOut.println("\nDistribution of tag lengths");
        for (int i = 0; i < tagLengthDistribution[0].length && i <= maxRealTagLength; i++) {
            fileOut.write(i + "\t" + tagLengthDistribution[0][i] + "\t" + tagLengthDistribution[1][i]);
            fileOut.newLine();
        }
        fileOut.close();
        fileOut = new BufferedWriter(new FileWriter(this.outputFolder + "/" + this.outputPrefix + ".linker_composition_distribution.txt", false));
        //fileOut.println("\nDistribution of linker compositions");
        String lin1 = "#"; // default value, currently only two half-linkers A and B are considered
        String lin2 = "#"; // default value, currently only two half-linkers A and B are considered
        for (int i = 0; i < nLinkers; i++) {
            for (int j = 0; j < nLinkers; j++) {
                if (i == 0) {
                    lin1 = "A";
                } else if (i == 1) {
                    lin1 = "B";
                }
                if (j == 0) {
                    lin2 = "A";
                } else if (j == 1) {
                    lin2 = "B";
                }

                fileOut.write(lin1 + "_" + lin2 + "\t");
            }
        }
        fileOut.write("Ambiguous\tTotal");
        fileOut.newLine();
        int total = 0;
        for (int i = 0; i < nLinkers; i++) {
            for (int j = 0; j < nLinkers; j++) {
                fileOut.write(linkerCompositionDistribution[i][j] + "\t");
                total += linkerCompositionDistribution[i][j];
            }
        }
        total += nAmbiguousLinkerComposition;
        fileOut.write(nAmbiguousLinkerComposition + "\t" + total);
        fileOut.newLine();
        for (int i = 0; i < nLinkers; i++) {
            for (int j = 0; j < nLinkers; j++) {
                fileOut.write(String.format("%.2f", (100.0 * linkerCompositionDistribution[i][j] / total)) + "%" + "\t");
            }
        }
        fileOut.write(String.format("%.2f", (100.0 * nAmbiguousLinkerComposition / total)) + "%" + "\t" + "100%");
        fileOut.close();
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
    // testBarcode
    public int getnOutput() {
    	return nOutput;
    }
    
    public int getiOutput() {
    	return iOutput;
    }
    
    public void setiOutput(int iOutput) {
    	this.iOutput = iOutput;
    }
}