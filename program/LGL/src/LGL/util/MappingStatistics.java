package LGL.util;

/**
 * @author update by sun
 */
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;

import multhread.BigFileProcess;

public class MappingStatistics {

    private int scorecutoff = 30;
    private int threadNum = 0;
    static int maxThreadNum = 16;
    private int[][] statistics;

    public MappingStatistics(String samFile1, String samFile2, String outputPrefix, String mappingScoreCutoff, String threadNum) throws IOException {
    	new File(outputPrefix + ".bedpe").delete();
    	new File(outputPrefix + ".mapping_statistics.txt").delete();
    	BufferedWriter outBedpe = new BufferedWriter(new FileWriter(outputPrefix + ".bedpe"));
    	BufferedWriter  outStatistics = new BufferedWriter(new FileWriter(outputPrefix + ".mapping_statistics.txt"));
    	this.scorecutoff = Integer.parseInt(mappingScoreCutoff);
    	this.threadNum = Integer.valueOf(threadNum);
        if (this.threadNum > maxThreadNum) {
        	System.out.println("Error: parameter error. The maximun numbers of threads must <= " + maxThreadNum);
            System.exit(0);
        }
        statistics = new int[3][3];
        Arrays.fill(statistics[0], 0);
        Arrays.fill(statistics[1], 0);
        Arrays.fill(statistics[2], 0);
        //System.out.println("XSCSCSCSCssssss");
    	BigFileProcess bfp = new BigFileProcess(samFile1, samFile2, scorecutoff, this.threadNum, outBedpe, this);
    	bfp.start_mapping();
        bfp.join_mapping();
        outStatistics.write("\tNon-mappable\tUniquely-mapped\tOthers");
        outStatistics.newLine();
        String[] label = new String[3];
        label[0] = "Non-mappable";
        label[1] = "Uniquely-mapped";
        label[2] = "Others";
        for (int i = 0; i < 3; i++) {
            outStatistics.write(label[i] + "\t");
            for (int j = 0; j < 2; j++) {
                outStatistics.write(statistics[i][j] + "\t");
            }
            outStatistics.write(statistics[i][2] + "");
            outStatistics.newLine();
        }
        outStatistics.close();
    }

	public void addStatistics(int[][] statistics) {
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				this.statistics[i][j] = this.statistics[i][j] + statistics[i][j];
			}
		}
	}

	public static void main(String[] args) throws IOException {
        if (args.length == 5) {
            new MappingStatistics(args[0], args[1], args[2], args[3], args[4]);
        } else {
            System.out.println("Usage: java MappingStatistics <samFile1> <samFile2> <outputPrefix> <mappingScoreCutoff> <threadNum>");
            System.out.println("       <samFile1>: String, short read mapping data in SAM format");
            System.out.println("       <samFile2>: String, short read mapping data in SAM format");
            System.out.println("       <outputPrefix>: String, prefix for output statistics and SAM files");
            System.out.println("       <mappingScoreCutoff>: Integer, mapping score cutoff for filtering the low quality reads");
            System.out.println("       <threadNum>: Integer, number of thread");
        }
    }
}
