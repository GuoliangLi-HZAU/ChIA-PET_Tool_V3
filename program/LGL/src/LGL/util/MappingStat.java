package LGL.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.Arrays;

/**
 *
 * @author changhd
 */
public class MappingStat {

    int classify;
    int scorecutoff = 20;
    int[][] statistics;
    String[] label;
    String header1 = "";
    String header2 = "";
    String sequences1 = "";
    String sequences2 = "";

    public MappingStat(String InputFile1, String InputFile2, String OutputPrefix, String cutoff) throws Exception {
        this.scorecutoff = Integer.parseInt(cutoff);
        BufferedReader fileIn1 = new BufferedReader(new InputStreamReader(new FileInputStream(InputFile1)));
        BufferedReader fileIn2 = new BufferedReader(new InputStreamReader(new FileInputStream(InputFile2)));
        new File(OutputPrefix + ".mapping_statistics.txt").delete();
        new File(OutputPrefix + ".bedpe.selected.uniq.1.sam").delete();
        new File(OutputPrefix + ".bedpe.selected.uniq.2.sam").delete();
        PrintWriter fileOut1 = new PrintWriter(new FileOutputStream(new File(OutputPrefix + ".mapping_statistics.txt")));
        PrintWriter fileOut2 = new PrintWriter(new FileOutputStream(new File(OutputPrefix + ".bedpe.selected.uniq.1.sam")));
        PrintWriter fileOut3 = new PrintWriter(new FileOutputStream(new File(OutputPrefix + ".bedpe.selected.uniq.2.sam")));
        String line1;
        String line2;
        int mapping1 = -1;
        int mapping2 = -1;
        statistics = new int[3][3];
        Arrays.fill(statistics[0], 0);
        Arrays.fill(statistics[1], 0);
        while ((line1 = fileIn1.readLine()) != null && (line2 = fileIn2.readLine()) != null) {
            if (!line1.startsWith("@")) {
                String fields1[] = line1.split("\t");
                String fields2[] = line2.split("\t");
                if (!(fields1[0].equalsIgnoreCase(fields2[0]))) {
                    System.out.println("Error in lines...");
                } else {
                    for (int i = 0; i < fields1.length; i++) {
                        if (fields1[i].startsWith("AS:i")) {
                            mapping1 = classify2(line1);
                            break;
                        } else {
                            mapping1 = -1;
                        }
                    }
                    for (int i = 0; i < fields2.length; i++) {
                        if (fields2[i].startsWith("AS:i")) {
                            mapping2 = classify2(line2);
                            break;
                        } else {
                            mapping2 = -1;
                        }
                    }
                    if (mapping1 == -1) {
                        mapping1 = classify1(line1);
                    }
                    if (mapping2 == -1) {
                        mapping2 = classify1(line2);
                    }
                }
                if (mapping1 == 0 && mapping2 == 0) {
                    statistics[0][0]++;
                } else if (mapping1 == 1 && mapping2 == 0) {
                    statistics[1][0]++;
                } else if (mapping1 == 2 && mapping2 == 0) {
                    statistics[2][0]++;
                } else if (mapping1 == 0 && mapping2 == 1) {
                    statistics[0][1]++;
                } else if (mapping1 == 1 && mapping2 == 1) {
                    statistics[1][1]++;
                    fileOut2.println(line1);
                    fileOut3.println(line2);
                } else if (mapping1 == 2 && mapping2 == 1) {
                    statistics[2][1]++;
                } else if (mapping1 == 0 && mapping2 == 2) {
                    statistics[0][2]++;
                } else if (mapping1 == 1 && mapping2 == 2) {
                    statistics[1][2]++;
                } else {
                    statistics[2][2]++;
                }
            } else {
                fileOut2.println(line1);
                fileOut3.println(line2);
            }
        }
        label = new String[3];
        label[0] = "Non-mappable";
        label[1] = "unique-mapped";
        label[2] = "multiple-mapped";
        fileOut1.println("\tNon-mappable\tunique-mapped\tmultiple-mapped");
        for (int i = 0; i < 3; i++) {
            fileOut1.print(label[i] + "\t");
            for (int j = 0; j < 2; j++) {
                fileOut1.print(statistics[i][j] + "\t");
            }
            fileOut1.print(statistics[i][2]);
            fileOut1.println();
        }
        fileIn1.close();
        fileIn2.close();
        fileOut1.close();
        fileOut2.close();
        fileOut3.close();
    }

    public int classify1(String line) {
        int ok1 = -1;
        int ok2 = -1;
        int ok3 = -1;
        String[] fields = line.split("\t");
        for (int i = 0; i < fields.length; i++) {
            if (fields[i].equalsIgnoreCase("XT:A:U")) {
                ok1 = 1;
            } else if (fields[i].equalsIgnoreCase("X0:i:1")) {
                ok2 = 1;
            } else if (fields[i].equalsIgnoreCase("X1:i:0")) {
                ok3 = 1;
            }
        }
        int classify1 = -1;
        if (fields.length <= 11 || ((Integer.parseInt(fields[1]) & 0x4) != 0)) {
            classify1 = 0;
        } else if (ok1 == 1 && (Integer.parseInt(fields[4]) >= scorecutoff) && ok2 == 1 && ok3 == 1 && ((Integer.parseInt(fields[1]) & 0x4) == 0)) {
            classify1 = 1;
        } else {
            classify1 = 2;
        }
        return classify1;
    }

    public int classify2(String line) {
        String[] fields = line.split("\t");
        int classify2 = -1;
        String[] ASs = new String[3];
        String[] XSs = new String[3];
        Arrays.fill(ASs, null);
        Arrays.fill(XSs, null);
        for (int i = 0; i < fields.length; i++) {
            if (fields[i].startsWith("AS:i")) {
                ASs = fields[i].split(":");
            }
            if (fields[i].startsWith("XS:i")) {
                XSs = fields[i].split(":");
            }
        }
        if ((Integer.parseInt(fields[4]) > scorecutoff) && (Integer.parseInt(ASs[2]) - Integer.parseInt(XSs[2]) > 20) && 
        		((Integer.parseInt(fields[1]) & 0x4) == 0)) {
            classify2 = 1;
        } else if (Integer.parseInt(fields[4]) == 0 || ((Integer.parseInt(fields[1]) & 0x4) != 0)) {
            classify2 = 0;
        } else {
            classify2 = 2;
        }
        return classify2;
    }

    public static void main(String[] args) throws Exception {
        if (args.length == 4) {
            new MappingStat(args[0], args[1], args[2], args[3]);
        } else {
            System.out.println("Usage: java MappingStat <inputFile1> <inputFile1> <outputPrefix> <mappingCutOff>");
        }
    }
}
