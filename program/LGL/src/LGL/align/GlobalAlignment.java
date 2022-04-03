package LGL.align;

import java.io.*;
import java.util.*;

public class GlobalAlignment {

    int debugLevel = 1;
    int[][] scoreMatrix;
    int score;
    int m;
    int n;
    int MatchScore = 2;
    int MismatchScore = -2;
    int IndelScore = -1;
    String str1;
    String str2;
    String alignedStr1 = "";
    String alignedStr2 = "";
    String alignedStatus = "";

    public GlobalAlignment(int m, int n) {
        this.m = m;
        this.n = n;
        scoreMatrix = new int[m + 1][n + 1];
        for (int i = 0; i < m; i++) {
            Arrays.fill(scoreMatrix[i], 0);
        }
    }

    public void initMatrix(int score) {
        for (int i = 0; i < m; i++) {
            Arrays.fill(scoreMatrix[i], score);
        }
    }

    public void align(String str1, String str2) {
        this.str1 = str1;
        this.str2 = str2;
        int length1 = str1.length();
        int length2 = str2.length();
        if (m < length1) {
            m = length1;
            scoreMatrix = new int[m + 1][n + 1];
        }
        if (n < length2) {
            n = length2;
            scoreMatrix = new int[m + 1][n + 1];
        }
        initMatrix(0);
        int index1, index2;
        for (index1 = 1; index1 <= length1; index1++) {//???从1开始
            scoreMatrix[index1][0] = index1 * MismatchScore;
        }
        for (index2 = 1; index2 <= length2; index2++) {
            scoreMatrix[0][index2] = index2 * MismatchScore;
        }

        for (index1 = 0; index1 < length1; index1++) {
            for (index2 = 0; index2 < length2; index2++) {
                if (str1.charAt(index1) == str2.charAt(index2)) {
                    scoreMatrix[index1 + 1][index2 + 1] = scoreMatrix[index1][index2] + MatchScore;
                } else {
                    scoreMatrix[index1 + 1][index2 + 1] = scoreMatrix[index1][index2] + MismatchScore;
                }
                if (scoreMatrix[index1 + 1][index2 + 1] < scoreMatrix[index1 + 1][index2] + IndelScore) {
                    scoreMatrix[index1 + 1][index2 + 1] = scoreMatrix[index1 + 1][index2] + IndelScore;
                }
                if (scoreMatrix[index1 + 1][index2 + 1] < scoreMatrix[index1][index2 + 1] + IndelScore) {
                    scoreMatrix[index1 + 1][index2 + 1] = scoreMatrix[index1][index2 + 1] + IndelScore;
                }
            }
        }
        score = scoreMatrix[length1][length2];
        System.out.println("length1 = " + length1);
        System.out.println("length2 = " + length2);

        System.out.print("i = 0:");
        for (index2 = 0; index2 <= length1; index2++) {
            System.out.print("\t" + index2);
        }
        System.out.println();
        for (index1 = 0; index1 <= length1; index1++) {
            System.out.print("i = " + index1 + ":");
            for (index2 = 0; index2 <= length2; index2++) {
                System.out.print("\t" + scoreMatrix[index1][index2]);
            }
            System.out.println();
        }

        alignedStr1 = "";
        alignedStr2 = "";
        alignedStatus = "";
        index1 = length1 - 1;
        index2 = length2 - 1;
        while ((index1 >= 0) || (index2 >= 0)) {
            if (index1 < 0) {
                alignedStr1 = alignedStr1 + "-";
                alignedStr2 = alignedStr2 + str2.charAt(index2);
                alignedStatus = alignedStatus + " ";
                index2--;
            } else if (index2 < 0) {
                alignedStr1 = alignedStr1 + str1.charAt(index1);
                alignedStr2 = alignedStr2 + "-";
                alignedStatus = alignedStatus + " ";
                index1--;
            } else if ((scoreMatrix[index1 + 1][index2 + 1] == scoreMatrix[index1][index2] + MatchScore) && (str1.charAt(index1) == str2.charAt(index2))) {
                alignedStr1 = alignedStr1 + str1.charAt(index1);
                alignedStr2 = alignedStr2 + str2.charAt(index2);
                alignedStatus = alignedStatus + "|";
                index1--;
                index2--;
            } else if (scoreMatrix[index1 + 1][index2 + 1] == scoreMatrix[index1][index2] + MismatchScore) {
                alignedStr1 = alignedStr1 + str1.charAt(index1);
                alignedStr2 = alignedStr2 + str2.charAt(index2);
                alignedStatus = alignedStatus + "X";
                index1--;
                index2--;
            } else if (scoreMatrix[index1 + 1][index2 + 1] == scoreMatrix[index1 + 1][index2] + IndelScore) {
                alignedStr1 = alignedStr1 + "-";
                alignedStr2 = alignedStr2 + str2.charAt(index2);
                alignedStatus = alignedStatus + " ";
                index2--;
            } else {
                //if (scoreMatrix[index1 + 1][index2 + 1] == scoreMatrix[index1][index2 + 1] + IndelScore)
                alignedStr1 = alignedStr1 + str1.charAt(index1);
                alignedStr2 = alignedStr2 + "-";
                alignedStatus = alignedStatus + " ";
                index1--;
            }
            System.out.println("index1 = " + (index1 + 1) + "; index2 = " + (index2 + 1));
        }
        //alignedStr1 = (new StringBuilder(alignedStr1)).reverse().toString();
        alignedStr1 = (new StringBuffer(alignedStr1)).reverse().toString();
        //alignedStr2 = (new StringBuilder(alignedStr2)).reverse().toString();
        alignedStr2 = (new StringBuffer(alignedStr2)).reverse().toString();
        //alignedStatus = (new StringBuilder(alignedStatus)).reverse().toString();
        alignedStatus = (new StringBuffer(alignedStatus)).reverse().toString();
        System.out.println("str1: " + str1);
        System.out.println("str2: " + str2);
        System.out.println("aligned score: " + score);
        System.out.println("aligned str1      : " + alignedStr1);
        System.out.println("aligned strStatus : " + alignedStatus);
        System.out.println("aligned str2      : " + alignedStr2);
    }

    public static void main(String[] args) throws IOException {
        /*File dir1 = new File (".");
        File dir2 = new File ("..");
        try {
        System.out.println ("Current dir : " + dir1.getCanonicalPath());
        System.out.println ("Parent  dir : " + dir2.getCanonicalPath());
        }
        catch(Exception e) {
        e.printStackTrace();
        }*/
        if (args.length == 2) {
            GlobalAlignment globalAligner = new GlobalAlignment(args[0].length(), args[1].length());
            //globalAligner.align("AACCGGTT", "ACCGTATT");
            globalAligner.align(args[0], args[1]);
        } else {
            System.out.println("Usage: java GlobalAlignment <sequence 1> <sequence 2>");
        }
    }
}
