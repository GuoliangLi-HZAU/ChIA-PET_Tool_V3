package LGL.align;

import java.util.*;

public class Align {

    int debugLevel = 1;
    private int[][] scoreMatrix;
    private int score;
    private int maxI;
    private int maxJ;
    private int m;
    private int n;
    private int MatchScore = 1;
    private int MismatchScore = -1;
    private int IndelScore = -1;
    private String str1 = "";
    private String str2 = "";
    private String alignedStr1 = "";
    private String alignedStr2 = "";
    private String alignedStatus = "";

    public Align(int m, int n) {
        this.m = m;
        this.n = n;
        scoreMatrix = new int[m + 1][n + 1];
        init();
    }

    public void init() {
        initMatrix(0);
        setScore(0);
        setMaxI(0);
        setMaxJ(0);
        this.setStr1("");
        this.setStr2("");
        this.setAlignedStatus("");
        this.setAlignedStr1("");
        this.setAlignedStr2("");
    }

    public void initMatrix(int score) {
        for (int i = 0; i < getM(); i++) {
            Arrays.fill(getScoreMatrix()[i], score);//填充数组
        }
    }

    public void outputScoreMatrix() {//没用到
        int index1, index2;
        int length1 = this.getStr1().length();
        int length2 = this.getStr2().length();
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
                System.out.print("\t" + getScoreMatrix()[index1][index2]);//???越界
            }
            System.out.println();
        }
    }

    /**
     * @return the scoreMatrix
     */
    public int[][] getScoreMatrix() {
        return scoreMatrix;
    }

    /**
     * @param scoreMatrix the scoreMatrix to set
     */
    public void setScoreMatrix(int[][] scoreMatrix) {
        this.scoreMatrix = scoreMatrix;
    }

    /**
     * @return the score
     */
    public int getScore() {
        return score;
    }

    /**
     * @param score the score to set
     */
    public void setScore(int score) {
        this.score = score;
    }

    /**
     * @return the maxI
     */
    public int getMaxI() {
        return maxI;
    }

    /**
     * @param maxI the maxI to set
     */
    public void setMaxI(int maxI) {
        this.maxI = maxI;
    }

    /**
     * @return the maxJ
     */
    public int getMaxJ() {
        return maxJ;
    }

    /**
     * @param maxJ the maxJ to set
     */
    public void setMaxJ(int maxJ) {
        this.maxJ = maxJ;
    }

    /**
     * @return the m
     */
    public int getM() {
        return m;
    }

    /**
     * @param m the m to set
     */
    public void setM(int m) {
        this.m = m;
    }

    /**
     * @return the n
     */
    public int getN() {
        return n;
    }

    /**
     * @param n the n to set
     */
    public void setN(int n) {
        this.n = n;
    }

    /**
     * @return the MatchScore
     */
    public int getMatchScore() {
        return MatchScore;
    }

    /**
     * @param MatchScore the MatchScore to set
     */
    public void setMatchScore(int MatchScore) {
        this.MatchScore = MatchScore;
    }

    /**
     * @return the MismatchScore
     */
    public int getMismatchScore() {
        return MismatchScore;
    }

    /**
     * @param MismatchScore the MismatchScore to set
     */
    public void setMismatchScore(int MismatchScore) {
        this.MismatchScore = MismatchScore;
    }

    /**
     * @return the IndelScore
     */
    public int getIndelScore() {
        return IndelScore;
    }

    /**
     * @param IndelScore the IndelScore to set
     */
    public void setIndelScore(int IndelScore) {
        this.IndelScore = IndelScore;
    }

    /**
     * @return the str1
     */
    public String getStr1() {
        return str1;
    }

    /**
     * @param str1 the str1 to set
     */
    public void setStr1(String str1) {
        this.str1 = str1;
    }

    /**
     * @return the str2
     */
    public String getStr2() {
        return str2;
    }

    /**
     * @param str2 the str2 to set
     */
    public void setStr2(String str2) {
        this.str2 = str2;
    }

    /**
     * @return the alignedStr1
     */
    public String getAlignedStr1() {
        return alignedStr1;
    }

    /**
     * @param alignedStr1 the alignedStr1 to set
     */
    public void setAlignedStr1(String alignedStr1) {
        this.alignedStr1 = alignedStr1;
    }

    /**
     * @return the alignedStr2
     */
    public String getAlignedStr2() {
        return alignedStr2;
    }

    /**
     * @param alignedStr2 the alignedStr2 to set
     */
    public void setAlignedStr2(String alignedStr2) {
        this.alignedStr2 = alignedStr2;
    }

    /**
     * @return the alignedStatus
     */
    public String getAlignedStatus() {
        return alignedStatus;
    }

    /**
     * @param alignedStatus the alignedStatus to set
     */
    public void setAlignedStatus(String alignedStatus) {
        this.alignedStatus = alignedStatus;
    }
}
