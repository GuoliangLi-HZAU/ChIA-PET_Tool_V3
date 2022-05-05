package LGL.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.Arrays;

/**
 *
 * @author changhd
 */
public class UniqueSam {
	

    public UniqueSam(String Inputfile, String OutputFile, String runmode) throws IOException {
        BufferedReader FileIn = new BufferedReader(new InputStreamReader(new FileInputStream(Inputfile)));
        new File(OutputFile).delete();
        PrintWriter fileOut = new PrintWriter(new FileOutputStream(new File(OutputFile)));
        String line;
        String[] lastlines = new String[5];
        Arrays.fill(lastlines, null);
        int AS = 0;
        int XS = 0;
        while ((line = FileIn.readLine()) != null) {
            if (!line.startsWith("@")) {
            	String[] fields = line.split("\t");
            	if(runmode.equalsIgnoreCase("HiChIP")) {
            		//0x100   256  SECONDARY      secondary alignment
            		if( (Integer.parseInt(fields[1]) & 0x100) != 0) {
            			continue;
            		}
            	}
                if (lastlines[0] == null) {
                    
                    lastlines[0] = fields[0];
                    lastlines[1] = fields[4];
                    lastlines[4] = line;
                    for (int i = 0; i < fields.length; i++) {
                        if (fields[i].startsWith("AS:i")) {
                            String[] ASs = fields[i].split(":");
                            lastlines[2] = ASs[2];
                        }
                        if (fields[i].startsWith("XS:i")) {
                            String[] XSs = fields[i].split(":");
                            lastlines[3] = XSs[2];
                        }
                    }
                } else {
                    //String[] fields = line.split("\t");
                    if (lastlines[0].equalsIgnoreCase(fields[0])) {
                        for (int i = 0; i < fields.length; i++) {
                            if (fields[i].startsWith("AS:i")) {
                                String[] ASs = fields[i].split(":");
                                AS = Integer.parseInt(ASs[2]);
                            }
                            if (fields[i].startsWith("XS:i")) {
                                String[] XSs = fields[i].split(":");
                                XS = Integer.parseInt(XSs[2]);
                            }
                        }
                        
                        if (AS - XS > Integer.parseInt(lastlines[2]) - Integer.parseInt(lastlines[3])) {
                            lastlines[2] = String.valueOf(AS);
                            lastlines[3] = String.valueOf(XS);
                            lastlines[0] = fields[0];
                            lastlines[1] = fields[4];
                            lastlines[4] = line;
                        } else if ((AS - XS == Integer.parseInt(lastlines[2]) - Integer.parseInt(lastlines[3])) && (Integer.parseInt(fields[4]) > Integer.parseInt(lastlines[1]))) {
                            lastlines[2] = String.valueOf(AS);
                            lastlines[3] = String.valueOf(XS);
                            lastlines[0] = fields[0];
                            lastlines[1] = fields[4];
                            lastlines[4] = line;
                        }
                    } else {
                        fileOut.println(lastlines[4]);
                        Arrays.fill(lastlines, null);
                        lastlines[0] = fields[0];
                        lastlines[1] = fields[4];
                        lastlines[4] = line;
                        for (int i = 0; i < fields.length; i++) {
                            if (fields[i].startsWith("AS:i")) {
                                String[] ASs = fields[i].split(":");
                                lastlines[2] = ASs[2];
                            }
                            if (fields[i].startsWith("XS:i")) {
                                String[] XSs = fields[i].split(":");
                                lastlines[3] = XSs[2];
                            }
                        }
                    }
                }
            } else {
                fileOut.println(line);
            }

        }
        fileOut.println(lastlines[4]);
        FileIn.close();
        fileOut.close();
    }
    
    public static void main(String[] args) throws IOException {
        if (args.length == 3) {
            new UniqueSam(args[0], args[1], args[2]);
        } else {
            System.out.println("Usage: java UniqueSam <inputFile> <OutputFile> <runmode>");
        }
    }

}
