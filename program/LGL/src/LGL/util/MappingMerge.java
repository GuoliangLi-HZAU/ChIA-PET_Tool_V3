package LGL.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Collections;
import java.util.Vector;


public class MappingMerge {
	

    public MappingMerge(String Inputfile, String OutputFile, int cutoff, String runmode) throws IOException {
        BufferedReader FileIn = new BufferedReader(new InputStreamReader(new FileInputStream(Inputfile)));
        new File(OutputFile+".pe.sam").delete();
        PrintWriter fileOut = new PrintWriter(new FileOutputStream(new File(OutputFile+".pe.sam")));
        new File(OutputFile+".bedpe").delete();
        PrintWriter bedpeoutf = new PrintWriter(new FileOutputStream(new File(OutputFile+".bedpe")));
        String line;
        //String[] lastlines = new String[5];
        //0 ID;
        // 1 2 3, flag mapQ readline
        //2 3 4 ... N, store map case
        Vector<ReadMap> lastlines = new Vector<ReadMap>();
        //Arrays.fill(lastlines, null);
        int Nmap = 0;
        int mapq = 0;
        String lastread = ""; int span = 0;
        String bedpeline = "";
        ReadMap readtmp1 = null, readtmp2 = null;
        int readlen = 0;
        Vector<ReadPair> storepepair = new Vector<ReadPair>();
        while ((line = FileIn.readLine()) != null) {
            if (!line.startsWith("@")) {
            	String[] fields = line.split("\t");
            	mapq = Integer.parseInt(fields[4]);
            	readlen = getreadlen(fields[5]);
            	ReadMap readtmp = new ReadMap(fields[0], Integer.parseInt(fields[1]), fields[2], Integer.parseInt(fields[3]), line, 
            			mapq, readlen);
            	
            	if(runmode.equalsIgnoreCase("HiChIP")) {
            		//0x100   256  SECONDARY      secondary alignment
            		//0x800	SUPPLEMENTARY .. supplementary alignment
            		if( (Integer.parseInt(fields[1]) & 0x100) != 0) {
            			//continue;
            		}
            	}
                if (lastlines.size() == 0) {
                	lastread = fields[0];
                	if(mapq>=cutoff) {
                		lastlines.add(readtmp);
                	}
                } else {
                    //String[] fields = line.split("\t");
                    if (lastread.equalsIgnoreCase(fields[0])) { // same read
                        // store
                    	if(mapq>=cutoff) {
                    		lastlines.add(readtmp);
                    	}
                    	lastread = fields[0];
                    } else { // new read, process last read
                        storepepair.clear();
                        if(lastlines.size()==1) {
                        	//single, print sam and not print bedpe
                        	fileOut.write(lastlines.get(0).getline()+"\n");
                        }else if(lastlines.size()==2) { //PE
                        	//print pe
                        	readtmp1 = lastlines.get(0);
                        	readtmp2 = lastlines.get(1);
                        	fileOut.write(readtmp1.getline()+"\n");
                        	fileOut.write(readtmp2.getline()+"\n");
                        	ReadPair readptmp = new ReadPair(readtmp1, readtmp2);
                        	if(runmode.equalsIgnoreCase("ChIA-PET")) {
                        		if(readptmp.getsup()==0) {
                        			bedpeoutf.write(readptmp.getbedpe()+"\n");
                        		}
                        	}else {
                        		bedpeoutf.write(readptmp.getbedpe()+"\n");
                        	}
                        }else if(lastlines.size()>2) {
                        	for(int i=0;i<lastlines.size();i++) {
                        		for(int j=i+1; j< lastlines.size(); j++) {
                        			ReadPair readptmp = new ReadPair(lastlines.get(i), lastlines.get(j));
                        			storepepair.add(readptmp);
                        		}
                        	}
                        	//sort by span
                        	Collections.sort(storepepair);
                        	int diffspan = 0, chrspan1 = 0, issamestrand = 0;
                        	ReadPair readptmp = storepepair.get(0);
                        	//print readptmp
                        	fileOut.write(readptmp.getline1()+"\n");
                        	fileOut.write(readptmp.getline2()+"\n");
                        	bedpeoutf.write(readptmp.getbedpe()+"\n");
                        	chrspan1 = readptmp.getspan();
                        	issamestrand = readptmp.issamestrand();
                        	for(int i=1;i<storepepair.size();i++) {
                        		diffspan = readptmp.getdiffspan(storepepair.get(i));
                        		readptmp = storepepair.get(i);
                        		if(diffspan < 1000 || readptmp.getspan() < 1000) {
                        			//here is correct or not??
                        			if(chrspan1<300 && issamestrand == 1 && readptmp.getspan() == -1) {
                        				//print readptmp
                            			fileOut.write(readptmp.getline1()+"\n");
                                    	fileOut.write(readptmp.getline2()+"\n");
                                    	bedpeoutf.write(readptmp.getbedpe()+"\n");
                                    	break;
                            		}
                        			continue;
                        		}else {
                        			//print readptmp
                        			fileOut.write(readptmp.getline1()+"\n");
                                	fileOut.write(readptmp.getline2()+"\n");
                                	bedpeoutf.write(readptmp.getbedpe()+"\n");
                        		}
                        	}
                        }
                        lastlines.clear();
                        if(mapq>=cutoff) {
                    		lastlines.add(readtmp);
                    	}
                        lastread = fields[0];
                    }
                }
            } else {
                fileOut.println(line);
                Nmap = 0;
            }

        }
        // last read
        storepepair.clear();
        if(lastlines.size()==1) {
        	//single, print sam and not print bedpe
        	fileOut.write(lastlines.get(0).getline()+"\n");
        }else if(lastlines.size()==2) { //PE
        	//print pe
        	readtmp1 = lastlines.get(0);
        	readtmp2 = lastlines.get(1);
        	fileOut.write(readtmp1.getline()+"\n");
        	fileOut.write(readtmp2.getline()+"\n");
        	ReadPair readptmp = new ReadPair(readtmp1, readtmp2);
        	bedpeoutf.write(readptmp.getbedpe()+"\n");
        }else if(lastlines.size()>2) {
        	for(int i=0;i<lastlines.size();i++) {
        		for(int j=i+1; j< lastlines.size(); j++) {
        			ReadPair readptmp = new ReadPair(lastlines.get(i), lastlines.get(j));
        			storepepair.add(readptmp);
        		}
        	}
        	//sort by span
        	Collections.sort(storepepair);
        	int diffspan = 0;
        	ReadPair readptmp = storepepair.get(0);
        	//print readptmp
        	fileOut.write(readptmp.getline1()+"\n");
        	fileOut.write(readptmp.getline2()+"\n");
        	bedpeoutf.write(readptmp.getbedpe()+"\n");
        	for(int i=1;i<storepepair.size();i++) {
        		diffspan = readptmp.getdiffspan(storepepair.get(i));
        		readptmp = storepepair.get(i);
        		if(diffspan < 1000) {
        			continue;
        		}else {
        			//print readptmp
        			fileOut.write(readptmp.getline1()+"\n");
                	fileOut.write(readptmp.getline2()+"\n");
                	bedpeoutf.write(readptmp.getbedpe()+"\n");
        		}
        	}
        }
        lastlines.clear();
        //end process last read
        FileIn.close();
        fileOut.close();
        bedpeoutf.close();
    }
    
    public int getreadlen(String cigar) {
    	if(cigar==null) return 0;
    	String readseq = "";
    	int readlen = 0;
    	for(int i=0;i<cigar.length();i++) {
    		if(cigar.charAt(i) >=48 && cigar.charAt(i) <=57) {
    			readseq += cigar.charAt(i);
    		}else if(cigar.charAt(i) == 'M') {
    			readlen += Integer.parseInt(readseq);
    			readseq = "";
    		}else if(cigar.charAt(i) == 'I') {
    			readlen += Integer.parseInt(readseq);
    			readseq = "";
    		}else {
    			readseq = "";
    		}
    	}
    	return readlen;
    }
    
    public static void main(String[] args) throws IOException {
        if (args.length == 4) {
            new MappingMerge(args[0], args[1], Integer.parseInt(args[2]), args[3]);
        } else {
            System.out.println("Usage: java MappingMerge <inputFile> <OutputFile> <cutoff> <runmode>");
        }
    }

}
