package process;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Calendar;
import java.util.Hashtable;
import java.util.Map;
import java.util.Vector;

/*
 * by momocoding
 * */

public class peak2interaction {

    int debugLevel = 4;
    Calendar rightNow = Calendar.getInstance();
    Hashtable<String, Integer> peaksmap1 = new Hashtable<String, Integer>();
    Hashtable<String, Integer> peaksmap2 = new Hashtable<String, Integer>();
    int minOverlapSize = 1;

    public peak2interaction(String regionFile, String outputFile) throws IOException {
        rightNow = Calendar.getInstance();
        System.out.println("[" + rightNow.getTime().toString() + "] start peak2interaction ... ");
        BufferedWriter outbuffer = new BufferedWriter(new FileWriter(outputFile, false));
        // load regions
        int intraloopN=0, interloopN=0;
        File f = new File(regionFile);
		if (f.exists()) {
		    try {
		    	BufferedReader reader = new BufferedReader(new FileReader(f));
		    	Vector<String> storelines = new Vector<String>();
		    	String line = reader.readLine();
		    	int ireads= 0, selfligationN = 0, samechrom = 0, distance = 0;
			    while (line != null) {
			    	line = line.trim();
			    	storelines.add(line);
			    	ireads++;
			    	line = reader.readLine();
			    	
			    	if (ireads==2)
			        {
			    		ireads = 0;
			            // split lines and determine bin
			    		String[] currEalla = storelines.get(0).split("[ \t]+");
			    		String[] currEallb = storelines.get(1).split("[ \t]+");
			    		String peakchrom1=null, peakstart1=null, peakend1=null, peakname1=null, peakdep1=null;
		                String peakchrom2=null, peakstart2=null, peakend2=null, peakname2=null, peakdep2=null;
		                if (currEalla.length>13){
		                	peakchrom1 = currEalla[10];
		                    peakstart1 = currEalla[11];
		                    peakend1   = currEalla[12];
		                    peakname1  = currEalla[13];
		                    //peakdep1  = currEalla[14];

		                    peakchrom2 = currEallb[10];
		                    peakstart2 = currEallb[11];
		                    peakend2   = currEallb[12];
		                    peakname2  = currEallb[13];
		                    //peakdep2  = currEallb[14];
		                    //System.out.println(peakchrom1 + " " + peakchrom2 + " " + peakstart2 +
		                    //		" " + peakname1 + " " + peakname2);
		                }else {
		                	System.out.println(currEalla.length);
		                }
		                if (!currEalla[1].equals(currEallb[1])){
	                        //cout << "Pair overlap ERROR!!!"<<endl;
	                        String tline = storelines.get(1);
	                        storelines.clear();
	                        storelines.add(tline);
	                        ireads = 1;
	                        continue;
		                }
		                if (peakname1.equals(peakname2)) {
	                        selfligationN++;
	                        storelines.clear();
	                        continue;
		                }
		                
		               // interaction name
		                String longname = peakchrom1 + "\t" + peakstart1 + "\t" + peakend1 +
		                		"\t0\t0\t0\t" + 
		                		peakchrom2 + "\t" + peakstart2 + "\t" + peakend2 + 
		                		"\t0\t0\t0"; // + peakname1+"\t"+peakname2+ "\t"+ peakdep1 + "\t"+ peakdep2
		                // add peak to map
		                if ( peakchrom1.equals(peakchrom2) && !peaksmap1.containsKey(longname)) {
		                    peaksmap1.put(longname, 1);
		                } else if ( peakchrom1.equals(peakchrom2)) {
		                    peaksmap1.put(longname, peaksmap1.get(longname).intValue()+1);
		                } else if ( !peaksmap2.containsKey(longname)) {
		                	peaksmap2.put(longname, 1);
		                } else {
		                	peaksmap2.put(longname, peaksmap2.get(longname).intValue()+1);
		                }// increment count
		                storelines.clear();
		             }
			    	
			    }
			    reader.close();
			    System.out.println("[size of pets] " + peaksmap1.size() + " " + peaksmap2.size());
			    // print out info, samechrom-0/1 distance score petindex
			    for(Map.Entry<String, Integer> entry: peaksmap1.entrySet()) {
			    	if(entry.getValue()>1) { //pet count >= 2
				    	outbuffer.write(entry.getKey()+"\t"+entry.getValue());
				    	outbuffer.newLine();
				    	intraloopN++;
			    	}
			    }
			    for(Map.Entry<String, Integer> entry: peaksmap2.entrySet()) {
			    	if(entry.getValue()>1) { //pet count >= 2
				    	outbuffer.write(entry.getKey()+"\t"+entry.getValue());
				    	outbuffer.newLine();
				    	interloopN++;
			    	}
			    }
			    outbuffer.close();
		    } catch (IOException e) {
			    e.printStackTrace();  
			}
		} else {
			System.out.println("Error: "+regionFile+" doesn't exist");
		}
    }




    public static void main(String[] args) throws IOException {
        if (args.length == 2) {
            new peak2interaction(args[0], args[1]);
        } else {
            System.out.println("Usage: java peak2interaction <input_aln_file> <out_file>");
            // huang 
            // System.out.println("       <extensionMode>: 1 - from 5' to 3';  2 - in both directions");
            System.out.println("       <extensionMode>: 1 - from 3' to 5';  2 - in both directions");
            System.exit(1);
        }
    }
}
