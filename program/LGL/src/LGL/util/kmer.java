package LGL.util;

import java.io.*;
import java.text.DecimalFormat;
import java.util.*;
import java.util.zip.GZIPInputStream;

import LGL.align.LocalAlignment;
import LGL.data.FastQMulthread;

public class kmer {
	static int mLengthSum = 0;
	static int mReads = 0;
	
	static int KMR = 0;// 0x3FC // 0011 1111 1100, 5, 0101, 32 byte
	static int Klen = 5;
	static int kit = 0x3; //0011
	static long NprocessRead = 1000000;
	static int headS = 0; //10
	static int tailE = 0;
	static double percent = 0.1;
	static String linkerfile = "";
	static int linkerScore = 14;
	static int minLinkerLen = 12;
	static int skipLine = 1000000;
	static int even = 0;
	static String kmerfileprefix = "kmer";
	static int maxKmer = 1000;
	
	public static void main(String[] args) throws FileNotFoundException, IOException, InterruptedException {
		String fastQFile1 = "";
		String fastQFile2 = "";
		
		for(int i=0;i<args.length; i++) {
			if(args[i].equals("--fq1")) {
				fastQFile1 = args[i+1];
			}else if(args[i].equals("--fq2")) {
				fastQFile2 = args[i+1];
			}else if(args[i].equals("-k")) {
				Klen = Integer.parseInt(args[i+1]);
			}else if(args[i].equals("-s")) {
				headS = Integer.parseInt(args[i+1]);
			}else if(args[i].equals("-e")) {
				tailE = Integer.parseInt(args[i+1]);
			}else if(args[i].equals("-N")) {
				NprocessRead = Integer.parseInt(args[i+1]);
			}else if(args[i].equals("-L")) {
				linkerfile = args[i+1];
			}else if(args[i].equals("--score")) {
				linkerScore = Integer.parseInt(args[i+1]);
			}else if(args[i].equals("-p")) {
				percent = Double.parseDouble(args[i+1]);
			}else if(args[i].equals("--minlinker")) {
				minLinkerLen = Integer.parseInt(args[i+1]);
			}else if(args[i].equals("--skip")) {
				skipLine = Integer.parseInt(args[i+1]);
			}else if(args[i].equals("--even")) {
				even = Integer.parseInt(args[i+1]);
			}else if(args[i].equals("--prefix")) {
				kmerfileprefix = args[i+1];
			}else if(args[i].equals("--maxkmer")) {
				maxKmer = Integer.parseInt(args[i+1]);
			}
		}
		if(Klen > 15) {
			System.out.println("Please provide a kmer length shoter than 15");
		}else if(Klen % 2 == 1 || even == 1) {
			for(int i=0;i<Klen-1;i++) {
				KMR = (KMR | kit) << 2;
			}
			
			System.out.println("Kmer unit is " + Integer.toBinaryString(KMR));
		}else {
			System.out.println("Please provide a odd value of kmer length!!!");
		}
		if(Klen<=5 && percent < 0.2) {
			System.out.println("Warining: Kmer too short and percent too small, we changed percent from " + percent + " to 0.2 !!!");
		}
		linkerfile = kmerfileprefix + ".linkerprocess.txt";
		//0x3fffffff 1073741823 0011 1111 1111 1111 1111 1111 1111 1111
		//max 15 kmer, too many memory
		// ex. 5, 10
		//1000 0000 0000 0000 0000 0000
		//0011 1111 1111 1111 1111 1100
		int mKmerBufLen = 2<<(Klen * 2 - 1);
		//even
		//0001 0000 0000
		//0000 1111 1100
		int[] mKmer1 = null;// = new long[1073741824];
		int[] mKmer2 = null;// = new long[1073741824];
		//
		//static ConcurrentHashMap<Integer, Integer> mKmer1 = new ConcurrentHashMap<Integer, Integer>();
		//static ConcurrentHashMap<Integer, Integer> mKmer2 = new ConcurrentHashMap<Integer, Integer>();
		System.out.println("mKmerBufLen " + mKmerBufLen);
		mKmer1 = new int[mKmerBufLen];
		mKmer2 = new int[mKmerBufLen];
		
		//int[] randomLine = RandomInt(1000000, 2000000, 10000);
		if(!fastQFile1.equals("") && !fastQFile2.equals("")) {
			String[] fastq1s = fastQFile1.split(",");
		    String[] fastq2s = fastQFile2.split(",");
		    for(int jk = 0; jk < fastq1s.length; jk++) {
		    	File fastq1 = new File(fastq1s[jk]);
		        File fastq2 = new File(fastq2s[jk]);
		        //process 
		        Calendar rightNow = Calendar.getInstance();
		        System.out.println("[" + rightNow.getTime().toString() +"]" + fastq1s.length + ": " + fastq1s[jk] + ", " + fastq2s[jk]);
		        getfastq(fastq1, fastq2, mKmer1, mKmer2, skipLine);
		    }
		}else if(!fastQFile1.equals("")) {
			String[] fastq1s = fastQFile1.split(",");
		    for(int jk = 0; jk < fastq1s.length; jk++) {
		    	File fastq1 = new File(fastq1s[jk]);
		        //process 
		        Calendar rightNow = Calendar.getInstance();
		        System.out.println("[" + rightNow.getTime().toString() +"]" + fastq1s.length + ": " + fastq1s[jk] );
		        getfastq(fastq1, mKmer1, skipLine);
		    }
		}else {
			System.out.println("please input --fastq1");
			System.exit(0);
		}
		if(mLengthSum<1000 && NprocessRead>1000) {
			skipLine = 100;
			if(!fastQFile1.equals("") && !fastQFile2.equals("")) {
				String[] fastq1s = fastQFile1.split(",");
			    String[] fastq2s = fastQFile2.split(",");
			    for(int jk = 0; jk < fastq1s.length; jk++) {
			    	File fastq1 = new File(fastq1s[jk]);
			        File fastq2 = new File(fastq2s[jk]);
			        //process 
			        Calendar rightNow = Calendar.getInstance();
			        System.out.println("[" + rightNow.getTime().toString() +"]" + fastq1s.length + ": " + fastq1s[jk] + ", " + fastq2s[jk]);
			        getfastq(fastq1, fastq2, mKmer1, mKmer2, skipLine);
			    }
			}else if(!fastQFile1.equals("")) {
				String[] fastq1s = fastQFile1.split(",");
			    for(int jk = 0; jk < fastq1s.length; jk++) {
			    	File fastq1 = new File(fastq1s[jk]);
			        //process 
			        Calendar rightNow = Calendar.getInstance();
			        System.out.println("[" + rightNow.getTime().toString() +"]" + fastq1s.length + ": " + fastq1s[jk] );
			        getfastq(fastq1, mKmer1, skipLine);
			    }
			}
		}
		float averlen = (float)mLengthSum/mReads;
	    
	    Map<String, Integer> printmKmer1 = new HashMap<String, Integer>();
	    Map<String, Integer> printmKmer2 = new HashMap<String, Integer>();
	    Map<String, Integer> printmKmer = new HashMap<String, Integer>();

	    BufferedWriter kmerfilep1 = new BufferedWriter(new FileWriter(kmerfileprefix+".kmer.1.txt"));
	    BufferedWriter kmerfilep2 = new BufferedWriter(new FileWriter(kmerfileprefix+".kmer.2.txt"));
	    BufferedWriter kmerfilep = new BufferedWriter(new FileWriter(kmerfileprefix+".kmer.txt"));
	    String Seq = "";
	    //only print kmer N > 100
	    for(int i=0; i<mKmerBufLen;i++) {
	    	Seq = getkmerseq(i);
	    	if(seqcomplex(Seq) == 1) {
	    		mKmer1[i] = 0;
	    		mKmer2[i] = 0;
	    	}
	    	if(mKmer1[i]+mKmer2[i]>100) {
	    		printmKmer.put(Seq, mKmer1[i]+mKmer2[i]);
	    	}
	    	if(mKmer1[i]>100){
	    		printmKmer1.put(Seq, mKmer1[i]);
	    	}
	    	if(mKmer2[i]>100){
	    		printmKmer2.put(Seq, mKmer2[i]);
	    	}
	    }
	    Map<String, Integer> printmKmer1_s = sortMap(printmKmer1);
	    Map<String, Integer> printmKmer2_s = sortMap(printmKmer2);
	    Map<String, Integer> printmKmer_s = sortMap(printmKmer);
	    int Nkmer = 0;
	    for(Map.Entry<String, Integer> entry: printmKmer1_s.entrySet()) {
	    	if(entry.getValue()<100) {
	    		break;
	    	}else {
	    		Nkmer++;
	    		kmerfilep1.write(entry.getKey() + "\t" + entry.getValue());
	    		kmerfilep1.newLine();
	    	}
	    	if(Nkmer>100) {
	    		break;
	    	}
	    }
	    Nkmer = 0;
	    for(Map.Entry<String, Integer> entry: printmKmer2_s.entrySet()) {
	    	if(entry.getValue()<100) {
	    		break;
	    	}else {
	    		Nkmer++;
	    		kmerfilep2.write(entry.getKey() + "\t" + entry.getValue());
	    		kmerfilep2.newLine();
	    	}
	    	if(Nkmer>100) {
	    		break;
	    	}
	    }
	    Nkmer = 0;
	    for(Map.Entry<String, Integer> entry: printmKmer_s.entrySet()) {
	    	if(entry.getValue()<100) {
	    		break;
	    	}else {
	    		Nkmer++;
	    		kmerfilep.write(entry.getKey() + "\t" + entry.getValue());
	    		kmerfilep.newLine();
	    	}
	    	if(Nkmer>100) {
	    		break;
	    	}
	    }
	    kmerfilep1.close();
	    kmerfilep2.close();
	    kmerfilep.close();
	    /*
	    String Seq = ""; int key = 0, value = 0; String binaryk="";
	    for (Map.Entry<Integer, Integer> entry : mKmer1.entrySet()) {
	    	key = entry.getKey();
	    	value = entry.getValue();
	    	if(value > NprocessRead*0.2) {
	    		Seq = getkmerseq(key);
	    		System.out.println("fq1\t" + Integer.toBinaryString(key) + "\t" + Seq + "\t" + value);	
	    	}
	    }
	    
	    for (Map.Entry<Integer, Integer> entry : mKmer2.entrySet()) {
	    	key = entry.getKey();
	    	value = entry.getValue();
	    	if(value > NprocessRead*0.2) {
	    		Seq = getkmerseq(key);
	    		System.out.println("fq2\t" + Integer.toBinaryString(key) + "\t" + Seq + "\t" + value);	
	    	}
	    }
	    */
	    BufferedWriter linkerfilep = new BufferedWriter(new FileWriter(linkerfile));
	    
	    ArrayList<KmerGraph> KmerGraphs_all1 = new ArrayList<KmerGraph>();
	    ArrayList<KmerGraph> KmerGraphs_all2 = new ArrayList<KmerGraph>();
	    Map<Integer, Integer> KmerG_idx1 = new HashMap<Integer, Integer>();
	    Map<Integer, Integer> KmerG_idx2 = new HashMap<Integer, Integer>();
	    int cutoff = (int) (NprocessRead*percent);
	    int[] up =new int[4];
	    //int[] down =new int[4];
	    int Nkm1 = 0, Nkm2 = 0;
	    for(int i=0; i<mKmerBufLen;i++) {
	    	if(mKmer1[i]>=NprocessRead*percent) {
	    		//Seq = getkmerseq(i);
	    		up = upnode(i, mKmer1);
			    //down = downnode(i);
			    int ishead = isvalid(up, mKmer1, cutoff);
			    KmerGraphs_all1.add(new KmerGraph(i, mKmer1[i], ishead));
			    KmerG_idx1.put(i, Nkm1);
			    Nkm1++;
	    	}else {
	    		mKmer1[i] = 0;
	    	}
	    	
	    	if(mKmer2[i]>=NprocessRead*percent) {
	    		//Seq = getkmerseq(i);
	    		up = upnode(i, mKmer2);
			    //down = downnode(i);
			    int ishead = isvalid(up, mKmer1, cutoff);
			    //System.out.println("Head " + new String(getkmerseq(i)));
			    // Here is the head
			    KmerGraphs_all2.add(new KmerGraph(i, mKmer2[i], ishead));
			    KmerG_idx2.put(i, Nkm2);
			    Nkm2++;
	    	}else {
	    		mKmer2[i] = 0;
	    	}
	    }
	    
	    /*
	    ArrayList<KmerGraph> KmerGraphs_all1_final = new ArrayList<KmerGraph>();
	    ArrayList<KmerGraph> KmerGraphs_all2_final = new ArrayList<KmerGraph>();
	    if(KmerGraphs_all1.size()>maxKmer) {
	    	KmerGraphs_all1.sort(new Comparator<KmerGraph>() {  
	            @Override  
	            public int compare(KmerGraph arg0, KmerGraph arg1) {  
	                return (int) (arg1.value - arg0.value);  // down sort
	            }  
	        });
	    	for(int i=0;i<100;i++) {
	    		KmerGraphs_all1_final.add(KmerGraphs_all1.get(i));
	    	}
	    }else {
	    	KmerGraphs_all1_final = KmerGraphs_all1;
	    }
	    
	    if(KmerGraphs_all2.size()>maxKmer) {
	    	KmerGraphs_all2.sort(new Comparator<KmerGraph>() {  
	            @Override  
	            public int compare(KmerGraph arg0, KmerGraph arg1) {  
	                return (int) (arg1.value - arg0.value);  // down sort
	            }  
	        });
	    	for(int i=0;i<100;i++) {
	    		KmerGraphs_all2_final.add(KmerGraphs_all2.get(i));
	    	}
	    }else {
	    	KmerGraphs_all2_final = KmerGraphs_all2;
	    }
	    
	    for(KmerGraph S:KmerGraphs_all1_final) {
	    	outnode(S, KmerGraphs_all1_final, mKmer1, cutoff, KmerG_idx1);
	    }
	    for(KmerGraph S:KmerGraphs_all2_final) {
	    	outnode(S, KmerGraphs_all2_final, mKmer2, cutoff, KmerG_idx2);
	    }
	    */
	    for(KmerGraph S:KmerGraphs_all1) {
	    	outnode(S, KmerGraphs_all1, mKmer1, cutoff, KmerG_idx1);
	    }
	    for(KmerGraph S:KmerGraphs_all2) {
	    	outnode(S, KmerGraphs_all2, mKmer2, cutoff, KmerG_idx2);
	    }

	    Calendar rightNow = Calendar.getInstance();
	    System.out.println("[" + rightNow.getTime().toString() +"] Get linker by kmers ... " + KmerGraphs_all1.size() + " " + KmerGraphs_all2.size());
	    ArrayList<KmerGraph> KmerGraphs_extend1 = new ArrayList<KmerGraph>();
		for(KmerGraph S:KmerGraphs_all1) {
			if(S.ishead == -1) {
				//System.out.println("HHHx " + getkmerseq(S.node()) );
				KmerGraphs_extend1.add(S);
				if(Klen>5) {
					//extend 1 time
					for(int i=0;i<S.next.size();i++) {
					    KmerGraphs_extend1.add(S.next.get(i));
					}
				}
			}
		}
		ArrayList<KmerGraph> KmerGraphs_extend2 = new ArrayList<KmerGraph>();
		for(KmerGraph S:KmerGraphs_all2) {
			if(S.ishead == -1) {
				KmerGraphs_extend2.add(S);
				if(Klen>5) {
					//extend 1 time
					for(int i=0;i<S.next.size();i++) {
					    KmerGraphs_extend2.add(S.next.get(i));
					}
				}
			}
		}
		ArrayList<DNAString> mylinkers1;
		if(KmerGraphs_extend1.size()>1) {
			rightNow = Calendar.getInstance();
		    System.out.println("[" + rightNow.getTime().toString() +"] Head node kmer: " + KmerGraphs_extend1.size());
			mylinkers1 = mergeKmer(KmerGraphs_extend1, (int) (NprocessRead*percent));
		}else {
			rightNow = Calendar.getInstance();
		    System.out.println("[" + rightNow.getTime().toString() +"] Head node kmer all: " + KmerGraphs_all1.size());
			mylinkers1 = mergeKmer(KmerGraphs_all1, (int) (NprocessRead*percent));
		}
		ArrayList<DNAString> mylinkers2;
		if(KmerGraphs_extend2.size()>1) {
			rightNow = Calendar.getInstance();
		    System.out.println("[" + rightNow.getTime().toString() +"] Head node kmer: " + KmerGraphs_extend2.size());
			mylinkers2 = mergeKmer(KmerGraphs_extend2, (int) (NprocessRead*percent));
		}else {
			rightNow = Calendar.getInstance();
		    System.out.println("[" + rightNow.getTime().toString() +"] Head node kmer all: " + KmerGraphs_all2.size());
			mylinkers2 = mergeKmer(KmerGraphs_all2, (int) (NprocessRead*percent));
		}

		linkerfilep.write("Kmer after filter: " + KmerGraphs_extend1.size() + " " + KmerGraphs_extend2.size());
	    linkerfilep.newLine();
	    
	    mylinkers1.addAll(mylinkers2);
	    //need remove contains short pa

	    // remove short reads
	    for (int i = 0; i < mylinkers1.size(); i++) {
            if(mylinkers1.get(i).length()<minLinkerLen) {
            	mylinkers1.remove(i);
                i--;
            }
        }
	    //sort with score
	    mylinkers1 = ADSsort(mylinkers1);
	    for(DNAString DS:mylinkers1) {
			System.out.println("Linker detected in fq1: " + DS + " " + DS.Svalue);
		}
	    System.out.println("==------------------------------------------------------==");
	    for(DNAString DS:mylinkers2) {
			System.out.println("Linker detected in fq2: " + DS + " " + DS.Svalue);
		}
	    /*
	    mylinkers1 = removeDuplicate(mylinkers1);
	    linkerfilep.write("Linker detected all: " + mylinkers1.size());
	    linkerfilep.newLine();
	    mylinkers1 = removeSimaliar(mylinkers1);
	    */

	    removeContains(mylinkers1);
	    linkerfilep.write("Linker detected after remove contaion: " + mylinkers1.size());
	    removeHeadSame(mylinkers1);
	    linkerfilep.write("\nLinker detected after remove similar: " + mylinkers1.size());
	    linkerfilep.newLine();
  	    //sw linker
        ArrayList<Linker> alllinkers = storeLiners();
	    for(int i=0;i<mylinkers1.size();i++) {
	    	DNAString DS = mylinkers1.get(i);
	    	int[] alignsw = swlinkers(alllinkers, DS.getseq());
	    	if(alignsw[0]>-1) {
	    	    System.out.println("Final linker: " + i + " " + DS + " " + DS.Svalue + " " + alignsw[0] + "," + alllinkers.get(alignsw[0]).getlinkerseq() + " " + alignsw[1]+","+alignsw[2]);
	    	}else {
	    		System.out.println("Final linker: " + i + " " + DS + " " + DS.Svalue + " -1 -1");
	    	}
	    	if(alignsw[0]>-1) {
	    		linkerfilep.write("Linker: " + i + " " + DS + " " + DS.Svalue + " " + alignsw[0] + "," + alllinkers.get(alignsw[0]).getlinkerseq() + " " + alignsw[1]+","+alignsw[2]);
	    	}else {
	    		linkerfilep.write("Linker: " + i + " " + DS + " " + DS.Svalue + " -1 -1");
	    	}
	    	if(alignsw[1]>=linkerScore) {
	    		alllinkers.get(alignsw[0]).setalignl(alignsw[1]);
	    		//alllinkers.get(alignsw[0]).findLinker.add("");
	    	}
	    	if(alignsw[2]>=linkerScore) {
	    		alllinkers.get(alignsw[0]).setalignr(alignsw[2]);
	    	}
	    	//linkerfilep.write("Linker: " + DS + " " + DS.Svalue);
	    	linkerfilep.newLine();
		}
	    
	    // get the last linker file~!!
	    ArrayList<Integer> linkersw = new ArrayList<Integer>();
	    ArrayList<Integer> linkersw_single = new ArrayList<Integer>();
	    String newlinkerfile = kmerfileprefix + ".linker.txt";
	    BufferedWriter Newlinkerfile = new BufferedWriter(new FileWriter(newlinkerfile));
	    for(int i=0; i<alllinkers.size(); i++) {
	    	if(alllinkers.get(i).getalignl()>=1 && alllinkers.get(i).getalignr()>=1) {
	    		linkersw.add(i);
	    	}else if(alllinkers.get(i).getalignl()>=1 || alllinkers.get(i).getalignr()>=1) {
	    		linkersw_single.add(i);
	    	}
	    }
	    //System.out.println("Size LL " + mylinkers1.size());

	    if(linkersw.size()==1) {// match in linker lib
	        int linkeridx = linkersw.get(0);
	        Newlinkerfile.write(alllinkers.get(linkeridx).printlinkerseq());//System.out.println("Size LL " + mylinkers1.size());

	        Newlinkerfile.newLine();
        }else if(linkersw_single.size()>=1) { // part linker in lib
        	if(linkersw_single.size()==1) {
	        	int linkeridx = linkersw_single.get(0);
	        	Newlinkerfile.write(alllinkers.get(linkeridx).printlinkerseq_half());
	        	Newlinkerfile.newLine();
        	}else {
        		int maxSidx = 0, maxScore =0;
        		for(int i=0;i<linkersw_single.size();i++) {
        			int midx = linkersw_single.get(i);
        			if(alllinkers.get(midx).alignl > maxScore) {
        				maxScore = alllinkers.get(midx).alignl;
        				maxSidx = midx;
        			}
        			if(alllinkers.get(midx).alignr > maxScore) {
        				maxScore = alllinkers.get(midx).alignr;
        				maxSidx = midx;
        			}
        		}
        		//max score
        		int linkeridx = linkersw_single.get(maxSidx);
        		Newlinkerfile.write(alllinkers.get(linkeridx).printlinkerseq_half());
	        	Newlinkerfile.newLine();
        	}
        }
	    else if(!fastQFile2.equals("")){ //PE
		    DecimalFormat df = new DecimalFormat("0.000");
		    for(int jk = 0; jk < 1; jk++) { //fastq1s.length
		    	File fastq1 = new File(fastQFile1.split(",")[jk]);
		        File fastq2 = new File(fastQFile2.split(",")[jk]);
		        //process 
		        Map<String, Integer> linkerComb = evaluateLinker(fastq1, fastq2, mylinkers1, skipLine);
		        Map<String, Integer> linkerComb_sort = sortMap(linkerComb);
		        
		        //print
		        String mpercent; int Np = 0;
		        for(Map.Entry<String, Integer> eLinker:linkerComb_sort.entrySet()) {
		        	String[] linkerID = eLinker.getKey().split("_");
		        	mpercent = df.format((float)eLinker.getValue()/NprocessRead);
		        	if(Np<=2 || Double.parseDouble(mpercent) > 0.05) {
		        	    linkerfilep.write(eLinker.getKey() + " " + eLinker.getValue() + " " + mpercent + " " + mylinkers1.get(Integer.parseInt(linkerID[0])) + " " + mylinkers1.get(Integer.parseInt(linkerID[1])));
		        	    linkerfilep.newLine();
		        	    Np++;
		        	}
		        }
		        
		        /*if(linkersw.size()>1) {// match in linker lib, but more than 1
			        int linkeridx = linkersw.get(0);	
			        Newlinkerfile.write(alllinkers.get(linkeridx).printlinkerseq());
			        Newlinkerfile.newLine();
		        }*/
		        
		        // process and make decision of last linker
		        int MODE = -1;
		        if(mylinkers1.size()==1) {
		        	selectLinker(1, linkerComb_sort, 0, 1, Newlinkerfile, mylinkers1.get(0), mylinkers1.get(0));
		        }else if(mylinkers1.size()==2){
		        	selectLinker(2, linkerComb_sort, 0, 1, Newlinkerfile, mylinkers1.get(0), mylinkers1.get(1));
		        }else { //more than 2 linkers
		        	//ArrayList<String> finallinkers = new ArrayList<String>();
		        	// process with same linker in library
		        	Np = 0;int prinID=-1;
			        for(Map.Entry<String, Integer> eLinker:linkerComb_sort.entrySet()) {
			        	String[] linkerID = eLinker.getKey().split("_");
			        	float lpercent = (float)eLinker.getValue()/NprocessRead;
			        	if(Np>=2) {
			        		break;
			        	}
			        	if(lpercent>0.01) {
			        		if(linkerID[0].equals(linkerID[1])) {
			        			MODE=0;
			        			if(Np==0) {
			        				Newlinkerfile.write("Linker_mode: 0");// 0, AA
			        				Newlinkerfile.newLine();
			        			}
			        			Newlinkerfile.write(mylinkers1.get(Integer.parseInt(linkerID[0])).getseq());
			        			Newlinkerfile.newLine();
				    	        prinID = Integer.parseInt(linkerID[0]);
				    	        Np++;
			        		}else if(MODE==-1){
			        			MODE = 1;
			        			if(Np==0) {
			        				Newlinkerfile.write("Linker_mode: 1");// 1, AB
			        				Newlinkerfile.newLine();
			        			}
			        			Newlinkerfile.write(mylinkers1.get(Integer.parseInt(linkerID[0])).getseq());
			        			Newlinkerfile.newLine();
			        			Newlinkerfile.write(mylinkers1.get(Integer.parseInt(linkerID[1])).getseq());
			        			Newlinkerfile.newLine();
			        		}
			        	}
			        	
			        }
			        if(Np==1) {
			        	if(prinID==0 && mylinkers1.size()>0) {
			        		Newlinkerfile.write(mylinkers1.get(1).getseq());
			        		Newlinkerfile.newLine();
			        	}else if(prinID>0) {
			        		Newlinkerfile.write(mylinkers1.get(0).getseq());
			        		Newlinkerfile.newLine();
			        	}
			        }
		        	//linkerfilep.write("ERR");
	    	        //linkerfilep.newLine();
		        }
		    }
        }else { //SE
        	if(mylinkers1.size()>0) {
	        	Newlinkerfile.write(mylinkers1.get(0).getseq());
	    	    Newlinkerfile.newLine();
        	}else {
	        	//Newlinkerfile.write("Single-end ERR");
	        	//Newlinkerfile.newLine();
        	}
        }
	    
	    Newlinkerfile.write("Seq_len: " + averlen);
	    Newlinkerfile.newLine();
	    
	    linkerfilep.close();
	    Newlinkerfile.close();
	}

	  //不用这个了
	  public static void getLCString(char[] str1, char[] str2) {
	    int len1, len2;
	    len1 = str1.length;
	    len2 = str2.length;
	    int maxLen = len1 > len2 ? len1 : len2;
	    int[] max = new int[maxLen];// 保存最长子串长度的数组
	    int[] maxIndex = new int[maxLen];// 保存最长子串长度最大索引的数组
	    int[] c = new int[maxLen];
	    int i, j;
	    for (i = 0; i < len2; i++) {
	      for (j = len1 - 1; j >= 0; j--) {
	        if (str2[i] == str1[j]) {
	          if ((i == 0) || (j == 0))
	            c[j] = 1;
	          else
	            c[j] = c[j - 1] + 1;//此时C[j-1]还是上次循环中的值，因为还没被重新赋值
	        } else {
	          c[j] = 0;
	        }
	        // 如果是大于那暂时只有一个是最长的,而且要把后面的清0;
	        if (c[j] > max[0]) {
	          max[0] = c[j];
	          maxIndex[0] = j;
	          for (int k = 1; k < maxLen; k++) {
	            max[k] = 0;
	            maxIndex[k] = 0;
	          }
	        }
	        // 有多个是相同长度的子串
	        else if (c[j] == max[0]) {
	          for (int k = 1; k < maxLen; k++) {
	            if (max[k] == 0) {
	              max[k] = c[j];
	              maxIndex[k] = j;
	              break; // 在后面加一个就要退出循环了
	            }
	          }
	        }
	      }
	      //for (int temp : c) {
	      //  System.out.print(temp);
	      //}
	      //System.out.println();
	    }
	    //打印最长子字符串
	    for (j = 0; j < maxLen; j++) {
	      if (max[j] > 0) {
	        System.out.println("第" + (j + 1) + "个公共子串:");
	        for (i = maxIndex[j] - max[j] + 1; i <= maxIndex[j]; i++)
	          System.out.print(str1[i]);
	        System.out.println(" ");
	      }
	    }
	  }
	
	public static void selectLinker(int linkersize, Map<String, Integer> linkerComb_sort, int linker1, int linker2, BufferedWriter Newlinkerfile,
			DNAString linkerseq1, DNAString linkerseq2) throws IOException {
		String AA = linker1+"_"+linker1;
		String BB = linker2+"_"+linker2;
		String AB = linker1+"_"+linker2;
		String BA = linker2+"_"+linker1;
		if(linkersize==1) {
        	if(linkerComb_sort.containsKey(AA)) {
	        	float lpercent = (float)linkerComb_sort.get(AA)/NprocessRead;
	        	if(lpercent>0.01) {
	        		Newlinkerfile.write(linkerseq1.getseq());
	        		Newlinkerfile.newLine();
	        	}else {
	        		System.out.println("Unexpected: no linker detected, please input --linker !!!");
	        		System.exit(0);
	        	}
        	}else {
        		System.out.println("Unexpected: no linker detected, please input --linker !!!");
        		System.exit(0);
        	}
        }else if(linkersize==2){
        	int[] newlinker = new int[4];
        	if(!linkerComb_sort.containsKey(AA)) {
        		linkerComb_sort.put(AA, 0);
        		newlinker[0]=-1;
        	}else if(!linkerComb_sort.containsKey(BB)) {
        		linkerComb_sort.put(BB, 0);
        		newlinker[1]=-1;
        	}else if(!linkerComb_sort.containsKey(AB)){
        		linkerComb_sort.put(AB, 0);
        		newlinker[2]=-1;
        	}else if(!linkerComb_sort.containsKey(BA)) {
        		linkerComb_sort.put(BA, 0);
        		newlinker[3]=-1;
        	}
        	float AA_BB = (float)(linkerComb_sort.get(AA)+linkerComb_sort.get(BB))/NprocessRead;
        	float AB_BA = (float)(linkerComb_sort.get(AB)+linkerComb_sort.get(BA))/NprocessRead;
        	if(AA_BB > AB_BA) {
        		Newlinkerfile.write("Linker_mode: 0"); //AA
				Newlinkerfile.newLine();
        		if(newlinker[0]!=-1) {
        			Newlinkerfile.write(linkerseq1.getseq());
        			Newlinkerfile.newLine();
        		}
        		if(newlinker[1]!=-1) {
        			Newlinkerfile.write(linkerseq2.getseq());
        			Newlinkerfile.newLine();
        		}
        	}else {
        		Newlinkerfile.write("Linker_mode: 1"); //AB
        		Newlinkerfile.newLine();
        		Newlinkerfile.write(linkerseq1.getseq());
        		Newlinkerfile.newLine();
        		Newlinkerfile.write(linkerseq2.getseq());
        		Newlinkerfile.newLine();
        	}
        }
	}
	
	public static int[] swlinkers(ArrayList<Linker> alllinker, String seq) {
		int lsize = alllinker.size();
		int[] max_score1 = new int[lsize];
		int[] max_score2 = new int[lsize];
		int[] minI1 = new int[lsize];
		int[] minI2 = new int[lsize];
		int MismatchScore = -100;
		int IndelScore = -100;
		for(int i=0;i<lsize;i++) {
			Linker S = alllinker.get(i);
			int maxlen = S.getmaxlen();
			if(seq.length()>maxlen) {
				maxlen = seq.length();
			}
			LocalAlignment localAligner = new LocalAlignment(maxlen,maxlen, MismatchScore, IndelScore);
	    	localAligner.align(S.linker[0], seq, 0);
	    	max_score1[i] = localAligner.getMaxScore();
	    	minI1[i] = localAligner.getMinI();
	    	
	    	localAligner.align(S.linker[1], seq, 0);
	    	max_score2[i] = localAligner.getMaxScore();
	    	minI2[i] = localAligner.getMinI();
	    	
	    	if(minI2[i]==minI1[i]) {
	    		if(max_score1[i]>max_score2[i]) {
	    			max_score2[i] = 0;
	    		}else {
	    			max_score1[i]=0;
	    		}
	    	}
		}
		int max_score = -1, max_index = -1;
		for(int i=0;i<lsize;i++) {
			if(max_score1[i]>max_score) {
				max_score = max_score1[i];
				max_index = i;
			}
			if(max_score2[i]>max_score) {
				max_score = max_score2[i];
				max_index = i;
			}
		}
		int[] alignsw = new int[3];
		alignsw[0] = max_index; //index
		alignsw[1] = max_score1[max_index]; //linker A score
		alignsw[2] = max_score2[max_index]; //B
		return alignsw;
	}
	
	public static ArrayList<Linker> storeLiners() {
		ArrayList<Linker> alllinkers = new ArrayList<Linker>();
		String[] L1s = new String[] {"GTTGGATAAGATATCGCGG", "GTTGGAATGTATATCGCGG"}; //AB
		Linker L1 = new Linker(L1s, 0, L1s[0].length());
		L1.findLinker = new ArrayList<String>();
		alllinkers.add(L1);
		String[] L2s = new String[] {"GTTGGATCCGATATCGCGG", "GTTGGATCATATATCGCGG"}; //AB
		Linker L2 = new Linker(L2s, 0, L2s[0].length());
		L2.findLinker = new ArrayList<String>();
		alllinkers.add(L2);
		String[] L3s = new String[] {"CTGCTGTCCGATATCGCGGCCGC", "CTGCTGTCATATATCGCGGCCGC"}; //AB
		Linker L3 = new Linker(L3s, 0, L3s[0].length());
		L3.findLinker = new ArrayList<String>();
		alllinkers.add(L3);
		String[] L4s = new String[] {"CGCGATATCTTATCTGACT", "GTCAGATAAGATATCGCGT"};
		Linker L4 = new Linker(L4s, 1, L4s[0].length());
		L4.findLinker = new ArrayList<String>();
		alllinkers.add(L4);
		String[] L5s = new String[] {"AGTTGGATACCTGCAGTACTAGTCAGTGGGCCC", "GGGCCCACTGACTAGTACTGCAGGTATCCAACT"};
		Linker L5 = new Linker(L5s, 1, L5s[0].length());
		L5.findLinker = new ArrayList<String>();
		alllinkers.add(L5);
		String[] L6s = new String[] {"CTGCTGATGTATATCGG", "CTGCTGTAAGATATCGG"}; //AB
		Linker L6 = new Linker(L6s, 0, L6s[0].length());
		L6.findLinker = new ArrayList<String>();
		alllinkers.add(L6);
		String[] L7s = new String[] {"ACGCGATATCTTATCTGACT", "AGTCAGATAAGATATCGCGT"};
		Linker L7 = new Linker(L7s, 1, L7s[0].length());
		L7.findLinker = new ArrayList<String>();
		alllinkers.add(L7);
		String[] L8s = new String[] {"AGATCGGAAGAGCGTCGTGTAG", "CTACACGACGCTCTTCCGATCA"};
		Linker L8 = new Linker(L8s, 1, L8s[0].length());
		L8.findLinker = new ArrayList<String>();
		alllinkers.add(L8);
		//String[] L9s = new String[] {"GTTGGATCCGATATCGCGG", "GTTGGATCATATATCGCGG"};
		//Linker L9 = new Linker(L9s);
        //alllinkers.add(L9);
        return alllinkers;
	}
	
	public static void removeHeadSame(ArrayList<DNAString> linkers) {
		int linkersize = linkers.size();
		for(int i=0;i<linkersize;i++) {
			String seq1 = linkers.get(i).getseq();
			int len1 = seq1.length();
			String seq1_h = seq1.substring(1);
			String seq1_t = seq1.substring(0, len1-1);
			String seq1_ht = seq1.substring(1, len1-1);
			for(int j=i+1;j<linkersize;j++) {
				String seq2 = linkers.get(j).getseq();
				int len2 = seq2.length();
				String seq2_h = seq2.substring(1);
				String seq2_t = seq2.substring(0, len2-1);
				String seq2_ht = seq2.substring(1, len2-1);
				if(seq1.length()>seq2.length()) {
					if(seq1_h.contains(seq2_h) || seq1_t.contains(seq2_t) || seq1_ht.contains(seq2_ht) ) {
						linkers.remove(j);
						linkersize--;
						j--;
						continue;
						//seq2rm
					}
				}else if(seq2_h.contains(seq1_h) || seq2_t.contains(seq1_t) || seq2_ht.contains(seq1_ht) ){
					linkers.remove(i);
					linkersize--;
					i--;
					break;
					//seq1rm
				}
			}
		}
		
	}
	public static void removeContains(ArrayList<DNAString> linkers) {
		int linkersize = linkers.size();
		//int[] linkersrm = new int[linkersize];
		for(int i=0;i<linkersize;i++) {
			String seq1 = linkers.get(i).getseq();
			for(int j=i+1;j<linkersize;j++) {
				String seq2 = linkers.get(j).getseq();
				if(seq1.length()>seq2.length()) {
					if(seq1.contains(seq2)) {
						linkers.remove(j);
						linkersize--;
						j--;
						//linkersrm[j] = 1;
						continue;
						//seq2rm
					}
				}else if(seq2.contains(seq1)){
					linkers.remove(i);
					linkersize--;
					i--;
					break;
					//seq1rm
				}
			}
		}
		
	}

	public static int seqcomplex(String kmer) {
		if(kmer.contains("AAAA") || kmer.contains("CCCC") || kmer.contains("TTTT") || kmer.contains("GGGG")) {
			return 1;
		}else {
			return 0;
		}
	}
	 
    /**
     * 首先生成一个不重复的数集（0~9999），然后每次从这个集合中随机的取出一个数字作为生成的随机数（并且从集合中移除它）
     */
    public static int[] RandomInt(int start, int end, int number) {
    	if(number>end-start) {
    		System.out.println("Error: Random Int function, end minus start is smaller than number");
    		System.exit(0);
    	}
        long startTime = System.currentTimeMillis(); //开始测试时间
        Random rd = new Random();
        int[] rds = new int[number];//随机数数组
        List<Integer> lst = new ArrayList<Integer>();//存放有序数字集合
        int index = 0;//随机索引
        for (int i = start; i <end; i++) {
            lst.add(i);
        }
        //int span = end - start - 1;
        for (int i = 0; i < number; i++) {
            index = rd.nextInt(number - i);
            rds[i] = lst.get(index);
            lst.remove(index);
        }
        long endTime = System.currentTimeMillis(); //获取结束时间
        System.out.println("testC运行时间： " + (endTime - startTime) + "ms");
        return rds;
    }
	
	public static ArrayList<DNAString> ADSsort(ArrayList<DNAString> oldlist) {  
        ArrayList<DNAString> list = oldlist;  
        Collections.sort(list, new Comparator<DNAString>() {  
            @Override  
            public int compare(DNAString arg0, DNAString arg1) {  
                return (int) (arg1.Svalue - arg0.Svalue);  // down sort
            }  
        });
        return list;  
    } 
	
	public static Map<String, Integer> sortMap(Map<String, Integer> oldMap) {  
        ArrayList<Map.Entry<String, Integer>> list = new ArrayList<Map.Entry<String, Integer>>(oldMap.entrySet());  
        Collections.sort(list, new Comparator<Map.Entry<String, Integer>>() {  
            @Override  
            public int compare(Map.Entry<java.lang.String, Integer> arg0,  
                    Map.Entry<java.lang.String, Integer> arg1) {  
                return (int) (arg1.getValue() - arg0.getValue());  // down sort
            }  
        });  
        Map<String, Integer> newMap = new LinkedHashMap<String, Integer>();  
        for (int i = 0; i < list.size(); i++) {  
            newMap.put(list.get(i).getKey(), list.get(i).getValue());  
        }  
        return newMap;  
    } 
	
	public static Map<String, Integer> evaluateLinker(File fastq1, File fastq2, ArrayList<DNAString> mylinkers, int skipLine) throws FileNotFoundException, IOException, InterruptedException {
		BufferedReader reader1;
		BufferedReader reader2;
		Calendar rightNow = Calendar.getInstance();
		//get max len of linker
		int maxlinkerlen = maxlinkerlen(mylinkers);
		Map<String, Integer> linkerComb = new HashMap<String,Integer>();
		if(isGZipped(fastq1)) {
			System.out.println("[" + rightNow.getTime().toString() +"] Evalute linker by paired-end reads! ...");
			reader1  = new BufferedReader(
	                new InputStreamReader(new GZIPInputStream(new FileInputStream(fastq1))));
			reader2  = new BufferedReader(
	                new InputStreamReader(new GZIPInputStream(new FileInputStream(fastq2))));
		}else {
			System.out.println("[" + rightNow.getTime().toString() +"] Evalute linker by paired-end reads! ...");
			reader1 = new BufferedReader(new FileReader(fastq1));
			reader2 = new BufferedReader(new FileReader(fastq2));
		}
		//check valid fastq file and store
		String line1 = "";
		String line2 = "";
		ArrayList<String> line = new ArrayList<String>();
		//ArrayList<FastQMulthread> fastqList = new ArrayList<FastQMulthread>();
		int lineNum = 0;
		int mark = 0;
		int TempN = 0;
		while ((line1 = reader1.readLine()) != null && (line2 = reader2.readLine()) != null) {
			lineNum++;
			if(lineNum <= skipLine*4) {
				continue;
			}
			line.add(line1);
			line.add(line2);
			mark = 1;
			if (lineNum % 4 == 0) {
				FastQMulthread fq1 = new FastQMulthread();
				FastQMulthread fq2 = new FastQMulthread();
				for (int i = 0; i < line.size(); i += 2) {
					fq1.setFastq(i / 2, line.get(i));
					fq2.setFastq(i / 2, line.get(i + 1));
				}
				line.clear();
				//System.out.println(" zx " + fq1.getFastq()[1] + " " + fq2.getFastq()[1]);
				//LocalAlign(String seq, int iRead, int maxlinkerlin, ArrayList<DNAString> mylinkers)
				int[] arrayOfInt1 = LocalAlign(fq1.getFastq()[1], maxlinkerlen, mylinkers);
				int[] arrayOfInt2 = LocalAlign(fq2.getFastq()[1], maxlinkerlen, mylinkers);
				
				if(arrayOfInt1[3]>=linkerScore && arrayOfInt2[3]>=linkerScore) {
					//System.out.println(arrayOfInt1[2] + " " + arrayOfInt1[3] + ", " + arrayOfInt2[2] + " " + arrayOfInt2[3]);
					String linkerType =  Integer.toString(arrayOfInt1[2]) + "_" + Integer.toString(arrayOfInt2[2]);
					if(linkerComb.containsKey(linkerType)) {
						int newvalue= linkerComb.get(linkerType)+1;
						linkerComb.put(linkerType, newvalue);
					}else {
						linkerComb.put(linkerType, 1);
					}
					TempN++;
				}
			}
			if(lineNum - skipLine*4 > NprocessRead*4) {
				System.out.println("NNx " + TempN + " " + (int)(lineNum/4));
				break;
			}
			if (lineNum % 2000 == 0) {
				//queue.put(fastqList);
				//fastqList = new ArrayList<FastQMulthread>();
				mark = 2;
			}
		}
		if (mark == 1) {// put <2000 into queue
			;//queue.put(fastqList);
		}
		//queue.put(END);

		reader1.close();
		reader2.close();
		return linkerComb;
	}
	
	public static int maxlinkerlen( ArrayList<DNAString> mylinkers) {
		int maxlen = 0;
		for(DNAString DS:mylinkers) {
			if(DS.length()>maxlen) {
				maxlen = DS.length();
			}
		}
		return maxlen;
	}
	
	public static int[] LocalAlign(String seq, int maxlinkerlin, ArrayList<DNAString> mylinkers) throws IOException {
		
		LocalAlignment localAligner = new LocalAlignment(maxlinkerlin, maxlinkerlin);
		
	    int bestScore = -1;
	    int secondBestScore = -1;
	    int bestLinkerIndex = -1;
	    int minJ = -1;
	    int minI = -1;
	    int maxJ = -1;
	    int maxI = -1;
	    int secondminJ = -1;
	    int secondminI = -1;
	    int secondmaxJ = -1;
	    int secondmaxI = -1;
	    for (int i = 0; i < mylinkers.size(); i++) {
	    	//System.out.println(mylinkers.get(i).getseq() + " zz " + seq );
	    	localAligner.align(mylinkers.get(i).getseq(), seq, 0); // 17
	    	int score = localAligner.getMaxScore();
	    	if (bestScore < score || (bestScore == score && maxJ > localAligner.getMaxJ() ) ) {
	    		secondBestScore = bestScore;
	    		secondminJ = minJ;
	    		secondminI = minI;
	    		secondmaxJ = maxJ;
	    		secondmaxI = maxI;
	    		bestScore = score;
		        bestLinkerIndex = i;
		        minJ = localAligner.getMinJ();
		        minI = localAligner.getMinI();
		        maxJ = localAligner.getMaxJ();
		        maxI = localAligner.getMaxI();
	    	}else if(bestScore == score && maxJ < localAligner.getMaxJ()) { 
	    		secondBestScore = bestScore;
	    		secondminJ = localAligner.getMinJ();
	    		secondminI = localAligner.getMinI();
	    		secondmaxJ = localAligner.getMaxJ();
	    		secondmaxI = localAligner.getMaxI();
	    		//bestScore = score;
		        //bestLinkerIndex = i;
	    	}else if (secondBestScore < score) {
	    		secondBestScore = score;
	    	}
	    }

	    int secondBestScoreDiff = bestScore - secondBestScore;
	    int tag_Start = 0;
	    int tag_End = minJ - minI;
	    if (tag_End < 0) {
	    	tag_End = 0;
	    }

	    int[] arrayOfInt = new int[14];
	    arrayOfInt[0] = tag_Start;
	    arrayOfInt[1] = tag_End;
	    arrayOfInt[2] = bestLinkerIndex; // index
	    arrayOfInt[3] = bestScore;
	    arrayOfInt[4] = secondBestScoreDiff;
	    arrayOfInt[5] = minJ;
	    arrayOfInt[6] = minI;
	    arrayOfInt[7] = maxJ;
	    arrayOfInt[8] = maxI;
	    arrayOfInt[9] = secondminJ;
	    arrayOfInt[10] = secondminI;
	    arrayOfInt[11] = secondmaxJ;
	    arrayOfInt[12] = secondmaxI;
	    return arrayOfInt;
    }
	/*
	 *  <-   00 1111 1111 1111 ??
	 *   s 0011 1111 1111 1111
	 *  -> 00?? 1111 1111 1111
	 *  s find in: s >> 2 | max?
	 *  s find out: (KMR & (s << 2) ) | min?
	 * */
	public static ArrayList<DNAString> mergeKmer(ArrayList<KmerGraph> KmerGraphs_extend, int cutoff) {
		int maxA = 0;
		int maxT = 1<<((Klen-1) * 2); // 0011 1111 1100, 11; even 1111 1100; 
		int maxC = 2<<((Klen-1) * 2);
		int maxG = 3<<((Klen-1) * 2);
		int minA=0, minT=1, minC=2,minG=3;
			
		//记住下一级开始的，做一次延伸。
		/*
		for(KmerGraph S:KmerGraphs_extend) {
			System.out.println("size " + new String(getkmerseq(S.node())) + " " + S.next.size() + " " + cutoff);
		}
		*/
		//System.out.println("size=--- " + KmerGraphs_extend.size());
		//需要提前记录节点，否则递归出现java.lang.StackOverflowError不好处理
		//for(KmerGraph S:KmerGraphs_extend) {
		//	S.next = outnode(S, mKmer, cutoff);
		//}
		/*
		for(KmerGraph S:KmerGraphs_extend) {
			System.out.println("size " + new String(getkmerseq(S.node())) + " " + S.next.size() + " " + cutoff);
		}
		*/
		//System.out.println("size " + KmerGraphs_extend.size());
		ArrayList<DNAString> mylinkers = FindPath(KmerGraphs_extend);
		return mylinkers;
		
		// for DB graph
		/*
		LinkedHashMap<DNAString, Integer> kmers = new LinkedHashMap<DNAString, Integer>(mKmerGood.size());
		for(int i:mKmerGood) {
		    byte[] Seq = getkmerseqArray(i);
		    DNAString kmerseq = new DNAString(Seq);
		    kmers.put(kmerseq, mKmer[i]);
		}
		DBGraph G = new DBGraph(kmers , Klen, true);
		DNAString[] contigs = G.getSequences(false, Klen);
		String t = new String(getkmerseqArray(4171932));
		System.out.println("Length " + contigs.length + ", " + mKmerGood.size() + ", " + t );
		for(DNAString dnaseq:contigs) {
			System.out.println("DNAS " + dnaseq.toString());
		}
		*/
	}

	public static ArrayList<DNAString> removeSimaliar(ArrayList<DNAString> list) throws IOException{  
		
		DNAString DS1, DS2;
		int[] arrayrm = new int[list.size()];
        for(int i=0;i<list.size();i++){
        	DS1 = list.get(i);
        	for(int j=i+1; j<list.size(); j++) {
        		DS2 = list.get(j);
        		int maxlinkerlen = DS2.length();
        		int minlinkerlen = DS1.length();
        		if(DS1.length()>DS2.length()) {
        			maxlinkerlen = DS1.length();
        			minlinkerlen = DS2.length();
        		}
        	    int[] arrayOfInt = LocalAlign(DS1.getseq(), DS2.getseq(), maxlinkerlen);
        	    //System.out.println(DS1.getseq() + " zxs " + DS2.getseq() + " ss " + arrayOfInt[3] + " " + minlinkerlen + " " +i + " " +j);
        	    if(arrayOfInt[3]>=(minlinkerlen*0.9)) {
        	    	if(DS1.Svalue > DS2.Svalue*1.1) {
        	    		arrayrm[j]=1;
        	    	}else if(DS2.Svalue > DS1.Svalue*1.1){
        	    		arrayrm[i]=1;
        	    	}else if(DS1.length() > DS2.length()) {
        	    		arrayrm[j]=1;
        	    	}else {
        	    		arrayrm[i]=1;
        	    	}
        	    	/*
        	    	if(DS1.length() > DS2.length()-3) {
        	    		arrayrm[j]=1;
        	    	}else if(DS1.length() < DS2.length()-3) {
        	    		arrayrm[i]=1;
        	    	}
        	    	*/
        	    }
        	}
        }
        ArrayList<DNAString> listTemp = new ArrayList<DNAString>();
        for(int i=0;i<list.size();i++) {
        	if(arrayrm[i]==0) {
        		listTemp.add(list.get(i));
        	}
        }
        return listTemp;
    }
	
	public static int[] LocalAlign(String seq, String seq2, int maxlinkerlin) throws IOException {
		
		LocalAlignment localAligner = new LocalAlignment(maxlinkerlin, maxlinkerlin);
		
	    int bestScore = -1;
	    int secondBestScore = -1;
	    int bestLinkerIndex = -1;
	    int minJ = -1;
	    int minI = -1;
	    int maxJ = -1;
	    int maxI = -1;
	    int secondminJ = -1;
	    int secondminI = -1;
	    int secondmaxJ = -1;
	    int secondmaxI = -1;
	    for (int i = 0; i < 1; i++) {
	    	//System.out.println(mylinkers.get(i).getseq() + " zz " + seq );
	    	localAligner.align(seq2, seq, 0); // 17
	    	int score = localAligner.getMaxScore();
	    	if (bestScore < score || (bestScore == score && maxJ > localAligner.getMaxJ() ) ) {
	    		secondBestScore = bestScore;
	    		secondminJ = minJ;
	    		secondminI = minI;
	    		secondmaxJ = maxJ;
	    		secondmaxI = maxI;
	    		bestScore = score;
		        bestLinkerIndex = i;
		        minJ = localAligner.getMinJ();
		        minI = localAligner.getMinI();
		        maxJ = localAligner.getMaxJ();
		        maxI = localAligner.getMaxI();
	    	}else if(bestScore == score && maxJ < localAligner.getMaxJ()) { 
	    		secondBestScore = bestScore;
	    		secondminJ = localAligner.getMinJ();
	    		secondminI = localAligner.getMinI();
	    		secondmaxJ = localAligner.getMaxJ();
	    		secondmaxI = localAligner.getMaxI();
	    		//bestScore = score;
		        //bestLinkerIndex = i;
	    	}else if (secondBestScore < score) {
	    		secondBestScore = score;
	    	}
	    }

	    int secondBestScoreDiff = bestScore - secondBestScore;
	    int tag_Start = 0;
	    int tag_End = minJ - minI;
	    if (tag_End < 0) {
	    	tag_End = 0;
	    }

	    int[] arrayOfInt = new int[14];
	    arrayOfInt[0] = tag_Start;
	    arrayOfInt[1] = tag_End;
	    arrayOfInt[2] = bestLinkerIndex; // index
	    arrayOfInt[3] = bestScore;
	    arrayOfInt[4] = secondBestScoreDiff;
	    arrayOfInt[5] = minJ;
	    arrayOfInt[6] = minI;
	    arrayOfInt[7] = maxJ;
	    arrayOfInt[8] = maxI;
	    arrayOfInt[9] = secondminJ;
	    arrayOfInt[10] = secondminI;
	    arrayOfInt[11] = secondmaxJ;
	    arrayOfInt[12] = secondmaxI;
	    return arrayOfInt;
    }
	
	public static ArrayList<DNAString> removeDuplicate(ArrayList<DNAString> list){  
		ArrayList<DNAString> listTemp = new ArrayList<DNAString>();  
        for(int i=0;i<list.size();i++){  
        	if(list.get(i).length()>minLinkerLen) {
	            if(!listTemp.contains(list.get(i))){
	                listTemp.add(list.get(i));
	            }
        	}
        }  
        return listTemp;  
    }
	
	public static void outnode(KmerGraph node, ArrayList<KmerGraph> KmerGraph_all, int[] mKmer, int cutoff, Map<Integer, Integer> KmerG_idx){
		int[] down = downnode(node.node()); //node.outedge();
		//int[] up = node.inedge();
		int index = 0;
		if(down[0] != node.node() && mKmer[down[0]] >= cutoff) {
			index = KmerG_idx.get(down[0]);
			if(index>=0 && index < KmerGraph_all.size()) {
		        node.next.add(KmerGraph_all.get(index));
			}
		}
		if(down[1] != node.node() && mKmer[down[1]] >= cutoff) {
			index = KmerG_idx.get(down[1]);
			if(index>=0 && index < KmerGraph_all.size()) {
			    node.next.add(KmerGraph_all.get(index));
			}
		}
		if(down[2] != node.node() && mKmer[down[2]] >= cutoff) {
			index = KmerG_idx.get(down[2]);
			if(index>=0 && index < KmerGraph_all.size()) {
			    node.next.add(KmerGraph_all.get(index));
			}
		}
		if(down[3] != node.node() && mKmer[down[3]] >= cutoff) {
			index = KmerG_idx.get(down[3]);
			if(index>=0 && index < KmerGraph_all.size()) {
			    node.next.add(KmerGraph_all.get(index));
			}
		}
	}
	
	private static ArrayList<DNAString> FindPath(ArrayList<KmerGraph> input) {
        ArrayList<DNAString> result = new ArrayList<>();
        if (input == null || input.size() == 0) {
        	//System.out.println("EEEEEE");
            result.add(new DNAString(""));
        } else {
            for (int i = 0; i < input.size(); i++) {
            	//System.out.println("111x " + new String(getkmerseqArray(input.get(i).node())) + " " + input.size() + " " + input.get(i).node() + " " + input.get(i).value );
                if (input.get(i).Visited) {
                	//System.out.println("QWWWWW");
                    result.add(new DNAString(""));
                    continue;
                }
                input.get(i).Visited = true;
                //ArrayList<KmerGraph> KmerGraph_next = outnode(input.get(i), mKmer, cutoff);
                //System.out.println("KKK " + new String(getkmerseqArray(input.get(i).node())) + " " + input.get(i).next.size());
                ArrayList<DNAString> next_seq = FindPath(input.get(i).next); //next节点的后面seq..., iterator
                
                for (int j = 0; j < next_seq.size(); j++) {
                	DNAString s = next_seq.get(j);
                    if (s.length() == 0) {
                    	//result.add(new DNASequence(input.get(i).Seq.getSeq(), '+', input.get(i).Seq.Value));
                    	result.add(new DNAString(input.get(i).node(), Klen, input.get(i).value));
                        //System.out.println("HHXX " + new String(getkmerseqArray(input.get(i).node())) + " " + input.size() + " " + input.get(i).node() + " " + input.get(i).value );
                        //result.add(new DNAString("HHXX "));
                    } else {
                    	//result.add(new DNAString("nnnn "));
                    	//System.out.println("nnxn " + new String(getkmerseqArray(input.get(i).node())) + " " + s);
                    	result.add(new DNAString(input.get(i).node(), s.subSequence(Klen-1, s.length()).toByteArray(), Klen, Math.min(input.get(i).value, s.Svalue)));
                        //result.add(new DNAString(input.get(i).toString() + s.getSeq().substring(input.get(i).Seq.getSeq().length() - 1), '+', Math.min(input.get(i).Seq.Value, s.Value)));
                    }
                }
                input.get(i).Visited = false;
            }
        }
        return result;
    }
	
	public static void findpath(LinkedHashMap<Integer, KmerGraph> KmerGraphs, int i, KmerGraph km, String seq, int cutoff) {
		/*int hasup = isvalid();
		if(km.hasup()) {
			// make edge = 0
			int up[] = km.inedge();
			if(up[0]>cutoff) {
				KmerGraph km_tmp = KmerGraphs.get(km)
			}
		}*/
	}
	
	public static int isvalid(int[] up, int[] mKmer, int cutoff) {
		if(mKmer[up[0]] > cutoff) {
			return up[0];
		}else if(mKmer[up[1]] > cutoff) {
			return up[1];
		}else if(mKmer[up[2]] > cutoff) {
			return up[2];
		}else if(mKmer[up[3]] > cutoff) {
			return up[3];
		}else {
			return -1;
		}
	}
	
	public static int[] upnode(int i) {
		int[] up = new int[4];
		int maxA = 0;
		int maxT = 1<<((Klen-1) * 2); // 0011 1111 1100, 11
		int maxC = 2<<((Klen-1) * 2);
		int maxG = 3<<((Klen-1) * 2);
		up[0] = i >> 2 | maxA;
	    up[1] = i >> 2 | maxT;
 		up[2] = i >> 2 | maxC;
	    up[3] = i >> 2 | maxG;
		return up;
	}
	
	public static int[] upnode(int i, int[] mKmer) {
		int[] up = new int[4];
		int maxA = 0;
		int maxT = 1<<((Klen-1) * 2); // 0011 1111 1100, 11
		int maxC = 2<<((Klen-1) * 2);
		int maxG = 3<<((Klen-1) * 2);
		up[0] = i >> 2 | maxA;
	    up[1] = i >> 2 | maxT;
 		up[2] = i >> 2 | maxC;
	    up[3] = i >> 2 | maxG;
		/*
		System.out.println("up-" + Integer.toBinaryString(i) + " " + mKmer[i]);
		System.out.println("up" + Integer.toBinaryString(up[0]) + " " + mKmer[up[0]]);
		System.out.println("up" + Integer.toBinaryString(up[1]) + " " + mKmer[up[1]]);
		System.out.println("up" + Integer.toBinaryString(up[2]) + " " + mKmer[up[2]]);
		System.out.println("up" + Integer.toBinaryString(up[3]) + " " + mKmer[up[3]]);
		*/
		return up;
	}
	
	public static int[] downnode(int i) {
		int[] down = new int[4];
		int minA=0, minT=1, minC=2,minG=3;
		down[0] =  (KMR & (i << 2) ) | minA;
		down[1] =  (KMR & (i << 2) ) | minT;
		down[2] =  (KMR & (i << 2) ) | minC;
		down[3] =  (KMR & (i << 2) ) | minG;
		return down;
	}
	
	public static byte[] getkmerseqArray(int kmer) {
		int base = 0;
		byte[] Seq = new byte[Klen];
		for(int i=0;i<Klen;i++) {
			base = kmer & kit;
			kmer = kmer >> 2;
		    switch(base){
                case 0:
                	Seq[Klen-i-1] = 'A';
                	break;
                case 1:
                	Seq[Klen-i-1] = 'T';
                	break;
                case 2:
                	Seq[Klen-i-1] = 'C';
                	break;
                case 3:
                	Seq[Klen-i-1] = 'G';
                	break;
                default:
                	System.out.println("E " + base + "\t" + kmer);
                	Seq[0] = 'X';
		    }
		}
		return Seq;
	}
	
	public static String getkmerseq(int kmer) {
		int base = 0;
		String Seq = "";
		for(int i=0;i<Klen;i++) {
			base = kmer & kit;
			kmer = kmer >> 2;
		    switch(base){
                case 0:
                	Seq = "A" + Seq;
                	break;
                case 1:
                	Seq = "T" + Seq;
                	break;
                case 2:
                	Seq = "C" + Seq;
                	break;
                case 3:
                	Seq = "G" + Seq;
                	break;
                default:
                	System.out.println("E " + base + "\t" + kmer);
                	Seq = "";
		    }
		}
		return Seq;
	}
	
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
	
	public static void getfastq(File fastq1, File fastq2, int[] mKmer1, int[] mKmer2, int skipLine) throws FileNotFoundException, IOException, InterruptedException {
		BufferedReader reader1;
		BufferedReader reader2;
		Calendar rightNow = Calendar.getInstance();
		if(isGZipped(fastq1)) {
			System.out.println("[" + rightNow.getTime().toString() +"] find kmer with gzip fastq file! ...");
			reader1  = new BufferedReader(
	                new InputStreamReader(new GZIPInputStream(new FileInputStream(fastq1))));
			reader2  = new BufferedReader(
	                new InputStreamReader(new GZIPInputStream(new FileInputStream(fastq2))));
		}else {
			System.out.println("[" + rightNow.getTime().toString() +"] find kmer with fastq file! ...");
			reader1 = new BufferedReader(new FileReader(fastq1));
			reader2 = new BufferedReader(new FileReader(fastq2));
		}
		//check valid fastq file and store
		String line1 = "";
		String line2 = "";
		ArrayList<String> line = new ArrayList<String>();
		//ArrayList<FastQMulthread> fastqList = new ArrayList<FastQMulthread>();
		int lineNum = 0;
		int mark = 0;
		while ((line1 = reader1.readLine()) != null && (line2 = reader2.readLine()) != null) {
			
			lineNum++;
			if(lineNum <= skipLine*4) {
				continue;
			}
			line.add(line1);
			line.add(line2);
			mark = 1;
			if (lineNum % 4 == 0) {
				FastQMulthread fq1 = new FastQMulthread();
				FastQMulthread fq2 = new FastQMulthread();
				for (int i = 0; i < line.size(); i += 2) {
					fq1.setFastq(i / 2, line.get(i));
					fq2.setFastq(i / 2, line.get(i + 1));
				}
				line.clear();
				ProcessRead(fq1, mKmer1);
				ProcessRead(fq2, mKmer2);
				//fastqList.add(fq1);
				//fastqList.add(fq2);
			}
			if(lineNum - skipLine*4 > NprocessRead*4) {
				break;
			}
			if (lineNum % 2000 == 0) {
				//queue.put(fastqList);
				//fastqList = new ArrayList<FastQMulthread>();
				mark = 2;
			}
		}
		if (mark == 1) {// put <2000 into queue
			;//queue.put(fastqList);
		}
		//queue.put(END);

		reader1.close();
		reader2.close();
	}
	
	public static void getfastq(File fastq1, int[] mKmer1, int skipLine) throws FileNotFoundException, IOException, InterruptedException {
		BufferedReader reader1;
		Calendar rightNow = Calendar.getInstance();
		if(isGZipped(fastq1)) {
			System.out.println("[" + rightNow.getTime().toString() +"] find kmer with gzip fastq file! ...");
			reader1  = new BufferedReader(
	                new InputStreamReader(new GZIPInputStream(new FileInputStream(fastq1))));
		}else {
			System.out.println("[" + rightNow.getTime().toString() +"] find kmer with fastq file! ...");
			reader1 = new BufferedReader(new FileReader(fastq1));
		}
		//check valid fastq file and store
		String line1 = "";
		String line2 = "";
		ArrayList<String> line = new ArrayList<String>();
		//ArrayList<FastQMulthread> fastqList = new ArrayList<FastQMulthread>();
		int lineNum = 0;
		int mark = 0;
		while ((line1 = reader1.readLine()) != null ) {
			
			lineNum++;
			if(lineNum <= skipLine*4) {
				continue;
			}
			line.add(line1);
			line.add(line2);
			mark = 1;
			if (lineNum % 4 == 0) {
				FastQMulthread fq1 = new FastQMulthread();
				FastQMulthread fq2 = new FastQMulthread();
				for (int i = 0; i < line.size(); i += 2) {
					fq1.setFastq(i / 2, line.get(i));
					fq2.setFastq(i / 2, line.get(i + 1));
				}
				line.clear();
				ProcessRead(fq1, mKmer1);
				//fastqList.add(fq1);
				//fastqList.add(fq2);
			}
			if(lineNum - skipLine*4 > NprocessRead*4) {
				break;
			}
			if (lineNum % 2000 == 0) {
				//queue.put(fastqList);
				//fastqList = new ArrayList<FastQMulthread>();
				mark = 2;
			}
		}
		if (mark == 1) {// put <2000 into queue
			;//queue.put(fastqList);
		}
		//queue.put(END);

		reader1.close();
	}
	
	//public static void ProcessRead(FastQMulthread fastQ1, ConcurrentHashMap<Integer, Integer> mKmer) {
	public static void ProcessRead(FastQMulthread fastQ1, int[] mKmer) {
	    int len = fastQ1.getFastq()[1].length();

	    mLengthSum += len;

	    String seqstr = fastQ1.getFastq()[1];
	    String qualstr = fastQ1.getFastq()[3];

	    int kmer = 0;
	    boolean needFullCompute = true;
	    for(int i=headS; i<len-tailE; i++) {
	        char base = seqstr.charAt(i);
	        //char qual = qualstr.charAt(i);
	        // get last 3 bits
	        
	        if(base == 'N'){
	            needFullCompute = true;
	            continue;
	        }

	        // Klen bases required for kmer computing
	        if(i<Klen-1)
	            continue;

	        // calc 5 KMER
	        // 0x3FC == 0011 1111 1100
	        // calc 11 kmer
	        // 0x3ffffc 0011 1111 1111 1111 1111 1100
	        // calC 15 kmer
	        // 0x3ffffffc 0011 1111 1111 1111 1111 1111 1111 1100
	        if(!needFullCompute){
	            int val = baseval(base);
	            if(val < 0){
	                needFullCompute = true;
	                continue;
	            } else {
	                kmer = ((kmer<<2) & KMR ) | val;
	                //System.out.println(seqstr.substring(i-Klen+1, i+1) + " -- " + Integer.toBinaryString(kmer));
	                /*
	                if(mKmer.containsKey(kmer)) {
	                    mKmer.put(kmer, mKmer.get(kmer)+1);
	                }else {
	                	mKmer.put(kmer, 1);
	                }
	                */
	                mKmer[kmer]++;
	            }
	        } else {
	            boolean valid = true;
	            kmer = 0;
	            for(int k=0; k<Klen; k++) {
	                int val = baseval(seqstr.charAt(i - Klen + 1 + k));
	                if(val < 0) {
	                    valid = false;
	                    break;
	                }
	                kmer = ((kmer<<2) & KMR ) | val;
	            }
	            if(!valid) {
	                needFullCompute = true;
	                continue;
	            } else {
	            	//System.out.println(seqstr.substring(i-Klen+1, i+1) + " " + Integer.toBinaryString(kmer));
	            	/*
	            	if(mKmer.containsKey(kmer)) {
	                    mKmer.put(kmer, mKmer.get(kmer)+1);
	                }else {
	                	mKmer.put(kmer, 1);
	                }
	                */
	            	mKmer[kmer]++;
	                needFullCompute = false;
	            }
	        }

	    }
	    mReads++;
	}
	
	public static int baseval(char base) {
	    switch(base){
	        case 'A':
	            return 0;
	        case 'T':
	            return 1;
	        case 'C':
	            return 2;
	        case 'G':
	            return 3;
	        default:
	            return -1;
	    }
	}
}


class Linker{
	String[] linker;
	//int swscore;
	int alignl;
	int alignr;
	int linkermode;
	int linkerlen;
	ArrayList<String> findLinker;

	public Linker(String[] linker, int linkermode, int linkerlen) {
		this.linker = linker;
		this.linkermode = linkermode; // 0 AB, 1 AA'
		this.linkerlen = linkerlen;
	}
	public void setalignl(int alignl) {
		this.alignl = alignl;
	}
	public void setalignr(int alignr) {
		this.alignr = alignr;
	}
	public int getalignl() {
		return this.alignl;
	}
	public int getalignr() {
		return this.alignr;
	}
	public int getmaxlen() {
		int len = 0;
		for(int i=0;i<linker.length; i++) {
			if(this.linker[i].length()>len) {
				len = this.linker[i].length();
			}
		}
		return len;
	}
	public String getlinkerseq() {
		String seq = "";
		for(int i=0;i<linker.length; i++) {
			if(i==0) {
				seq = this.linker[i];
			}else {
				seq = seq + "," + this.linker[i];
			}
		}
		return seq;
	}
	public String printlinkerseq() {
		String seq = "";
		for(int i=0;i<linker.length; i++) {
			if(i==0) {
				seq = "Linker_mode: " + this.linkermode + "\n" + this.linker[i];
			}else {
				seq = seq + "\n" + this.linker[i];
			}
		}
		return seq;
	}
	public String printlinkerseq_half() {
		String seq = "";
		
		if(this.alignl>=1) {
			seq = "Linker_mode: 0\n" + this.linker[0];
		}else if(this.alignr>=1) {
			seq = "Linker_mode: 0\n" + this.linker[1];
		}
		
		return seq;
	}
}


