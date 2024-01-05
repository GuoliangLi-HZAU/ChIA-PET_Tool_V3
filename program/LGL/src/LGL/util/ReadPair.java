package LGL.util;

public class ReadPair implements Comparable<ReadPair> {
	private String readid;
    private int flag1;
    private String chrom1;
    private int pos1;
    private int readlen1;
    private String line1;
    private int flag2;
    private String chrom2;
    private int pos2;
    private int readlen2;
    private String line2;
    private int paired;
    private int mapQ1;
    private int mapQ2;
    private int supp;
    
    public ReadPair(ReadMap read1, ReadMap read2){
    	int chrspan = read1.getchrom().compareTo(read2.getchrom());
    	this.readid = read1.getname();
    	if(chrspan>0) {
    		this.flag1 = read1.getflag();
            this.chrom1 = read1.getchrom();
            this.pos1 = read1.getpos();
            this.line1 = read1.getline();
            this.readlen1 = read1.getreadlen();
            this.mapQ1 = read1.getmapq();
            this.flag2 = read2.getflag();
            this.chrom2 = read2.getchrom();
            this.pos2 = read2.getpos();
            this.line2 = read2.getline();
            this.mapQ2 = read2.getmapq();
            this.readlen2 = read2.getreadlen();
            this.paired = -1;
    	}else if(chrspan==0){
    		if(read1.getpos()<read2.getpos()) {
    			this.flag1 = read1.getflag();
                this.chrom1 = read1.getchrom();
                this.pos1 = read1.getpos();
                this.line1 = read1.getline();
                this.readlen1 = read1.getreadlen();
                this.mapQ1 = read1.getmapq();
                this.flag2 = read2.getflag();
                this.chrom2 = read2.getchrom();
                this.pos2 = read2.getpos();
                this.line2 = read2.getline();
                this.mapQ2 = read2.getmapq();
                this.readlen2 = read2.getreadlen();
                this.paired = 1;
    		}else {
    			this.flag1 = read2.getflag();
                this.chrom1 = read2.getchrom();
                this.pos1 = read2.getpos();
                this.line1 = read2.getline();
                this.readlen1 = read2.getreadlen();
                this.mapQ1 = read1.getmapq();
                this.flag2 = read1.getflag();
                this.chrom2 = read1.getchrom();
                this.pos2 = read1.getpos();
                this.line2 = read1.getline();
                this.mapQ2 = read2.getmapq();
                this.readlen2 = read1.getreadlen();
                this.paired = 1;
    		}
    	}else {
    		this.flag1 = read2.getflag();
            this.chrom1 = read2.getchrom();
            this.pos1 = read2.getpos();
            this.line1 = read2.getline();
            this.readlen1 = read2.getreadlen();
            this.mapQ1 = read1.getmapq();
            this.flag2 = read1.getflag();
            this.chrom2 = read1.getchrom();
            this.pos2 = read1.getpos();
            this.line2 = read1.getline();
            this.mapQ2 = read2.getmapq();
            this.readlen2 = read1.getreadlen();
            this.paired = -1;
    	}
    }
    
    public int getsup() {
    	this.supp = 0;
    	if((this.flag1 & 0x800)!=0) {
    		this.supp+=1;
    	}
    	if((this.flag2 & 0x800)!=0) {
    		this.supp+=1;
    	}
    	return this.supp;
    }
    
    public ReadPair(int flag1, String chrom1, int pos1, String line1,
    		int flag2, String chrom2, int pos2, String line2, int paired){
    	int chrspan = chrom1.compareTo(chrom2);
    	if(chrspan>0) {
    		this.flag1 = flag1;
            this.chrom1 = chrom1;
            this.pos1 = pos1;
            this.line1 = line1;
            this.flag2 = flag2;
            this.chrom2 = chrom2;
            this.pos2 = pos2;
            this.line2 = line2;
            this.paired = paired;
    	}else if(chrspan==0){
    		if(pos1<pos2) {
		        this.flag1 = flag1;
		        this.chrom1 = chrom1;
		        this.pos1 = pos1;
		        this.line1 = line1;
		        this.flag2 = flag2;
		        this.chrom2 = chrom2;
		        this.pos2 = pos2;
		        this.line2 = line2;
		        this.paired = paired;
    		}else {
    			this.flag1 = flag2;
		        this.chrom1 = chrom2;
		        this.pos1 = pos2;
		        this.line1 = line2;
		        this.flag2 = flag1;
		        this.chrom2 = chrom1;
		        this.pos2 = pos1;
		        this.line2 = line1;
		        this.paired = paired;
    		}
    	}else {
    		this.flag1 = flag2;
	        this.chrom1 = chrom2;
	        this.pos1 = pos2;
	        this.line1 = line2;
	        this.flag2 = flag1;
	        this.chrom2 = chrom1;
	        this.pos2 = pos1;
	        this.line2 = line1;
	        this.paired = paired;
    	}
    }
    
    public String getreadid() {
    	return readid;
    }
    public int getpos1() {
    	return pos1;
    }
    public String getchrom1() {
    	return chrom1;
    }
    public int getflag1() {
    	return flag1;
    }
    public String getline1() {
    	return line1;
    }
    public int getpos2() {
    	return pos2;
    }
    public String getchrom2() {
    	return chrom2;
    }
    public int getflag2() {
    	return flag2;
    }
    public String getline2() {
    	return line2;
    }
    public int getpaired() {
    	return paired;
    }
    
    public int getmapq1() {
    	return mapQ1;
    }
    public int getmapq2() {
    	return mapQ2;
    }
    public int getmapq() {
    	return (int)((mapQ1+mapQ2)/2);
    }

    public int getspan(){
        if(this.chrom1.equalsIgnoreCase(this.chrom2)) {
        	return Math.abs(this.pos1 - this.pos2);
        }else {
        	return -1;
        }
    }
    
    public int issamestrand() {
    	String strand1 = "+";
    	String strand2 = "+";
    	if((this.flag1 & 0x10) !=0) {
    		strand1 = "-";
    	}else {
    		strand1 = "+";
    	}
    	if((this.flag2 & 0x10) !=0 ) {
    		strand2 = "-";
    	}else {
    		strand2 = "+";
    	}
    	if(strand1.equalsIgnoreCase(strand2)) {
    		return 1;
    	}else {
    		return 0;
    	}
    }
    
    public String getbedpe() {
    	String strand1 = "+";
    	String strand2 = "+";
    	if((this.flag1 & 0x10) !=0) {
    		strand1 = "-";
    	}else {
    		strand1 = "+";
    	}
    	if((this.flag2 & 0x10) !=0 ) {
    		strand2 = "-";
    	}else {
    		strand2 = "+";
    	}

    	return this.chrom1 + "\t" + this.pos1 + "\t" + (this.pos1 + this.readlen1) + 
    			"\t" + this.chrom2 + "\t" + this.pos2 + "\t" + (this.pos2 + this.readlen2) + 
    			"\t" + this.readid + "\t" + this.getmapq() + "\t" + strand1 + "\t" + strand2;
    }
    
    public int getdiffspan(Object st){
        if(this == st) return 0;
        if(st ==null) return -2;
        if(this.getClass() != st.getClass()) return -2;
 
        ReadPair another = (ReadPair) st;  // convert ReadMap class
        if(another.getchrom1().equalsIgnoreCase(this.chrom1) && 
        		another.getchrom2().equalsIgnoreCase(this.chrom2) ) {
        	return Math.abs(another.getpos1() - this.pos1) + Math.abs(another.getpos2() - this.pos2);
        }else {
        	return -1;
        }
    }

    @Override
    public int compareTo(ReadPair compareRP) {
    	//System.out.println(this.readid + " " + this.getpaired() + " " + compareRP.getpaired());
        if (this == null || compareRP == null) {
            return 0;
        }

        if (this.getpaired() != compareRP.getpaired()) { // single vs paired
            if (this.getpaired() == 1) {
                return -1;
            } else if (compareRP.getpaired() == 1) {
                return 1;
            } else {
                return Integer.compare(compareRP.getmapq(), this.getmapq());
            }
        } else if (this.getpaired() == 1) { // this paired and compareRP paired
        	int span1 = Math.abs(compareRP.pos1 - compareRP.pos2);
            int span2 = Math.abs(this.pos1 - this.pos2);
            return Integer.compare(span1, span2);
            //现在这里逻辑并不好，但是新版本报错Exception in thread "main" java.lang.IllegalArgumentException: Comparison method violates its general contract
            //先用着
            /*
        	if(compareRP.getchrom1().equalsIgnoreCase(this.chrom1) && 
        			compareRP.getchrom2().equalsIgnoreCase(this.chrom2) && this.chrom1.equalsIgnoreCase(this.chrom2)) {
                int span1 = Math.abs(compareRP.pos1 - compareRP.pos2);
                int span2 = Math.abs(this.pos1 - this.pos2);
                System.out.println(this.chrom1 + " " + this.chrom2 + " " + span1 + " " + span2);
                if(span1 == span2) {
                	return Integer.compare(compareRP.getmapq(), this.getmapq());
                }else {
        		    return Integer.compare(span1, span2);
                }
            }else {
            	System.out.println(compareRP.chrom1 + " " + compareRP.chrom2 + " " + this.chrom1 + " " + this.chrom2 + 
            			" " + compareRP.getmapq() + " " + this.getmapq());
                return Integer.compare(compareRP.getmapq(), this.getmapq());
            }
            */
        } else {
            return Integer.compare(compareRP.getmapq(), this.getmapq());
        }
    }
	
	public int compareTo_test(ReadPair compareRP) {
		// TODO Auto-generated method stub
		if(this != null && compareRP != null) {
			if(compareRP.getpaired() != this.getpaired()) {
				if(this.getpaired() == 1) {
					return -1;
				}else if(compareRP.getpaired() == 1) {
					return 1;
				}else {
					return compareRP.getmapq() - this.getmapq();
				}
			}else { // l and r is same : single or paired
				if(this.getpaired() == 1)
			    {
				    int asbspan = this.getdiffspan(compareRP);
				    if(asbspan<100) {
					    return compareRP.getmapq() - this.getmapq();
				    }else {
					    return compareRP.getspan() - this.getspan();
				    }
			    } else {
				    return compareRP.getmapq() - this.getmapq();
			    }
			}
		}else {
			return 0;
		}

	}

}