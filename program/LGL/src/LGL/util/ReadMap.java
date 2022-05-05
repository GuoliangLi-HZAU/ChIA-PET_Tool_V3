package LGL.util;

public class ReadMap {
	private String name;
    private int flag;
    private String chrom;
    private int pos;
    private String line;
    private int mapQ;
    private int readlen;
    
    public ReadMap(String name, int flag, String chrom, int pos, String line, int mapQ, int readlen){
        this.name = name;
        this.flag = flag;
        this.chrom = chrom;
        this.pos = pos;
        this.line = line;
        this.mapQ = mapQ;
        this.readlen = readlen;
    }
    
    public String getname() {
    	return name;
    }
    public int getpos() {
    	return pos;
    }
    public String getchrom() {
    	return chrom;
    }
    public int getflag() {
    	return flag;
    }
    public String getline() {
    	return line;
    }
    public int getmapq() {
    	return mapQ;
    }
    public int getreadlen() {
    	return readlen;
    }

    public int getspan(Object st){
        if(this == st) return 0;
        if(st ==null) return -2;
        if(this.getClass() != st.getClass()) return -2;
 
        ReadMap another = (ReadMap) st;  // convert ReadMap class
        if(another.getchrom().equalsIgnoreCase(this.chrom)) {
        	return Math.abs(another.getpos() - this.pos);
        }else {
        	return -1;
        }
    }

}