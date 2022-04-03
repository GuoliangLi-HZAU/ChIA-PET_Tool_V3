package LGL.util;

import java.util.ArrayList;

public class KmerGraph{
	
	//private int[] inedge; // ATCG
	//private int[] outedge;
	private int node;
	public int value;
    public boolean Visited = false;
    public ArrayList<KmerGraph> next = new ArrayList<>();
    public int ishead;
    
	public KmerGraph(int node, int value) { //int[] inedge, int[] outedge, 
		//this.inedge = inedge;
		//this.outedge = outedge;
		this.node = node;
		this.value = value;
		this.ishead = 0;
	}
	
	public KmerGraph(int node, int value, int head) { //int[] inedge, int[] outedge,
		//this.inedge = inedge;
		//this.outedge = outedge;
		this.node = node;
		this.value = value;
		this.ishead = head;
	}
	
	public void setnext(ArrayList<KmerGraph> next) {
		this.next = next;
	}
	
	/*
	public int[] inedge() {
		return inedge;
	}
	
	public int[] outedge() {
		return outedge;
	}
	*/
	
	public int node() {
		return node;
	}
	
	public int getcase() {
		return node & 3;
	}
	
}