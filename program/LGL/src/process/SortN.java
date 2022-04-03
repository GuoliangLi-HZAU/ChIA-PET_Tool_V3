package process;

@SuppressWarnings("rawtypes")
public class SortN implements Comparable {
	
	private String content;// one row
    private int column = 0;// the column to sort
    
    public SortN() {
    	
    }
    
    public void setContent(String content) {
    	this.content = content;
    }
    
    public void setColumn(int value) {
    	column = value;
    }
    
    public String getContent() {
    	return content;
    }
    
    public int getColumn() {
    	return column;
    }
    
	public int compareTo(Object o) {
		return column - ((SortN)o).getColumn();
	}
}
