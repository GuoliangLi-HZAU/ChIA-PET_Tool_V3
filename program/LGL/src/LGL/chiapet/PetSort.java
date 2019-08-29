package LGL.chiapet;

import LGL.data.CLUSTER3;
import LGL.data.REGION;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Collections;

/**
 *
 * @author ligl
 */
public class PetSort {
	
    int debugLevel = 4;
    int sortingLabel = 1;//ascending order [default]
    Calendar rightNow = Calendar.getInstance();
    ArrayList<CLUSTER3> cluster3 = new ArrayList<CLUSTER3>();

    public PetSort(String preclusterFile, String clusterFileName, int sortingLabel) throws IOException {
        rightNow = Calendar.getInstance();
        System.out.println("[" + rightNow.getTime().toString() + "] start PetCluster2 ... ");
        this.sortingLabel = sortingLabel;

        PrintWriter clusterFileOut = new PrintWriter(new BufferedWriter(new FileWriter(clusterFileName)));
        BufferedReader fileIn = new BufferedReader(new InputStreamReader(new FileInputStream
        		(new File(preclusterFile))));
        String line = null;
        int nPETs = 0;
        while ((line = fileIn.readLine()) != null) {
            if (line.length() <= 0)//skip the short lines
            {
                continue;
            }
            String[] fields = line.split("\t");
            REGION head = new REGION(fields[0], Integer.parseInt(fields[1]), Integer.parseInt(fields[2]));
            REGION tail = new REGION(fields[3], Integer.parseInt(fields[4]), Integer.parseInt(fields[5]));
            double weight = 1.0;
            if (fields.length >= 7) {
                weight = Double.parseDouble(fields[6]);
            }
            CLUSTER3 cluster = new CLUSTER3(head, tail, weight, this.sortingLabel);
            cluster3.add(cluster);

            /*
             * 1. combine the chromosomes from head and tail as chromPair
             * 2. get the vector with chromPair from the hashtable; if not such a vector, create a vector this
             *  chromPair
             * 3. add the current cluster to the vector
             */

            nPETs++;
            if (nPETs % 1000000 == 0) {
                rightNow = Calendar.getInstance();
                System.out.println("[" + rightNow.getTime().toString() + "] " + (nPETs / 1000000) + 
                		"M PETs read from " + preclusterFile);
            }
        }
        fileIn.close();
        
        /*
         * 1. get all keys from the hashtable
         * 2. for each key, get the vector
         * 3. for each vector, output all the clusters inside the vector
         */

        Collections.sort(cluster3);
        outputClusters(clusterFileOut);
        clusterFileOut.close();
        System.out.println("number of clusters = " + cluster3.size());
    }
    
    public void outputClusters(PrintWriter clusterFileOut) throws IOException {
        for (int i = 0; i < cluster3.size(); i++) {
        	CLUSTER3 anchorCluster = ((ArrayList<CLUSTER3>)cluster3).get(i);
            clusterFileOut.println(anchorCluster.toString());
            if (i % 100000 == 0) {
                clusterFileOut.flush();
            }
        }
        clusterFileOut.flush();
    }

    public static void main(String[] args) throws IOException {
        if (args.length == 2) {
            int sortingLabel = 1;//ascending order [default]
            new PetSort(args[0], args[1], sortingLabel);
        } else if (args.length == 3) {
            int sortingLabel = Integer.parseInt(args[2]);
            new PetSort(args[0], args[1], sortingLabel);
        } else {
            System.out.println("Usage: java PetSort <pre_cluster_file> <sorted_pre_cluster_file> [<sortingLabel>]");
            System.out.println("[<sortingLabel>] is optional. The default is sorting the head and tail in ascending order");
            System.out.println("sortingLabel:  1: ascending order [default]");
            System.out.println("              -1: descending order");
            System.out.println("              otherwise: no sorting");
            System.exit(1);
        }
    }
}
