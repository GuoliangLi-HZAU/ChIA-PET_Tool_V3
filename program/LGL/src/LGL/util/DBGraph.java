package LGL.util;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.*;

/*
 * de Bruijn graph
 * */

public class DBGraph {

    //TODO choose different data structure to avoid problems with duplicate nodes (node.seq = node.twin.seq) created by collapse
    /**
     * Stores all nodes in the graph.
     */
    private LinkedHashMap<DNAString, DBGNode> nodesMap;

    /**
     * k-mer size
     */
    private int k;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * Nodes in the de Bruijn graph.
     * <p>
     * Every node corresponds to a k-mer or a simplified linear stretch.
     * . The duality of nodes edges
     * leads to a bipartite graph.
     */
    private class DBGNode {
        /**
         * outgoing edges (incoming edges are parallel to the outgoing edges of the xxxtwin node)
         */
        private LinkedList<DBGEdge> outgoing;
        /**
         * node sequence
         */
        private DNAString seq;
        /**
         * coverage
         */
        private int count;

        /**
         * Create new DBGNode from a sequence and its read count.
         *
         * @param seq   sequence to store in the new node.
         * @param count number of reads containing any of the k-mers represented in this node.
         */
        public DBGNode(DNAString seq, int count) {
            outgoing = new LinkedList<>();
            this.seq = seq;
            this.count = count;
        }

        /**
         * Get the in-degree of this node (which is equal to the out-degree of its twin).
         *
         * @return the number of incoming edges.
         */
        private int getIndegree() {
            int indegree = 0;
            return indegree;
        }

        /**
         * Checks whether or not this node is a tip.
         *
         * @return true if node is a tip.
         */
        private boolean isTip() {
            if (seq.length() > 2 * k) return false;
            else {
                if (getOutdegree() == 1 && getIndegree() == 0) return true;
                if (getIndegree() == 1 && getOutdegree() == 0) return true;
                return false;
            }
        }

        /**
         * Get the out-degree of this node.
         *
         * @return number of outgoing edges.
         */
        private int getOutdegree() {
            return outgoing.size();
        }

        /**
         * Compute the size of this DBGNode, i.e. the number of k-mers represented by this node.
         *
         * @return the number of k-mers represented by this node.
         */
        private int getSize() {
            return seq.length() - k + 1;
        }

        /**
         * Adds a new edge to the specified node or updates an existing one.
         * <p>
         * As all edges are mirrored by a twin edge, this will also create the reverse edge from n.twin
         * to this nodes's twin. If the edges already exists, this method will either increase or
         * overwrite their multiplicity depending on the value specified in update.
         *
         * @param targetNode   DBGNode to create an edge to.
         * @param multiplicity multiplicity of edge to create.
         * @param update       in case the edge already exists, add to edge multiplicity (true) or overwrite (false)
         */
        private void addEdgeTo(DBGNode targetNode, int multiplicity, boolean update) {
            boolean edgeExists = false;
            for (DBGEdge edge : outgoing) {
                if (edge.target == targetNode) {
                    if (update) {
                        edge.multiplicity += multiplicity;
                    } else {
                        edge.multiplicity = multiplicity;
                    }
                    edgeExists = true;
                    break;
                }
            }
            if (!edgeExists) {
                DBGEdge forward = new DBGEdge();
                forward.target = targetNode;
                forward.multiplicity = multiplicity;
                this.outgoing.add(forward);
            }
        }

        /**
         * Removes an edge and its twin reverse-edge.
         *
         * @param targetNode target node.
         * @return true if the edge existed, false otherwise.
         */
        private boolean rmvEdgeTo(DBGNode targetNode) {
            DBGEdge forward = null;
            for (DBGEdge edge : outgoing) {
                if (edge.target == targetNode) {
                    forward = edge;
                    break;
                }
            }
            if (forward == null) {
                return false;
            } else {
                outgoing.remove(forward);
                return true;
            }
        }

        /**
         * Removes all incoming and outgoing edges from this node and its twin.
         */
        private void isolate() {
            while (outgoing.size() > 0) {
                // remove edge
                DBGEdge forward = outgoing.remove();
            }
        }

        //TODO find better way of computing extended sequences
        private DNAString getShortSeq() {
            return seq.subSequence(k - 1, seq.length());
        }


        /**
         * Generates a String representation of this node.
         *
         * @return String representation of this node.
         */
        @Override
        public String toString() {
            String s = "[Node: " + seq + "]";
            s = s + "\nOut:[";
            for (DBGEdge edge : outgoing) {
                s = s + edge.target.seq.toString() + " ";
            }
            s = s + "]";
            return s;
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * Edges in the de Bruijn graph.
     * <p>
     * DBGEdges are directed edges which connect overlapping sequences. Every edge stores a pointer to its twin
     * edge which connects the target twin node to the source twin. The duality of nodes edges
     * leads to a bipartite graph.
     */
    private class DBGEdge {
        /**
         * Target node of the directed DBGEdge.
         */
        private DBGNode target;
        /**
         * Number of reads using this edge.
         */
        private int multiplicity;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * Create new de Bruijn Graph (DBGraph) object from an array of DNAString reads for a given k-mer length k.
     * <p>
     * The graph is created in 3 steps:
     * (1) unique k-mers are counted in the set of reads and turned into new DBGNodes
     * (2) nodes are connected based on k-1 overlap of their sequences
     * (3) node multiplicities are computed by iterating over the read array
     *
     * @param reads   array of DNAString reads.
     * @param kmers   LinkedHashSet mapping k-mers to counts.
     * @param k       k-mer size for creating this new graph.
     * @param verbose be verbose.
     */
    public DBGraph(Map<DNAString, Integer> kmers, int k, boolean verbose) {

        // track time needed to create graph
        long startTime = System.currentTimeMillis();

        // initialize instance variables
        this.k = k;
        this.nodesMap = new LinkedHashMap<>(kmers.size());

        if (verbose) {
            System.out.println(" adding nodes ... ");
        }

        // add nodesMap
        for (DNAString kmer : kmers.keySet()) {
            // (if kmer is in nodesMap, so is its complement)
            if (!nodesMap.containsKey(kmer)) {
                int count = kmers.get(kmer);

                // create new node
                // no need store reverse complement, canse linker in ChIA-PET is sensitive 
                DBGNode newNode = new DBGNode(kmer, count);
                nodesMap.put(kmer, newNode);
            }
        }

        if (verbose) {
            System.out.println(" adding edges ... ");
        }

        // add edges
        for (DNAString nodeSeq : nodesMap.keySet()) {
            DBGNode node = nodesMap.get(nodeSeq);
            System.out.println("nodeS " + nodeSeq.toString());
            for (DNAString s : nodeSeq.subSequence(1, k).allVariations(k)) {
            	System.out.println("S " + s.toString());
                DBGNode n = nodesMap.get(s);
                // do not add edges to twin node or itself
                if (n != null && n != node) {
                	System.out.println(n.seq.toString());
                    node.addEdgeTo(n, 1, false);
                }
            }
        }

        if (verbose) {
            System.out.println(" computing edge multiplicities ... ");
        }

        if (verbose) {
            long endTime = System.currentTimeMillis();
            System.out.println("done (" + (endTime - startTime) / 1000 +
                    " seconds). " + getSize() + " nodes.\n");
        }
    }

    /**
     * Get number of nodes in the graph.
     *
     * @return number of nodes in the graph.
     */
    public int getSize() {
        return nodesMap.size();
    }

    /**
     * Check if the graph contains a node corresponding to the given sequence.
     *
     * @param s DNAString to look for.
     * @return true if the graph contains a node with the argument sequence.
     */
    public boolean contains(DNAString s) {
        return nodesMap.containsKey(s);
    }

    /**
     * Get all node sequences from the graph.
     *
     * @param includeRC         include the reverse complement of every sequence.
     * @param minSequenceLength minimum sequence length.
     * @return DNAString array containing all node sequences.
     */
    public DNAString[] getSequences(boolean includeRC, int minSequenceLength) {
        LinkedHashSet<DNAString> visited = new LinkedHashSet<>(nodesMap.size());
        LinkedList<DNAString> contigs = new LinkedList<>();
        for (DNAString nodeSeq : nodesMap.keySet()) {
            DBGNode n = nodesMap.get(nodeSeq);
            if (!visited.contains(n.seq) && n.seq.length() > minSequenceLength) {
                contigs.add(n.seq);
                visited.add(n.seq);
            }
        }
        return contigs.toArray(new DNAString[contigs.size()]);
    }

    /**
     * Get all node sequences from the graph.
     *
     * @return DNAString array containing all node sequences.
     */
    public DNAString[] getSequences() {
        return getSequences(true, 0);
    }

    /**
     * Find the length of the largest sequence.
     *
     * @return maximum sequence length.
     */
    public int getMaxSequenceLength() {
        int max = 0;
        for (DNAString nodeSeq : nodesMap.keySet()) {
            if (nodeSeq.length() > max) {
                max = nodeSeq.length();
            }
        }
        return max;
    }


    /**
     * Generate String representation of the graph (for debugging purposes).     *
     *
     * @return String representation of this DBGraph.
     */
    @Override
    public String toString() {
        String s = "";
        for (DNAString nodeSeq : nodesMap.keySet()) {
            DBGNode n = nodesMap.get(nodeSeq);
            s = s + n.toString() + "\n";
        }
        return s;
    }


    public static void main(String[] args) {
        String inputFile = "/home/tohei/Data/ASAData/reads_complex.fasta";
        String outputFile = "/home/tohei/Data/ASAData/out_new.fasta";

        //DBGraph G = new DBGraph(mKer , k, true);
        //DNAString[] contigs = G.getSequences(false, k);

        /*
        System.out.println();
        int max = G.getMaxSequenceLength();
        System.out.println("Max contig length: " + max);
        System.out.println();

        // save
        //contigs
        System.out.println("done.");
        */
    }
}