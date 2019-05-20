import java.io.*; 
import java.util.ArrayList;


public class ContigWriter {
    public static void writeFASTA (AlignedQ[] contigs, PrintStream out) throws Exception {
	writeFASTX(contigs, "FASTA", out);
    }
    public static void writeFASTQ (AlignedQ[] contigs, PrintStream out) throws Exception {
        writeFASTX(contigs, "FASTQ", out);
    }
    private static void writeFASTX ( AlignedQ[] contigs, String format, PrintStream out) throws Exception {
	if (contigs != null) {
	    for (AlignedQ seq: contigs) {
    		if (format.equals("FASTQ")) out.print('@');
		else out.print('>');
		out.println(seq.getName());
		AlignedQ.Base base;
		int length = seq.length();
		String qualities = "";
		for (int i = 0; i < length; ++i) {
		    base = seq.getConBase(i);
		    if (base.nuc != '*') {
			out.print(base.nuc);
			if (format.equals("FASTQ")) {
			    if (base.qual > 93) qualities += "~";
			    else qualities += (char) (base.qual+33);
			}
		    }
		}
		out.println("");
		if (format.equals("FASTQ")) {
		    out.println("+");
		    out.println(qualities);
		}
	    }
	}
    }
    public static void writeACE (AlignedQ[] contigs, PrintStream out) throws Exception {
	if (contigs == null) return;
	int n_seq = 0;
	for (AlignedQ i : contigs) {
	    n_seq += i.nSeq();
	}
	out.println("AS " + contigs.length + " " + n_seq);
	for (AlignedQ i : contigs) {
	    int [] numb = i.getFirsts();
	    /*for (int k=0; k < numb.length; ++k) {
		System.err.print(numb[k] + " ");
	    }
	    System.err.println();*/
	    out.println("");
	    out.print("CO " + i.getName() + " " + i.length() + " " + i.nSeq() + " ");
	    ArrayList<Integer> qual = new ArrayList<Integer>();
	    String seq = "";
	    //ArrayList<String> BS = new ArrayList<String>();
	    String oneline = "";
	    int prevseq = -1;
	    for (int j=0; j < i.length(); ++j) {
		AlignedQ.Base temp = i.getConBase(j);
		seq += temp.nuc;
		if (temp.nuc != '*') qual.add(temp.qual);
		
	    }
	    if (oneline.length() > 0) oneline += seq.length() + " " + i.getSequence(prevseq).getFileName();
	    out.println(0 /*BS.size()*/ + " U");
	    for (int j=0; j < seq.length(); ++j) {
		out.print(seq.charAt(j));
		if ((j+1)%50 == 0) out.println("");
	    }
	    out.println("");
	    out.println("");
	    out.println("BQ");
	    for (int j=0; j < qual.size(); ++j) {
		out.print(" " + qual.get(j));
		if ((j+1)%50 == 0) out.println("");
	    }
	    out.println(""); // start printing AF
	    out.println("");
	    for (int j=0; j < i.nSeq(); ++j) {
		out.print("AF " + i.getSequence(j).getID() + " ");
		if (i.isRevComp(j)) out.print("C");
		else out.print("U");
		out.println(" " + (i.getStartInAlignment(j)+1));
	    }
	    out.println("");
	    for (int j=0; j < i.nSeq(); ++j) {
		out.println("RD " + i.getSequence(j).getID() + " " + i.getPaddedLength(j) + " 0 0");
		String paddedSeq = i.getPaddedSequence(j);
		int pos;
		for (pos =0; pos < paddedSeq.length(); ++pos) {
		    System.out.print(paddedSeq.charAt(pos));
		    if ((pos+1)%50 == 0) System.out.println("");
		}
		if ((pos+1)%50 > 0) System.out.println("");
		out.println("");
		out.println("QA " + (i.getStartInSeq(j)+1) + " " + (i.getEndInAlignment(j)-i.getStartInAlignment(j)+i.getStartInSeq(j)) + " " + (i.getStartInSeq(j)+1) + " " + (i.getEndInAlignment(j)-i.getStartInAlignment(j)+i.getStartInSeq(j)));
		out.println("DS FASTQ_FILE: " + i.getSequence(j).getFileName());
		out.println("");
	    }
	}
    }
}
