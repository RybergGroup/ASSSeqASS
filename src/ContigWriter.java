import java.io.*; 
import java.util.ArrayList;


public class ContigWriter {
    public static void writeFASTA (AlignedQ[] contigs, String [] contigNames, PrintStream out) throws Exception {
	writeFASTX(contigs, contigNames, "FASTA", out);
    }
    public static void writeFASTA (AlignedQ[] contigs, PrintStream out) throws Exception {
	String [] contigNames = new String[contigs.length];
	for (int i=0; i < contigs.length; ++i) contigNames[i] = "Contig_" + i;
	writeFASTX(contigs, contigNames, "FASTA", out);
    }
    public static void writeFASTQ (AlignedQ[] contigs, String [] contigNames, PrintStream out) throws Exception {
        writeFASTX(contigs, contigNames, "FASTQ", out);
    }
    public static void writeFASTQ (AlignedQ[] contigs, PrintStream out) throws Exception {
	String [] contigNames = new String[contigs.length];
	for (int i=0; i < contigs.length; ++i) contigNames[i] = "Contig_" + i;
        writeFASTX(contigs, contigNames, "FASTQ", out);
    }
    public static void writeMFASTA (AlignedQ[] contigs, String [] contigNames, PrintStream out) throws Exception {
        writeFASTX(contigs, contigNames, "MFASTA", out);
    }
    public static void writeMFASTA (AlignedQ[] contigs, PrintStream out) throws Exception {
	String [] contigNames = new String[contigs.length];
	for (int i=0; i < contigs.length; ++i) contigNames[i] = "Contig_" + i;
        writeFASTX(contigs, contigNames, "MFASTA", out);
    }
    public static void writeMFASTQ (AlignedQ[] contigs, String [] contigNames, PrintStream out) throws Exception {
        writeFASTX(contigs, contigNames, "MFASTQ", out);
    }
    public static void writeMFASTQ (AlignedQ[] contigs, PrintStream out) throws Exception {
	String [] contigNames = new String[contigs.length];
	for (int i=0; i < contigs.length; ++i) contigNames[i] = "Contig_" + i;
        writeFASTX(contigs, contigNames, "MFASTQ", out);
    }
    private static void writeFASTX ( AlignedQ[] contigs, String[] contigNames, String format, PrintStream out) throws Exception {
	if (contigs != null && contigNames!= null) {
	    int contigNo = 0;
	    for (AlignedQ seq: contigs) {
		// print contig	
    		if (format.charAt(format.length()-1)=='Q') out.print('@');
		else out.print('>');
		if (contigNo < contigNames.length) out.println(contigNames[contigNo]);
		else out.println("null");
		++contigNo;
		//out.println("Contig" + contigNo);
		AlignedQ.Base base;
		int length = seq.length();
		String qualities = "";
		int Npreseeding=0;
		int Nposterior=0;
		if (format.charAt(0) == 'M') {
		    for (int read = 0; read < seq.getNseq(); ++read) {
			if (seq.getStartInSeq(read)-seq.getStartInAlignment(read) > Npreseeding) Npreseeding=seq.getStartInSeq(read)-seq.getStartInAlignment(read);
		    }
		    for (int i=Npreseeding; i>0; --i) { out.print('-'); qualities += (char) (33); }
		}
		for (int i = 0; i < length; ++i) {
		    base = seq.getConBase(i);
		    if (format.charAt(0) == 'M') {
			if (base.nuc != '*') out.print(base.nuc);
			else out.print('-');
                        if (format.charAt(format.length()-1)=='Q') {
                            if (base.qual > 93) qualities += "~";
                            else qualities += (char) (base.qual+33);
                        }
		    }
		    else if (base.nuc != '*') {
			out.print(base.nuc);
			if (format.charAt(format.length()-1)=='Q') {
			    if (base.qual > 93) qualities += "~";
			    else qualities += (char) (base.qual+33);
			}
		    }
		}
		if (format.charAt(0) == 'M') {
		    for (int read = 0; read < seq.getNseq(); ++read) {
			if (seq.getSubSequence(read,seq.getEndInSeq(read)).length()-1+seq.getEndInAlignment(read)-seq.length() > Nposterior) Nposterior=seq.getSubSequence(read,seq.getEndInSeq(read)).length()-1+seq.getEndInAlignment(read)-seq.length();
		    }
		    for (int i=Nposterior; i>0; --i) { out.print('-'); qualities += (char) (33); }
		}
		out.println("");
		if (format.charAt(format.length()-1)=='Q') {
		    out.println("+");
		    out.println(qualities);
		}
		if (format.charAt(0) == 'M') {
		    String[] names = seq.getSequenceNames();
		    for (int read = 0; read < seq.getNseq(); ++read) {
			qualities = "";
			if (format.charAt(format.length()-1)=='Q') out.print('@');
			else out.print('>');
			out.println(names[read]);
			for (int i = Npreseeding+seq.getStartInAlignment(read)-seq.getStartInSeq(read); i>0; --i) out.print('-');
			out.print(seq.getSubSequence(read,0,seq.getStartInSeq(read)).toLowerCase());
			for (int i = seq.getStartInAlignment(read); i <= seq.getEndInAlignment(read); ++i) {
			    char nuc = seq.getBase(read,i);
			    int qual = seq.getBaseQual(read,i);
			    if (nuc == '*') out.print('-');
			    else out.print(nuc);
			    if (format.charAt(format.length()-1)=='Q') {
				if (qual > 93) qualities += "~";
				else qualities += (char) (qual+33);
			    }
			}
			if (seq.getEndInSeq(read)+1 < seq.getSequence(read).length()) out.print(seq.getSubSequence(read,seq.getEndInSeq(read)+1).toLowerCase());
			for (int i = length+Nposterior-(seq.getEndInAlignment(read)+seq.getSubSequence(read,seq.getEndInSeq(read)).length()); i>0; --i) out.print('-');
			out.println("");
			if (format.charAt(format.length()-1)=='Q') {
	 		    out.println("+");
			    out.println(qualities);
			}
		    }
		}
	    }
	}
    }
    public static void writeACE (AlignedQ[] contigs, String[] contigNames, boolean writeReadFileAsABI, PrintStream out) throws Exception {
	if (contigs == null) return;
	if (contigNames == null) {
	    contigNames = new String[contigs.length];
	    for (int i=0; i < contigs.length; ++i) contigNames[i] = "Contig_" + i;
	}
	int n_seq = 0;
	int contigNo = 0;
	for (AlignedQ i : contigs) {
	    n_seq += i.nSeq();
	}
	out.println("AS " + contigs.length + " " + n_seq);
	for (AlignedQ i : contigs) {
	    int [] numb = i.getFirsts();
	    out.println("");
	    if (contigNo < contigNames.length) out.print("CO " + contigNames[contigNo]);
	    else out.print("CO null");
	    ++contigNo;
	    out.print(" " + i.length() + " " + i.nSeq() + " ");
	    ArrayList<Integer> qual = new ArrayList<Integer>();
	    String seq = "";
	    String oneline = "";
	    int prevseq = -1;
	    for (int j=0; j < i.length(); ++j) {
		AlignedQ.Base temp = i.getConBase(j);
		seq += temp.nuc;
		if (temp.nuc != '*') {
		    if (temp.qual > 97) qual.add(97);
		    else qual.add(temp.qual);
		}
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
		out.println(" " + (i.getStartInAlignment(j)+1-i.getStartInSeq(j)));
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
		if (writeReadFileAsABI) {
		    String ABIname = i.getSequence(j).getFileName();
		    int index = ABIname.length()-1;
		    while (Character.isDigit(ABIname.charAt(index))) --index;
		    while (ABIname.charAt(index) == '.') --index;
		    while (ABIname.charAt(index) != '.') --index;
		    out.println("DS FASTQ_FILE: " + ABIname.substring(0,index) + ".ab1");
		}
	       	else out.println("DS FASTQ_FILE: " + i.getSequence(j).getFileName());
		out.println("");
	    }
	}
    }
}
