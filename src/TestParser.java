import java.io.*;

public class TestParser {
    public static void main (String[] args) throws Exception {
	ACEparser contigs = new ACEparser(new FileInputStream(args[0]));
	System.out.println(contigs.nContigs());
	for (int i=0; i < contigs.nContigs(); ++i) { 
	    System.out.println(contigs.getContigName(0));
	    System.out.println(contigs.getContigSeq(0));
	    for (int j=0; j < contigs.getNbasesInContig(i); ++j) {
		if (j >0) System.out.print(" ");
		System.out.print(contigs.getQualScore(i, j));
	    }
	    System.out.println();
	    for (int j=0; j < contigs.getNreadsForContig(i); ++j) {
		System.out.print(contigs.getReadNameForContig(i,j) + " ");
		System.out.println(contigs.getReadSeqForContig(i,j));
		System.out.println("Start in contig: " + contigs.getStartInContigForRead(i,j));
		System.out.println("Quality start: " + contigs.getReadQualStartForContig(i,j) + " end: " + contigs.getReadQualEndForContig(i,j));
		System.out.println(contigs.getDescriptionOfReadForContig(i,j));
	    }
	}
    }
}
