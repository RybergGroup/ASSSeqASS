import java.io.*;
import java.util.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

public class Test {
    private static void help () {
	System.out.println("ASSSeq assemble assembels sequences given in FastQ or PHD.1 format.");
	System.out.println();
	System.out.println("--batch / -b  will assemble contigs from different batches of sequences, given");
	System.out.println("              by separate -f/--file, separately. Alternatively, a batch file");
	System.out.println("              where the name of the contig is given in the first collumn and");
	System.out.println("              consecative columns name files with reads");
	System.out.println("--files / -f  give batch of sequence files, each as separate argument");
	System.out.println("--format / -o give ouput format as ACE, FASTA, or FASTQ");
	System.out.println("--help / -h   print this help");
	//System.out.println("--lowestalignmentscore / -a lowest score for reporting contig, default 100");
	System.out.println("--names       give names for different batches given by -f/--files as a space");
	System.out.println("              separated list");
	System.out.println("--sep         give character for separation of columns in files, default '\\t'");
	System.out.println("--trimqual / -q trim edges to given qual value, default 15");
	System.out.println();
    }
    public static void main (String[] args) throws Exception {
	//int lowestAlScore = 100;
	//int lowestDip = 75;
	int trim_qual = 15;
	String output_format = "ACE";
	int n_seq = 0;
	ArrayList<ArrayList<String>> inputFiles = new ArrayList<ArrayList<String>>();
	boolean batch = false;
	String[] batchNames = {"contig"};
	String batchFileName = "";
	char separator = '\t';
	////////////////////
	// read arguments //
	////////////////////
	for (int i=0; i < args.length; ++i ) {
	    //System.err.println(args[i]);
	    if (args[i].equals("--batch") || args[i].equals("-b")) {
		if (i+1 < args.length && args[i+1].charAt(0) != '-') {
		    batchFileName = args[++i];
		}
		batch = true;
	    }
	    else if (args[i].equals("--files") || args[i].equals("-f")) {
		ArrayList<String> files = new ArrayList<String>();
		while (i+1 < args.length && args[i+1].charAt(0) != '-') {
		    ++i;
		    files.add(args[i]);
		}
		inputFiles.add(files);
	    }
	    else if (args[i].equals("--format") || args[i].equals("-o")) {
		if (i+1 < args.length && args[++i].charAt(0) != '-') {
		    if (args[i].toUpperCase().equals("ACE")) {
			output_format = "ACE";
		    }
		    if (args[i].toUpperCase().equals("ACEABI")) {
                        output_format = "ACEABI";
                    }
		    else if (args[i].toUpperCase().equals("FASTA")) {
			output_format = "FASTA";
		    }
		    else if (args[i].toUpperCase().equals("FASTQ")) {
                        output_format = "FASTQ";
                    }
		    else if (args[i].toUpperCase().equals("MFASTA")) {
                        output_format = "MFASTA";
                    }
		    else if (args[i].toUpperCase().equals("MFASTQ")) {
                        output_format = "MFASTQ";
                    }
		    else {
			System.err.println(args[i] + " is not a supported output file format.");
			System.exit(1);
		    }
		}
		else {
		    System.err.println("--format / -o require a supported file format as additional argument. ");
		    System.exit(1);
		}
	    }
	    /*else if (args[i].equals("--lowestalignmentscore") || args[i].equals("-a")) {
                if (i+1 < args.length && args[++i].charAt(0) != '-') {
                    lowestAlScore = Integer.parseInt(args[i]);
                }
                else {
                    System.err.println("--lowestalignmentscore / -a require an integer as next argument.");
                    System.exit(1);
                }
            }*/
	    else if (args[i].equals("--names")) {
		ArrayList<String> names = new ArrayList<String>();
		while (i+1 < args.length && args[i+1].charAt(0) != '-') {
		    ++i;
		    names.add(args[i]);
		}
		if (!names.isEmpty()) batchNames = names.toArray(new String[names.size()]);
	    }
	    else if (args[i].equals("--sep")) {
		if (i+1 < args.length && args[i+1].charAt(0) != '-') {
                    separator = args[++i].charAt(0);
                }
                else {
                    System.err.println("-- sep require aaracter as next argument.");
                    System.exit(1);
                }
	    }
	    else if (args[i].equals("--trimqual") || args[i].equals("-q")) {
		if (i+1 < args.length && args[i+1].charAt(0) != '-') {
		    trim_qual = Integer.parseInt(args[++i]);
		}
		else {
		    System.err.println("--trimqual / -q require an integer as next argument.");
		    System.exit(1);
		}
	    }
	    else if (args[i].equals("--help") || args[i].equals("-h")) {
		help();
		System.exit(0);
	    }
	    else if (i == args.length -1) {
		inputFiles.add(new ArrayList<String>(Arrays.asList(args[i])));
	    }
	    else {
		System.err.println("Unrecognized argument '" + args[i] + "'.");
		System.exit(1);
	    }
	}
	//////////////////////////////////////////////
	// Read sequences and order them in batches //
	//////////////////////////////////////////////
	AlignedQ[] contigs = null;
	SequenceQ[][] sequenceBatches;
	if (!batchFileName.isEmpty()) {
	    HashMap<String,ArrayList<String>> sequenceBatchesFromFile = new HashMap<String,ArrayList<String>>();
	    ArrayList<String> sequenceFiles = new ArrayList<String>();
	    InputStream fr = new FileInputStream(batchFileName);
	    int c;
	    char mode = 'b';
	    boolean escape = false;
	    String word = "";
	    String name = "";
	    while ((c=fr.read()) !=-1) {
		if (escape) {
		    word += (char) c;
		    escape = false;
		}
		else {
		    if ((char) c == '\n' || (char) c == '\r') {
			if (!name.isEmpty()) {
			    if (!word.isEmpty()) sequenceFiles.add(word);
			    if (!sequenceFiles.isEmpty()) {
				if (sequenceBatchesFromFile.containsKey(name)) {
				    sequenceBatchesFromFile.get(name).addAll(sequenceFiles);
				}
				else {
				    sequenceBatchesFromFile.put(name, new ArrayList<String>(sequenceFiles));
				}
			    }
			}
			mode = 'b';
			word = "";
			name = "";
			sequenceFiles.clear();
		    }
		    else if (mode != 'i') {
			if ((char) c == '"' || (char) c == '\'') {
			    char quote = 'c';
			    while ((c=fr.read()) !=-1) {
				if ((char) c == '\\' && !escape) escape = true;
				else if (c == quote && !escape) {
				    quote = 'e';
				    break;
				}
				else {
				    word += c;
				    escape = false;
				}
			    }
			    if (quote != 'e') {
				System.err.println("Error when reading batch file " + batchFileName + ". Reach end of file while in quote.");
				System.exit(1);
			    }
			    escape = false;
			}
			else if ((char) c == separator && mode != 'i') {
			    if (mode == 'b') {
				name = word.trim();
				mode = 'f';
			    }
			    else { sequenceFiles.add(word); }
			    word = "";
			}
			else if ((char) c == '#') mode = 'i';
			else if ((char) c == '\\') escape = true;
			else {
			    word += (char) c;
			}
		    }
		}
	    }
	    String[] new_names = new String[sequenceBatchesFromFile.size()+batchNames.length];
	    ArrayList<ArrayList<String>> tempList = new ArrayList<ArrayList<String>>();
	    int counter = 0;
	    for (String key : sequenceBatchesFromFile.keySet()) {
		new_names[counter] = key;
		tempList.add(sequenceBatchesFromFile.get(key));
		++counter;
	    }
	    tempList.addAll(inputFiles);
	    inputFiles = tempList;
	    System.arraycopy(batchNames,0,new_names,sequenceBatchesFromFile.size(),batchNames.length);
	    batchNames = new_names;
	}
	/*{
    	for (ArrayList<String> files : inputFiles) {
	    System.err.println("#################");
	    for (String temp : files) {
		System.err.print(temp + ",");
	    }
	    System.err.println("");
	    System.err.println("#################");
	}
	System.exit(0);
	}*/
       	if (batch) sequenceBatches = new SequenceQ[inputFiles.size()][];
	else sequenceBatches = new SequenceQ[1][];
	for (int i=0; i < inputFiles.size(); ++i) {
	    SequenceQ[] sequences = new SequenceQ[inputFiles.get(i).size()];
	    int added = 0;
	    for (String fileName : inputFiles.get(i)) {
		InputStream fr = new FileInputStream(fileName);
		char inputFormat = 'U';
		int c;
		while ((c=fr.read()) !=-1) {
		    if ((char) c == '@') { // FastQ format
			inputFormat = 'Q';
			break;
		    }
		    else if ((char) c == 'B' || (char) c == 'b') { //PHD.1 format
			inputFormat = 'p';
			break;
		    }
		}
		fr.close();
		if (inputFormat == 'Q') {
		    fr = new FileInputStream(fileName);
		    int j = 1;
		    while ((c=fr.read()) !=-1) {
			SequenceQ seq;
			String name = fileName;
			if ((char) c == '@') {
			    seq = new SequenceQ(fr, name, 'i');
			    ++j;
			}
			else seq = new SequenceQ(fr, name);
			//System.err.println(added);
			if (seq.length()>0) {
			    if (added >= sequences.length) {
				SequenceQ[] temp = new SequenceQ[sequences.length+1];
				for (int k=0; k < sequences.length; ++k) {
				    temp[k] = sequences[k];
				}
				sequences = temp;
			    }
			    sequences[added] = seq;
			    ++added;
			}
		    }
		}
		else if (inputFormat == 'p') {
		    SequenceQ seq = new SequenceQ(fileName);
		    if (seq.length()>0) sequences[added++] = seq;
		}
		else {
		    System.err.println("Unrecognized file format of file " + fileName + ".");
		    System.exit(1);
		}
	    }
	    if (batch) {
		sequenceBatches[i] = sequences;
	    }
	    else {
		if (sequenceBatches[0] == null) {
		    sequenceBatches[0] = sequences;
		}
		else {
		    SequenceQ[] temp = new SequenceQ[sequenceBatches[0].length+sequences.length];
		    System.arraycopy(sequenceBatches[0],0,temp,0,sequenceBatches[0].length);
		    System.arraycopy(sequences,0,temp,sequenceBatches[0].length,sequences.length);
		    sequenceBatches[0] = temp;
		}
	    }
	}
	//////////////////////////////////
	// Align sequences into contigs //
	//////////////////////////////////
	String[] contigNames = null; 
	for (int b=0; b < sequenceBatches.length; ++b) {
	    AlignedQ[] addContigs = MSAQ.align(sequenceBatches[b],trim_qual); // Do Multiple Sequence Alignment considering quality scores
	    if (contigs == null) {
		contigs = addContigs;
		contigNames = new String[addContigs.length];
		for (int i=0; i < addContigs.length; ++ i) {
		    contigNames[i] = batchNames[b];
		    if (addContigs.length > 1) contigNames[i] += "_" + i;
		    //System.err.println(contigNames[i]);
		}
	    }
	    else {
		AlignedQ[] newContigs = new AlignedQ[contigs.length + addContigs.length];
		String[] newContigNames = new String[contigNames.length + addContigs.length];
		System.arraycopy(contigs,0,newContigs,0,contigs.length);
		System.arraycopy(contigNames,0,newContigNames,0,contigNames.length);
		for (int i=0; i < addContigs.length; ++i) {
		    newContigs[contigs.length+i] = addContigs[i];
		    if (b > batchNames.length) newContigNames[contigNames.length+i] = batchNames[batchNames.length-1] + "_" + (b-batchNames.length);
		    else newContigNames[contigNames.length+i] = batchNames[b];
		    if (addContigs.length > 1) newContigNames[contigNames.length+i] += "_" + i;
		}
		//System.arraycopy(addContigs,0,newContigs,contigs.length,addContigs.length);
		contigs = newContigs;
		contigNames = newContigNames;
	    }
	}
	//for (String name : contigNames) { System.err.println(name); }
	//////////////////
	// Print output //
	//////////////////
     	if (output_format.equals("ACE")) { ContigWriter.writeACE(contigs, contigNames, false, System.out); }
	else if (output_format.equals("ACEABI")) { ContigWriter.writeACE(contigs, contigNames, true, System.out); }
	else if (output_format.equals("FASTA")) { ContigWriter.writeFASTA(contigs, contigNames, System.out); }
	else if (output_format.equals("FASTQ")) { ContigWriter.writeFASTQ(contigs, contigNames, System.out); }
	else if (output_format.equals("MFASTA")) { ContigWriter.writeMFASTA(contigs, contigNames, System.out); }
	else if (output_format.equals("MFASTQ")) { ContigWriter.writeMFASTQ(contigs, contigNames, System.out); }
    }
}
