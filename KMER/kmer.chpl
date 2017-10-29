module KMER {
  param defaultfilename="none";
  config const infile=defaultfilename;
  config const k=4;

  // use stdin if filename == defaultfilename
  proc input_channel(filename: string, defaultfilename: string) throws {
    var channel = stdin;

    if infile != defaultfilename {
      var file = open(infile, iomode.r);
      channel = file.reader();
    }
    return channel;
  }

  iter fasta_iterator(inchannel) throws {
    var sequence = "", seqlabel = "";
    var started = false;
    var line = "";

    while (inchannel.readline(line)) {
      if (line[1] == '>') {
        if (started) then
          yield (seqlabel, sequence);

        started = true;
        seqlabel = line(2..);
        sequence = "";
      } else {
        sequence += line.strip();
      }
    }

    if (started) then
      yield (seqlabel, sequence);
  }

  proc composition(sequence, k) {
    const basevalue  = ['A' => 0, 'C' => 1, 'G' => 2, 'T' => 3];
    const n = sequence.length;
    const nkmers = 4**k;
    const lastplace = 4**(k-1);
    var intseq : [1..#n] int;

    for i in 1..#n do
      intseq[i] = basevalue[sequence(i)];

    var kmer_value = 0;
    for i in 1..k do
      kmer_value = kmer_value * 4 + intseq[i];

    var kmer_counts : [0..#nkmers] int;
    kmer_counts[kmer_value] += 1;
    for i in k+1..n {
        kmer_value = (kmer_value-intseq[i-k]*lastplace)*4 + intseq[i];
        kmer_counts[kmer_value] += 1;
    }

    return kmer_counts;
  }

  proc main() {
    var inchannel = try! input_channel(infile, defaultfilename);
    try! {
        var sequences : [1..0] string;
        for (seqlabel, sequence) in fasta_iterator(inchannel) do
          writeln(composition(sequence, k));
    }
  }
}
