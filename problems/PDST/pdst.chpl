module SSEQ {
  param defaultfilename="none";
  config const infile=defaultfilename;

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
        seqlabel = line(2..).strip();
        sequence = "";
      } else {
        sequence += line.strip();
      }
    }

    if (started) then
      yield (seqlabel, sequence);
  }

  proc pdist(seq1, seq2) {
    var n = seq1.length;
    if seq2.length != n then
      return -1;

    var ndiff = 0 : uint;
    for (b1, b2) in zip(seq1, seq2) do
      if b1 != b2 then
        ndiff += 1;

    return 1.0 * ndiff / n;
  }

  proc dist_matrix(sequences) {
    const nseqs = sequences.size;
    var dmatrix : [1..#nseqs, 1..#nseqs] real;

    for i in 1..nseqs {
      for j in 1..i-1 {
        dmatrix[i, j] = pdist(sequences[i], sequences[j]);
        dmatrix[j, i] = dmatrix(i, j);
      }
    }

    return dmatrix;
  }

  proc main() {
    var inchannel = try! input_channel(infile, defaultfilename);
    try! {
        var sequences : [1..0] string;
        for (seqlabel, sequence) in fasta_iterator(inchannel) do
          sequences.push_back(sequence);

        var dmatrix = dist_matrix(sequences);
        writeln(dmatrix);
    }
  }
}
