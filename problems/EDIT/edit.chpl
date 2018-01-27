module EDIT {
  const defaultfilename="none";
  config const infile=defaultfilename;

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
    if started then
      yield (seqlabel, sequence);
  }

  proc edit_dist(seq1, seq2) {
    const n = seq1.length;
    const m = seq2.length;
    var score : [0..n, 0..m] int;
    var i, j : int;

    score[0, 1..m] = 1..m;
    score[1..n, 0] = 1..n;

    // calculate subsequence lengths
    for i in 1..n {
      for j in 1..m {
        const match = score[i-1, j-1];
        const mismatch = score[i-1, j-1] + 1;
        const insertion = score[i-1, j] + 1;
        const deletion = score[i, j-1] + 1;
        if ascii(seq1[i]) == ascii(seq2[j]) then
          score[i, j] = min(match, insertion, deletion);
        else 
          score[i, j] = min(mismatch, insertion, deletion);
      }
    }

    return score[n, m];
  }

  proc main() {
    var inchannel = try! input_channel(infile, defaultfilename);
    try! {
      var sequences : [1..0] string;
      for (seqlabel, sequence) in fasta_iterator(inchannel) {
        sequences.push_back(sequence);
      }
      writeln(edit_dist(sequences[1], sequences[2]));
    }
  }
}
