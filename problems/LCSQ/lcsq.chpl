module LCSQ {
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

  proc lcs(seq1, seq2) {
    const n = seq1.length;
    const m = seq2.length;
    var length : [0..n, 0..m] int;
    var i, j : int;

    // calculate subsequence lengths
    length = 0;
    for i in 1..n {
      for j in 1..m {
        if seq1[i] == seq2[j] then
          length[i, j] = length[i-1, j-1] + 1;
        else 
          length[i, j] = max(length[i-1, j], length[i, j-1]);
      }
    }

    // backtrack
    const l = length[n, m];
    var subsequence : string;
    var position = l;
    i = n;
    j = m;

    while (i > 0 && j > 0) {
        var up = length[i-1, j];
        var left = length[i, j-1];
        if up < position && left < position {
          subsequence = seq1[i] + subsequence;
          position -= 1;
          i -= 1;
          j -= 1;
        } else {
          if left > up then
            j -= 1;
          else
            i -= 1;
        }
    }
    return subsequence;
  }

  proc main() {
    var inchannel = try! input_channel(infile, defaultfilename);
    try! {
      var sequences : [1..0] string;
      for (seqlabel, sequence) in fasta_iterator(inchannel) {
        sequences.push_back(sequence);
      }
      writeln(lcs(sequences[1], sequences[2]));
    }
  }
}
