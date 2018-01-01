module CTEA {
  const defaultfilename="none";
  config const infile=defaultfilename;
  config const modulus=134217727;

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

  proc num_alignments(seq1, seq2, modulus) {
    const n = seq1.length;
    const m = seq2.length;
    var score : [0..n, 0..m] int;
    var paths : [0..n, 0..m] int;

    score[0, 1..m] = 1..m;
    score[1..n, 0] = 1..n;
    paths[0, 0..m] = 1;
    paths[1..n, 0] = 1;

    var i, j: int;
    for i in 1..n {
      for j in 1..m {
        var match = score[i-1, j-1];
        const ins = score[i-1, j] + 1;
        const del = score[i, j-1] + 1;

        if seq1[i] != seq2[j] then
          match = score[i-1, j-1] + 1;

        const minscore = min(match, ins, del);
        score[i, j] = minscore;

        if minscore == match {
          paths[i, j] += paths[i-1, j-1];
          paths[i, j] %= modulus;
        }

        if minscore == ins {
          paths[i, j] += paths[i-1, j];
          paths[i, j] %= modulus;
        }

        if minscore == del {
          paths[i, j] += paths[i, j-1];
          paths[i, j] %= modulus;
        }
      }
    }

    return paths[n, m];
  }

  proc main() {
    try! {
      var inchannel = input_channel(infile, defaultfilename);
      var seqs : [1..0] string;
      for (seqlabel, seq) in fasta_iterator(inchannel) do
        seqs.push_back(seq.strip());

      var n = num_alignments(seqs[1], seqs[2], modulus);
      writeln(n);
    }
  }
}
