module GLOB {
  config const defaultfilename = "none";
  config const infile = defaultfilename;
  config const scorefile = "blosum.txt";
  config const gap_penalty = -5;

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

  proc read_table(infile) throws {
    try! {
      var file = open(infile, iomode.r);
      var channel = file.reader();

      var line: string;
      channel.readline(line);

      var bases: domain(string);
      var baseorder: [bases] int;
      for (i, base) in zip(1.., line.split()) {
        bases.add(base);
        baseorder[base] = i;
      }

      const n = bases.size;
      var table : [1..n, 1..n] int;

      while (channel.readline(line)) {
        const items = line.split();
        const base = items[1];
        const i = baseorder[base];
        for (j, valstr) in zip(1.., items[2..]) {
          const val = valstr:int;
          table[i, j] = val;
        }
      }
      return (baseorder, table);
    }
  }

  proc alignment(seq1, seq2, gap_penalty, bases, blosum) {
    const n = seq1.length;
    const m = seq2.length;
    var score : [0..n, 0..m] int;

    for i in 1..n do
      score[i, 0] = i * gap_penalty;

    for j in 1..m do
      score[0, j] = j * gap_penalty;

    for i in 1..n {
      const basei = seq1(i);
      const basei_idx = bases[basei];
      for j in 1..m {
        const basej = seq2(j);
        const basej_idx = bases[basej];

        const match = score[i-1, j-1] + blosum[basei_idx, basej_idx];
        const ins = score[i-1, j] + gap_penalty;
        const del = score[i, j-1] + gap_penalty;

        score[i, j] = max(ins, del, match);
      }
    }

    return score[n, m];
  }

  proc main() {
    try! {
      const (bases, blosum) = read_table(scorefile);
      var inchannel = input_channel(infile, defaultfilename);
      var seqs : [1..0] string;
      for (seqlabel, seq) in fasta_iterator(inchannel) do
        seqs.push_back(seq.strip());

      var score = alignment(seqs[1], seqs[2], gap_penalty, bases, blosum);
      writeln(score);
    }
  }
}
