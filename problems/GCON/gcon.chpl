module GCON {
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

    var match_score : [0..n, 0..m] int;
    var ins_score : [0..n, 0..m] int;
    var del_score : [0..n, 0..m] int;
    const minint = min(int(32));

    for i in 1..n {
      match_score[i, 0] = minint;
      ins_score[i, 0] = gap_penalty;
      del_score[i, 0] = minint;
    }

    for j in 1..m {
      match_score[0, j] = minint;
      ins_score[0, j] = minint;
      del_score[0, j] = gap_penalty;
    }

    for i in 1..n {
      const basei = seq1(i);
      const basei_idx = bases[basei];
      for j in 1..m {
        const basej = seq2(j);
        const basej_idx = bases[basej];

        match_score[i, j] = max(match_score[i-1, j-1], ins_score[i-1, j-1], del_score[i-1, j-1]);
        match_score[i, j] += blosum[basei_idx, basej_idx];

        ins_score[i, j] = max(match_score[i-1, j] + gap_penalty,
                              ins_score[i-1, j],
                              del_score[i-1, j] + gap_penalty);

        del_score[i, j] = max(match_score[i, j-1] + gap_penalty,
                              ins_score[i, j-1] + gap_penalty,
                              del_score[i, j-1]);
      }
    }

    return max(match_score[n, m], ins_score[n, m], del_score[n, m]);
  }

  proc main() {
    try! {
      const (bases, blosum) = read_table(scorefile);
      var inchannel = input_channel(infile, defaultfilename);
      var seqs : [1..0] string;
      for (seqlabel, seq) in fasta_iterator(inchannel) do
        seqs.push_back(seq.strip());

      const score = alignment(seqs[1], seqs[2], gap_penalty, bases, blosum);
      writeln(score);
    }
  }
}
