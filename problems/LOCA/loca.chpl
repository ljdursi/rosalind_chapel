module LOCA {
  config const defaultfilename = "none";
  config const infile = defaultfilename;
  config const scorefile = "pam250.txt";
  config const gap_penalty = -5;

  enum move {ins, del, diag, start};

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
    var moves : [0..n, 0..m] move;

    score = 0;
    moves = move.start;

    for i in 1..n {
      const basei = seq1(i);
      const basei_idx = bases[basei];
      for j in 1..m {
        const basej = seq2(j);
        const basej_idx = bases[basej];

        const match = score[i-1, j-1] + blosum[basei_idx, basej_idx];
        const ins = score[i-1, j] + gap_penalty;
        const del = score[i, j-1] + gap_penalty;
        const start = 0;

        (score[i, j], moves[i,j]) = max((ins, move.ins), (del, move.del),
                                        (match, move.diag), (start, move.start));
      }
    }

    // backtrack
    var (dist, (i, j)) = maxloc reduce zip(score, score.domain);
    var aug1, aug2: string;

    while true {
      if i <= 0 || j <= 0 then
        break;

      if moves[i, j] == move.start then
        break;

      if moves[i, j] == move.diag {
        aug1 = seq1[i] + aug1;
        aug2 = seq2[j] + aug2;
        i -= 1;
        j -= 1;
      } else if moves[i, j] == move.del {
        aug2 = seq2[j] + aug2;
        j -= 1;
      } else {
        aug1 = seq1[i] + aug1;
        i -= 1;
      }
    }

    return (dist, aug1, aug2);
  }

  proc main() {
    try! {
      const (bases, blosum) = read_table(scorefile);
      var inchannel = input_channel(infile, defaultfilename);
      var seqs : [1..0] string;
      for (seqlabel, seq) in fasta_iterator(inchannel) do
        seqs.push_back(seq.strip());

      var score = alignment(seqs[1], seqs[2], gap_penalty, bases, blosum);
      writeln(score[1]);
      writeln(score[2]);
      writeln(score[3]);
    }
  }
}
