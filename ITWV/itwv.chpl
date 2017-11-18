module ITWV {
  const defaultfilename="none";
  config const infile=defaultfilename;

  const baseidx = ['A' => 1, 'C' => 2, 'G' => 3, 'T' => 4];

  // use stdin if filename == defaultfilename
  iter lines_from_file(filename: string, defaultfilename: string) throws {
    var channel = stdin;

    if infile != defaultfilename {
      var file = open(infile, iomode.r);
      channel = file.reader();
    }

    var line: string;
    while (channel.readline(line)) do
      yield line.strip();
  }

  proc is_interleaving(superseq, seq1, seq2) {
    const n = seq1.length;
    const m = seq2.length;
    const l = superseq.length;

    if n+m != l then
      return false;

    var scores : [0..n, 0..m] int;

    for i in 1..n do
      if seq1(i) == superseq(i) then
        scores[i, 0] = scores[i-1, 0] + 1;

    for j in 1..m do
      if seq2(j) == superseq(j) then
        scores[0, j] = scores[0, j-1] + 1;

    for i in 1..n do
      for j in 1..m {
        const idx = i + j;
        var left = scores[i, j-1];
        var down = scores[i-1, j];

        if seq2(j) == superseq(idx) then
          left += 1;

        if seq1(i) == superseq(idx) then
          down += 1;

        scores[i, j] = max(left, down);
      }

    return scores[n, m] == l;
  }

  proc calc_counts(sequence) {
    const n = sequence.length;
    const nb = baseidx.size;

    var counts : [1..nb, 0..n] int;
    for (i, base) in zip(1.., sequence) {
      counts[1..nb, i] = counts[1..nb, i-1];
      counts[baseidx[base], i] = counts[baseidx[base], i-1] + 1;
    }

    return counts;
  }

  proc find_windows(sseq_counts, sseq_len, seq1, seq2) {
    const totlen = seq1.length + seq2.length;
    const nb = baseidx.size;
    var basecounts : [1..nb] int;

    for base in seq1 + seq2 do
        basecounts[baseidx[base]] += 1;

    var windows : [1..0] int;
    for idx in 1..(sseq_len-totlen+1) {
        const delta = sseq_counts[1..nb, idx+totlen-1] - sseq_counts[1..nb, idx-1];
        if basecounts.equals(delta) then
            windows.push_back(idx);
    }

    return windows;
  }

  proc contains_interleaving(sequence, counts, seq1, seq2) {
    const n = sequence.length;
    const totlen = seq1.length + seq2.length;
    const windows = find_windows(counts, n, seq1, seq2);

    for window in windows {
      const superseq = sequence(window..window+totlen-1);
      if is_interleaving(superseq, seq1, seq2) then
        return true;
    }
    return false;
  }

  proc main() {
    try! {
      var lines : [1..0] string;
      for sequence in lines_from_file(infile, defaultfilename) do
        lines.push_back(sequence);
      
      const superseq = lines[1];
      const counts = calc_counts(superseq);

      const n = lines.size - 1;
      var mat : [1..n, 1..n] int;
      for i in 1..n do 
        for j in 1..i do 
          if contains_interleaving(superseq, counts, lines[i+1], lines[j+1]) {
            mat[i, j] = 1;
            mat[j, i] = 1;
          }

      writeln(mat);
    }
  }
}
