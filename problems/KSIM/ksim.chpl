module KSIM {
  // Finding All Similar Motifs
  // http://rosalind.info/problems/ksim/
  const defaultfilename="none";
  config const infile=defaultfilename;

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

  proc banded_alignment_mtx(ref score, const seq1, const seq2, const bandwidth) {
    // Performs a banded alignment of given bandwidth between
    // seq1 + seq2
    const n = seq1.size;
    const m = seq2.size;

    for i in 1..n {
      const basei = seq1[i];

      for j in max(1, i-bandwidth)..min(m, i+bandwidth) {
        const basej = seq2[j];
        const match = if basei == basej then score[i-1, j-1] else score[i-1, j-1]+1;
        const ins = score[i-1, j] + 1;
        const del = score[i, j-1] + 1;

        score[i, j] = min(match, ins, del);
      }
    }
  }


  proc main() {
    try! {
      var lines = lines_from_file(infile, defaultfilename);
      const k = lines[1] : int;
      const (s, t) = (lines[2], lines[3]);
      const (m, n) = (s.length, t.length);
      const motif : [1..m] int = [char in s] ascii(char);
      const sequence : [1..n] int = [char in t] ascii(char);

      // set up boundary conditions for banded alignment mtx once
      var score : [0..m, 0..m+k] int = 10*k;
      score[0, 0..k] = 0..k;
      score[1..k, 0] = 1..k;

      for i in 1..(n-m+k+1) {
        const last = min(n, i+m+k-1);
        const len = last-i+1;

        // perform a banded alignment starting at sequence(i)
        banded_alignment_mtx(score, motif, sequence[i..last].reindex(1..len), k);

        // output any full motif alignments with score
        // at or under the given threshold
        for col in max(1, len-2*k)..len do
          if score[m, col] <= k then
            writeln(i, ' ', col);
      }
    }
  }
}
