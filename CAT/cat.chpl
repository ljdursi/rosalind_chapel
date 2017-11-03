module CAT {
  const defaultfilename="none";
  config const infile=defaultfilename;
  config const modulus=1000000;

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

  proc inc(base, base1, base2) {
    if base == base1 then
      return +1;
    else if base == base2 then
      return -1;
    return 0;
  }

  proc delta_base_counts(seq, base1, base2) {
    const n=seq.length;
    var deltas : [1..n] int;

    deltas[1] = inc(seq(1), base1, base2);
    for i in 2..n do
      deltas[i] = deltas[i-1] + inc(seq(i), base1, base2);
   
    return deltas;
  }

  iter possible_partitions(seq, a_minus_u, c_minus_g) {
    const matchings = ['A' => 'U', 'C' => 'G', 'G' => 'C', 'U' => 'A'];
    const base1 = seq(1);
    const base2 = matchings[base1];

    var partitions : [1..0] int;
    for i in 2..seq.length {
        if seq(i) == base2 then
            if (a_minus_u[i] == 0) && (c_minus_g[i] == 0) then
                yield i;
    }
  }

  proc non_crossing_matchings(seq) {
    const n = seq.length;
    var nmatchings : [0..n, 1..n+1] int;
    nmatchings[0, 1..n] = 1;
    nmatchings[0..n, n+1] = 1;

    var a_minus_u = delta_base_counts(seq, "A", "U");
    var c_minus_g = delta_base_counts(seq, "C", "G");

    for length in 2..n by 2 {
      for start in 1..n-length+1 {
        const part_seq = seq(start..#length);
        var part_amu : [1..length] int = a_minus_u[start..#length];
        var part_cmg : [1..length] int = c_minus_g[start..#length];
        if start > 1 {
          part_amu -= a_minus_u[start-1];
          part_cmg -= c_minus_g[start-1];
        }
 
        var nmatches = 0;
        for partition in possible_partitions(part_seq, part_amu, part_cmg) {
           const (left_start, left_len) = (start + 1, partition - 2);
           const (right_start, right_len) = (start + partition, length - partition);
                
           const nleft = nmatchings[left_len, left_start];
           const nright = nmatchings[right_len, right_start];
           nmatches += (nleft * nright) % modulus;
        }
        nmatchings[length, start] = nmatches % modulus;
      }
    }
    return nmatchings[n, 1];
  }

  proc main() {
    var inchannel = try! input_channel(infile, defaultfilename);
    try! {
      for (seqlabel, sequence) in fasta_iterator(inchannel) do
          writeln(non_crossing_matchings(sequence));
    }
  }
}
