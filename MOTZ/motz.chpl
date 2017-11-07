module MOTZ {
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

  iter possible_partitions(seq) {
    const matchings = ['A' => 'U', 'C' => 'G', 'G' => 'C', 'U' => 'A'];
    const n = seq.length;
    const base2 = matchings[seq(n)];

    for i in 1..seq.length-1 do
      if seq(i) == base2 then
        yield i;
  }

  proc non_crossing_matchings(seq) {
    const n = seq.length;
    var nmatchings : [0..n, 0..n] int;
    nmatchings[0, 0..n] = 1;
    nmatchings[1, 0..n] = 1;

    for length in 2..n {
      for start in 1..n-length+1 {
        const part_seq = seq(start..#length);

        var nmatches = nmatchings[length-1, start];
        for partition in possible_partitions(part_seq) {
           const (left_start, left_len) = (start, partition - 1);
           const (right_start, right_len) = (start + partition, length - partition - 1);
                
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
