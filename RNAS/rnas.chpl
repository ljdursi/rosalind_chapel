module MOTZ {
  use BigInteger;
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

  iter lines_iterator(inchannel) throws {
    var sequence = "";
    var line = "";

    while (inchannel.readline(line)) do
      yield line.strip();
  }

  iter possible_partitions(seq) {
    const matchings = ['A' => 'U', 'C' => 'G', 'G' => 'C', 'U' => 'A'];
    const wobble_matchings = ['G' => 'U', 'U' => 'G'];
    const n = seq.length;
    const base2 = matchings[seq(n)];
    const can_wobble = {'G', 'U'}.member(seq(n));
    var base3 = '';
    if can_wobble then
      base3 = wobble_matchings[seq(n)];

    for i in 1..seq.length-4 do
      if (seq(i) == base2) || (can_wobble && seq(i) == base3) then
        yield i;
  }

  proc non_crossing_wobble_matchings(seq) {
    const n = seq.length;
    var nmatchings : [0..n, 0..n] bigint;
    nmatchings[0, 0..n] = new bigint(1);
    nmatchings[1, 0..n] = new bigint(1);

    for length in 2..n {
      for start in 1..n-length+1 {
        const part_seq = seq(start..#length);

        var nmatches = nmatchings[length-1, start];
        for partition in possible_partitions(part_seq) {
           const (left_start, left_len) = (start, partition - 1);
           const (right_start, right_len) = (start + partition, length - partition - 1);
                
           const nleft = nmatchings[left_len, left_start];
           const nright = nmatchings[right_len, right_start];
           nmatches += (nleft * nright);
        }
        nmatchings[length, start] = nmatches;
      }
    }
    return nmatchings[n, 1];
  }

  proc main() {
    var inchannel = try! input_channel(infile, defaultfilename);
    try! {
      for sequence in lines_iterator(inchannel) do
          writeln(non_crossing_wobble_matchings(sequence));
    }
  }
}
