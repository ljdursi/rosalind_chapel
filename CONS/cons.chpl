module CONS {
  const defaultfilename="none";
  config const infile=defaultfilename;

  const base_to_idx = ["A" => 0, "C" => 1, "G" => 2, "T" => 3];
  const idx_to_base = [0 => "A", 1 => "C", 2 => "G", 3 => "T"];

  // use stdin if filename == defaultfilename
  proc input_channel(filename: string, defaultfilename: string) throws {
    var channel = stdin;

    if infile != defaultfilename {
      try! {
        var file = open(infile, iomode.r);
        channel = file.reader();
      }
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

    if (started) then
      yield (seqlabel, sequence);
  }

  proc profile(sequences) {
    const l = sequences[1].len;
    var counts : [0..3, 1..l] uint;
    for sequence in sequences do
      for i in 1..l {
        const base = sequence[i];
        counts[base_to_idx[base], i] += 1;
      }
    return counts;
  }

  proc consensus_from_profile(profile) : string {
    var consensus = "";
    const len = profile[0, 1..].size;
    for col in 1..len {
       const (val, idx) = maxloc reduce zip(profile(0..3, col), 0..3);
       consensus += idx_to_base[idx];
    }
    return consensus;
  }

  proc main() {
    try! {
      var inchannel = input_channel(infile, defaultfilename);
      var maxlabel = "none", maxgc = -1.0;

      var sequences : [1..0] string;
      for (seqlabel, sequence) in fasta_iterator(inchannel) {
        sequences.push_back(sequence);
      }
      var counts = profile(sequences);
      var consensus = consensus_from_profile(counts);
      writeln(consensus);
      for i in 0..3 do
        writeln(idx_to_base[i], ": ", counts[i, 1..]);
    }
  }
}
