module PMCH {
  use BigInteger;
  use Assert;
  const defaultfilename="none";
  config const infile=defaultfilename;

  // use stdin if filename == defaultfilename
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

    while (inchannel.read(line)) {
      if (line[1] == '>') {
        if (started) then
          yield (seqlabel, sequence);

        started = true;
        seqlabel = line(2..);
        sequence = "";
      } else {
        sequence += line;
      }
    }

    if (started) then
      yield (seqlabel, sequence);
  }

  proc basecounts(seq: string) {
    var counts = ['A' => 0, 'C' => 0, 'G' => 0, 'T' => 0];
    for base in seq do
        counts[base] += 1;

    return counts;
  }


  proc factorial(n: int) {
    var nfact = new bigint(1);
    for i in 2..n do
        nfact *= i;

    return nfact;
  }


  proc main() {
    var inchannel = try! input_channel(infile, defaultfilename);
    try! {
        for (seqlabel, sequence) in fasta_iterator(inchannel) {
          var counts = basecounts(sequence);
          assert(counts['A'] == counts['U']);
          assert(counts['C'] == counts['G']);
          var result = factorial(counts['A'])*factorial(counts['C']);
          writeln(result);
        }
    }
  }
}
