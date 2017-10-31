module REVP {
  const defaultfilename="none";
  config const infile=defaultfilename;
  config const minlen=4;
  config const maxlen=12;

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
        seqlabel = line(2..).strip();
        sequence = "";
      } else {
        sequence += line.strip();
      }
    }

    if (started) then
      yield (seqlabel, sequence);
  }

  proc complement(base: string) : string {
    const complements = ['A' => 'T', 'C' => 'G', 'G' => 'C', 'T' => 'A'];
    return complements[base];
  }

  proc reversecomplement(line: string) : string {
    var rc:string = "";
    for i in 1..line.length by -1 do
      rc += complement(line[i]);

    return rc;
  }

  proc findallminpalindromes(seq, minlen) {
    const halflen : int = minlen / 2;
    const seqlen : int = seq.length;

    var minpalindromeposns: [1..0] int;
    for i in 1..seqlen-minlen+1 {
        const kmer = seq(i..#halflen);
        const nextkmer = seq(i+halflen..#halflen);
        if reversecomplement(nextkmer) == kmer then
            minpalindromeposns.push_back(i);
    }

    return minpalindromeposns;
  }

  proc growpalindrome(seq, pos, knownlen, maxlen) {
    var palindromelen = knownlen;
    const halflen = knownlen/2;

    for newhalflen in (knownlen/2+1)..maxlen/2 {
        const first = pos + halflen - newhalflen;
        const last = pos + halflen - 1 + newhalflen;
        if first < 1 then
            break;
        if last > seq.length then
            break;
        if seq[first] != complement(seq[last]) then
            break;
        palindromelen += 2;
    }
    
    return palindromelen;
  }

  proc main() {
    var inchannel = try! input_channel(infile, defaultfilename);
    try! {
        for (seqlabel, sequence) in fasta_iterator(inchannel) {
          var minpalindromes = findallminpalindromes(sequence, minlen);
          for posn in minpalindromes {
              const palindromelength = growpalindrome(sequence, posn, minlen, maxlen);
              for len in minlen..palindromelength by -2 {
                 writeln(posn+(minlen-len)/2, " ", len);
              }
          }
        }
    }
  }
}
