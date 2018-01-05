module LONG {
  use List;
  param defaultfilename="none";
  config const infile=defaultfilename;

  record overlap_edge {
    var seqno: int;
    var overlap: int;
  }

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

  proc overlap_graph(sequences) {
    const nsequences = sequences.size;
    const mink = (sequences[1].length+1)/2 : int;

    var graph : [1..#nsequences] list(overlap_edge);
    var nins : [1..#nsequences] int;
    var suffixdomain : domain(string);
    var suffixes : [suffixdomain] list(int);

    for (seqno, sequence) in zip(1.., sequences) {
        for k in 1..mink-1 {
            const suffix = sequence[k..];
            suffixes[suffix].append(seqno);
        }
    }

    for (seqno, sequence) in zip(1.., sequences) {
        for k in mink..sequence.length {
            const prefix = sequence[1..#k];
            if suffixdomain.member(prefix) then
                for loc in suffixes[prefix] {
                    if loc != seqno {
                        graph[loc].append(new overlap_edge(seqno, k));
                        nins[seqno] += 1;
                    }
                }
        }
    }

    var starts: [1..0] int;
    for i in 1..nsequences do
        if nins[i] == 0 then
            starts.push_back(i);

    return (graph, starts);
  }

  proc reconstruct_chromosome(graph, starts, sequences) {
    // we are assured that there is one possible reconstruction
    var current = starts[1];
    var chromosome = sequences[current];

    while graph[current].length > 0 {
       const edge = graph[current].pop_front();
       current = edge.seqno;
       chromosome += sequences[current][edge.overlap+1..];
    }

    return chromosome;
  }

  proc main() {
    try! {
      var inchannel = input_channel(infile, defaultfilename);
      var sequences : [1..0] string;
  
      for (seqlabel, sequence) in fasta_iterator(inchannel) do
        sequences.push_back(sequence);
     
      var (graph, starts) = overlap_graph(sequences);
      const chromosome = reconstruct_chromosome(graph, starts, sequences);
      writeln(chromosome);
    }
  }
}
